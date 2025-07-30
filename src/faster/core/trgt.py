import json
import subprocess
import logging
import os
import stat
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

logger = logging.getLogger(__name__)

class TRGTAnalyzer:
    """Class for running TRGT analysis and parsing VCF output."""
    
    def __init__(self):
        """Initialize the TRGTAnalyzer class."""
        pass
    
    def run_trgt_analysis(self, 
                          input_bam: str,
                          reference: str,
                          repeat_annotation_bed: str,
                          output_prefix: str,
                          threads: int = 4) -> str:
        """Run TRGT genotype analysis.
        
        Args:
            input_bam: Path to the input BAM file
            reference: Path to the reference genome fasta file
            repeat_annotation_bed: Path to the repeat annotation BED file
            output_prefix: Output prefix for TRGT results
            threads: Number of threads to use (default: 4)
            
        Returns:
            Path to the generated VCF file (gzipped)
            
        Raises:
            FileNotFoundError: If TRGT binary or input files not found
            RuntimeError: If TRGT analysis fails
        """
        try:
            # Set paths for TRGT binary
            package_dir = Path(__file__).parent  # src/faster/core/
            trgt_binary = package_dir.parent / 'bin' / 'trgt'
            
            # Ensure TRGT binary is executable
            if trgt_binary.exists():
                os.chmod(trgt_binary, os.stat(trgt_binary).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            
            # Verify binary and input files exist
            if not trgt_binary.exists():
                raise FileNotFoundError(f"TRGT binary not found at: {trgt_binary}")
            if not Path(input_bam).exists():
                raise FileNotFoundError(f"Input BAM file not found: {input_bam}")
            if not Path(reference).exists():
                raise FileNotFoundError(f"Reference genome not found: {reference}")
            if not Path(repeat_annotation_bed).exists():
                raise FileNotFoundError(f"Repeat annotation BED file not found: {repeat_annotation_bed}")
            
            # Extract directory from output_prefix if it contains a path
            output_prefix_path = Path(output_prefix)
            if '/' in str(output_prefix_path) or '\\' in str(output_prefix_path):
                output_dir = output_prefix_path.parent
                output_dir.mkdir(parents=True, exist_ok=True)
            
            # Prepare command
            cmd = [
                str(trgt_binary),
                'genotype',
                '-g', str(reference),
                '-r', str(input_bam),
                '-b', str(repeat_annotation_bed),
                '-t', str(threads),
                '-o', str(output_prefix)
            ]
            
            # Run TRGT
            logger.info("Running TRGT genotype analysis...")
            logger.debug(f"Command: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            logger.info(f"TRGT stdout:\n{result.stdout}")
            logger.info(f"TRGT stderr:\n{result.stderr}")
            
            if result.returncode != 0:
                logger.error(f"TRGT failed with exit code {result.returncode}")
                raise RuntimeError(f"TRGT analysis failed. See logs above for details.")
            
            # Check if VCF file was generated (TRGT generates gzipped VCF by default)
            vcf_file = f"{output_prefix}.vcf.gz"
            if not Path(vcf_file).exists():
                # Try uncompressed VCF as fallback
                vcf_file = f"{output_prefix}.vcf"
                if not Path(vcf_file).exists():
                    raise RuntimeError(f"TRGT VCF file not generated: {vcf_file}")
            
            logger.info("TRGT analysis completed successfully")
            logger.info(f"Results saved to: {vcf_file}")
            
            return vcf_file
            
        except FileNotFoundError as e:
            logger.error(f"Error: {str(e)}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error: {str(e)}")
            raise

    def convert_vcf_to_exhunter_json(self, vcf_file: str, sample_id: str = None) -> Dict:
        """Convert TRGT VCF file to ExpansionHunter JSON format.
        
        Args:
            vcf_file: Path to the TRGT VCF file
            sample_id: Sample ID to use in the output (default: extracted from filename)
            
        Returns:
            Dictionary in ExpansionHunter JSON format
        """
        try:
            # Parse VCF file
            parsed_data = self.parse_trgt_vcf(vcf_file)
            
            # Extract sample ID if not provided
            if sample_id is None:
                sample_id = self._extract_sample_id_from_vcf(vcf_file, parsed_data)
            
            # Group variants by TRID (marker)
            markers_data = self._group_variants_by_marker(parsed_data['variants'])
            
            # Convert to ExpansionHunter format
            exhunter_data = {
                "LocusResults": {},
                "SampleParameters": {
                    "SampleId": sample_id,
                    "Sex": ""  # TRGT doesn't provide sex information
                }
            }
            
            # Process each marker
            for marker_name, variants in markers_data.items():
                exhunter_data["LocusResults"][marker_name] = self._convert_marker_to_exhunter(
                    marker_name, variants
                )
            
            logger.info(f"Successfully converted {len(markers_data)} markers to ExpansionHunter format")
            return exhunter_data
            
        except Exception as e:
            logger.error(f"Error converting VCF to ExpansionHunter JSON: {str(e)}")
            raise

    def _extract_sample_id_from_vcf(self, vcf_file: str, parsed_data: Dict) -> str:
        """Extract sample ID from VCF file or filename.
        
        Args:
            vcf_file: Path to the VCF file
            parsed_data: Parsed VCF data
            
        Returns:
            Sample ID string
        """
        # Try to get sample ID from VCF headers
        if 'headers' in parsed_data and len(parsed_data['headers']) > 9:
            sample_id = parsed_data['headers'][9]  # First sample column
            if sample_id and sample_id != '.':
                return sample_id
        
        # Fallback to filename
        filename = Path(vcf_file).stem
        if filename.endswith('.vcf'):
            filename = filename[:-4]
        return filename

    def _group_variants_by_marker(self, variants: List[Dict]) -> Dict[str, List[Dict]]:
        """Group variants by TRID (marker name).
        
        Args:
            variants: List of variant dictionaries
            
        Returns:
            Dictionary mapping marker names to lists of variants
        """
        markers = {}
        for variant in variants:
            trid = variant['INFO'].get('TRID')
            if trid:
                if trid not in markers:
                    markers[trid] = []
                markers[trid].append(variant)
        
        return markers

    def _convert_marker_to_exhunter(self, marker_name: str, variants: List[Dict]) -> Dict:
        """Convert a single marker's variants to ExpansionHunter format.
        
        Args:
            marker_name: Name of the marker
            variants: List of variants for this marker
            
        Returns:
            Dictionary in ExpansionHunter marker format
        """
        # Calculate coverage from SD (spanning reads) field
        total_coverage = 0
        for variant in variants:
            sample_data = variant.get('SAMPLE', {})
            sd_values = sample_data.get('SD', [])
            
            # Handle tuple format from pysam
            if isinstance(sd_values, tuple):
                sd_values = list(sd_values)
            
            if isinstance(sd_values, list):
                total_coverage += sum(int(x) for x in sd_values if str(x).isdigit())
            elif isinstance(sd_values, str):
                # Handle comma-separated string
                for val in sd_values.split(','):
                    if val.strip().isdigit():
                        total_coverage += int(val.strip())
        
        # Get first variant for basic info
        first_variant = variants[0]
        
        # Create ExpansionHunter format marker
        exhunter_marker = {
            "AlleleCount": 2,  # Default for diploid
            "Coverage": total_coverage,
            "FragmentLength": "",  # Not available in TRGT
            "LocusId": marker_name,
            "ReadLength": "",  # Not available in TRGT
            "Variants": {}
        }
        
        # Convert each variant
        for i, variant in enumerate(variants):
            variant_id = self._generate_variant_id(variant)
            exhunter_marker["Variants"][variant_id] = self._convert_variant_to_exhunter(variant)
        
        return exhunter_marker

    def _generate_variant_id(self, variant: Dict) -> str:
        """Generate variant ID in ExpansionHunter format.
        
        Args:
            variant: Variant dictionary
            
        Returns:
            Variant ID string
        """
        chrom = variant['CHROM']
        pos = variant['POS']
        return f"{chrom}_{pos}_{variant['INFO'].get('END', pos)}"

    def _convert_variant_to_exhunter(self, variant: Dict) -> Dict:
        """Convert a single variant to ExpansionHunter format.
        
        Args:
            variant: Variant dictionary
            
        Returns:
            Dictionary in ExpansionHunter variant format
        """
        sample_data = variant.get('SAMPLE', {})
        
        # Extract genotype
        genotype = self._extract_genotype(sample_data, variant)
        
        # Extract repeat unit from MOTIFS
        motifs = variant['INFO'].get('MOTIFS', '')
        # Handle both string and tuple formats
        if isinstance(motifs, tuple):
            repeat_unit = motifs[0] if motifs else ""
        elif isinstance(motifs, str):
            repeat_unit = motifs.split(',')[0] if motifs else ""
        else:
            repeat_unit = ""
        
        # Generate reference region
        chrom = variant['CHROM']
        pos = variant['POS']
        end = variant['INFO'].get('END', pos)
        reference_region = f"{chrom}:{pos}-{end}"
        
        # Create ExpansionHunter variant format
        exhunter_variant = {
            "CountsOfFlankingReads": "()",  # Not available in TRGT
            "CountsOfInrepeatReads": "()",  # Not available in TRGT
            "CountsOfSpanningReads": self._format_spanning_reads(sample_data),
            "Genotype": genotype,
            "GenotypeConfidenceInterval": self._generate_confidence_interval(sample_data),
            "ReferenceRegion": reference_region,
            "RepeatUnit": repeat_unit,
            "VariantId": self._generate_variant_id(variant),
            "VariantSubtype": "Repeat",
            "VariantType": "Repeat"
        }
        
        return exhunter_variant

    def _extract_genotype(self, sample_data: Dict, variant_info: Dict = None) -> str:
        """Extract genotype from TRGT sample data.
        
        Args:
            sample_data: Sample data dictionary
            
        Returns:
            Genotype string in format "allele1/allele2"
        """
        gt = sample_data.get('GT', '')
        al = sample_data.get('AL', [])
        
        if gt == '.' or not gt:
            return "INVALID_GT"  # Changed from "N/A" to be more explicit
        
        # Handle tuple format from pysam
        if isinstance(al, tuple):
            al = list(al)
        
        # If we have allele lengths, convert to repeat counts
        if al and isinstance(al, list) and len(al) >= 2:
            try:
                allele1 = self._convert_allele_length_to_repeat_count(al[0], sample_data, variant_info, allele_index=0)
                allele2 = self._convert_allele_length_to_repeat_count(al[1], sample_data, variant_info, allele_index=1)
                return f"{allele1:.1f}/{allele2:.1f}"
            except (ValueError, TypeError) as e:
                logger.warning(f"Error converting allele lengths to repeat counts: {str(e)}")
                # Fallback to original allele lengths
                return f"{al[0]}/{al[1]}"
        
        # Fallback to GT field
        if '/' in gt:
            return gt
        else:
            return f"{gt}/{gt}"  # Homozygous

    def _convert_allele_length_to_repeat_count(self, allele_length: int, sample_data: Dict, variant_info: Dict = None, allele_index: int = 0) -> float:
        """Convert allele length to repeat count.
        
        Args:
            allele_length: Allele length from AL field
            sample_data: Sample data dictionary containing MC and other fields
            variant_info: Variant info dictionary containing STRUC field
            allele_index: Index of the allele (0 for first, 1 for second)
            
        Returns:
            Repeat count as float
        """
        # First try to use MC (Motif Counts) field if available
        mc = sample_data.get('MC', [])
        
        # Handle tuple format from pysam
        if isinstance(mc, tuple):
            mc = list(mc)
        
        if mc:
            try:
                return self._parse_motif_counts(mc, allele_length, allele_index)
            except Exception as e:
                logger.warning(f"Error parsing MC field: {str(e)}")
        
        # Try to use STRUC field if available
        if variant_info:
            struc = variant_info.get('INFO', {}).get('STRUC', '')
            if struc:
                try:
                    return self._parse_structure_to_repeat_count(struc, allele_length)
                except Exception as e:
                    logger.warning(f"Error parsing STRUC field: {str(e)}")
        
        # Fallback to calculation based on MOTIFS field
        if variant_info:
            motifs = variant_info.get('INFO', {}).get('MOTIFS', '')
            if motifs:
                # Handle both string and tuple formats
                if isinstance(motifs, tuple):
                    motif = motifs[0] if motifs else ""
                elif isinstance(motifs, str):
                    motif = motifs.split(',')[0] if motifs else ""
                else:
                    motif = ""
                
                if motif and motif != '.':
                    motif_length = len(motif)
                    if motif_length > 0:
                        return float(allele_length) / motif_length
        
        # Final fallback to simple calculation based on average motif length
        return float(allele_length) / 4.0  # Assume average motif length of 4

    def _parse_motif_counts(self, mc_data, allele_length: int, allele_index: int = 0) -> float:
        """Parse MC (Motif Counts) field to calculate total repeat count.
        
        Args:
            mc_data: MC field data (can be list or string)
            allele_length: Allele length for validation
            allele_index: Index of the allele (0 for first, 1 for second)
            
        Returns:
            Total repeat count as float
        """
        if isinstance(mc_data, str):
            # Handle comma-separated string (multiple alleles)
            mc_parts = mc_data.split(',')
            # Use specified allele index
            mc_str = mc_parts[allele_index] if allele_index < len(mc_parts) else mc_parts[0] if mc_parts else ""
        elif isinstance(mc_data, list):
            # Use specified allele index
            mc_str = str(mc_data[allele_index]) if allele_index < len(mc_data) else str(mc_data[0]) if mc_data else ""
        else:
            mc_str = str(mc_data)
        
        if not mc_str or mc_str == '.':
            return float(allele_length) / 4.0
        
        # Parse underscore-separated motif counts
        # Example: "3_0_0_0_0_16" -> [3, 0, 0, 0, 0, 16]
        try:
            motif_counts = [float(x) for x in mc_str.split('_')]
            total_repeats = sum(motif_counts)
            return total_repeats
        except (ValueError, TypeError):
            logger.warning(f"Could not parse MC field: {mc_str}")
            return float(allele_length) / 4.0

    def _parse_structure_to_repeat_count(self, struc: str, allele_length: int) -> float:
        """Parse STRUC field to calculate repeat count.
        
        Args:
            struc: STRUC field string (e.g., "(GGAA)n(GGAG)n(AAAG)n(AGAA)n(AAAA)n(GAAA)n")
            allele_length: Allele length for validation
            
        Returns:
            Total repeat count as float
        """
        if not struc or struc == '.':
            return float(allele_length) / 4.0
        
        try:
            # Parse structure like "(GGAA)n(GGAG)n(AAAG)n(AGAA)n(AAAA)n(GAAA)n"
            # Extract motif patterns and their repeat counts
            import re
            
            # Find all patterns like (MOTIF)n
            pattern = r'\(([A-Z]+)\)n'
            matches = re.findall(pattern, struc)
            
            if matches:
                # For now, use a simplified calculation
                # In a more sophisticated implementation, we would parse the exact repeat counts
                total_motifs = len(matches)
                avg_motif_length = sum(len(motif) for motif in matches) / total_motifs
                estimated_repeats = allele_length / avg_motif_length
                return estimated_repeats
            else:
                # Fallback to simple calculation
                return float(allele_length) / 4.0
                
        except Exception as e:
            logger.warning(f"Error parsing STRUC field: {str(e)}")
            return float(allele_length) / 4.0

    def _format_spanning_reads(self, sample_data: Dict) -> str:
        """Format spanning reads information.
        
        Args:
            sample_data: Sample data dictionary
            
        Returns:
            Formatted spanning reads string
        """
        sd = sample_data.get('SD', [])
        al = sample_data.get('AL', [])
        
        if not sd or not al:
            return "()"
        
        try:
            # Handle tuple format from pysam
            if isinstance(sd, tuple):
                sd = list(sd)
            if isinstance(al, tuple):
                al = list(al)
            
            if isinstance(sd, list) and isinstance(al, list):
                pairs = []
                for i, (span_count, allele_len) in enumerate(zip(sd, al)):
                    if i < 2:  # Only first two alleles
                        pairs.append(f"({allele_len}, {span_count})")
                return ", ".join(pairs)
            elif isinstance(sd, str) and isinstance(al, str):
                sd_vals = sd.split(',')
                al_vals = al.split(',')
                pairs = []
                for i, (span_count, allele_len) in enumerate(zip(sd_vals, al_vals)):
                    if i < 2:  # Only first two alleles
                        pairs.append(f"({allele_len.strip()}, {span_count.strip()})")
                return ", ".join(pairs)
        except (ValueError, IndexError):
            pass
        
        return "()"

    def _generate_confidence_interval(self, sample_data: Dict) -> str:
        """Generate genotype confidence interval.
        
        Args:
            sample_data: Sample data dictionary
            
        Returns:
            Confidence interval string
        """
        al = sample_data.get('AL', [])
        allr = sample_data.get('ALLR', [])
        
        if not al:
            return "N/A-N/A/N/A-N/A"
        
        try:
            # Handle tuple format from pysam
            if isinstance(al, tuple):
                al = list(al)
            if isinstance(allr, tuple):
                allr = list(allr)
            
            if isinstance(al, list) and len(al) >= 2:
                allele1 = al[0]
                allele2 = al[1]
                
                # Try to get ranges from ALLR
                if allr and isinstance(allr, list) and len(allr) >= 2:
                    range1 = allr[0].split('-')[0] if '-' in str(allr[0]) else str(allele1)
                    range2 = allr[1].split('-')[0] if '-' in str(allr[1]) else str(allele2)
                    return f"{range1}-{range1}/{range2}-{range2}"
                else:
                    return f"{allele1}-{allele1}/{allele2}-{allele2}"
        except (ValueError, IndexError):
            pass
        
        return "N/A-N/A/N/A-N/A"

    def _is_gzipped(self, file_path: str) -> bool:
        """Check if a file is gzip compressed.
        
        Args:
            file_path: Path to the file
            
        Returns:
            True if the file is gzip compressed, False otherwise
        """
        try:
            with open(file_path, 'rb') as f:
                return f.read(2) == b'\x1f\x8b'
        except Exception:
            return False
    
    def _open_file(self, file_path: str):
        """Open a file, handling both regular and gzipped files.
        
        Args:
            file_path: Path to the file
            
        Returns:
            File object (regular or gzipped)
        """
        if self._is_gzipped(file_path):
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')
    
    def parse_trgt_vcf(self, vcf_file: str) -> Dict:
        """Parse TRGT VCF file and extract basic information.
        
        Args:
            vcf_file: Path to the TRGT VCF file
            
        Returns:
            Dictionary containing parsed VCF information
        """
        try:
            if not Path(vcf_file).exists():
                raise FileNotFoundError(f"VCF file not found: {vcf_file}")
            
            # Use manual parsing
            return self._parse_vcf_manual(vcf_file)
                
        except Exception as e:
            logger.error(f"Error parsing VCF file {vcf_file}: {str(e)}")
            raise
    
    def _parse_vcf_manual(self, vcf_file: str) -> Dict:
        """Parse VCF file manually (fallback method).
        
        Args:
            vcf_file: Path to the VCF file
            
        Returns:
            Dictionary containing parsed VCF information
        """
        parsed_data = {
            'header': {},
            'variants': []
        }
        
        with self._open_file(vcf_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Parse header lines
                if line.startswith('##'):
                    if line.startswith('##INFO='):
                        # Parse INFO field definitions
                        info_match = re.search(r'##INFO=<ID=(\w+),.*Description="([^"]*)"', line)
                        if info_match:
                            field_name = info_match.group(1)
                            description = info_match.group(2)
                            parsed_data['header'][f'INFO_{field_name}'] = description
                    elif line.startswith('##FORMAT='):
                        # Parse FORMAT field definitions
                        format_match = re.search(r'##FORMAT=<ID=(\w+),.*Description="([^"]*)"', line)
                        if format_match:
                            field_name = format_match.group(1)
                            description = format_match.group(2)
                            parsed_data['header'][f'FORMAT_{field_name}'] = description
                    continue
                
                # Parse column headers
                if line.startswith('#CHROM'):
                    headers = line.split('\t')
                    parsed_data['headers'] = headers
                    continue
                
                # Parse variant lines
                if not line.startswith('#'):
                    fields = line.split('\t')
                    if len(fields) >= 9:  # Minimum required fields
                        variant = {
                            'CHROM': fields[0],
                            'POS': int(fields[1]),
                            'ID': fields[2],
                            'REF': fields[3],
                            'ALT': fields[4].split(',') if fields[4] != '.' else [],
                            'QUAL': fields[5],
                            'FILTER': fields[6].split(';') if fields[6] != '.' else [],
                            'INFO': self._parse_info_field(fields[7]),
                            'FORMAT': fields[8].split(':'),
                            'SAMPLE': self._parse_sample_field(fields[8], fields[9]) if len(fields) > 9 else {}
                        }
                        parsed_data['variants'].append(variant)
        
        logger.info(f"Successfully parsed {len(parsed_data['variants'])} variants from {vcf_file} using manual parsing")
        return parsed_data
    
    def _parse_info_field(self, info_str: str) -> Dict:
        """Parse INFO field string into dictionary.
        
        Args:
            info_str: INFO field string from VCF
            
        Returns:
            Dictionary of INFO field key-value pairs
        """
        info_dict = {}
        if info_str == '.':
            return info_dict
        
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
            else:
                info_dict[item] = True
        
        return info_dict
    
    def _parse_sample_field(self, format_str: str, sample_str: str) -> Dict:
        """Parse SAMPLE field string into dictionary.
        
        Args:
            format_str: FORMAT field string from VCF
            sample_str: SAMPLE field string from VCF
            
        Returns:
            Dictionary of sample field key-value pairs
        """
        sample_dict = {}
        if sample_str == '.' or format_str == '.':
            return sample_dict
        
        format_fields = format_str.split(':')
        sample_values = sample_str.split(':')
        
        for i, field in enumerate(format_fields):
            if i < len(sample_values):
                value = sample_values[i]
                if value == '.':
                    sample_dict[field] = None
                else:
                    # Try to convert to appropriate type
                    try:
                        if ',' in value:
                            # Multiple values (e.g., for alleles)
                            sample_dict[field] = value.split(',')
                        elif value.replace('.', '').replace('-', '').isdigit():
                            # Numeric value
                            sample_dict[field] = float(value)
                        else:
                            # String value
                            sample_dict[field] = value
                    except ValueError:
                        sample_dict[field] = value
        
        return sample_dict
    
    def get_basic_vcf_info(self, vcf_file: str) -> Dict:
        """Get basic information from TRGT VCF file.
        
        Args:
            vcf_file: Path to the TRGT VCF file
            
        Returns:
            Dictionary containing basic VCF information
        """
        try:
            parsed_data = self.parse_trgt_vcf(vcf_file)
            
            basic_info = {
                'total_variants': len(parsed_data['variants']),
                'markers': [],
                'sample_id': None
            }
            
            # Extract unique markers
            markers = set()
            for variant in parsed_data['variants']:
                trid = variant['INFO'].get('TRID')
                if trid:
                    markers.add(trid)
            
            basic_info['markers'] = sorted(list(markers))
            
            # Try to extract sample ID from first variant
            if parsed_data['variants']:
                sample_data = parsed_data['variants'][0]['SAMPLE']
                # Sample ID might be in the filename or other fields
                # For now, we'll use a placeholder
                basic_info['sample_id'] = 'sample'
            
            return basic_info
            
        except Exception as e:
            logger.error(f"Error getting basic VCF info: {str(e)}")
            raise 