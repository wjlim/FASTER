import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
import math
import os

logger = logging.getLogger(__name__)

class ResultComparator:
    def __init__(self):
        """Initialize the ResultComparator class."""
        pass

    def _parse_genotype(self, genotype: str) -> Tuple[float, float]:
        """Parse genotype string into two allele numbers.
        
        Args:
            genotype: Genotype string in format "allele1/allele2" or single "allele"
            
        Returns:
            Tuple of (allele1, allele2) as floats
        """
        if '/' in genotype:
            alleles = genotype.split('/')
            allele1 = float(alleles[0])
            allele2 = float(alleles[1])
        else:
            # Single allele case - homozygous
            allele1 = allele2 = float(genotype)
            
        return tuple(sorted([allele1, allele2]))

    def _combine_exhunter_genotypes(self, variants: Dict) -> Tuple[float, float]:
        """Combine multiple ExpansionHunter variant genotypes into a single genotype.
        
        Args:
            variants: Dictionary of variants from ExpansionHunter
            
        Returns:
            Combined genotype as (allele1, allele2)
        """
        total_allele1 = 0.0
        total_allele2 = 0.0
        
        for variant in variants.values():
            if 'Genotype' not in variant:
                continue
                
            allele1, allele2 = self._parse_genotype(variant['Genotype'])
            total_allele1 += allele1
            total_allele2 += allele2
            
        return tuple(sorted([total_allele1, total_allele2]))

    def _match_alleles(self, str_alleles: Tuple[float, float], 
                      eh_alleles: Tuple[float, float], 
                      tolerance: float = 0.5) -> bool:
        """Check if alleles match within tolerance, considering rounding.
        
        Args:
            str_alleles: Tuple of alleles from STR analysis
            eh_alleles: Tuple of alleles from ExpansionHunter
            tolerance: Maximum difference allowed for match
            
        Returns:
            True if alleles match, False otherwise
        """
        # Try exact match
        if str_alleles == eh_alleles:
            return True
            
        # Try rounding
        str_rounded = tuple(sorted([round(x) for x in str_alleles]))
        eh_rounded = tuple(sorted([round(x) for x in eh_alleles]))
        if str_rounded == eh_rounded:
            return True
            
        # Try ceiling
        str_ceiling = tuple(sorted([math.ceil(x) for x in str_alleles]))
        eh_ceiling = tuple(sorted([math.ceil(x) for x in eh_alleles]))
        if str_ceiling == eh_ceiling:
            return True
            
        # Check within tolerance
        return (abs(str_alleles[0] - eh_alleles[0]) <= tolerance and 
                abs(str_alleles[1] - eh_alleles[1]) <= tolerance)

    def compare_results(self, str_json_path: str, eh_json_path: str) -> Dict:
        """Compare STR analysis and ExpansionHunter results using only compare_markers from marker_info.json."""
        # Load JSON files
        with open(str_json_path) as f:
            str_data = json.load(f)
        with open(eh_json_path) as f:
            eh_data = json.load(f)

        # Get sample ID from STR data
        sample_id = str_data.get("SampleParameters", {}).get("SampleId")
        if not sample_id:
            # Try to get sample ID from filename if not found in data
            sample_id = Path(str_json_path).stem.split('.')[0]
            logger.warning(f"SampleId not found in STR data, using filename: {sample_id}")

        # Load compare_markers from marker_info.json
        config_path = os.path.join(os.path.dirname(__file__), '../config/marker_info.json')
        with open(config_path) as f:
            config = json.load(f)
            compare_markers = set(config.get('compare_markers', []))

        comparison_results = {
            "sample_id": sample_id,
            "matching_markers": [],
            "mismatching_markers": [],
            "missing_markers": [],
            "summary": {
                "total_markers": 0,
                "matching": 0,
                "mismatching": 0,
                "missing": 0
            }
        }

        # Compare only markers in compare_markers
        for marker, str_info in str_data["LocusResults"].items():
            if marker not in compare_markers:
                continue
            comparison_results["summary"]["total_markers"] += 1
            # Check if marker exists in ExpansionHunter results
            if marker not in eh_data["LocusResults"]:
                str_variant = list(str_info["variants"].values())[0]
                str_genotype = str_variant.get("genotype", "N/A")
                comparison_results["missing_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype
                })
                comparison_results["summary"]["missing"] += 1
                continue

            eh_info = eh_data["LocusResults"][marker]
            # Get STR genotype
            str_variant = list(str_info["variants"].values())[0]
            str_genotype = str_variant.get("genotype", "N/A")
            try:
                str_alleles = self._parse_genotype(str_genotype)
            except (ValueError, TypeError):
                logger.warning(f"Could not parse STR genotype for {marker}: {str_genotype}")
                comparison_results["mismatching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": "N/A",
                    "error": "Invalid STR genotype format"
                })
                comparison_results["summary"]["mismatching"] += 1
                continue

            # Get ExpansionHunter combined genotype
            try:
                eh_alleles = self._combine_exhunter_genotypes(eh_info["Variants"])
            except (KeyError, ValueError) as e:
                logger.warning(f"Could not parse ExpansionHunter genotype for {marker}: {str(e)}")
                comparison_results["mismatching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": "N/A",
                    "error": f"Invalid ExpansionHunter data: {str(e)}"
                })
                comparison_results["summary"]["mismatching"] += 1
                continue

            # Compare genotypes
            if self._match_alleles(str_alleles, eh_alleles):
                comparison_results["matching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": f"{eh_alleles[0]}/{eh_alleles[1]}"
                })
                comparison_results["summary"]["matching"] += 1
            else:
                comparison_results["mismatching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": f"{eh_alleles[0]}/{eh_alleles[1]}"
                })
                comparison_results["summary"]["mismatching"] += 1

        return comparison_results

    def save_results(self, results: Dict, output_prefix: str):
        """Save comparison results to JSON file.
        
        Args:
            results: Comparison results dictionary
            output_prefix: Output file prefix
        """
        # Calculate match ratio excluding missing markers
        total_compared = results['summary']['matching'] + results['summary']['mismatching']
        match_ratio = results['summary']['matching'] / total_compared if total_compared > 0 else 0
        results['summary']['match_ratio'] = round(match_ratio * 100, 2)  # as percentage
        
        output_path = Path(f"{output_prefix}.comparison.json")
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Comparison results saved to: {output_path}")
        
        # Generate summary report
        summary_path = Path(f"{output_prefix}.comparison.txt")
        with open(summary_path, 'w') as f:
            f.write(f"Comparison Results for {results['sample_id']}\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("Summary:\n")
            f.write(f"Total markers analyzed: {total_compared}\n")
            f.write(f"Matching: {results['summary']['matching']}\n")
            f.write(f"Mismatching: {results['summary']['mismatching']}\n")
            f.write(f"Match ratio: {results['summary']['match_ratio']}%\n")
            f.write(f"Missing in ExpansionHunter: {results['summary']['missing']}\n\n")
            
            if results['mismatching_markers']:
                f.write("Mismatching Markers:\n")
                for mm in results['mismatching_markers']:
                    f.write(f"{mm['marker']}: STR={mm['str_genotype']}, EH={mm['eh_genotype']}\n")
                f.write("\n")
                
            if results['missing_markers']:
                f.write("Markers not available in ExpansionHunter:\n")
                for mm in results['missing_markers']:
                    f.write(f"{mm['marker']}: STR={mm['str_genotype']}\n")
                    
        logger.info(f"Summary report saved to: {summary_path}") 