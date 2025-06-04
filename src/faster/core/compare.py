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
                      tolerance: float = 0.5) -> float:
        """Check if alleles match within tolerance, considering rounding,
        and also check for single allele exact matches.
        
        Args:
            str_alleles: Tuple of alleles from STR analysis (sorted)
            eh_alleles: Tuple of alleles from ExpansionHunter (sorted)
            tolerance: Maximum difference allowed for full pair match
            
        Returns:
            1.0 if alleles form a full pair match (exact, rounded, ceiling, or tolerance).
            0.5 if one allele from STR exactly matches one from ExpansionHunter.
            0.0 otherwise.
        """
        # Check for full pair match (score 1.0)
        # Exact match
        if str_alleles == eh_alleles:
            return 1.0
            
        # Rounded match
        str_rounded = tuple(sorted([round(x) for x in str_alleles]))
        eh_rounded = tuple(sorted([round(x) for x in eh_alleles]))
        if str_rounded == eh_rounded:
            return 1.0
            
        # Ceiling match
        str_ceiling = tuple(sorted([math.ceil(x) for x in str_alleles]))
        eh_ceiling = tuple(sorted([math.ceil(x) for x in eh_alleles]))
        if str_ceiling == eh_ceiling:
            return 1.0
            
        # Tolerance match for the pair
        if (abs(str_alleles[0] - eh_alleles[0]) <= tolerance and 
            abs(str_alleles[1] - eh_alleles[1]) <= tolerance):
            return 1.0

        # Check for single exact allele match (score 0.5)
        # Alleles within tuples are already sorted by _parse_genotype and _combine_exhunter_genotypes
        if (str_alleles[0] == eh_alleles[0] or
            str_alleles[0] == eh_alleles[1] or
            str_alleles[1] == eh_alleles[0] or
            str_alleles[1] == eh_alleles[1]):
            return 0.5
            
        # No match
        return 0.0

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
                "total_markers": 0, # Markers from compare_markers list found in STR data
                "matching": 0.0,    # Sum of match scores (0.5 or 1.0)
                "mismatching": 0, # Count of 0-score markers
                "missing": 0      # Markers from compare_markers not in EH data
            }
        }

        # Compare only markers in compare_markers
        for marker, str_info in str_data["LocusResults"].items():
            if marker not in compare_markers:
                continue
            
            # Increment total_markers only if the marker from STR data is in our comparison list
            comparison_results["summary"]["total_markers"] += 1

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
            str_variant = list(str_info["variants"].values())[0]
            str_genotype = str_variant.get("genotype", "N/A")
            
            try:
                str_alleles = self._parse_genotype(str_genotype)
            except (ValueError, TypeError):
                logger.warning(f"Could not parse STR genotype for {marker}: {str_genotype}")
                # This counts as a mismatch as comparison cannot proceed reliably
                comparison_results["mismatching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": "N/A", # EH genotype not parsed yet / irrelevant if STR is bad
                    "error": "Invalid STR genotype format"
                })
                comparison_results["summary"]["mismatching"] += 1
                continue

            try:
                eh_alleles = self._combine_exhunter_genotypes(eh_info["Variants"])
            except (KeyError, ValueError) as e:
                logger.warning(f"Could not parse ExpansionHunter genotype for {marker}: {str(e)}")
                # This counts as a mismatch
                comparison_results["mismatching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": "N/A",
                    "error": f"Invalid ExpansionHunter data: {str(e)}"
                })
                comparison_results["summary"]["mismatching"] += 1
                continue

            # Compare genotypes and get score
            match_score = self._match_alleles(str_alleles, eh_alleles)

            if match_score > 0:
                comparison_results["matching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": f"{eh_alleles[0]}/{eh_alleles[1]}", # Store combined EH genotype
                    "match_score": match_score 
                })
                comparison_results["summary"]["matching"] += match_score
            else: # match_score is 0
                comparison_results["mismatching_markers"].append({
                    "marker": marker,
                    "str_genotype": str_genotype,
                    "eh_genotype": f"{eh_alleles[0]}/{eh_alleles[1]}" # Store combined EH genotype
                })
                comparison_results["summary"]["mismatching"] += 1
                
        return comparison_results

    def save_results(self, results: Dict, output_prefix: str):
        """Save comparison results to JSON file and generate a text summary.
        
        Args:
            results: Comparison results dictionary
            output_prefix: Output file prefix
        """
        markers_compared_for_score = results['summary']['total_markers'] - results['summary']['missing']
        
        sum_of_match_scores = results['summary']['matching']

        if markers_compared_for_score > 0:
            match_ratio = (sum_of_match_scores / markers_compared_for_score) * 100
        else:
            match_ratio = 0.0
            
        results['summary']['match_ratio'] = round(match_ratio, 2)  # as percentage
        results['summary']['markers_compared_for_score'] = markers_compared_for_score # For clarity in JSON
        
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
            f.write(f"Total markers from input list (in STR data): {results['summary']['total_markers']}\n")
            f.write(f"Markers missing in ExpansionHunter: {results['summary']['missing']}\n")
            f.write(f"Markers effectively compared: {markers_compared_for_score}\n") 
            f.write(f"Sum of match scores (0.5 partial, 1.0 full): {results['summary']['matching']}\n") 
            f.write(f"Completely mismatching markers (0 score): {results['summary']['mismatching']}\n") 
            f.write(f"Overall match ratio: {results['summary']['match_ratio']}%\n\n")
            
            if results['matching_markers']:
                f.write("Matching/Partially Matching Markers (see .json for scores):\n")
                for m in results['matching_markers']:
                    score = m.get('match_score', 'N/A') 
                    f.write(f"{m['marker']}: STR={m['str_genotype']}, EH={m['eh_genotype']}, Score={score}\n")
                f.write("\n")

            if results['mismatching_markers']:
                f.write("Completely Mismatching Markers (0 score):\n") 
                for mm in results['mismatching_markers']:
                    error_info = f", Error: {mm['error']}" if 'error' in mm else ""
                    f.write(f"{mm['marker']}: STR={mm['str_genotype']}, EH={mm['eh_genotype']}{error_info}\n")
                f.write("\n")
                
            if results['missing_markers']:
                f.write("Markers not available in ExpansionHunter (not scored):\n") 
                for mm in results['missing_markers']:
                    f.write(f"{mm['marker']}: STR={mm['str_genotype']}\n")
                    
        logger.info(f"Summary report saved to: {summary_path}") 