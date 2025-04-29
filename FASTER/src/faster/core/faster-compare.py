import json
from typing import Dict, Tuple
import csv
import sys

def parse_genotype(genotype: str) -> Tuple[float, float]:
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

def match_alleles(str_alleles: Tuple[float, float], 
                  eh_alleles: Tuple[float, float], 
                  tolerance: float = 0.5) -> bool:
    """Check if alleles match within tolerance.
    
    Args:
        str_alleles: Tuple of alleles from STR analysis
        eh_alleles: Tuple of alleles from ExpansionHunter
        tolerance: Maximum difference allowed for match
        
    Returns:
        True if alleles match, False otherwise
    """
    return (abs(str_alleles[0] - eh_alleles[0]) <= tolerance and 
            abs(str_alleles[1] - eh_alleles[1]) <= tolerance)

def fast_compare(str_json_path: str, eh_json_path: str) -> Dict:
    """Fast comparison of STR analysis and ExpansionHunter results.
    Only uses markers that exist in both files.
    
    Args:
        str_json_path: Path to STR analysis JSON file
        eh_json_path: Path to ExpansionHunter JSON file
        
    Returns:
        Dictionary with comparison statistics
    """
    # Load JSON files
    with open(str_json_path) as f:
        str_data = json.load(f)
    with open(eh_json_path) as f:
        eh_data = json.load(f)
        
    total_markers = 0
    matched = 0
    mismatched = 0
    
    # Get common markers
    str_markers = set(str_data["LocusResults"].keys())
    eh_markers = set(eh_data["LocusResults"].keys())
    common_markers = str_markers.intersection(eh_markers)
    
    # Compare only common markers
    for marker in common_markers:
        str_info = str_data["LocusResults"][marker]
        eh_info = eh_data["LocusResults"][marker]
        
        # Get STR genotype
        str_variant = list(str_info["variants"].values())[0]
        str_genotype = str_variant["genotype"]
        try:
            str_alleles = parse_genotype(str_genotype)
        except ValueError:
            # Skip if genotype cannot be parsed
            continue
        
        # Get ExpansionHunter genotype
        eh_variants = eh_info.get("Variants", {})
        if not eh_variants:
            continue
            
        total_eh_allele1 = 0.0
        total_eh_allele2 = 0.0
        valid_eh_genotype = False
        
        for variant in eh_variants.values():
            if 'Genotype' not in variant:
                continue
            try:
                allele1, allele2 = parse_genotype(variant['Genotype'])
                total_eh_allele1 += allele1
                total_eh_allele2 += allele2
                valid_eh_genotype = True
            except ValueError:
                continue
                
        if not valid_eh_genotype:
            continue
            
        eh_alleles = tuple(sorted([total_eh_allele1, total_eh_allele2]))
        
        # Compare genotypes
        total_markers += 1
        if match_alleles(str_alleles, eh_alleles):
            matched += 1
        else:
            mismatched += 1
    
    # Calculate percentages only if we have valid markers
    if total_markers > 0:
        mismatch_percent = (mismatched / total_markers * 100)
        match_percent = (matched / total_markers * 100)
    else:
        mismatch_percent = 0
        match_percent = 0
    
    return {
        "total_markers": total_markers,
        "matched": matched,
        "mismatched": mismatched,
        "mismatch_percent": round(mismatch_percent, 2),
        "match_percent": round(match_percent, 2)
    }

def main():
    if len(sys.argv) != 4:
        print("Usage: python faster-compare.py <str_json> <eh_json> <output_csv>")
        sys.exit(1)
        
    str_json_path = sys.argv[1]
    eh_json_path = sys.argv[2]
    output_csv = sys.argv[3]
    
    # Run comparison
    results = fast_compare(str_json_path, eh_json_path)
    
    # Write CSV
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['total_markers', 'matched', 'mismatched', 'mismatch_percent', 'match_percent'])
        writer.writerow([
            results['total_markers'],
            results['matched'],
            results['mismatched'],
            results['mismatch_percent'],
            results['match_percent']
        ])

if __name__ == "__main__":
    main() 