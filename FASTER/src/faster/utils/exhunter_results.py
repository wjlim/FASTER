import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)

class ExhunterResultGenerator:
    """Class for generating standardized results from ExpansionHunter output."""
    
    def __init__(self, marker_config_file: Optional[str] = None):
        """Initialize the ExhunterResultGenerator.
        
        Args:
            marker_config_file: Optional path to marker configuration JSON file
        """
        self.marker_config = {}
        if marker_config_file:
            with open(marker_config_file, 'r') as f:
                self.marker_config = json.load(f)

    def _get_motif(self, marker: str) -> str:
        """Get the repeat motif for a marker.
        
        Args:
            marker: Marker name
            
        Returns:
            Repeat motif string
        """
        if marker in self.marker_config:
            return self.marker_config[marker].get('motif', '')
        return ''

    def _convert_to_csv_data(self, results: Dict) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
        """Convert results dictionary to DataFrames for CSV output.
        
        Args:
            results: Results dictionary
            
        Returns:
            Tuple of (marker_results_df, None) - contamination not applicable for EH
        """
        # Extract sample info
        sample_id = results['SampleParameters']['sample_id']
        
        # Prepare data for marker results DataFrame
        marker_data = []
        
        for marker, info in results['LocusResults'].items():
            if 'Variants' not in info:
                continue
                
            for variant_id, variant in info['Variants'].items():
                marker_data.append({
                    'Sample_ID': sample_id,
                    'Marker': marker,
                    'Variant_ID': variant_id,
                    'Genotype': variant.get('Genotype', ''),
                    'RepeatUnit': variant.get('RepeatUnit', ''),
                    'RepeatSize': variant.get('RepeatSize', ''),
                    'Motif': self._get_motif(marker)
                })
        
        marker_results_df = pd.DataFrame(marker_data)
        
        return marker_results_df, None

    def generate_results(self, eh_json: Dict, sample_id: str) -> Dict:
        """Generate standardized results dictionary from ExpansionHunter output.
        
        Args:
            eh_json: ExpansionHunter JSON results
            sample_id: Sample identifier
            
        Returns:
            Standardized results dictionary
        """
        results = {
            'LocusResults': {},
            'SampleParameters': {
                'sample_id': sample_id,
                'analysis_date': eh_json.get('SampleParameters', {}).get('AnalysisDate', ''),
                'gender': eh_json.get('SampleParameters', {}).get('Gender', '')
            }
        }
        
        # Process each locus
        for locus in eh_json.get('LocusResults', {}):
            variants = eh_json['LocusResults'][locus].get('Variants', {})
            if not variants:
                continue
                
            results['LocusResults'][locus] = {
                'Variants': variants,
                'motif': self._get_motif(locus)
            }
        
        return results

    def save_results(self, results: Dict, output_prefix: str):
        """Save results in both JSON and CSV formats.
        
        Args:
            results: Results dictionary
            output_prefix: Output file prefix
        """
        # Save JSON results
        json_path = f"{output_prefix}.STR_analysis.json"
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Saved JSON results to: {json_path}")
        
        # Convert to CSV format
        marker_results_df, _ = self._convert_to_csv_data(results)
        
        # Save marker results CSV
        csv_path = f"{output_prefix}.marker_results.csv"
        marker_results_df.to_csv(csv_path, index=False)
        logger.info(f"Saved marker results CSV to: {csv_path}") 