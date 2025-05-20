from pathlib import Path
from typing import List, Optional, Dict, Tuple
import glob
import json
import pandas as pd
import numpy as np
from ..utils.plotting import VectorPlotter

class VectorReport:
    """Generates reports for STR vector analysis."""
    
    def __init__(self, vector_dir: str, output_dir: str):
        """Initialize the vector report generator.
        
        Args:
            vector_dir: Directory containing json files
            output_dir: Directory to save report outputs
        """
        self.vector_dir = Path(vector_dir)
        self.output_dir = Path(output_dir)
        self.plotter = VectorPlotter()
        
    def _calculate_distance(self, coord1: Dict[str, float], coord2: Dict[str, float]) -> float:
        """Calculate Euclidean distance between two coordinate points.
        
        Args:
            coord1: First coordinate point (x, y)
            coord2: Second coordinate point (x, y)
            
        Returns:
            Euclidean distance between points
        """
        return np.sqrt(
            (coord1['x'] - coord2['x'])**2 + 
            (coord1['y'] - coord2['y'])**2
        )
        
    def _find_closest_pairs(self, vector_data: List[Dict]) -> List[Dict]:
        """Find closest pairs between STR and EH samples.
        
        Args:
            vector_data: List of vector data dictionaries
            
        Returns:
            List of dictionaries containing closest pair information
        """
        # Separate STR and EH data
        str_data = [d for d in vector_data if d.get('source_type') == 'str']
        eh_data = [d for d in vector_data if d.get('source_type') == 'eh']
        
        closest_pairs = []
        
        # Calculate distances between all STR-EH pairs
        for str_sample in str_data:
            str_coords = str_sample['cartesian_coordinates']
            str_id = str_sample['sample_id']
            
            # Find closest EH sample
            min_distance = float('inf')
            closest_eh = None
            
            for eh_sample in eh_data:
                eh_coords = eh_sample['cartesian_coordinates']
                eh_id = eh_sample['sample_id']
                
                distance = self._calculate_distance(str_coords, eh_coords)
                
                if distance < min_distance:
                    min_distance = distance
                    closest_eh = {
                        'sample_id': eh_id,
                        'distance': distance,
                        'str_coords': str_coords,
                        'eh_coords': eh_coords
                    }
            
            if closest_eh:
                closest_pairs.append({
                    'str_sample': str_id,
                    'eh_sample': closest_eh['sample_id'],
                    'distance': closest_eh['distance'],
                    'str_coords': closest_eh['str_coords'],
                    'eh_coords': closest_eh['eh_coords']
                })
        
        # Sort by distance
        closest_pairs.sort(key=lambda x: x['distance'])
        return closest_pairs
        
    def _save_closest_pairs(self, closest_pairs: List[Dict], output_dir: Path):
        """Save closest pairs information to CSV and JSON files.
        
        Args:
            closest_pairs: List of closest pair information
            output_dir: Output directory
        """
        # Save as CSV
        df = pd.DataFrame(closest_pairs)
        df = df[['str_sample', 'eh_sample', 'distance']]  # Select relevant columns
        df.columns = ['STR Sample', 'EH Sample', 'Distance']
        csv_path = output_dir / 'closest_pairs.csv'
        df.to_csv(csv_path, index=False)
        
        # Save as JSON
        json_path = output_dir / 'closest_pairs.json'
        with open(json_path, 'w') as f:
            json.dump(closest_pairs, f, indent=2)
            
        return csv_path, json_path
        
    def generate_report(self) -> str:
        """Generate vector analysis report.
        
        Returns:
            Path to the generated report HTML file
        """
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find all json files
        vector_files = list(self.vector_dir.glob("*.json"))
        
        if not vector_files:
            raise ValueError(f"No json files found in {self.vector_dir}")
            
        # Load vector data
        vector_data = []
        for file_path in vector_files:
            with open(file_path) as f:
                data = json.load(f)
                data['filename'] = file_path.stem
                vector_data.append(data)
                
        # Find closest pairs
        closest_pairs = self._find_closest_pairs(vector_data)
        
        # Save closest pairs information
        csv_path, json_path = self._save_closest_pairs(closest_pairs, self.output_dir)
        
        # Generate vector plot
        plot_path = self.plotter.plot_vectors(
            vector_files=[str(f) for f in vector_files],
            output_dir=str(self.output_dir)
        )
        
        return plot_path 