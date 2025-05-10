from pathlib import Path
from typing import List, Optional
import glob
from ..utils.plotting import VectorPlotter

class VectorReport:
    """Generates reports for STR vector analysis."""
    
    def __init__(self, vector_dir: str, output_dir: str):
        """Initialize the vector report generator.
        
        Args:
            vector_dir: Directory containing vector.json files
            output_dir: Directory to save report outputs
        """
        self.vector_dir = Path(vector_dir)
        self.output_dir = Path(output_dir)
        self.plotter = VectorPlotter()
        
    def generate_report(self) -> str:
        """Generate vector analysis report.
        
        Returns:
            Path to the generated report HTML file
        """
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find all vector.json files
        vector_files = list(self.vector_dir.glob("*.vector.json"))
        
        if not vector_files:
            raise ValueError(f"No vector.json files found in {self.vector_dir}")
        
        # Generate vector plot
        plot_path = self.plotter.plot_vectors(
            vector_files=[str(f) for f in vector_files],
            output_dir=str(self.output_dir)
        )
        
        return plot_path 