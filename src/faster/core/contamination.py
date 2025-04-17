import pandas as pd
from typing import Optional, Dict, List
from ..models.data_models import ContaminationInfo, PeakInfo

class ContaminationDetector:
    def __init__(self):
        """Initialize contamination detector."""
        pass
    
    def detect_contamination(self, peaks: pd.DataFrame) -> Optional[ContaminationInfo]:
        """
        Detect contamination in peaks.
        
        Args:
            peaks (pd.DataFrame): DataFrame containing peak information
            
        Returns:
            Optional[ContaminationInfo]: Contamination information if detected
        """
        if len(peaks) <= 2:
            return None
        
        # Sort peaks by height in descending order and reset index
        peaks = peaks.sort_values('height', ascending=False).reset_index(drop=True)
        
        # Calculate relative heights
        max_height = peaks['height'].iloc[0]
        peaks['relative_height'] = round((peaks['height'] / max_height) * 100, 2)
        
        # Split into main and contamination clusters
        main_cluster = peaks.iloc[:2]
        contamination_cluster = peaks.iloc[2:]
        
        # Calculate relative distance
        relative_distance = round(
            contamination_cluster['height'].mean() / main_cluster['height'].mean(),
            2
        )
        
        return self._create_contamination_info(main_cluster, contamination_cluster, relative_distance)
    
    def _create_contamination_info(self, main_cluster: pd.DataFrame, contamination_cluster: pd.DataFrame,
                                 relative_distance: float) -> ContaminationInfo:
        """Create ContaminationInfo object from cluster data."""
        # Create peak info objects
        main_peaks = [
            PeakInfo(
                allele=str(row['allele']),
                height=float(row['height']),
                size=float(row['size']),
                relative_height=float(row['relative_height'])
            )
            for _, row in main_cluster.iterrows()
        ]
        
        contamination_peaks = [
            PeakInfo(
                allele=str(row['allele']),
                height=float(row['height']),
                size=float(row['size']),
                relative_height=float(row['relative_height'])
            )
            for _, row in contamination_cluster.iterrows()
        ]
        
        return ContaminationInfo(
            is_contaminated=True,
            main_profile_peaks=main_peaks,
            contamination_peaks=contamination_peaks,
            relative_distance=relative_distance
        )

    def _check_contamination(self, peaks: pd.DataFrame, primary_peaks: pd.DataFrame) -> Dict:
        """
        Simple contamination detection based on peak count.
        If there are 3 or more peaks, consider it as contamination.
        """
        if peaks.empty or len(peaks) < 3:  # No contamination if less than 3 peaks
            return None
        
        # Sort peaks by height in descending order
        sorted_peaks = peaks.sort_values('height', ascending=False)
        max_height = sorted_peaks['height'].iloc[0]
        
        # First two peaks are considered main profile
        main_cluster = sorted_peaks.iloc[:2]
        # Remaining peaks are considered contamination
        contamination_cluster = sorted_peaks.iloc[2:]
        
        # Calculate relative distance using mean heights
        relative_distance = round(
            contamination_cluster['height'].mean() / main_cluster['height'].mean(),
            2
        )
        
        return {
            'is_contaminated': True,
            'main_profile_peaks': self._format_peaks(main_cluster),
            'contamination_peaks': self._format_peaks(contamination_cluster),
            'relative_distance': relative_distance
        }

    def _format_peaks(self, peaks: pd.DataFrame) -> List[Dict]:
        """Format peaks for output."""
        max_height = peaks['height'].max()
        return [
            {
                'allele': str(row['allele']),
                'height': float(row['height']),
                'size': float(row['size']),
                'relative_height': round(row['height'] / max_height * 100, 2)
            }
            for _, row in peaks.iterrows()
        ] 