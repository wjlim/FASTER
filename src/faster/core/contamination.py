import pandas as pd
import numpy as np
from typing import Optional, Dict, List, Tuple, NamedTuple
from ..models.data_models import ContaminationInfo, PeakInfo

class JoinPoint(NamedTuple):
    """Represents a joining point in neighbor joining."""
    size: float
    height: float
    joined_peaks: List[int]  # indices of peaks that were joined

class ContaminationDetector:
    def __init__(self):
        """Initialize contamination detector with new parameters."""
        self.MAX_MAIN_PEAKS = 4  # Maximum number of peaks in main profile
        self.RELATIVE_DISTANCE_THRESHOLD = 0.3  # Maximum relative distance for main profile peaks
        
    def _calculate_distance_matrix(self, peaks: pd.DataFrame) -> np.ndarray:
        """Calculate distance matrix between peaks based on relative height differences."""
        n = len(peaks)
        D = np.zeros((n, n))
        max_height = peaks['height'].iloc[0]
        
        for i in range(n):
            for j in range(i+1, n):
                # Calculate relative distance based on height differences
                rel_dist = abs(1 - (peaks['height'].iloc[j] / peaks['height'].iloc[i]))
                D[i,j] = D[j,i] = rel_dist
                
        return D
        
    def _find_joining_point(self, peak1: pd.Series, peak2: pd.Series) -> Tuple[float, float]:
        """Calculate the joining point between two peaks."""
        # Calculate joining point coordinates
        x1, y1 = peak1['size'], peak1['height']
        x2, y2 = peak2['size'], peak2['height']
        
        # Join at the weighted average point
        weight1 = y1 / (y1 + y2)
        weight2 = y2 / (y1 + y2)
        
        join_x = x1 * weight1 + x2 * weight2
        join_y = y1 * weight1 + y2 * weight2
        
        return join_x, join_y

    def _neighbor_joining(self, peaks: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, List[JoinPoint]]:
        """
        Implement neighbor joining for peaks.
        
        Args:
            peaks: DataFrame of peaks sorted by height
            
        Returns:
            Tuple of (main_peaks, contamination_peaks, join_points)
        """
        if len(peaks) <= 2:
            return peaks, pd.DataFrame(), []
            
        join_points = []
        current_peaks = peaks.copy()
        main_peaks = peaks.iloc[[0]].copy()  # Start with highest peak
        remaining_idx = list(range(1, len(peaks)))
        
        # Process up to 4 peaks
        while len(main_peaks) < self.MAX_MAIN_PEAKS and remaining_idx:
            # Get next peak
            next_idx = remaining_idx[0]
            next_peak = peaks.iloc[next_idx]
            
            # Calculate relative distance to reference peak
            ref_peak = peaks.iloc[0]  # Always use highest peak as reference
            relative_distance = 1 - (next_peak['height'] / ref_peak['height'])
            
            # If distance is within threshold, add to main peaks
            if relative_distance <= self.RELATIVE_DISTANCE_THRESHOLD:
                # Calculate joining point
                join_x, join_y = self._find_joining_point(ref_peak, next_peak)
                join_points.append(JoinPoint(
                    size=join_x,
                    height=join_y,
                    joined_peaks=[0, next_idx]
                ))
                
                main_peaks = pd.concat([main_peaks, next_peak.to_frame().T])
                remaining_idx.pop(0)
            else:
                # If distance exceeds threshold, stop joining
                break
        
        # All remaining peaks are contamination
        contamination_peaks = peaks.iloc[remaining_idx] if remaining_idx else pd.DataFrame()
        
        return main_peaks, contamination_peaks, join_points

    def detect_contamination(self, peaks: pd.DataFrame) -> Tuple[Optional[ContaminationInfo], List[JoinPoint]]:
        """
        Detect contamination using neighbor joining approach.
        
        Args:
            peaks: DataFrame containing peak information
            
        Returns:
            Tuple of (ContaminationInfo, join_points)
        """
        if len(peaks) <= 2:
            return self._create_contamination_info(peaks, pd.DataFrame()), []
            
        # Sort peaks by height in descending order
        peaks = peaks.sort_values('height', ascending=False).reset_index(drop=True)
        
        # Calculate relative heights
        max_height = peaks['height'].iloc[0]
        peaks['relative_height'] = (peaks['height'] / max_height)
        
        # Perform neighbor joining
        main_peaks, contamination_peaks, join_points = self._neighbor_joining(peaks)
        
        # Create contamination info
        contamination_info = self._create_contamination_info(main_peaks, contamination_peaks)
        
        return contamination_info, join_points
    
    def _create_contamination_info(self, main_cluster: pd.DataFrame, contamination_cluster: pd.DataFrame) -> ContaminationInfo:
        """Create ContaminationInfo object from peak data."""
        if main_cluster.empty:
            return None
            
        max_height = main_cluster['height'].iloc[0]
        
        # Create peak info objects for main profile
        main_peaks = [
            PeakInfo(
                allele=str(row['allele']),
                height=float(row['height']),
                size=float(row['size']),
                relative_height=round(float(row['height']) / max_height * 100, 2)
            )
            for _, row in main_cluster.iterrows()
        ]
        
        # Create peak info objects for contamination
        contamination_peaks = [
            PeakInfo(
                allele=str(row['allele']),
                height=float(row['height']),
                size=float(row['size']),
                relative_height=round(float(row['height']) / max_height * 100, 2)
            )
            for _, row in contamination_cluster.iterrows()
        ]
        
        # Calculate relative distance as mean height ratio of contamination to main peaks
        if not contamination_cluster.empty and not main_cluster.empty:
            relative_distance = round(
                contamination_cluster['height'].mean() / main_cluster['height'].mean(),
                2
            )
        else:
            relative_distance = 0
        
        return ContaminationInfo(
            is_contaminated=len(contamination_peaks) > 0,
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