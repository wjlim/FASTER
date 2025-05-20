import pandas as pd
import numpy as np
import json
import os
from typing import Optional, Dict, List, Tuple, NamedTuple
from ..models.data_models import ContaminationInfo, PeakInfo

class JoinPoint(NamedTuple):
    """Represents a joining point in neighbor joining."""
    size: float
    height: float
    joined_peaks: List[int]  # indices of peaks that were joined

class ContaminationDetector:
    def __init__(self, config_path=None, threshold: Optional[float] = None):
        """Initialize contamination detector with new parameters and dye cutoffs."""
        self.MAX_MAIN_PEAKS = 4  # Maximum number of peaks in main profile
        self.RELATIVE_DISTANCE_THRESHOLD = 0.3  # Maximum relative distance for main profile peaks
        self.MAX_relative_peak_height = 0.5 # Maximum relative height for peaks in main profile
        # Load dye cutoffs from marker_info.json
        if config_path is None:
            config_path = os.path.join(os.path.dirname(__file__), '../config/marker_info.json')
        with open(config_path) as f:
            config = json.load(f)
            self.dye_cutoffs = config['dye_cutoffs']
        self.threshold = threshold

    def _get_min_cutoff(self, dye):
        return self.dye_cutoffs.get(dye, {}).get('min', 0)

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

    def _neighbor_joining(self, peaks: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, List[JoinPoint], Optional[float]]:
        """
        Implement sequential neighbor joining for peaks.
        Reference peak updates sequentially.
        Ensures that the top two peaks are always included in main_peaks_df.
        
        Args:
            peaks: DataFrame of peaks sorted by height
            
        Returns:
            Tuple of (main_peaks, contamination_peaks, join_points, triggering_distance)
            triggering_distance is the relative distance that caused the first contamination, or None.
        """
        peaks['min_cutoff'] = peaks['dye'].apply(self._get_min_cutoff)
        if self.threshold is not None:
            peaks['min_cutoff'] = peaks['height'].max() * self.threshold
        peaks = peaks[peaks['height'] > peaks['min_cutoff']]
        if len(peaks) <= 1:
            return peaks, pd.DataFrame(), [], None
        
        join_points = []
        max_peak = peaks.iloc[0]
        main_peaks_list = [max_peak]  # Start with highest peak
        current_ref_idx = 0
        triggering_distance = None
        
        # Always include the second peak in main_peaks_list if it exists
        if len(peaks) > 1:
            main_peaks_list.append(peaks.iloc[1])
            # Calculate joining point between first and second peak
            join_x, join_y = self._find_joining_point(peaks.iloc[0], peaks.iloc[1])
            join_points.append(JoinPoint(
                size=join_x,
                height=join_y,
                joined_peaks=[0, 1]
            ))
            current_ref_idx = 1
        
        # Process remaining peaks sequentially, starting from the third
        for i in range(2, len(peaks)):
            current_ref_peak = peaks.iloc[current_ref_idx]
            next_peak = peaks.iloc[i]
            dye = next_peak['dye'] if 'dye' in next_peak else 'B'
            min_cutoff = next_peak['min_cutoff']  # Get the min_cutoff for the current peak
            relative_height_diff = (next_peak['height'] - min_cutoff) / (max_peak['height'] - min_cutoff)
            # Calculate relative distance to the *current reference* peak
            current_ref_height = current_ref_peak['height']
            if current_ref_height == 0:
                relative_distance = float('inf')
            else:
                relative_distance = 1 - (next_peak['height'] / current_ref_height)
            
            # Check if peak should be added to main cluster
            if relative_distance <= self.RELATIVE_DISTANCE_THRESHOLD and len(main_peaks_list) < self.MAX_MAIN_PEAKS and relative_height_diff >= self.MAX_relative_peak_height:
                main_peaks_list.append(next_peak)
                join_x, join_y = self._find_joining_point(current_ref_peak, next_peak)
                join_points.append(JoinPoint(
                    size=join_x,
                    height=join_y,
                    joined_peaks=[current_ref_idx, i]
                ))
                current_ref_idx = i
            else:
                triggering_distance = relative_distance
                break
        
        main_peaks_df = pd.DataFrame(main_peaks_list)
        contamination_start_idx = len(main_peaks_list)
        contamination_peaks_df = peaks.iloc[contamination_start_idx:].copy()
        
        return main_peaks_df, contamination_peaks_df, join_points, triggering_distance

    def detect_contamination(self, peaks: pd.DataFrame) -> Tuple[Optional[ContaminationInfo], List[JoinPoint]]:
        """
        Detect contamination using sequential neighbor joining approach.
        
        Args:
            peaks: DataFrame containing peak information
            
        Returns:
            Tuple of (ContaminationInfo, join_points)
        """
        if len(peaks) <= 2:
             # No contamination if 2 or fewer peaks, pass empty contamination df
            contamination_info = self._create_contamination_info(peaks, pd.DataFrame(), None)
            return contamination_info, []
            
        # Sort peaks by height in descending order
        peaks = peaks.sort_values('height', ascending=False).reset_index(drop=True)
        
        # Calculate relative heights (used in _create_contamination_info)
        max_height = peaks['height'].iloc[0]
        if max_height > 0:
            peaks['relative_height'] = (peaks['height'] / max_height)
        else:
            peaks['relative_height'] = 0.0
        
        # Perform neighbor joining
        main_peaks, contamination_peaks, join_points, triggering_distance = self._neighbor_joining(peaks)
        
        # Create contamination info, passing the triggering distance
        contamination_info = self._create_contamination_info(main_peaks, contamination_peaks, triggering_distance)
        
        return contamination_info, join_points
    
    def _create_contamination_info(self, main_cluster: pd.DataFrame, 
                                 contamination_cluster: pd.DataFrame, 
                                 triggering_distance: Optional[float]) -> ContaminationInfo:
        """Create ContaminationInfo object from peak data, including triggering distance."""
        is_contaminated = not contamination_cluster.empty
        
        # Handle case where main_cluster might be empty (though unlikely with current logic)
        if main_cluster.empty:
             # If main cluster is somehow empty but contamination exists, treat all as contamination?
             # Or return None? Returning None for now.
             if is_contaminated: return None
             # If both are empty, return non-contaminated info
             else: return ContaminationInfo(is_contaminated=False, main_profile_peaks=[], contamination_peaks=[], relative_distance=0)

            
        max_height = main_cluster['height'].iloc[0]
        # Avoid division by zero if max_height is 0
        if max_height == 0: max_height = 1.0
        
        # Create peak info objects for main profile
        main_peaks = [
            PeakInfo(
                allele=str(row['allele']),
                height=float(row['height']),
                size=float(row['size']),
                # Use pre-calculated relative_height if available, else calculate
                relative_height=float(row.get('relative_height', row['height'] / max_height)) * 100
            )
            for _, row in main_cluster.iterrows()
        ]
        
        # Create peak info objects for contamination
        contamination_peaks = [
            PeakInfo(
                allele=str(row['allele']),
                height=float(row['height']),
                size=float(row['size']),
                relative_height=float(row.get('relative_height', row['height'] / max_height)) * 100
            )
            for _, row in contamination_cluster.iterrows()
        ]
        
        # Use triggering_distance if contamination occurred, otherwise 0
        effective_relative_distance = round(triggering_distance, 2) if is_contaminated and triggering_distance is not None else 0
        
        return ContaminationInfo(
            is_contaminated=is_contaminated,
            main_profile_peaks=main_peaks,
            contamination_peaks=contamination_peaks,
            # Store the distance that *triggered* contamination
            relative_distance=effective_relative_distance 
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
        """Format peaks for output, omitting relative_height from JSON serialization."""
        return [
            {
                'allele': str(row['allele']),
                'height': float(row['height']),
                'size': float(row['size'])
                # 'relative_height' intentionally omitted for JSON output
            }
            for _, row in peaks.iterrows()
        ] 