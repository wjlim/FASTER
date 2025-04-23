import pandas as pd
import numpy as np
from typing import Optional, Dict, List, Tuple
from ..models.data_models import ContaminationInfo, PeakInfo
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

class ContaminationDetector:
    def __init__(self):
        """Initialize contamination detector."""
        self.MAX_MAIN_PEAKS = 4  # Maximum number of peaks in main profile
        self.HEIGHT_RATIO_THRESHOLD = 0.2  # Minimum height ratio for significant peaks
        self.CLUSTER_DISTANCE_THRESHOLD = 0.2  # Threshold for determining separate clusters
    
    def detect_contamination(self, peaks: pd.DataFrame) -> Optional[ContaminationInfo]:
        """
        Detect contamination in peaks using enhanced clustering.
        
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
        
        # If more than 2 peaks, perform clustering analysis
        if len(peaks) > 2:
            return self._analyze_peaks(peaks)
        
        return None
    
    def _analyze_peaks(self, peaks: pd.DataFrame) -> Optional[ContaminationInfo]:
        """Analyze peaks to determine contamination."""
        max_height = peaks['height'].iloc[0]
        
        # If we have 3-4 peaks, check if they're close enough to be main profile
        if len(peaks) <= self.MAX_MAIN_PEAKS:
            height_ratios = peaks['height'] / max_height
            
            # If all peaks are relatively close in height (>20% of max)
            if all(ratio >= self.HEIGHT_RATIO_THRESHOLD for ratio in height_ratios):
                # Check if peaks form a tight cluster
                if self._check_cluster_tightness(peaks):
                    # All peaks are main profile
                    return self._create_contamination_info(peaks, pd.DataFrame())
        
        # For more peaks or if peaks don't form tight cluster
        # First identify main profile peaks (up to 4)
        X = peaks[['height']].values  # Using only height for clustering
        X_scaled = StandardScaler().fit_transform(X)
        
        # Try clustering with 2 clusters
        kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X_scaled)
        
        # Calculate cluster centers and distances
        cluster_centers = kmeans.cluster_centers_
        cluster_distance = abs(cluster_centers[0] - cluster_centers[1])[0]
        
        # If clusters are well separated
        if cluster_distance > self.CLUSTER_DISTANCE_THRESHOLD:
            # Identify main cluster (higher peaks)
            clusters = pd.DataFrame({'height': peaks['height'], 'cluster': labels})
            cluster_means = clusters.groupby('cluster')['height'].mean()
            main_cluster_id = cluster_means.idxmax()
            
            main_peaks = peaks[labels == main_cluster_id]
            contamination_peaks = peaks[labels != main_cluster_id]
            
            # Ensure main peaks don't exceed MAX_MAIN_PEAKS
            if len(main_peaks) > self.MAX_MAIN_PEAKS:
                # Move excess peaks to contamination
                contamination_peaks = pd.concat([
                    contamination_peaks,
                    main_peaks.iloc[self.MAX_MAIN_PEAKS:]
                ])
                main_peaks = main_peaks.iloc[:self.MAX_MAIN_PEAKS]
        else:
            # If clusters aren't well separated, take top 4 as main peaks
            main_peaks = peaks.iloc[:self.MAX_MAIN_PEAKS]
            contamination_peaks = peaks.iloc[self.MAX_MAIN_PEAKS:]
            
            # Only keep contamination peaks if they're significantly lower
            contamination_peaks = contamination_peaks[
                contamination_peaks['height'] < main_peaks['height'].min() * 0.5
            ]
            
            if contamination_peaks.empty:
                return self._create_contamination_info(main_peaks, pd.DataFrame())
        
        return self._create_contamination_info(main_peaks, contamination_peaks)
    
    def _check_cluster_tightness(self, peaks: pd.DataFrame) -> bool:
        """Check if peaks form a tight cluster based on heights."""
        heights = peaks['height'].values
        height_diffs = np.diff(heights)
        max_height = heights[0]
        
        # Check if consecutive peaks have similar heights
        return all(abs(diff) < max_height * 0.4 for diff in height_diffs)
    
    def _create_contamination_info(self, main_cluster: pd.DataFrame, contamination_cluster: pd.DataFrame) -> ContaminationInfo:
        """Create ContaminationInfo object from cluster data."""
        if main_cluster.empty:
            return None
            
        max_height = main_cluster['height'].iloc[0]
        
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
                relative_height=round(float(row['height']) / max_height * 100, 2)
            )
            for _, row in contamination_cluster.iterrows()
        ]
        
        # Calculate relative distance only if there are contamination peaks
        if not contamination_cluster.empty:
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