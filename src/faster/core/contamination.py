import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from typing import Optional, Dict, List, Tuple
from ..models.data_models import ContaminationInfo, PeakInfo

class ContaminationDetector:
    def __init__(self):
        """Initialize contamination detector."""
        pass
    
    def detect_contamination(self, peaks: pd.DataFrame) -> Optional[ContaminationInfo]:
        """
        Detect contamination in peaks using clustering.
        
        Args:
            peaks (pd.DataFrame): DataFrame containing peak information
            
        Returns:
            Optional[ContaminationInfo]: Contamination information if detected
        """
        # Return None if there are 2 or fewer peaks
        if len(peaks) <= 2:
            return None
        
        # Sort peaks by height in descending order and reset index
        peaks = peaks.sort_values('height', ascending=False).reset_index(drop=True)
        
        # Calculate relative heights
        max_height = peaks['height'].iloc[0]
        peaks['relative_height'] = round((peaks['height'] / max_height) * 100, 2)
        
        # Prepare data for clustering
        X = peaks[['height', 'size']].values
        
        # Normalize features for clustering
        X_scaled = StandardScaler().fit_transform(X)
        
        # For small number of peaks, use a simpler approach
        if len(peaks) < 4:
            return self._simple_detection(peaks)
        else:
            return self._clustering_detection(peaks, X_scaled)
    
    def _simple_detection(self, peaks: pd.DataFrame) -> Optional[ContaminationInfo]:
        """Simple detection for small number of peaks."""
        # Calculate height differences between consecutive peaks
        height_diffs = peaks['height'].diff()
        significant_drops = height_diffs[height_diffs < -peaks['height'].iloc[0] * 0.3]
        
        if not significant_drops.empty:
            # Get the first significant drop index
            split_idx = significant_drops.index[0]
            
            # Split peaks into two groups using integer indexing
            cluster1 = peaks.iloc[:split_idx]
            cluster2 = peaks.iloc[split_idx:]
            
            if not cluster1.empty and not cluster2.empty:
                relative_distance = round(cluster2['height'].mean() / cluster1['height'].mean(), 2)
                
                # Check if the ratio indicates potential contamination
                if 0.1 <= relative_distance <= 0.8:
                    return self._create_contamination_info(cluster1, cluster2, relative_distance)
        
        return None
    
    def _clustering_detection(self, peaks: pd.DataFrame, X_scaled: np.ndarray) -> Optional[ContaminationInfo]:
        """Detection using KMeans clustering."""
        # Try different numbers of clusters
        best_inertia_ratio = float('inf')
        best_clusters = None
        
        for n_clusters in range(2, min(4, len(peaks))):
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            kmeans.fit(X_scaled)
            
            # Calculate the ratio of inertias between consecutive clusters
            if n_clusters == 2:
                inertia_ratio = kmeans.inertia_
            else:
                inertia_ratio = kmeans.inertia_ / prev_inertia
            
            if inertia_ratio < best_inertia_ratio:
                best_inertia_ratio = inertia_ratio
                best_clusters = kmeans.labels_
            
            prev_inertia = kmeans.inertia_
        
        if best_clusters is not None:
            peaks['cluster'] = best_clusters
            
            # Sort clusters by mean height
            cluster_stats = peaks.groupby('cluster').agg({
                'height': 'mean'
            }).sort_values('height', ascending=False)
            
            # Get the two highest clusters
            if len(cluster_stats) >= 2:
                cluster1_height = cluster_stats['height'].iloc[0]
                cluster2_height = cluster_stats['height'].iloc[1]
                relative_distance = round(cluster2_height / cluster1_height, 2)
                
                # Check if the ratio indicates potential contamination
                if 0.1 <= relative_distance <= 0.8:
                    main_cluster = peaks[peaks['cluster'] == cluster_stats.index[0]]
                    contamination_cluster = peaks[peaks['cluster'] == cluster_stats.index[1]]
                    
                    return self._create_contamination_info(main_cluster, contamination_cluster, relative_distance)
        
        return None
    
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