import pandas as pd
import numpy as np
from pathlib import Path
import json
from typing import Dict, List, Optional, Tuple
from ..models.data_models import PeakInfo, VariantInfo, LocusResult

class PeakCaller:
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize peak caller with configuration.
        
        Args:
            config_path (str, optional): Path to configuration file
        """
        self.config = self._load_config(config_path)
        self.dye_cutoffs = self.config['dye_cutoffs']
        self.marker_info = self.config['markers']
        self.marker_order = self.config['marker_order']
    
    def _load_config(self, config_path: Optional[str] = None) -> Dict:
        """Load configuration from file."""
        if config_path is None:
            config_path = Path(__file__).parent.parent / 'config' / 'marker_info.json'
        
        with open(config_path) as f:
            return json.load(f)
    
    def _preprocess_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Preprocess the input data.
        
        Args:
            data (pd.DataFrame): Raw input data
            
        Returns:
            pd.DataFrame: Preprocessed data
        """
        # Create a copy of the data to avoid SettingWithCopyWarning
        data = data.copy()
        
        # Extract dye from Dye/Sample Peak column
        data.loc[:, 'dye'] = data['Dye/Sample Peak'].str.split(',').str[0]
        data.loc[:, 'dye'] = data['dye'].str.replace('"O', 'O')
        
        # Use full Sample File Name as sample identifier
        data.loc[:, 'sample'] = data['Sample File Name']
        
        # Filter out rows without marker information but keep OL alleles
        valid_data = data[data['Marker'].notna()].copy()
        
        # Convert numeric columns
        valid_data.loc[:, 'Size'] = pd.to_numeric(valid_data['Size'], errors='coerce')
        valid_data.loc[:, 'Height'] = pd.to_numeric(valid_data['Height'], errors='coerce')
        
        # Rename columns to match expected format
        valid_data = valid_data.rename(columns={
            'Marker': 'marker',
            'Allele': 'allele',
            'Size': 'size',
            'Height': 'height'
        })
        
        return valid_data
    
    def _identify_primary_peaks(self, peaks_df: pd.DataFrame) -> pd.DataFrame:
        """
        Identify primary peaks using dye-specific cutoffs.
        
        Args:
            peaks_df (pd.DataFrame): DataFrame containing peak information
            
        Returns:
            pd.DataFrame: Selected peaks meeting the criteria
        """
        if peaks_df.empty:
            return pd.DataFrame(columns=peaks_df.columns)
        
        # Get sample name and check if it's a negative control
        sample_name = peaks_df['sample'].iloc[0]
        is_neg_control = 'NEG' in sample_name
        
        # Get dye information
        dye = peaks_df['dye'].iloc[0] if not peaks_df['dye'].empty else None
        if dye is None:
            return pd.DataFrame(columns=peaks_df.columns)
        
        # Get dye-specific cutoffs
        dye_key = dye.strip('"')
        default_cutoffs = {'min': 1000, 'max': float('inf')}
        cutoffs = self.dye_cutoffs.get(dye_key, default_cutoffs)
        
        # Create a copy of the DataFrame for filtering
        valid_peaks = peaks_df[peaks_df['allele'] != 'OL'].copy()
        
        if valid_peaks.empty:
            return pd.DataFrame(columns=peaks_df.columns)
        
        # First apply only minimum cutoff
        min_filtered_peaks = valid_peaks[valid_peaks['height'] >= cutoffs['min']].copy()
        
        if min_filtered_peaks.empty:
            return pd.DataFrame(columns=peaks_df.columns)
        
        # Sort peaks by height in descending order
        min_filtered_peaks = min_filtered_peaks.sort_values('height', ascending=False)
        
        # Check if there's only one significant peak after applying minimum cutoff
        if len(min_filtered_peaks) == 1:
            filtered_peaks = min_filtered_peaks
        else:
            if is_neg_control:
                filtered_peaks = min_filtered_peaks
            else:
                filtered_peaks = min_filtered_peaks[min_filtered_peaks['height'] <= cutoffs['max']].copy()
        
        if filtered_peaks.empty:
            return pd.DataFrame(columns=peaks_df.columns)
        
        # Get the highest peak
        max_height = filtered_peaks['height'].iloc[0]
        
        # Calculate relative heights
        filtered_peaks.loc[:, 'relative_height'] = round((filtered_peaks['height'] / max_height) * 100, 2)
        
        # Select peaks that are at least 10% of the highest peak
        significant_peaks = filtered_peaks[filtered_peaks['relative_height'] >= 10].copy()
        
        return significant_peaks
    
    def call_peaks(self, data: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Call peaks for all markers in the data.
        
        Args:
            data (pd.DataFrame): Raw input data
            
        Returns:
            Dict[str, pd.DataFrame]: Dictionary of primary peaks for each marker
        """
        preprocessed_data = self._preprocess_data(data)
        peaks_by_marker = {}
        
        for marker in self.marker_order:
            marker_data = preprocessed_data[preprocessed_data['marker'] == marker]
            if not marker_data.empty:
                primary_peaks = self._identify_primary_peaks(marker_data)
                if not primary_peaks.empty:
                    peaks_by_marker[marker] = primary_peaks
        
        return peaks_by_marker 