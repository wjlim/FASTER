import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import json
from pathlib import Path
from typing import Dict, Optional, Tuple
from ..models.data_models import ContaminationInfo

class PeakPlotter:
    """Plots peak data and contamination information."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the plotter with style settings.
        
        Args:
            config_path: Path to marker configuration file (JSON)
        """
        plt.style.use('default')
        
        # Load dye cutoffs from config
        if config_path is None:
            config_path = str(Path(__file__).parent.parent / 'config' / 'marker_info.json')
        
        with open(config_path) as f:
            config = json.load(f)
            self.dye_cutoffs = config['dye_cutoffs']
    
    def _create_hover_text(self, peak: pd.Series, dye: str) -> str:
        """Create hover text for plotly plots."""
        dye_limits = self.dye_cutoffs.get(dye, {"min": 1000, "max": 50000})
        return (f"Allele: {peak['allele']}<br>"
                f"Height: {int(peak['height'])}<br>"
                f"Size: {peak['size']:.2f}<br>"
                f"Dye: {dye}<br>"
                f"Dye Limits: {dye_limits['min']}-{dye_limits['max']}")

    def _get_peak_categories(self, peaks_df: pd.DataFrame, contamination_info: Optional[ContaminationInfo]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Categorize peaks into main profile, contamination, and other peaks.
        
        Args:
            peaks_df: DataFrame containing peak information
            contamination_info: Contamination detection results
            
        Returns:
            Tuple of (main_peaks, contamination_peaks, other_peaks)
        """
        # Get main profile peaks (top 2)
        main_peaks = peaks_df.iloc[:2].copy() if len(peaks_df) > 0 else pd.DataFrame(columns=peaks_df.columns)
        
        if contamination_info and contamination_info.is_contaminated:
            # Get contamination peaks
            contamination_alleles = [p.allele for p in contamination_info.contamination_peaks]
            contamination_mask = peaks_df['allele'].isin(contamination_alleles)
            contamination_peaks = peaks_df[contamination_mask].copy()
            
            # Get other peaks (excluding main and contamination)
            other_mask = ~contamination_mask & ~peaks_df.index.isin([0, 1])
            other_peaks = peaks_df[other_mask].copy()
        else:
            # If no contamination, all non-main peaks are other peaks
            contamination_peaks = pd.DataFrame(columns=peaks_df.columns)
            other_peaks = peaks_df.iloc[2:].copy() if len(peaks_df) > 2 else pd.DataFrame(columns=peaks_df.columns)
        
        return main_peaks, contamination_peaks, other_peaks

    def plot_peaks(self,
                  peaks_df: pd.DataFrame,
                  contamination_info: Optional[ContaminationInfo],
                  output_path: Path) -> Tuple[Path, str]:
        """
        Plot peaks and contamination information using both matplotlib and plotly.
        
        Args:
            peaks_df: DataFrame containing peak information
            contamination_info: Contamination detection results
            output_path: Path to save the plot
            
        Returns:
            Tuple of (static plot path, plotly HTML string)
        """
        # Sort peaks by height in descending order
        peaks_df = peaks_df.sort_values('height', ascending=False).reset_index(drop=True)
        
        # Categorize peaks
        main_peaks, contamination_peaks, other_peaks = self._get_peak_categories(peaks_df, contamination_info)
        
        # Matplotlib static plot
        plt.figure(figsize=(12, 6))
        
        # Plot connecting line for all peaks
        x = peaks_df['size']
        y = peaks_df['height']
        plt.plot(x, y, 'b-', linewidth=1, alpha=0.5)
        
        if len(peaks_df) <= 2:
            # Only main profile
            plt.scatter(x, y, c='green', alpha=0.7, label='Main Profile')
            for _, peak in peaks_df.iterrows():
                plt.annotate(f"{peak['allele']}\n({int(peak['height'])})",
                           (peak['size'], peak['height']),
                           xytext=(0, 10), textcoords='offset points',
                           ha='center', va='bottom')
        else:
            # Plot other peaks first (if any)
            if not other_peaks.empty:
                plt.scatter(other_peaks['size'], other_peaks['height'],
                          c='blue', alpha=0.5, label='All Peaks')
                for _, peak in other_peaks.iterrows():
                    plt.annotate(f"{peak['allele']}\n({int(peak['height'])})",
                               (peak['size'], peak['height']),
                               xytext=(0, -20), textcoords='offset points',
                               ha='center', va='top')
            
            # Plot contamination peaks (if any)
            if not contamination_peaks.empty:
                plt.scatter(contamination_peaks['size'], contamination_peaks['height'],
                          c='red', alpha=0.7, label='Contamination')
                for _, peak in contamination_peaks.iterrows():
                    plt.annotate(f"{peak['allele']}\n({int(peak['height'])})",
                               (peak['size'], peak['height']),
                               xytext=(0, -20), textcoords='offset points',
                               ha='center', va='top')
            
            # Plot main profile peaks
            plt.scatter(main_peaks['size'], main_peaks['height'],
                      c='green', s=100, label='Main Profile')
            for _, peak in main_peaks.iterrows():
                plt.annotate(f"{peak['allele']}\n({int(peak['height'])})",
                           (peak['size'], peak['height']),
                           xytext=(0, 10), textcoords='offset points',
                           ha='center', va='bottom')
        
        plt.xlabel('Size (bp)')
        plt.ylabel('Height (RFU)')
        title = 'Peak Analysis Results'
        if contamination_info and contamination_info.is_contaminated:
            title += f'\nCluster Distance: {contamination_info.relative_distance:.2f}'
        plt.title(title)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()

        # Plotly interactive plot
        fig = go.Figure()
        
        # Get dye color for hover text
        dye = peaks_df['dye'].iloc[0] if not peaks_df.empty else 'B'
        
        # Add line trace
        fig.add_trace(go.Scatter(
            x=x, y=y,
            mode='lines',
            line=dict(color='blue', width=1),
            opacity=0.5,
            showlegend=False,
            hoverinfo='skip'
        ))
        
        if len(peaks_df) <= 2:
            # Only main profile
            fig.add_trace(go.Scatter(
                x=x, y=y,
                mode='markers+text',
                marker=dict(color='green', size=10),
                name='Main Profile',
                text=[peak['allele'] for _, peak in peaks_df.iterrows()],
                textposition='top center',
                hovertemplate=[self._create_hover_text(peak, dye) for _, peak in peaks_df.iterrows()]
            ))
        else:
            # Plot other peaks first (if any)
            if not other_peaks.empty:
                fig.add_trace(go.Scatter(
                    x=other_peaks['size'], y=other_peaks['height'],
                    mode='markers+text',
                    marker=dict(color='blue', size=8),
                    name='All Peaks',
                    text=[peak['allele'] for _, peak in other_peaks.iterrows()],
                    textposition='bottom center',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in other_peaks.iterrows()]
                ))
            
            # Plot contamination peaks (if any)
            if not contamination_peaks.empty:
                fig.add_trace(go.Scatter(
                    x=contamination_peaks['size'], y=contamination_peaks['height'],
                    mode='markers+text',
                    marker=dict(color='red', size=8),
                    name='Contamination',
                    text=[peak['allele'] for _, peak in contamination_peaks.iterrows()],
                    textposition='bottom center',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in contamination_peaks.iterrows()]
                ))
            
            # Plot main profile peaks
            if not main_peaks.empty:
                fig.add_trace(go.Scatter(
                    x=main_peaks['size'], y=main_peaks['height'],
                    mode='markers+text',
                    marker=dict(color='green', size=12),
                    name='Main Profile',
                    text=[peak['allele'] for _, peak in main_peaks.iterrows()],
                    textposition='top center',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in main_peaks.iterrows()]
                ))
        
        title = 'Peak Analysis Results (Interactive)'
        if contamination_info and contamination_info.is_contaminated:
            title += f'<br>Cluster Distance: {contamination_info.relative_distance:.2f}'
        
        fig.update_layout(
            title=title,
            xaxis_title='Size (bp)',
            yaxis_title='Height (RFU)',
            hovermode='closest',
            showlegend=True,
            width=800,
            height=500,
            margin=dict(t=50, b=50, l=50, r=50)
        )
        
        plotly_html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        return output_path, plotly_html
    
    def plot_sample_summary(self, peaks_by_marker: Dict[str, pd.DataFrame],
                          contamination_by_marker: Dict[str, ContaminationInfo],
                          sample_name: str, output_dir: str) -> Dict[str, str]:
        """
        Create summary plots for all markers in a sample.
        
        Args:
            peaks_by_marker: Dictionary of peaks for each marker
            contamination_by_marker: Dictionary of contamination info
            sample_name: Sample name
            output_dir: Output directory for plots
            
        Returns:
            Dictionary mapping marker names to plotly HTML strings
        """
        plotly_plots = {}
        for marker, peaks in peaks_by_marker.items():
            contamination_info = contamination_by_marker.get(marker)
            save_path = f"{output_dir}/{sample_name}_{marker}_peaks.png"
            _, plotly_html = self.plot_peaks(peaks, contamination_info, save_path)
            plotly_plots[marker] = plotly_html
        return plotly_plots

    def generate_plotly_plots(self,
                          peaks_by_marker: Dict[str, pd.DataFrame],
                          contamination_by_marker: Dict[str, ContaminationInfo]) -> Dict[str, str]:
        """
        Generate plotly plots for each marker without saving static images.
        
        Args:
            peaks_by_marker: Dictionary of peaks for each marker
            contamination_by_marker: Dictionary of contamination info
            
        Returns:
            Dictionary mapping marker names to plotly HTML strings
        """
        plotly_plots = {}
        
        for marker, peaks_df in peaks_by_marker.items():
            if peaks_df.empty:
                continue
                
            # Sort peaks by height in descending order
            peaks_df = peaks_df.sort_values('height', ascending=False).reset_index(drop=True)
            
            # Categorize peaks
            contamination_info = contamination_by_marker.get(marker)
            main_peaks, contamination_peaks, other_peaks = self._get_peak_categories(peaks_df, contamination_info)
            
            # Create plotly figure
            fig = go.Figure()
            
            # Get dye color for hover text
            dye = peaks_df['dye'].iloc[0] if not peaks_df.empty else 'B'
            
            # Add line trace
            x = peaks_df['size']
            y = peaks_df['height']
            fig.add_trace(go.Scatter(
                x=x, y=y,
                mode='lines',
                line=dict(color='blue', width=1),
                opacity=0.5,
                showlegend=False,
                hoverinfo='skip'
            ))
            
            if len(peaks_df) <= 2:
                # Main profile only
                fig.add_trace(go.Scatter(
                    x=x, y=y,
                    mode='markers+text',
                    marker=dict(color='green', size=10),
                    name='Main Profile',
                    text=[peak['allele'] for _, peak in peaks_df.iterrows()],
                    textposition='top center',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in peaks_df.iterrows()]
                ))
            else:
                # Plot other peaks first (if any)
                if not other_peaks.empty:
                    fig.add_trace(go.Scatter(
                        x=other_peaks['size'], y=other_peaks['height'],
                        mode='markers+text',
                        marker=dict(color='blue', size=8),
                        name='All Peaks',
                        text=[peak['allele'] for _, peak in other_peaks.iterrows()],
                        textposition='bottom center',
                        hovertemplate=[self._create_hover_text(peak, dye) for _, peak in other_peaks.iterrows()]
                    ))
                
                # Plot contamination peaks (if any)
                if not contamination_peaks.empty:
                    fig.add_trace(go.Scatter(
                        x=contamination_peaks['size'], y=contamination_peaks['height'],
                        mode='markers+text',
                        marker=dict(color='red', size=8),
                        name='Contamination',
                        text=[peak['allele'] for _, peak in contamination_peaks.iterrows()],
                        textposition='bottom center',
                        hovertemplate=[self._create_hover_text(peak, dye) for _, peak in contamination_peaks.iterrows()]
                    ))
                
                # Plot main profile peaks
                if not main_peaks.empty:
                    fig.add_trace(go.Scatter(
                        x=main_peaks['size'], y=main_peaks['height'],
                        mode='markers+text',
                        marker=dict(color='green', size=12),
                        name='Main Profile',
                        text=[peak['allele'] for _, peak in main_peaks.iterrows()],
                        textposition='top center',
                        hovertemplate=[self._create_hover_text(peak, dye) for _, peak in main_peaks.iterrows()]
                    ))
            
            title = f'{marker} Peak Analysis'
            if contamination_info and contamination_info.is_contaminated:
                title += f'<br>Cluster Distance: {contamination_info.relative_distance:.2f}'
            
            fig.update_layout(
                title=title,
                xaxis_title='Size (bp)',
                yaxis_title='Height (RFU)',
                hovermode='closest',
                showlegend=True,
                width=800,
                height=500,
                margin=dict(t=50, b=50, l=50, r=50)
            )
            
            plotly_plots[marker] = fig.to_html(full_html=False, include_plotlyjs='cdn')
        
        return plotly_plots 