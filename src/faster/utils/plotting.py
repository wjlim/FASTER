import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from ..models.data_models import ContaminationInfo

class PeakPlotter:
    """Plots peak data and contamination information."""
    
    def __init__(self, max_height: int = 50000):
        """
        Initialize the plotter with style settings.
        
        Args:
            max_height: Maximum peak height cutoff (default: 50000)
        """
        plt.style.use('default')
        self.dye_cutoffs = {
            "B": {"min": 2500, "max": max_height},  # Blue
            "G": {"min": 5000, "max": max_height},  # Green
            "Y": {"min": 9000, "max": max_height},  # Yellow
            "R": {"min": 1000, "max": max_height},  # Red
            "P": {"min": 1000, "max": max_height},  # Purple
        }
    
    def _create_hover_text(self, peak: pd.Series, dye: str) -> str:
        """Create hover text for plotly plots."""
        dye_limits = self.dye_cutoffs.get(dye, {"min": 1000, "max": 50000})
        return (f"Allele: {peak['allele']}<br>"
                f"Height: {int(peak['height'])}<br>"
                f"Size: {peak['size']:.2f}<br>"
                f"Dye: {dye}<br>"
                f"Dye Limits: {dye_limits['min']}-{dye_limits['max']}")

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
        # Matplotlib static plot
        plt.figure(figsize=(12, 6))
        
        x = peaks_df['size']
        y = peaks_df['height']
        plt.plot(x, y, 'b-', linewidth=1, alpha=0.5)
        
        if len(peaks_df) <= 2:
            plt.scatter(x, y, c='green', alpha=0.7, label='Main Profile')
            for _, peak in peaks_df.iterrows():
                plt.annotate(f"{peak['allele']}\n({int(peak['height'])})",
                           (peak['size'], peak['height']),
                           xytext=(0, 10), textcoords='offset points',
                           ha='center', va='bottom')
        else:
            plt.scatter(x, y, c='blue', alpha=0.5, label='All Peaks')
            
            if contamination_info and contamination_info.is_contaminated:
                main_x = [p.size for p in contamination_info.main_profile_peaks]
                main_y = [p.height for p in contamination_info.main_profile_peaks]
                plt.scatter(main_x, main_y, c='green', s=100, label='Main Profile')
                
                cont_x = [p.size for p in contamination_info.contamination_peaks]
                cont_y = [p.height for p in contamination_info.contamination_peaks]
                plt.scatter(cont_x, cont_y, c='red', s=100, label='Contamination')
                
                for p in contamination_info.main_profile_peaks:
                    plt.annotate(f"{p.allele}\n({int(p.height)})",
                               (p.size, p.height),
                               xytext=(0, 10), textcoords='offset points',
                               ha='center', va='bottom')
                
                for p in contamination_info.contamination_peaks:
                    plt.annotate(f"{p.allele}\n({int(p.height)})",
                               (p.size, p.height),
                               xytext=(0, -20), textcoords='offset points',
                               ha='center', va='top')
        
        plt.xlabel('Size (bp)')
        plt.ylabel('Height (RFU)')
        plt.title('Peak Analysis Results')
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
            # All peaks
            fig.add_trace(go.Scatter(
                x=x, y=y,
                mode='markers',
                marker=dict(color='blue', size=8),
                name='All Peaks',
                hovertemplate=[self._create_hover_text(peak, dye) for _, peak in peaks_df.iterrows()]
            ))
            
            if contamination_info and contamination_info.is_contaminated:
                # Main profile peaks
                main_x = [p.size for p in contamination_info.main_profile_peaks]
                main_y = [p.height for p in contamination_info.main_profile_peaks]
                main_text = [p.allele for p in contamination_info.main_profile_peaks]
                fig.add_trace(go.Scatter(
                    x=main_x, y=main_y,
                    mode='markers+text',
                    marker=dict(color='green', size=12),
                    name='Main Profile',
                    text=main_text,
                    textposition='top center',
                    hovertemplate=[
                        self._create_hover_text(
                            pd.Series({
                                'allele': p.allele,
                                'height': p.height,
                                'size': p.size
                            }), dye
                        ) for p in contamination_info.main_profile_peaks
                    ]
                ))
                
                # Contamination peaks
                cont_x = [p.size for p in contamination_info.contamination_peaks]
                cont_y = [p.height for p in contamination_info.contamination_peaks]
                cont_text = [p.allele for p in contamination_info.contamination_peaks]
                fig.add_trace(go.Scatter(
                    x=cont_x, y=cont_y,
                    mode='markers+text',
                    marker=dict(color='red', size=12),
                    name='Contamination',
                    text=cont_text,
                    textposition='bottom center',
                    hovertemplate=[
                        self._create_hover_text(
                            pd.Series({
                                'allele': p.allele,
                                'height': p.height,
                                'size': p.size
                            }), dye
                        ) for p in contamination_info.contamination_peaks
                    ]
                ))
        
        fig.update_layout(
            title='Peak Analysis Results (Interactive)',
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
            
            contamination_info = contamination_by_marker.get(marker)
            
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
                # All peaks
                fig.add_trace(go.Scatter(
                    x=x, y=y,
                    mode='markers',
                    marker=dict(color='blue', size=8),
                    name='All Peaks',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in peaks_df.iterrows()]
                ))
                
                if contamination_info and contamination_info.is_contaminated:
                    # Main profile peaks
                    main_x = [p.size for p in contamination_info.main_profile_peaks]
                    main_y = [p.height for p in contamination_info.main_profile_peaks]
                    main_text = [p.allele for p in contamination_info.main_profile_peaks]
                    fig.add_trace(go.Scatter(
                        x=main_x, y=main_y,
                        mode='markers+text',
                        marker=dict(color='green', size=12),
                        name='Main Profile',
                        text=main_text,
                        textposition='top center',
                        hovertemplate=[
                            self._create_hover_text(
                                pd.Series({
                                    'allele': p.allele,
                                    'height': p.height,
                                    'size': p.size
                                }), dye
                            ) for p in contamination_info.main_profile_peaks
                        ]
                    ))
                    
                    # Contamination peaks
                    cont_x = [p.size for p in contamination_info.contamination_peaks]
                    cont_y = [p.height for p in contamination_info.contamination_peaks]
                    cont_text = [p.allele for p in contamination_info.contamination_peaks]
                    fig.add_trace(go.Scatter(
                        x=cont_x, y=cont_y,
                        mode='markers+text',
                        marker=dict(color='red', size=12),
                        name='Contamination',
                        text=cont_text,
                        textposition='bottom center',
                        hovertemplate=[
                            self._create_hover_text(
                                pd.Series({
                                    'allele': p.allele,
                                    'height': p.height,
                                    'size': p.size
                                }), dye
                            ) for p in contamination_info.contamination_peaks
                        ]
                    ))
            
            fig.update_layout(
                title=f'{marker} Peak Analysis',
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