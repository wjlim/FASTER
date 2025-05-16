import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import json
from pathlib import Path
from typing import Dict, Optional, Tuple, List, Any
from ..models.data_models import ContaminationInfo
from ..core.contamination import JoinPoint

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
                  output_path: Path,
                  join_points: Optional[List[JoinPoint]] = None) -> Tuple[Path, str]:
        """
        Plot peaks and contamination information using plotly.
        Handles visualization of dynamic joining points and connections.

        Args:
            peaks_df: DataFrame containing peak information
            contamination_info: Contamination detection results including triggering distance
            output_path: Path to save the plot
            join_points: Optional list of joining points from neighbor joining

        Returns:
            Tuple of (plot path, plotly HTML string)
        """
        if peaks_df.empty:
            return output_path, ""

        # Sort peaks by height
        peaks_df = peaks_df.sort_values('height', ascending=False).reset_index(drop=True)

        fig = go.Figure()
        dye = peaks_df['dye'].iloc[0]

        # 1. Plot Main profile peaks as green markers (no lines)
        main_peaks = pd.DataFrame()
        if contamination_info and contamination_info.main_profile_peaks:
            main_alleles = [p.allele for p in contamination_info.main_profile_peaks]
            main_peaks = peaks_df[peaks_df['allele'].astype(str).isin(main_alleles)].copy()
            if not main_peaks.empty:
                fig.add_trace(go.Scatter(
                    x=main_peaks['size'],
                    y=main_peaks['height'],
                    mode='markers+text', # Only markers and text
                    marker=dict(color='green', size=12),
                    name='Main Profile',
                    text=[str(x) for x in main_peaks['allele']],
                    textposition='top center',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in main_peaks.iterrows()]
                ))

        # 2. Plot Contamination peaks as red markers (no lines)
        contam_peaks = pd.DataFrame()
        if contamination_info and contamination_info.contamination_peaks:
            contam_alleles = [p.allele for p in contamination_info.contamination_peaks]
            contam_peaks = peaks_df[peaks_df['allele'].astype(str).isin(contam_alleles)].copy()
            if not contam_peaks.empty:
                fig.add_trace(go.Scatter(
                    x=contam_peaks['size'],
                    y=contam_peaks['height'],
                    mode='markers+text', # Only markers and text
                    marker=dict(color='red', size=8),
                    name='Contamination',
                    text=[str(x) for x in contam_peaks['allele']],
                    textposition='bottom center',
                    hovertemplate=[self._create_hover_text(peak, dye) for _, peak in contam_peaks.iterrows()]
                ))

        # 3. Plot all joining points as purple stars (no lines)
        if join_points:
            fig.add_trace(go.Scatter(
                x=[jp.size for jp in join_points],
                y=[jp.height for jp in join_points],
                mode='markers+text', # Only markers and text
                marker=dict(
                    symbol='star',
                    size=15,
                    color='purple',
                    line=dict(color='white', width=1)
                ),
                name='Joining Points',
                text=['*' for _ in join_points],
                textposition='top center',
                hovertemplate=[
                    f'Joining Point<br>Size: {jp.size:.2f}<br>Height: {jp.height:.2f}<br>Joined Peaks: {jp.joined_peaks}'
                    for jp in join_points
                ]
            ))

        # 4. If dynamic joining occurred, draw connecting lines from the last joining point
        if join_points and len(join_points) > 0:
            last_join = join_points[-1]
            # Draw black solid lines from last joining point to each main profile peak
            if not main_peaks.empty:
                for _, peak in main_peaks.iterrows():
                    fig.add_trace(go.Scatter(
                        x=[last_join.size, peak['size']],
                        y=[last_join.height, peak['height']],
                        mode='lines',
                        line=dict(color='black', width=2, dash='solid'),
                        showlegend=False,
                        hoverinfo='skip'
                    ))
            # Draw red dotted lines from last joining point to each contamination peak
            if not contam_peaks.empty:
                # Special handling for the *first* contamination peak to show triggering distance
                first_contam_peak = contam_peaks.iloc[0]
                triggering_dist_text = f"Triggering Dist: {contamination_info.relative_distance:.2f}" if contamination_info and contamination_info.is_contaminated else ""
                
                fig.add_trace(go.Scatter(
                    x=[last_join.size, first_contam_peak['size']],
                    y=[last_join.height, first_contam_peak['height']],
                    mode='lines',
                    line=dict(color='red', width=2, dash='dot'),
                    name=f"Contamination Trigger ({triggering_dist_text})", # Add info to legend/hover
                    showlegend=True, # Show this specific line legend
                    hoverinfo='name' # Show the name on hover
                ))
                
                # Draw lines for the rest of the contamination peaks without legend/hover info
                for i in range(1, len(contam_peaks)):
                    peak = contam_peaks.iloc[i]
                    fig.add_trace(go.Scatter(
                        x=[last_join.size, peak['size']],
                        y=[last_join.height, peak['height']],
                        mode='lines',
                        line=dict(color='red', width=2, dash='dot'),
                        showlegend=False,
                        hoverinfo='skip'
                    ))
        # Ensure no other lines are drawn between peaks

        # Update layout and save
        title = f'Peak Analysis Results'
        # Use the stored triggering distance if available
        if contamination_info and contamination_info.is_contaminated:
            title += f'<br>Triggering Relative Distance: {contamination_info.relative_distance:.2f}'
        
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
        
        fig.write_html(str(output_path))
        plotly_html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        return output_path, plotly_html
    
    def plot_sample_summary(self, peaks_by_marker: Dict[str, pd.DataFrame],
                          contamination_by_marker: Dict[str, Optional[ContaminationInfo]],
                          sample_name: str, 
                          output_dir: str,
                          join_points_by_marker: Dict[str, List[JoinPoint]]) -> Dict[str, str]:
        """
        Create summary plots for all markers in a sample.
        
        Args:
            peaks_by_marker: Dictionary of peaks for each marker
            contamination_by_marker: Dictionary of contamination info
            sample_name: Sample name
            output_dir: Output directory for plots
            join_points_by_marker: Dictionary mapping marker names to join points
            
        Returns:
            Dictionary mapping marker names to plotly HTML strings
        """
        plotly_plots = {}
        for marker, peaks in peaks_by_marker.items():
            contamination_info = contamination_by_marker.get(marker)
            join_points = join_points_by_marker.get(marker)
                
            save_path = Path(f"{output_dir}/{sample_name}_{marker}_peaks.html")
            _, plotly_html = self.plot_peaks(peaks, contamination_info, save_path, join_points)
            plotly_plots[marker] = plotly_html
        return plotly_plots

    def generate_plotly_plots(self,
                          peaks_by_marker: Dict[str, pd.DataFrame],
                          contamination_by_marker: Dict[str, Optional[ContaminationInfo]],
                          join_points_by_marker: Dict[str, List[JoinPoint]]) -> Dict[str, str]:
        """
        Generate plotly plots for each marker without saving static images.
        
        Args:
            peaks_by_marker: Dictionary of peaks for each marker
            contamination_by_marker: Dictionary of contamination info
            join_points_by_marker: Dictionary mapping marker names to join points
            
        Returns:
            Dictionary mapping marker names to plotly HTML strings
        """
        plotly_plots = {}
        
        for marker, peaks_df in peaks_by_marker.items():
            if peaks_df.empty:
                continue
                
            # Sort peaks by height in descending order
            peaks_df = peaks_df.sort_values('height', ascending=False).reset_index(drop=True)
            
            # Get contamination info and join points
            contamination_info = contamination_by_marker.get(marker)
            join_points = join_points_by_marker.get(marker)
            
            # Create plotly figure
            fig = go.Figure()
            
            # Get dye color for hover text
            dye = peaks_df['dye'].iloc[0] if not peaks_df.empty else 'B'
            
            # Plot Main profile peaks as green markers (no lines)
            main_peaks = pd.DataFrame()
            if contamination_info and contamination_info.main_profile_peaks:
                main_alleles = [p.allele for p in contamination_info.main_profile_peaks]
                main_peaks = peaks_df[peaks_df['allele'].astype(str).isin(main_alleles)].copy()
                if not main_peaks.empty:
                    fig.add_trace(go.Scatter(
                        x=main_peaks['size'],
                        y=main_peaks['height'],
                        mode='markers+text', # Only markers and text
                        marker=dict(color='green', size=12),
                        name='Main Profile',
                        text=[str(x) for x in main_peaks['allele']],
                        textposition='top center',
                        hovertemplate=[self._create_hover_text(peak, dye) for _, peak in main_peaks.iterrows()]
                    ))
    
            # Plot Contamination peaks as red markers (no lines)
            contam_peaks = pd.DataFrame()
            if contamination_info and contamination_info.contamination_peaks:
                contam_alleles = [p.allele for p in contamination_info.contamination_peaks]
                contam_peaks = peaks_df[peaks_df['allele'].astype(str).isin(contam_alleles)].copy()
                if not contam_peaks.empty:
                    fig.add_trace(go.Scatter(
                        x=contam_peaks['size'],
                        y=contam_peaks['height'],
                        mode='markers+text', # Only markers and text
                        marker=dict(color='red', size=8),
                        name='Contamination',
                        text=[str(x) for x in contam_peaks['allele']],
                        textposition='bottom center',
                        hovertemplate=[self._create_hover_text(peak, dye) for _, peak in contam_peaks.iterrows()]
                    ))
    
            # Plot all joining points as purple stars (no lines)
            if join_points:
                fig.add_trace(go.Scatter(
                    x=[jp.size for jp in join_points],
                    y=[jp.height for jp in join_points],
                    mode='markers+text', # Only markers and text
                    marker=dict(
                        symbol='star',
                        size=15,
                        color='purple',
                        line=dict(color='white', width=1)
                    ),
                    name='Joining Points',
                    text=['*' for _ in join_points],
                    textposition='top center',
                    hovertemplate=[
                        f'Joining Point<br>Size: {jp.size:.2f}<br>Height: {jp.height:.2f}<br>Joined Peaks: {jp.joined_peaks}'
                        for jp in join_points
                    ]
                ))
    
            # If dynamic joining occurred, draw connecting lines from the last joining point
            if join_points and len(join_points) > 0:
                last_join = join_points[-1]
                # Draw black solid lines from last joining point to each main profile peak
                if not main_peaks.empty:
                    for _, peak in main_peaks.iterrows():
                        fig.add_trace(go.Scatter(
                            x=[last_join.size, peak['size']],
                            y=[last_join.height, peak['height']],
                            mode='lines',
                            line=dict(color='black', width=2, dash='solid'),
                            showlegend=False,
                            hoverinfo='skip'
                        ))
                # Draw red dotted lines from last joining point to each contamination peak
                if not contam_peaks.empty:
                    # Special handling for the *first* contamination peak to show triggering distance
                    first_contam_peak = contam_peaks.iloc[0]
                    triggering_dist_text = f"Triggering Dist: {contamination_info.relative_distance:.2f}" if contamination_info and contamination_info.is_contaminated else ""
                    
                    fig.add_trace(go.Scatter(
                        x=[last_join.size, first_contam_peak['size']],
                        y=[last_join.height, first_contam_peak['height']],
                        mode='lines',
                        line=dict(color='red', width=2, dash='dot'),
                        name=f"Contamination Trigger ({triggering_dist_text})", # Add info to legend/hover
                        showlegend=True, # Show this specific line legend
                        hoverinfo='name' # Show the name on hover
                    ))
                    
                    # Draw lines for the rest of the contamination peaks without legend/hover info
                    for i in range(1, len(contam_peaks)):
                        peak = contam_peaks.iloc[i]
                        fig.add_trace(go.Scatter(
                            x=[last_join.size, peak['size']],
                            y=[last_join.height, peak['height']],
                            mode='lines',
                            line=dict(color='red', width=2, dash='dot'),
                            showlegend=False,
                            hoverinfo='skip'
                        ))
            # Ensure no other lines are drawn between peaks
    
            # Update layout and save
            title = f'{marker} Peak Analysis'
            if contamination_info and contamination_info.is_contaminated:
                title += f'<br>Relative Distance: {contamination_info.relative_distance:.2f}'
            
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

class VectorPlotter:
    """Plots STR vector data from vector.json files."""
    
    def __init__(self):
        """Initialize the vector plotter."""
        pass
        
    def _load_vector_data(self, vector_files: List[str]) -> List[Dict]:
        """Load vector data from JSON files.
        
        Args:
            vector_files: List of paths to vector.json files
            
        Returns:
            List of dictionaries containing vector data
        """
        vector_data = []
        for file_path in vector_files:
            with open(file_path) as f:
                data = json.load(f)
                # Add filename to data for hover info
                data['filename'] = Path(file_path).stem
                vector_data.append(data)
        return vector_data
    
    def _create_hover_text(self, data: Dict) -> str:
        """Create hover text for vector plot.
        
        Args:
            data: Vector data dictionary
            
        Returns:
            Formatted hover text
        """
        hover_text = [
            f"Sample: {data['sample_id']}",
            f"Magnitude: {data['vector_properties']['magnitude']:.3f}",
            f"Angle: {data['vector_properties']['angle_degrees']:.2f}°",
            f"X: {data['cartesian_coordinates']['x']:.3f}",
            f"Y: {data['cartesian_coordinates']['y']:.3f}",
            "\nMarkers:"
        ]
        
        # Add marker information
        for marker in data['markers']:
            hover_text.append(
                f"{marker['marker']}: {marker['allele1']}/{marker['allele2']}"
            )
            
        return "<br>".join(hover_text)
    
    def plot_vectors(self, vector_files: List[str], output_dir: str) -> str:
        """Create interactive vector plot using plotly.
        
        Args:
            vector_files: List of paths to vector.json files
            output_dir: Output directory for the plot
            
        Returns:
            Path to the generated HTML file
        """
        # Load vector data
        vector_data = self._load_vector_data(vector_files)
        
        # Create figure
        fig = go.Figure()
        
        # Separate data by source type
        str_data = [d for d in vector_data if d.get('source_type') == 'str']
        eh_data = [d for d in vector_data if d.get('source_type') == 'eh']
        
        # Plot STR data
        if str_data:
            x_coords = [d['cartesian_coordinates']['x'] for d in str_data]
            y_coords = [d['cartesian_coordinates']['y'] for d in str_data]
            hover_texts = [self._create_hover_text(d) for d in str_data]
            sample_ids = [d['sample_id'] for d in str_data]
            
            fig.add_trace(go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='markers+text',
                marker=dict(
                    size=12,
                    color='blue',
                    symbol='circle',
                    line=dict(color='black', width=1)
                ),
                text=sample_ids,
                textposition="top center",
                hovertext=hover_texts,
                hoverinfo='text',
                name='STR Analysis'
            ))
        
        # Plot EH data
        if eh_data:
            x_coords = [d['cartesian_coordinates']['x'] for d in eh_data]
            y_coords = [d['cartesian_coordinates']['y'] for d in eh_data]
            hover_texts = [self._create_hover_text(d) for d in eh_data]
            sample_ids = [d['sample_id'] for d in eh_data]
            
            fig.add_trace(go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='markers+text',
                marker=dict(
                    size=12,
                    color='red',
                    symbol='diamond',
                    line=dict(color='black', width=1)
                ),
                text=sample_ids,
                textposition="top center",
                hovertext=hover_texts,
                hoverinfo='text',
                name='ExpansionHunter'
            ))
        
        # Calculate axis ranges with padding
        all_x = [d['cartesian_coordinates']['x'] for d in vector_data]
        all_y = [d['cartesian_coordinates']['y'] for d in vector_data]
        x_range = [min(all_x), max(all_x)]
        y_range = [min(all_y), max(all_y)]
        
        # Add 10% padding to ranges
        x_padding = (x_range[1] - x_range[0]) * 0.1
        y_padding = (y_range[1] - y_range[0]) * 0.1
        x_range = [x_range[0] - x_padding, x_range[1] + x_padding]
        y_range = [y_range[0] - y_padding, y_range[1] + y_padding]
        
        # Update layout with improved styling
        fig.update_layout(
            title='STR Vector Analysis',
            xaxis=dict(
                title='X Coordinate',
                range=x_range,
                showgrid=True,
                gridcolor='lightgray',
                zeroline=True,
                zerolinecolor='black',
                zerolinewidth=2,
                showline=True,
                linewidth=2,
                linecolor='black',
                mirror=True
            ),
            yaxis=dict(
                title='Y Coordinate',
                range=y_range,
                showgrid=True,
                gridcolor='lightgray',
                zeroline=True,
                zerolinecolor='black',
                zerolinewidth=2,
                showline=True,
                linewidth=2,
                linecolor='black',
                mirror=True
            ),
            plot_bgcolor='white',
            hovermode='closest',
            width=1000,
            height=800,
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99,
                bgcolor='rgba(255, 255, 255, 0.8)',
                bordercolor='black',
                borderwidth=1
            ),
            margin=dict(t=50, b=50, l=50, r=50)
        )
        
        # Save plot
        output_path = Path(output_dir) / 'str_vectors.html'
        fig.write_html(str(output_path))
        
        return str(output_path) 