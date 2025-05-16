import os
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime

class ReportGenerator:
    """Generates HTML report with interactive navigation for STR analysis results."""
    
    def __init__(self):
        """Initialize report generator."""
        self.html_template = '''
        <!DOCTYPE html>
        <html>
        <head>
            <title>STR Analysis Report</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                body {
                    font-family: Arial, sans-serif;
                    margin: 20px;
                }
                .nav {
                    position: fixed;
                    top: 0;
                    left: 0;
                    width: 250px;
                    height: 100%;
                    overflow-y: auto;
                    background: #f8f9fa;
                    padding: 20px;
                }
                .content {
                    margin-left: 270px;
                }
                .marker-section {
                    margin-bottom: 30px;
                }
                .sample-header {
                    background: #e9ecef;
                    padding: 10px;
                    margin: 20px 0;
                }
                .contamination-summary {
                    background: #fff3cd;
                    padding: 15px;
                    margin: 10px 0;
                }
                .plot-container {
                    margin: 10px 0;
                }
                .nav-link {
                    color: #007bff;
                    text-decoration: none;
                    display: block;
                    margin: 5px 0;
                }
                .nav-link:hover {
                    text-decoration: underline;
                }
                .contamination-details {
                    background: #f8d7da;
                    padding: 10px;
                    margin: 5px 0;
                }
                .sample-section {
                    border: 1px solid #dee2e6;
                    padding: 20px;
                    margin: 20px 0;
                    border-radius: 5px;
                }
                .nav-sample {
                    font-weight: bold;
                    color: #495057;
                    margin: 15px 0 5px 0;
                }
                .nav-marker {
                    padding-left: 15px;
                }
                .sample-link {
                    color: #007bff;
                    text-decoration: none;
                }
                .sample-link:hover {
                    text-decoration: underline;
                }
                .plot-tabs {
                    margin-top: 10px;
                }
                .plot-tab {
                    display: inline-block;
                    padding: 8px 15px;
                    cursor: pointer;
                    background: #e9ecef;
                    border: 1px solid #dee2e6;
                    border-radius: 4px 4px 0 0;
                    margin-right: 5px;
                }
                .plot-tab.active {
                    background: #fff;
                    border-bottom: 1px solid #fff;
                }
                .plot-content {
                    border: 1px solid #dee2e6;
                    padding: 15px;
                    margin-top: -1px;
                }
                .plot-panel {
                    display: none;
                }
                .plot-panel.active {
                    display: block;
                }
                .nav-markers {
                    margin-left: 15px;
                }
                .nav-marker-header {
                    color: #666;
                    font-size: 0.9em;
                    margin: 5px 0;
                }
                .nav-sample {
                    margin: 15px 0;
                }
                .sample-link {
                    font-weight: bold;
                    font-size: 1.1em;
                    color: #007bff;
                }
                .markers-container {
                    margin-top: 20px;
                }
            </style>
            <script>
                function switchPlot(markerId, plotType) {
                    // Hide all plot panels for this marker
                    document.querySelectorAll(`#${markerId} .plot-panel`).forEach(panel => {
                        panel.classList.remove('active');
                    });
                    // Show selected plot panel
                    document.querySelector(`#${markerId} .${plotType}-plot`).classList.add('active');
                    
                    // Update tab styling
                    document.querySelectorAll(`#${markerId} .plot-tab`).forEach(tab => {
                        tab.classList.remove('active');
                    });
                    document.querySelector(`#${markerId} .${plotType}-tab`).classList.add('active');
                }
            </script>
        </head>
        <body>
            <div class="nav">
                <h3>Navigation</h3>
                <div id="sample-nav"></div>
            </div>
            <div class="content">
                <h1>STR Analysis Report</h1>
                <p>Generated: {generation_time}</p>
                <div id="report-content"></div>
            </div>
        </body>
        </html>
        '''

    def _create_navigation(self, samples_data: List[Dict[str, Any]]) -> str:
        """Create navigation HTML."""
        nav_html = []
        
        # Add each sample and its markers
        for results in samples_data:
            sample_id = results['SampleParameters']['SampleId']
            markers = list(results['LocusResults'].keys())
            
            # Add sample link with contamination info
            nav_html.extend([
                f'<div class="nav-sample">',
                f'<a class="nav-link sample-link" href="#sample_{sample_id}">{sample_id}</a>',
                '<div class="nav-markers">',
                '<div class="nav-marker-header">Markers:</div>'
            ])
            
            # Add marker links
            for marker in markers:
                nav_html.append(
                    f'<div class="nav-marker">'
                    f'<a class="nav-link" href="#{sample_id}_{marker}">{marker}</a>'
                    f'</div>'
                )
            
            nav_html.append('</div></div>')  # Close nav-markers and nav-sample
        
        return '\n'.join(nav_html)

    def _create_contamination_summary(self, results: Dict[str, Any]) -> str:
        """Create contamination summary HTML."""
        if 'SampleContamination' not in results:
            return ''
        
        contamination = results['SampleContamination']
        summary_html = [
            '<div class="contamination-summary">',
            '<h3>Sample Contamination Summary</h3>',
            f'<p>Contamination Rate: {contamination["contamination_rate"]}% ({contamination["total_contaminated_markers"]} out of {contamination["total_valid_markers"]} markers)</p>'
        ]
        
        if contamination['contaminated_markers']:
            summary_html.append('<h4>Contaminated Markers:</h4>')
            for marker_info in contamination['contaminated_markers']:
                summary_html.append(
                    f'<div class="contamination-details">'
                    f'<p><strong>{marker_info["marker"]}</strong></p>'
                    f'<p>Main Profile: {marker_info["main_profile"]}</p>'
                    f'<p>Contamination: {marker_info["contamination_peaks"]}</p>'
                    f'<p>Relative Distance: {marker_info["relative_distance"]}</p>'
                    f'</div>'
                )
        else:
            summary_html.append('<p>No contamination detected</p>')
        
        summary_html.append('</div>')
        return '\n'.join(summary_html)

    def _create_marker_info(self, marker_data: Dict[str, Any]) -> str:
        """Create HTML for marker information."""
        variant_info = marker_data["variants"][next(iter(marker_data["variants"]))]
        
        info_html = [
            '<div class="marker-info">',
            f'<p><strong>Genotype:</strong> {variant_info["genotype"]}</p>',
            f'<p><strong>Allele Count:</strong> {marker_data["allele_count"]}</p>',
            f'<p><strong>Median Height:</strong> {int(marker_data["median_height"])}</p>',
            f'<p><strong>Dye:</strong> {marker_data["dye"]} (Limits: {marker_data["height_limits"]["min"]}-{marker_data["height_limits"]["max"]})</p>'
        ]
        
        if marker_data["std_height"] is not None:
            info_html.append(f'<p><strong>Height StdDev:</strong> {round(marker_data["std_height"], 2)}</p>')
        
        if variant_info.get("contamination"):
            contamination = variant_info["contamination"]
            info_html.extend([
                '<div class="marker-contamination">',
                '<p><strong>Contamination Details:</strong></p>',
                f'<p>Main Profile: {"/".join(p["allele"] for p in contamination["main_profile_peaks"])}</p>',
                '<p>Contamination Peaks:</p>',
                '<ul>'
            ])
            
            for peak in contamination["contamination_peaks"]:
                info_html.append(
                    f'<li>Allele {peak["allele"]}: {int(peak["height"])} RFU ({peak["relative_height"]:.3f})</li>'
                )
            
            info_html.extend([
                '</ul>',
                f'<p>Relative Distance: {contamination["relative_distance"]}</p>',
                '</div>'
            ])
        
        info_html.append('</div>')
        return '\n'.join(info_html)

    def _create_marker_section(self, marker: str, marker_data: Dict[str, Any], sample_id: str, plot_dir: Optional[Path], plotly_plot: Optional[str] = None) -> str:
        """Create HTML for a marker section."""
        marker_id = f"{sample_id}_{marker}"
        section_html = [
            f'<div id="{marker_id}" class="marker-section">',
            f'<h3>{marker}</h3>',
            self._create_marker_info(marker_data)
        ]
        
        # Add plot tabs if we have plots
        if plot_dir or plotly_plot:
            section_html.extend([
                '<div class="plot-tabs">',
                '<div class="plot-tab static-tab active" onclick="switchPlot(\'' + marker_id + '\', \'static\')">Static Plot</div>'
            ])
            
            if plotly_plot:
                section_html.append('<div class="plot-tab interactive-tab" onclick="switchPlot(\'' + marker_id + '\', \'interactive\')">Interactive Plot</div>')
            
            section_html.append('</div>')
            
            # Add plot panels
            section_html.append('<div class="plot-content">')
            
            # Static plot panel
            if plot_dir:
                static_plot_path = plot_dir / f"{sample_id}_{marker}.png"
                if static_plot_path.exists():
                    section_html.extend([
                        '<div class="plot-panel static-plot active">',
                        f'<img src="{static_plot_path}" alt="{marker} plot" style="max-width: 100%;">',
                        '</div>'
                    ])
            
            # Interactive plot panel
            if plotly_plot:
                section_html.extend([
                    '<div class="plot-panel interactive-plot">',
                    f'<div id="plotly_{marker_id}"></div>',
                    '</div>'
                ])
            
            section_html.append('</div>')
        
        section_html.append('</div>')
        
        # Add plotly plot initialization if available
        if plotly_plot:
            section_html.append(f'<script>Plotly.newPlot("plotly_{marker_id}", {plotly_plot});</script>')
        
        return '\n'.join(section_html)

    def generate_combined_report(self, 
                               all_results: List[Dict[str, Any]],
                               plot_dir: Optional[Path],
                               output_dir: Path,
                               plotly_plots: Dict[str, Dict[str, str]]) -> None:
        """Generate a combined HTML report for all samples."""
        # Create navigation
        nav_html = self._create_navigation(all_results)
        
        # Create content
        content_html = []
        for results in all_results:
            sample_id = results['SampleParameters']['SampleId']
            
            # Add sample section
            content_html.extend([
                f'<div class="sample-section" id="sample_{sample_id}">',
                f'<h2>Sample: {sample_id}</h2>',
                self._create_contamination_summary(results)
            ])
            
            # Add marker sections
            for marker, marker_data in results['LocusResults'].items():
                content_html.append(
                    self._create_marker_section(
                        marker=marker,
                        marker_data=marker_data,
                        sample_id=sample_id,
                        plot_dir=plot_dir,
                        plotly_plot=plotly_plots.get(sample_id, {}).get(marker)
                    )
                )
            
            content_html.append('</div>')  # Close sample-section
        
        # Generate final HTML
        html_content = self.html_template.format(
            generation_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        )
        
        # Insert navigation and content
        html_content = html_content.replace('<div id="sample-nav"></div>', nav_html)
        html_content = html_content.replace('<div id="report-content"></div>', '\n'.join(content_html))
        
        # Save report
        output_path = output_dir / 'STR_analysis_report.html'
        with open(output_path, 'w') as f:
            f.write(html_content) 