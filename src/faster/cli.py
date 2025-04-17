import argparse
import pandas as pd
from pathlib import Path
from typing import Optional
from .core.peak_caller import PeakCaller
from .core.contamination import ContaminationDetector
from .utils.plotting import PeakPlotter
from .utils.results import ResultGenerator
from .utils.report_generator import ReportGenerator

def main():
    """Main entry point for the FASTER package."""
    parser = argparse.ArgumentParser(
        description='Forensic Analysis of STRs with Thermofisher Electrophoresis Result'
    )
    
    parser.add_argument('-i', '--input',
                       required=True,
                       help='Input data file (tab-separated)')
    
    parser.add_argument('-o', '--output',
                       required=True,
                       help='Output directory')
    
    parser.add_argument('--config',
                       help='Path to marker configuration file (JSON)',
                       default=None)
    
    parser.add_argument('--max-height',
                       type=int,
                       default=50000,
                       help='Maximum peak height cutoff (default: 50000)')
    
    parser.add_argument('--plot',
                       action='store_true',
                       help='Generate static PNG plots for each marker')
    
    parser.add_argument('--plotly',
                       action='store_true',
                       default=True,
                       help='Generate interactive Plotly plots in HTML report (default: True)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create plot directory if needed
    if args.plot:
        plot_dir = output_dir / 'plots'
        plot_dir.mkdir(exist_ok=True)
    
    try:
        # Initialize components
        peak_caller = PeakCaller(config_path=args.config)
        contamination_detector = ContaminationDetector()
        result_generator = ResultGenerator(config_path=args.config)
        plotter = PeakPlotter(max_height=args.max_height)
        report_generator = ReportGenerator()
        
        # Read input data
        data = pd.read_csv(args.input, sep='\t')
        
        # Store all results and plotly plots
        all_results = []
        all_plotly_plots = {}
        
        # Process each sample
        for sample_name in data['Sample File Name'].unique():
            sample_data = data[data['Sample File Name'] == sample_name]
            
            # Call peaks
            peaks_by_marker = peak_caller.call_peaks(sample_data)
            
            # Detect contamination
            contamination_by_marker = {}
            for marker, peaks in peaks_by_marker.items():
                contamination_info = contamination_detector.detect_contamination(peaks)
                if contamination_info:
                    contamination_by_marker[marker] = contamination_info
            
            # Generate and save results
            results = result_generator.generate_results(
                peaks_by_marker=peaks_by_marker,
                contamination_by_marker=contamination_by_marker,
                sample_name=sample_name
            )
            result_generator.save_results(results, output_dir)
            
            # Generate plots
            if args.plot:
                plotter.plot_sample_summary(
                    peaks_by_marker,
                    contamination_by_marker,
                    sample_name,
                    str(plot_dir)
                )
            
            # Generate plotly plots for HTML report
            if args.plotly:
                plotly_plots = plotter.generate_plotly_plots(
                    peaks_by_marker,
                    contamination_by_marker
                )
                all_plotly_plots[results['SampleParameters']['SampleId']] = plotly_plots
            
            # Store results
            all_results.append(results)
            
            print(f"Processed sample: {sample_name}")
            print(f"Results saved to: {output_dir}")
            if args.plot:
                print(f"Static plots saved to: {plot_dir}")
            print("---")
        
        # Generate combined HTML report with plotly plots
        report_generator.generate_combined_report(
            all_results=all_results,
            plot_dir=plot_dir if args.plot else None,
            output_dir=output_dir,
            plotly_plots=all_plotly_plots if args.plotly else {}
        )
        print(f"Combined HTML report generated: {output_dir / 'STR_analysis_report.html'}")
            
    except Exception as e:
        print(f"Error: {str(e)}")
        raise

if __name__ == '__main__':
    main() 