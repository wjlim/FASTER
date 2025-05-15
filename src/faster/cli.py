import sys
import argparse
import json
import subprocess
import logging
import pandas as pd
from pathlib import Path
from .core.peak_caller import PeakCaller
from .core.contamination import ContaminationDetector
from .utils.plotting import PeakPlotter
from .utils.results import ResultGenerator
from .utils.report_generator import ReportGenerator
from .core.compare import ResultComparator
from .core.vectorizer import GenotypeVectorizer
import os
import stat

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def main():
    """Main entry point for the FASTER package."""
    parser = argparse.ArgumentParser(
        description='Forensic Analysis of STRs with Thermofisher Electrophoresis Results and ExpansionHunter'
    )

    subparsers = parser.add_subparsers(title='subcommands', help='sub-command help', dest='command', required=True)
    
    # STR analysis subcommand
    str_parser = subparsers.add_parser('str', help='STR analysis help')
    str_parser.add_argument('-i', '--input',
                        required=True,
                        help='Input data file (tab-separated)')
    str_parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory')
    str_parser.add_argument('--config',
                        help='Path to marker configuration file (JSON)',
                        default=f"{Path(__file__).parent}/config/marker_info.json")
    str_parser.add_argument('--plot',
                        action='store_true',
                        help='Generate static PNG plots for each marker')
    str_parser.add_argument('--plotly',
                        default=False,
                        help='Generate interactive Plotly plots in HTML report (default: False)')

    # ExpansionHunter analysis subcommand
    exhunter_parser = subparsers.add_parser('exhunter', help='ExpansionHunter analysis help')
    exhunter_parser.add_argument('-i', '--input_bam',
                        required=True,
                        help='Path to the input BAM file')
    exhunter_parser.add_argument('-r', '--reference',
                        required=True,
                        help='Path to the reference fasta file')
    exhunter_parser.add_argument('-o', '--output_prefix',
                        required=True,
                        help='Output prefix for ExpansionHunter results')

    # Compare results subcommand
    compare_parser = subparsers.add_parser('compare', help='Compare STR analysis and ExpansionHunter results')
    compare_parser.add_argument('-i', '--str_json',
                        required=True,
                        help='Path to STR analysis JSON file')
    compare_parser.add_argument('-j', '--eh_json',
                        required=True,
                        help='Path to ExpansionHunter JSON file')
    compare_parser.add_argument('-o', '--output_prefix',
                        required=True,
                        help='Output prefix for comparison results')

    # Vectorize subcommand
    vectorize_parser = subparsers.add_parser(
        'vectorize',
        help='Vectorize genotype data from STR or ExpansionHunter analysis'
    )
    vectorize_parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input JSON file from STR or ExpansionHunter analysis'
    )
    vectorize_parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output text file to save vectorized results'
    )
    vectorize_parser.add_argument(
        '-t', '--type',
        choices=['str', 'eh'],
        required=True,
        help='Input type: str (STR analysis) or eh (ExpansionHunter)'
    )
    
    # Compare vectors subcommand
    compare_vectors_parser = subparsers.add_parser(
        'compare-vectors',
        help='Compare two compact vector files'
    )
    compare_vectors_parser.add_argument(
        '-i', '--vector1',
        required=True,
        help='First vector file'
    )
    compare_vectors_parser.add_argument(
        '-j', '--vector2',
        required=True,
        help='Second vector file'
    )
    compare_vectors_parser.add_argument(
        '-o', '--output',
        help='Output JSON file for comparison results'
    )
    
    args = parser.parse_args()
    
    if args.command == 'str':
        process_str_analysis(args)
    elif args.command == 'exhunter':
        process_exhunter_analysis(args)
    elif args.command == 'compare':
        process_compare(args)
    elif args.command == 'vectorize':
        run_vectorize(args)
    elif args.command == 'compare-vectors':
        vectorizer = GenotypeVectorizer()
        result = vectorizer.compare_vectors(args.vector1, args.vector2)
        
        # Print comparison to console
        vectorizer.print_comparison(args.vector1, args.vector2)
        
        # Save comparison results to JSON if output path is provided
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result.dict(), f, indent=2)
            logger.info(f"Comparison results saved to {args.output}")

def process_str_analysis(args):
    """Process STR analysis command."""
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
        plotter = PeakPlotter(config_path=args.config)
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
            
            logger.info(f"Processed sample: {sample_name}")
            logger.info(f"Results saved to: {output_dir}")
            if args.plot:
                logger.info(f"Static plots saved to: {plot_dir}")
            logger.info("---")
        
        # Generate combined HTML report with plotly plots
        report_generator.generate_combined_report(
            all_results=all_results,
            plot_dir=plot_dir if args.plot else None,
            output_dir=output_dir,
            plotly_plots=all_plotly_plots if args.plotly else {}
        )
        logger.info(f"Combined HTML report generated: {output_dir / 'STR_analysis_report.html'}")
            
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        raise

def process_exhunter_analysis(args):
    """Process ExpansionHunter analysis command."""
    try:
        # Set paths for ExpansionHunter binary and variant catalog
        package_dir = Path(__file__).parent  # src/faster/
        exhunter_binary = package_dir / 'bin' / 'ExpansionHunter'
        variant_catalog = package_dir / 'config' / 'variant_catalog.thermofisher_24markers.json'

        # Ensure ExpansionHunter binary is executable
        if exhunter_binary.exists():
            os.chmod(exhunter_binary, os.stat(exhunter_binary).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        # Verify binary and variant catalog exist
        if not exhunter_binary.exists():
            raise FileNotFoundError(f"ExpansionHunter binary not found at: {exhunter_binary}")
        if not variant_catalog.exists():
            raise FileNotFoundError(f"Variant catalog not found at: {variant_catalog}")

        # Extract directory from output_prefix if it contains a path
        output_prefix = Path(args.output_prefix)
        if '/' in str(output_prefix):
            output_dir = output_prefix.parent
            output_dir.mkdir(parents=True, exist_ok=True)

        # Prepare command
        cmd = [
            str(exhunter_binary),
            '--reads', str(args.input_bam),
            '--reference', str(args.reference),
            '--variant-catalog', str(variant_catalog),
            '--output-prefix', str(args.output_prefix)
        ]

        # Run ExpansionHunter
        logger.info("Running ExpansionHunter...")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        logger.info("ExpansionHunter analysis completed successfully")
        logger.info(f"Results saved to: {args.output_prefix}.*")
        
    except FileNotFoundError as e:
        logger.error(f"Error: {str(e)}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running ExpansionHunter: {str(e)}")
        logger.error(f"ExpansionHunter stderr output:\n{e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        raise

def process_compare(args):
    """Process comparison command."""
    try:
        comparator = ResultComparator()
        results = comparator.compare_results(args.str_json, args.eh_json)
        comparator.save_results(results, args.output_prefix)
        
    except Exception as e:
        logger.error(f"Error during comparison: {str(e)}")
        raise

def run_vectorize(args):
    """Run vectorization of genotype data."""
    try:
        # Load input data
        with open(args.input) as f:
            data = json.load(f)
            
        # Get sample_id from data or filename
        if args.type == 'eh':
            sample_id = data.get('SampleParameters', {}).get('SampleId', '')
        else:  # str
            sample_id = data.get('SampleParameters', {}).get('sample_id', '')
            
        # If no sample_id found in data, use input filename
        if not sample_id:
            sample_id = Path(args.input).stem
            
        # Initialize vectorizer
        vectorizer = GenotypeVectorizer()
        
        # Vectorize based on input type
        if args.type == 'str':
            result = vectorizer.vectorize_str(data, sample_id)
        else:  # eh
            result = vectorizer.vectorize_eh(data, sample_id)
            
        # Save results
        vectorizer.save_vector(result, args.output)
        
        logging.info(f"Successfully vectorized data and saved to {args.output}")
        
    except Exception as e:
        logging.error(f"Error during vectorization: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main() 