import sys
import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from pathlib import Path
import pandas as pd
import logging
import json

# Import FASTER pipeline components
from faster.core.peak_caller import PeakCaller
from faster.core.contamination import ContaminationDetector
from faster.utils.plotting import PeakPlotter
from faster.utils.results import ResultGenerator
from faster.utils.report_generator import ReportGenerator

# Configure logging to print to GUI
logger = logging.getLogger("faster_gui")
logger.setLevel(logging.INFO)

class GuiLogHandler(logging.Handler):
    def __init__(self, text_widget):
        super().__init__()
        self.text_widget = text_widget
    def emit(self, record):
        msg = self.format(record)
        self.text_widget.after(0, self.text_widget.insert, tk.END, msg + '\n')
        self.text_widget.after(0, self.text_widget.see, tk.END)

def get_resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.dirname(__file__), relative_path)

def run_str_analysis(input_path, output_dir, log_widget, generate_plots=True, generate_plotly=True):
    # Set up logger to GUI
    handler = GuiLogHandler(log_widget)
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.handlers = [handler]
    try:
        logger.info("Starting STR analysis...")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create plot directory if needed
        plot_dir = None
        if generate_plots:
            plot_dir = output_dir / 'plots'
            plot_dir.mkdir(exist_ok=True)
            
        # Use default config
        config_path = get_resource_path('config/marker_info.json')
        
        # Initialize components
        peak_caller = PeakCaller(config_path=config_path)
        contamination_detector = ContaminationDetector()
        result_generator = ResultGenerator(config_path=config_path)
        plotter = PeakPlotter(config_path=config_path)
        report_generator = ReportGenerator()
        
        # Read input data
        data = pd.read_csv(input_path, sep='\t')
        all_results = []
        all_plotly_plots = {}
        
        for sample_name in data['Sample File Name'].unique():
            sample_data = data[data['Sample File Name'] == sample_name]
            
            # Call peaks
            peaks_by_marker = peak_caller.call_peaks(sample_data)
            
            # Detect contamination
            contamination_by_marker = {}
            join_points_by_marker = {}
            for marker, peaks in peaks_by_marker.items():
                contamination_result = contamination_detector.detect_contamination(peaks)
                if contamination_result:
                    contamination_info, join_points = contamination_result
                    if contamination_info:
                        contamination_by_marker[marker] = contamination_info
                        join_points_by_marker[marker] = join_points
            
            # Generate and save results
            results = result_generator.generate_results(
                peaks_by_marker=peaks_by_marker,
                contamination_by_marker=contamination_by_marker,
                sample_name=sample_name
            )
            result_generator.save_results(results, output_dir)
            
            # Generate static PNG plots if requested
            if generate_plots:
                plotter.plot_sample_summary(
                    peaks_by_marker,
                    contamination_by_marker,
                    sample_name,
                    str(plot_dir),
                    join_points_by_marker
                )
            
            # Generate plotly plots if requested
            if generate_plotly:
                plotly_plots = plotter.generate_plotly_plots(
                    peaks_by_marker,
                    contamination_by_marker,
                    join_points_by_marker
                )
                all_plotly_plots[results['SampleParameters']['SampleId']] = plotly_plots
            
            all_results.append(results)
            logger.info(f"Processed sample: {sample_name}")
        
        # Generate HTML report if plotly is enabled
        if generate_plotly:
            report_generator.generate_combined_report(
                all_results=all_results,
                plot_dir=plot_dir if generate_plots else None,
                output_dir=output_dir,
                plotly_plots=all_plotly_plots
            )
            logger.info(f"Analysis complete. Report: {output_dir / 'STR_analysis_report.html'}")
        else:
            logger.info("Analysis complete.")

        with open(config_path) as f:
            config = json.load(f)
            marker_order = config.get('marker_order', list(config['markers'].keys()))

        genotype_rows = []
        main_profile_rows = []
        contamination_rows = []
        sample_names = []

        for results in all_results:
            sample_id = results['SampleParameters']['SampleId']
            sample_names.append(sample_id)
            locus = results['LocusResults']
            genotype_row = {}
            main_profile_row = {}
            contamination_row = {}
            for marker in marker_order:
                marker_data = locus.get(marker)
                if marker_data is None:
                    genotype_row[marker] = ''
                    main_profile_row[marker] = ''
                    contamination_row[marker] = ''
                    continue
                variant_info = list(marker_data['variants'].values())[0]
                # Genotype
                genotype_row[marker] = variant_info.get('genotype', '')
                # Main Profile
                if variant_info.get('contamination') and variant_info['contamination'].get('main_profile_peaks'):
                    main_profile = '/'.join(p['allele'] for p in variant_info['contamination']['main_profile_peaks'])
                else:
                    peaks = variant_info.get('peaks', [])
                    sorted_peaks = sorted(peaks, key=lambda x: x['height'], reverse=True)
                    main_profile = '/'.join(p['allele'] for p in sorted_peaks[:marker_data['allele_count']])
                main_profile_row[marker] = main_profile
                # Contamination
                contamination_row[marker] = (
                    '1' if variant_info.get('contamination') and variant_info['contamination'].get('is_contaminated') else ''
                )
            genotype_rows.append(genotype_row)
            main_profile_rows.append(main_profile_row)
            contamination_rows.append(contamination_row)

        # Write CSVs
        genotype_df = pd.DataFrame(genotype_rows, index=sample_names)[marker_order]
        genotype_df = genotype_df.astype(str)
        genotype_df.to_csv(output_dir / 'STR_analysis.genotype.csv', index_label='sample_name')
        genotype_df.to_excel(output_dir / 'STR_analysis.genotype.xlsx', index_label='sample_name')
        pd.DataFrame(main_profile_rows, index=sample_names)[marker_order].to_csv(output_dir / 'STR_analysis.main_profile.csv', index_label='sample_name')
        pd.DataFrame(main_profile_rows, index=sample_names)[marker_order].to_excel(output_dir / 'STR_analysis.main_profile.xlsx', index_label='sample_name')
        pd.DataFrame(contamination_rows, index=sample_names)[marker_order].to_csv(output_dir / 'STR_analysis.contamination.csv', index_label='sample_name')
        pd.DataFrame(contamination_rows, index=sample_names)[marker_order].to_excel(output_dir / 'STR_analysis.contamination.xlsx', index_label='sample_name')
            
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        messagebox.showerror("Error", str(e))

def threaded_run(input_path, output_dir, log_widget, generate_plots, generate_plotly):
    thread = threading.Thread(target=run_str_analysis, args=(input_path, output_dir, log_widget, generate_plots, generate_plotly))
    thread.start()

def main():
    root = tk.Tk()
    root.title("FASTER STR Analysis GUI")
    root.geometry("800x600")
    
    # Create main frame
    main_frame = tk.Frame(root, padx=20, pady=20)
    main_frame.pack(fill='both', expand=True)
    
    # Input file selection
    input_frame = tk.LabelFrame(main_frame, text="Input", padx=10, pady=10)
    input_frame.pack(fill='x', pady=(0, 10))
    
    input_var = tk.StringVar()
    input_entry = tk.Entry(input_frame, textvariable=input_var, width=70)
    input_entry.pack(side='left', padx=(0, 10), fill='x', expand=True)
    
    def browse_input():
        path = filedialog.askopenfilename(filetypes=[("All files", "*.*")])
        if path:
            input_var.set(path)
    tk.Button(input_frame, text="Browse", command=browse_input).pack(side='right')

    # Output directory selection
    output_frame = tk.LabelFrame(main_frame, text="Output", padx=10, pady=10)
    output_frame.pack(fill='x', pady=(0, 10))
    
    output_var = tk.StringVar()
    output_entry = tk.Entry(output_frame, textvariable=output_var, width=70)
    output_entry.pack(side='left', padx=(0, 10), fill='x', expand=True)
    
    def browse_output():
        path = filedialog.askdirectory()
        if path:
            output_var.set(path)
    tk.Button(output_frame, text="Browse", command=browse_output).pack(side='right')

    # Options frame
    options_frame = tk.LabelFrame(main_frame, text="Options", padx=10, pady=10)
    options_frame.pack(fill='x', pady=(0, 10))
    
    plot_var = tk.BooleanVar(value=True)
    plotly_var = tk.BooleanVar(value=True)
    
    tk.Checkbutton(options_frame, text="Generate static PNG plots", variable=plot_var).pack(anchor='w')
    tk.Checkbutton(options_frame, text="Generate interactive Plotly plots", variable=plotly_var).pack(anchor='w')

    # Log/output box
    log_frame = tk.LabelFrame(main_frame, text="Log", padx=10, pady=10)
    log_frame.pack(fill='both', expand=True)
    
    log_box = scrolledtext.ScrolledText(log_frame, height=10, width=80, state='normal')
    log_box.pack(fill='both', expand=True)

    # Run button
    def on_run():
        input_path = input_var.get()
        output_dir = output_var.get()
        if not input_path or not output_dir:
            messagebox.showwarning("Missing input", "Please select both input file and output directory.")
            return
        log_box.delete(1.0, tk.END)
        threaded_run(input_path, output_dir, log_box, plot_var.get(), plotly_var.get())
    
    button_frame = tk.Frame(main_frame)
    button_frame.pack(fill='x', pady=(10, 0))
    
    run_button = tk.Button(button_frame, text="Run Analysis", command=on_run, 
                          bg="#4CAF50", fg="white", height=2, width=20)
    run_button.pack(side='right')

    root.mainloop()

if __name__ == "__main__":
    main() 