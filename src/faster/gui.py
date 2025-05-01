import sys
import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from pathlib import Path
import pandas as pd
import logging

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
    # PyInstaller로 패키징된 경우
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    # 개발 환경(소스 실행)일 때
    return os.path.join(os.path.dirname(__file__), relative_path)

def run_str_analysis(input_path, output_dir, log_widget):
    # Set up logger to GUI
    handler = GuiLogHandler(log_widget)
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.handlers = [handler]
    try:
        logger.info("Starting STR analysis...")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
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
            peaks_by_marker = peak_caller.call_peaks(sample_data)
            contamination_by_marker = {}
            for marker, peaks in peaks_by_marker.items():
                contamination_result = contamination_detector.detect_contamination(peaks)
                if contamination_result:
                    contamination_info = contamination_result
                    contamination_by_marker[marker] = contamination_info
            results = result_generator.generate_results(
                peaks_by_marker=peaks_by_marker,
                contamination_by_marker=contamination_by_marker,
                sample_name=sample_name
            )
            result_generator.save_results(results, output_dir)
            # Generate plots
            plotter.plot_sample_summary(
                peaks_by_marker,
                contamination_by_marker,
                sample_name,
                str(plot_dir)
            )
            # Generate plotly plots for HTML report
            plotly_plots = plotter.generate_plotly_plots(
                peaks_by_marker,
                contamination_by_marker
            )
            all_plotly_plots[results['SampleParameters']['SampleId']] = plotly_plots
            all_results.append(results)
            logger.info(f"Processed sample: {sample_name}")
        report_generator.generate_combined_report(
            all_results=all_results,
            plot_dir=plot_dir,
            output_dir=output_dir,
            plotly_plots=all_plotly_plots
        )
        logger.info(f"Analysis complete. Report: {output_dir / 'STR_analysis_report.html'}")
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        messagebox.showerror("Error", str(e))

def threaded_run(input_path, output_dir, log_widget):
    thread = threading.Thread(target=run_str_analysis, args=(input_path, output_dir, log_widget))
    thread.start()

def main():
    root = tk.Tk()
    root.title("FASTER STR Analysis GUI")
    root.geometry("600x400")

    # Input file selection
    tk.Label(root, text="Input TSV file:").pack(anchor='w', padx=10, pady=(10,0))
    input_var = tk.StringVar()
    input_entry = tk.Entry(root, textvariable=input_var, width=60)
    input_entry.pack(padx=10, pady=2)
    def browse_input():
        path = filedialog.askopenfilename(filetypes=[("All files", "*.*")])
        if path:
            input_var.set(path)
    tk.Button(root, text="Browse", command=browse_input).pack(padx=10, pady=2, anchor='w')

    # Output directory selection
    tk.Label(root, text="Output directory:").pack(anchor='w', padx=10, pady=(10,0))
    output_var = tk.StringVar()
    output_entry = tk.Entry(root, textvariable=output_var, width=60)
    output_entry.pack(padx=10, pady=2)
    def browse_output():
        path = filedialog.askdirectory()
        if path:
            output_var.set(path)
    tk.Button(root, text="Browse", command=browse_output).pack(padx=10, pady=2, anchor='w')

    # Log/output box
    tk.Label(root, text="Log:").pack(anchor='w', padx=10, pady=(10,0))
    log_box = scrolledtext.ScrolledText(root, height=10, width=80, state='normal')
    log_box.pack(padx=10, pady=2, fill='both', expand=True)

    # Run button
    def on_run():
        input_path = input_var.get()
        output_dir = output_var.get()
        if not input_path or not output_dir:
            messagebox.showwarning("Missing input", "Please select both input file and output directory.")
            return
        log_box.delete(1.0, tk.END)
        threaded_run(input_path, output_dir, log_box)
    tk.Button(root, text="Run Analysis", command=on_run, bg="#4CAF50", fg="white", height=2).pack(pady=10)

    root.mainloop()

if __name__ == "__main__":
    main() 