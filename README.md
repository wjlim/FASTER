# FASTER (Forensic Analysis of STRs with Thermofisher Electrophoresis Result)

A robust tool for analyzing Short Tandem Repeat (STR) data from Thermofisher electrophoresis results, featuring advanced contamination detection and interactive visualization.

## Features

- **Advanced Peak Analysis**
  - Automatic peak detection with dye-specific height thresholds
  - Sophisticated contamination detection using machine learning
  - Support for all standard STR markers

- **Comprehensive Visualization**
  - Interactive HTML reports with plotly graphs
  - Peak height and size visualization
  - Contamination analysis visualization
  - Easy navigation between samples and markers

- **Quality Control**
  - Dye-specific height thresholds
  - Contamination detection and reporting
  - Detailed marker-level analysis

## Installation

### Prerequisites

- Python 3.10 or higher
- pip package manager

### Installation from Source

```bash
git clone https://github.com/wjlim/FASTER.git
cd FASTER
pip install .
```

## Usage

### Basic Command

```bash
faster -i <input_file> -o <output_directory>
```

### Command Line Options

- `-i, --input`: Input data file (tab-separated)
- `-o, --output`: Output directory
- `--config`: Path to marker configuration file (optional)
- `--max-height`: Maximum peak height cutoff (default: 50000)
- `--plot`: Generate static PNG plots for each marker
- `--plotly`: Generate interactive plots in HTML report (default: True)

### Example

```bash
# Basic usage
faster -i example/input.txt -o example_out/

# With custom maximum height
faster -i example/input.txt -o example_out/ --max-height 40000

# Generate both static and interactive plots
faster -i example/input.txt -o example_out/ --plot
```

## Output Files

The tool generates the following outputs:

```
output_directory/
├── {sample_name}.STR_analysis.json    # Analysis results in JSON format
├── STR_analysis_report.html           # Interactive HTML report
└── plots/                             # (Optional) Static plot images
    └── {sample_name}_{marker}_peaks.png
```

### HTML Report Features

- Sample-level navigation
- Marker-specific analysis views
- Interactive plotly graphs with hover information
- Contamination summary and details
- Easy switching between static and interactive plots

### JSON Results Structure

```json
{
  "LocusResults": {
    "marker_name": {
      "allele_count": int,
      "median_height": float,
      "dye": str,
      "std_height": float,
      "height_limits": {
        "min": int,
        "max": int
      },
      "variants": {
        "position": {
          "genotype": str,
          "allele_count": int,
          "motif": str,
          "contamination": {...},
          "peaks": [...]
        }
      }
    }
  },
  "SampleParameters": {...},
  "SampleContamination": {
    "contamination_rate": float,
    "contaminated_markers": [...],
    "total_valid_markers": int,
    "total_contaminated_markers": int
  }
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License

## Support

For support and questions, please create an issue in the GitHub repository.
