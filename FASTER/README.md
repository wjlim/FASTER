# FASTER (Forensic Analysis of STRs with Thermofisher Electrophoresis Result)

A robust tool for analyzing Short Tandem Repeat (STR) data from Thermofisher electrophoresis results, featuring advanced contamination detection and interactive visualization.

## Features

- **Peak Analysis**
  - Height-based peak detection with dye-specific thresholds
  - Main profile selection (top 2 peaks by height)
  - Maximum 4 peaks per marker consideration
  - Support for multiple dye channels (B, G, Y, R, P)

- **Contamination Detection**
  - Height-based clustering analysis
  - Relative distance calculation between clusters
  - Automatic contamination peak identification
  - Cluster distance reporting for contaminated markers

- **Interactive Visualization**
  - Dynamic HTML reports with plotly graphs
  - Color-coded peak display:
    - Main Profile peaks (green)
    - Contamination peaks (red)
    - Other peaks (blue)
  - Detailed hover information including:
    - Allele values
    - Peak heights
    - Size information
    - Dye channel and limits

- **Quality Control**
  - Dye-specific height thresholds from configuration
  - Height standard deviation analysis
  - Comprehensive contamination reporting
  - Marker-level quality metrics

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
- `--plot`: Generate static PNG plots for each marker
- `--plotly`: Generate interactive plots in HTML report (default: True)

### Example

```bash
# Basic usage
faster -i example/input.txt -o example_out/

# With custom configuration
faster -i example/input.txt -o example_out/ --config path/to/marker_info.json

# Generate both static and interactive plots
faster -i example/input.txt -o example_out/ --plot
```

## Output Files

The tool generates the following outputs:

```
output_directory/
├── {sample_name}.STR_analysis.json    # Analysis results in JSON format
├── {sample_name}.STR_report.html      # Interactive HTML report
└── {sample_name}_peaks/               # Static plot images (if --plot is used)
    └── {sample_name}_{marker}_peaks.png
```

### HTML Report Features

The interactive HTML report includes:

- Sample-level summary
- Marker-specific analysis views
- Interactive plotly graphs showing:
  - Main Profile peaks (green dots)
  - Contamination peaks (red dots)
  - Other peaks (blue dots)
  - Connecting lines between peaks
- Hover information for each peak:
  - Allele value
  - Peak height
  - Size (bp)
  - Dye channel and limits
- Cluster distance display for contaminated markers
- Easy navigation between markers

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

## Algorithm Details

### Peak Analysis

1. **Height-based Filtering**
   - Apply dye-specific thresholds from marker_info.json
   - Filter peaks based on minimum and maximum height limits
   - Normalize peak heights within valid range

2. **Main Profile Selection**
   - Select top 2 peaks by height as main profile
   - These peaks represent the primary genotype
   - Annotated in green in visualizations

3. **Contamination Detection**
   - Calculate height standard deviation
   - Perform height-based clustering analysis
   - Compute relative distance between clusters:
     ```
     relative_distance = cluster_distance / height_std_dev
     ```
   - Identify contamination peaks based on cluster analysis
   - Mark contamination peaks in red in visualizations

4. **Quality Control**
   - Monitor peak height distribution
   - Track contamination rates
   - Calculate height standard deviation
   - Apply dye-specific quality thresholds

## Support

For support and questions, please create an issue in the repository.

## License

This project is licensed under the MIT License - see the LICENSE file for details.