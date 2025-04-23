# FASTER

**F**orensic **A**nalysis of **ST**Rs with Th**ER**mofisher Electrophoresis Result

A robust tool for analyzing Short Tandem Repeat (STR) data from Thermofisher electrophoresis results.

## Features

- **Peak Analysis**
  - Height-based peak detection with dye-specific thresholds
  - Main profile selection (top 2 peaks by height)
  - Support for multiple dye channels (B, G, Y, R, P)

- **Contamination Detection**
  - Height-based clustering analysis
  - Relative distance calculation between clusters
  - Automatic contamination peak identification

- **Visualization**
  - Interactive HTML reports with plotly graphs
  - Static PNG plots (optional)
  - Detailed peak information display

## Installation

### Prerequisites
- Python 3.10 or higher
- pip package manager

### Setup
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

### Options
- `-i, --input`: Input data file (tab-separated)
- `-o, --output`: Output directory
- `--config`: Path to marker configuration file (optional)
- `--plot`: Generate static PNG plots
- `--plotly`: Generate interactive plots (default: True)

### Example
```bash
# Basic usage
faster -i example/input.txt -o example_out/

# With custom configuration
faster -i example/input.txt -o example_out/ --config path/to/marker_info.json
```

## Output Files

```
output_directory/
├── {sample_name}.STR_analysis.json    # Analysis results
├── {sample_name}.STR_report.html      # Interactive report
└── {sample_name}_peaks/               # Static plots (optional)
    └── {sample_name}_{marker}_peaks.png
```

## Algorithm

### Peak Analysis
1. **Height-based Filtering**
   - Apply dye-specific thresholds
   - Filter peaks based on height limits

2. **Main Profile Selection**
   - Select top 2 peaks by height
   - Represent primary genotype

3. **Contamination Detection**
   - Calculate height standard deviation
   - Perform clustering analysis
   - Compute relative distance:
     ```
     relative_distance = cluster_distance / height_std_dev
     ```

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

## Support
For support and questions, please create an issue in the repository.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
