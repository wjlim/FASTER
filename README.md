# FASTER

## **F**orensic **A**nalysis of **S**hort **T**and**E**m **R**epeats

A robust tool for analyzing Short Tandem Repeat (STR) data from Thermofisher electrophoresis and Illumina NGS results.

## Features

- **Peak Analysis**
  - Height-based peak detection with dye-specific thresholds
  - Main profile selection (top 2 peaks by height or all clustered alleles)
  - Support for multiple dye channels (B, G, Y, R, P)

- **Contamination Detection**
  - Height-based clustering analysis
  - Relative distance calculation between clusters
  - Automatic contamination peak identification

- **Visualization**
  - Interactive HTML reports with plotly graphs
  - Static PNG plots (optional)
  - Detailed peak information display

- **ExpansionHunter Integration**
  - Analysis of BAM files using ExpansionHunter
  - Comparison with STR analysis results
  - Combined result reporting
  - **Marker filtering:** Only markers listed in `compare_markers` in `marker_info.json` are used for comparison.

- **Genotype Vectorization**
  - Convert genotypes to compact vector representation
  - Support for both STR and ExpansionHunter results
  - Polar and Cartesian coordinate systems
  - Vector comparison and similarity scoring

- **Tabular Output (CSV/Excel)**
  - Automatic export of genotype, main profile, and contamination tables in CSV format after STR analysis
  - Genotype table also exported as Excel file (`.xlsx`)
  - Genotype values are quoted to prevent Excel from interpreting them as dates

## Installation

### Prerequisites
- Python 3.10 or higher
- pip package manager
- For Excel output: `openpyxl` (automatically installed via requirements)

### Setup
```bash
git clone https://github.com/wjlim/FASTER.git
cd FASTER
pip install .
```

## Usage

### STR Analysis
```bash
# Basic STR analysis
faster str -i <input_file> -o <output_directory>

# With custom configuration
faster str -i <input_file> -o <output_directory> --config path/to/marker_info.json
```

### ExpansionHunter Analysis
```bash
# Run ExpansionHunter analysis
faster exhunter -i <input_bam> -r <reference_fasta> -o <output_prefix>
```

### Compare Results
```bash
# Compare STR and ExpansionHunter results
faster compare -i <str_json> -j <eh_json> -o <output_prefix>
```

### Vectorization
```bash
# Vectorize STR results
faster vectorize -i <input_json> -o <output_file> -t str

# Vectorize ExpansionHunter results
faster vectorize -i <input_json> -o <output_file> -t eh
```

### Vector Comparison
```bash
# Compare two vectors
faster compare-vectors -i <vector1.json> -j <vector2.json>

# Compare and save results to JSON
faster compare-vectors -i <vector1.json> -j <vector2.json> -o <comparison.json>
```

## Output Formats by Submodule

### 1. STR Analysis (`faster str`)
Output files:
- `{sample_name}.STR_analysis.json`: Full analysis results (see below)
- `STR_analysis.genotype.csv`: Table of genotypes (rows: sample_name, columns: marker name, values: genotype, quoted to prevent Excel date conversion)
- `STR_analysis.main_profile.csv`: Table of main profiles (rows: sample_name, columns: marker name, values: main profile alleles)
- `STR_analysis.contamn.csv`: Table of contamination status (rows: sample_name, columns: marker name, values: 1 for contaminated, blank otherwise)
- `STR_analysis.genotype.xlsx`: Excel version of the genotype table (values quoted for Excel safety)

Example of `{sample_name}.STR_analysis.json`:
```json
{
  "LocusResults": {
    "marker_name": {
      "allele_count": int,
      "median_height": float,
      "dye": str,
      "std_height": float or null,
      "height_limits": {
        "min": int,
        "max": int
      },
      "variants": {
        "position": {
          "genotype": str,
          "allele_count": int,
          "motif": str,
          "peaks": [
            {
              "allele": str,
              "height": float,
              "size": float
            },
            // ...
          ],
          "contamination": {
            "is_contaminated": bool,
            "main_profile_peaks": [
              {
                "allele": str,
                "height": float,
                "size": float,
                "relative_height": float
              },
              // ...
            ],
            "contamination_peaks": [
              {
                "allele": str,
                "height": float,
                "size": float,
                "relative_height": float
              },
              // ...
            ],
            "relative_distance": float
          } or null
        }
      }
    }
  },
  "SampleParameters": {
    "SampleId": str,
    "analysis_date": str
  },
  "SampleContamination": {
    "contamination_rate": float,
    "contaminated_markers": [
      {
        "marker": str,
        "main_profile": str,
        "contamination_peaks": str,
        "relative_distance": float
      },
      // ...
    ],
    "total_valid_markers": int,
    "total_contaminated_markers": int
  }
}
```

### 2. ExpansionHunter Analysis (`faster exhunter`)
Output files:
- `{output_prefix}.vcf`: Variant calls in VCF format
- `{output_prefix}.json`: Detailed analysis results
```json
{
  "LocusResults": {
    "marker_name": {
      "Variants": {
        "variant_id": {
          "Genotype": str,
          "RepeatUnit": str,
          "RepeatSize": int
        }
      }
    }
  },
  "SampleParameters": {
    "SampleId": str,
    "Gender": str
  }
}
```

### 3. Results Comparison (`faster compare`)
- **Now only compares markers listed in `compare_markers` in `marker_info.json`.**
Output file: `{output_prefix}.comparison.json`
```json
{
  "sample_id": str,
  "matching_markers": [...],
  "mismatching_markers": [...],
  "missing_markers": [...],
  "summary": {
    "total_markers": int,
    "matching": int,
    "mismatching": int,
    "missing": int,
    "match_ratio": float
  }
}
```

### 4. Genotype Vectorization (`faster vectorize`)
Output file: `{output_file}.json`
```json
{
  "sample_id": str,
  "source_type": str,
  "markers": [
    {
      "marker": str,
      "allele1": float,
      "allele2": float
    }
  ],
  "vector_properties": {
    "magnitude": float,
    "angle_radians": float,
    "angle_degrees": float
  },
  "markers_used": int,
  "compact_vector": str,
  "cartesian_coordinates": {
    "x": float,
    "y": float
  }
}
```

### 5. Vector Comparison (`faster compare-vectors`)
Output file: `{output_file}.json`
```json
{
  "vector1": {
    "sample_id": str,
    "source_type": str,
    "markers": [...],
    "vector_properties": {...},
    "cartesian_coordinates": {...}
  },
  "vector2": {
    "sample_id": str,
    "source_type": str,
    "markers": [...],
    "vector_properties": {...},
    "cartesian_coordinates": {...}
  },
  "comparison_metrics": {
    "euclidean_distance": float,
    "magnitude_difference": float,
    "angle_difference_radians": float,
    "angle_difference_degrees": float,
    "similarity_score": float
  }
}
```

Note: All floating-point values in the output are rounded to 6 decimal places for consistency.

## Support
For support and questions, please create an issue in the repository.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

Note: This project uses ExpansionHunter binary which is licensed under the PolyForm Strict License 1.0.0. For more details about ExpansionHunter's license, please visit [ExpansionHunter's repository](https://github.com/Illumina/ExpansionHunter).
