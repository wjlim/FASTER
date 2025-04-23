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

- **ExpansionHunter Integration**
  - Analysis of BAM files using ExpansionHunter
  - Comparison with STR analysis results
  - Combined result reporting

- **Genotype Vectorization**
  - Convert genotypes to compact vector representation
  - Support for both STR and ExpansionHunter results
  - Polar and Cartesian coordinate systems
  - Vector comparison and similarity scoring

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
Output file: `{sample_name}.STR_analysis.json`
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
  "SampleParameters": {
    "sample_id": str,
    "analysis_date": str,
    "sample_name": str,
    "contamination_summary": {
      "contaminated_markers": [...],
      "total_markers": int,
      "contamination_percentage": float,
      "mean_contamination_rate": float
    }
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
Output file: `{output_prefix}.comparison.json`
```json
{
  "sample_id": str,
  "str_results": {
    "marker_name": {
      "genotype": str,
      "allele_count": int
    }
  },
  "eh_results": {
    "marker_name": {
      "genotype": str,
      "allele_count": int
    }
  },
  "comparison_summary": {
    "matching_markers": [...],
    "mismatched_markers": [...],
    "concordance_rate": float
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