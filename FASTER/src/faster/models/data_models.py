from dataclasses import dataclass
from typing import Dict, List, Optional
from pydantic import BaseModel, Field

@dataclass
class PeakInfo:
    allele: str
    height: float
    size: float
    relative_height: float

@dataclass
class ContaminationInfo:
    is_contaminated: bool
    main_profile_peaks: List[PeakInfo]
    contamination_peaks: List[PeakInfo]
    relative_distance: float

@dataclass
class VariantInfo:
    genotype: str
    median_height: float
    std_height: float
    dye: str
    allele_count: int
    contamination: Optional[Dict] = None
    peaks: List[Dict] = None

@dataclass
class LocusResult:
    allele_count: int
    median_height: float
    dye: str
    std_height: float
    variants: Dict[str, VariantInfo]

@dataclass
class ContaminationSummary:
    contaminated_markers: List[str]
    total_markers: int
    contamination_percentage: float
    mean_contamination_rate: float

@dataclass
class SampleParameters:
    sample_id: str
    analysis_date: str
    sample_name: str
    contamination_summary: Optional[ContaminationSummary] = None

@dataclass
class AnalysisResult:
    sample_name: str
    locus_results: Dict[str, LocusResult]
    sample_parameters: SampleParameters

class MarkerGenotype(BaseModel):
    """Model for a single marker's genotype."""
    marker: str = Field(..., description="Marker name")
    allele1: float = Field(..., description="First allele (smaller number)")
    allele2: float = Field(..., description="Second allele (larger number)")

class VectorProperties(BaseModel):
    """Model for vector properties in polar coordinates."""
    magnitude: float = Field(..., description="Magnitude (distance from origin)")
    angle_radians: float = Field(..., description="Angle in radians")
    angle_degrees: float = Field(..., description="Angle in degrees")

class CartesianCoordinates(BaseModel):
    """Model for Cartesian (x,y) coordinates."""
    x: float = Field(..., description="X coordinate")
    y: float = Field(..., description="Y coordinate")

class VectorizationResult(BaseModel):
    """Model for vectorization results."""
    sample_id: str = Field(..., description="Sample identifier")
    source_type: str = Field(..., description="Source type (str or eh)")
    markers: List[MarkerGenotype] = Field(..., description="List of marker genotypes")
    vector_properties: VectorProperties = Field(..., description="Vector properties")
    markers_used: int = Field(..., description="Number of markers used in vectorization")
    compact_vector: str = Field(..., description="Compact vector representation (magnitude,angle)")
    cartesian_coordinates: CartesianCoordinates = Field(None, description="Cartesian coordinates (x,y)")
    
    class Config:
        json_schema_extra = {
            "example": {
                "sample_id": "NA12878",
                "source_type": "str",
                "markers": [
                    {"marker": "D3S1358", "allele1": 15.0, "allele2": 16.0},
                    {"marker": "vWA", "allele1": 17.0, "allele2": 18.0}
                ],
                "vector_properties": {
                    "magnitude": 101.019800,
                    "angle_radians": 0.847139,
                    "angle_degrees": 48.537494
                },
                "markers_used": 21,
                "compact_vector": "101.019800,0.847139",
                "cartesian_coordinates": {
                    "x": 66.789123,
                    "y": 75.912345
                }
            }
        }

class VectorComparisonResult(BaseModel):
    """Model for vector comparison results."""
    vector1: VectorizationResult = Field(..., description="First vector result")
    vector2: VectorizationResult = Field(..., description="Second vector result")
    euclidean_distance: float = Field(..., description="Euclidean distance between vectors")
    magnitude_difference: float = Field(..., description="Absolute difference in magnitudes")
    angle_difference_radians: float = Field(..., description="Angular difference in radians")
    angle_difference_degrees: float = Field(..., description="Angular difference in degrees")
    similarity_score: float = Field(..., description="Normalized similarity score (0-1)")
    
    class Config:
        json_schema_extra = {
            "example": {
                "vector1": {
                    "sample_id": "NA12878-1",
                    "source_type": "str",
                    "markers": [
                        {"marker": "D3S1358", "allele1": 15.0, "allele2": 16.0}
                    ],
                    "vector_properties": {
                        "magnitude": 101.019800,
                        "angle_radians": 0.847139,
                        "angle_degrees": 48.537494
                    },
                    "markers_used": 21,
                    "compact_vector": "101.019800,0.847139",
                    "cartesian_coordinates": {
                        "x": 66.789123,
                        "y": 75.912345
                    }
                },
                "vector2": {
                    "sample_id": "NA12878-2",
                    "source_type": "eh",
                    "markers": [
                        {"marker": "D3S1358", "allele1": 15.0, "allele2": 16.0}
                    ],
                    "vector_properties": {
                        "magnitude": 98.123456,
                        "angle_radians": 0.785398,
                        "angle_degrees": 45.000000
                    },
                    "markers_used": 21,
                    "compact_vector": "98.123456,0.785398",
                    "cartesian_coordinates": {
                        "x": 69.456789,
                        "y": 69.456789
                    }
                },
                "euclidean_distance": 5.123456,
                "magnitude_difference": 2.896344,
                "angle_difference_radians": 0.061741,
                "angle_difference_degrees": 3.537494,
                "similarity_score": 0.987654
            }
        } 