from dataclasses import dataclass
from typing import Dict, List, Optional

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