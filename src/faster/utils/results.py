import json
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np
from datetime import datetime

class ResultGenerator:
    """Generates standardized JSON results for STR analysis."""
    
    def __init__(self):
        """Initialize the result generator."""
        self.marker_info = {
            "CSF1PO": {"chr": "chr5", "start": 150076323, "end": 150076375, "motif": "[ATCT]*"},
            "D10S1248": {"chr": "chr10", "start": 129294243, "end": 129294295, "motif": "[GGAA]*"},
            "D12S391": {"chr": "chr12", "start": 12297019, "end": 12297095, "motif": "[AGAT]+[AGAC]+AGAT"},
            "D13S317": {"chr": "chr13", "start": 82148024, "end": 82148068, "motif": "[TATC]*"},
            "D16S539": {"chr": "chr16", "start": 86352701, "end": 86352745, "motif": "[GATA]*"},
            "D18S51": {"chr": "chr18", "start": 63281666, "end": 63281738, "motif": "[AGAA]*"},
            "D19S433": {"chr": "chr19", "start": 29926234, "end": 29926298, "motif": "[CCTT]*cctaCCTTctttCCTT"},
            "D1S1656": {"chr": "chr1", "start": 230769615, "end": 230769683, "motif": "CCTA[TCTA]*"},
            "D21S11": {"chr": "chr21", "start": 19181972, "end": 19182099, "motif": "[TCTA]+[TCTG]+[TCTA]+ta[TCTA]+tca[TCTA]+tccata[TCTA]+"},
            "D22S1045": {"chr": "chr22", "start": 37140286, "end": 37140337, "motif": "[ATT]+ACT[ATT]+"},
            "D2S1338": {"chr": "chr2", "start": 218014858, "end": 218014950, "motif": "[GGAA]+GGAC[GGAA]+[GGCA]+"},
            "D2S441": {"chr": "chr2", "start": 68011947, "end": 68011994, "motif": "[TCTA]*"},
            "D3S1358": {"chr": "chr3", "start": 45540738, "end": 45540802, "motif": "TCTATCTG[TCTA]*"},
            "D5S818": {"chr": "chr5", "start": 123775555, "end": 123775599, "motif": "[ATCT]*"},
            "D7S820": {"chr": "chr7", "start": 84160225, "end": 84160277, "motif": "[TATC]*"},
            "D8S1179": {"chr": "chr8", "start": 124894864, "end": 124894916, "motif": "TCTATCTG[TCTA]*"},
            "FGA": {"chr": "chr4", "start": 154587735, "end": 154587823, "motif": "[GGAA]+GGAG[AAAG]+AGAAAAAA[GAAA]+"},
            "SE33": {"chr": "chr6", "start": 88277143, "end": 88277245, "motif": "[CTTT]+TT[CTTT]+"},
            "TH01": {"chr": "chr11", "start": 2171087, "end": 2171115, "motif": "[AATG]*"},
            "TPOX": {"chr": "chr2", "start": 1489652, "end": 1489684, "motif": "[AATG]*"},
            "vWA": {"chr": "chr12", "start": 5983976, "end": 5984044, "motif": "[TAGA]*[CAGA]*TAGA"},
            "AMEL": {"chr": "chrX", "start": 11293412, "end": 11300761, "motif": "null"}
        }
        
        # Dye-specific height cutoffs
        self.dye_cutoffs = {
            "B": {"min": 2500, "max": 50000},  # Blue
            "G": {"min": 5000, "max": 50000},  # Green
            "Y": {"min": 9000, "max": 50000},  # Yellow
            "R": {"min": 1000, "max": 50000},  # Red
            "P": {"min": 1000, "max": 50000},  # Purple
        }

    def _calculate_stats(self, peaks: List[Dict]) -> Dict[str, Any]:
        """Calculate statistics for a group of peaks."""
        heights = [p["height"] for p in peaks]
        return {
            "allele_count": len(peaks),
            "median_height": float(np.median(heights)) if heights else None,
            "std_height": float(np.std(heights)) if len(heights) > 1 else None
        }

    def _get_height_limits(self, dye: str) -> Dict[str, int]:
        """Get height limits based on dye color."""
        return self.dye_cutoffs.get(dye, {"min": 1000, "max": 50000})

    def _get_variant_key(self, marker: str) -> str:
        """Get the chromosome position key for a marker."""
        info = self.marker_info.get(marker, {})
        if info:
            return f"{info['chr']}_{info['start']}_{info['end']}"
        return "unknown_position"

    def _get_motif(self, marker: str) -> str:
        """Get the repeat motif for a marker."""
        return self.marker_info.get(marker, {}).get('motif', "[ATCT]*")

    def generate_results(self,
                        peaks_by_marker: Dict[str, pd.DataFrame],
                        contamination_by_marker: Dict[str, Any],
                        sample_name: str) -> Dict[str, Any]:
        """
        Generate a results dictionary for a sample.
        
        Args:
            peaks_by_marker: Dictionary mapping marker names to peak DataFrames
            contamination_by_marker: Dictionary mapping marker names to contamination info
            sample_name: Name of the sample
            
        Returns:
            Dictionary containing the analysis results
        """
        # Clean sample name - remove instrument ID and file extension
        clean_sample_name = sample_name.split('_AC')[0].split('.')[0]
        
        # Track contamination information
        contaminated_markers = []
        valid_markers = []
        
        results = {
            "LocusResults": {},
            "SampleParameters": {
                "SampleId": clean_sample_name,
                "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        }
        
        for marker, peaks_df in peaks_by_marker.items():
            if peaks_df.empty:
                continue
                
            valid_markers.append(marker)
            peaks = []
            for _, peak in peaks_df.iterrows():
                peak_info = {
                    "allele": peak["allele"],
                    "height": float(peak["height"]),
                    "size": round(float(peak["size"]), 2)
                }
                peaks.append(peak_info)
            
            # Calculate statistics
            stats = self._calculate_stats(peaks)
            dye = peaks_df["dye"].iloc[0] if not peaks_df.empty else None
            
            # Create variant info
            variant_key = self._get_variant_key(marker)
            
            # Sort peaks by height and get top 2 alleles for genotype
            sorted_peaks = sorted(peaks, key=lambda x: x["height"], reverse=True)
            top_alleles = [p["allele"] for p in sorted_peaks[:2]]
            genotype = "/".join(sorted(top_alleles)) if len(top_alleles) > 1 else top_alleles[0]
            
            variant_info = {
                "genotype": genotype,
                "allele_count": stats["allele_count"],
                "motif": self._get_motif(marker),
                "peaks": peaks
            }
            
            # Add contamination information if present
            if marker in contamination_by_marker:
                contamination = contamination_by_marker[marker]
                if contamination and contamination.is_contaminated:
                    variant_info["contamination"] = {
                        "is_contaminated": True,
                        "main_profile_peaks": [
                            {
                                "allele": p.allele,
                                "height": float(p.height),
                                "size": float(p.size),
                                "relative_height": float(p.relative_height)
                            }
                            for p in contamination.main_profile_peaks
                        ],
                        "contamination_peaks": [
                            {
                                "allele": p.allele,
                                "height": float(p.height),
                                "size": float(p.size),
                                "relative_height": float(p.relative_height)
                            }
                            for p in contamination.contamination_peaks
                        ],
                        "relative_distance": contamination.relative_distance
                    }
                    
                    # Add to contaminated markers list
                    contaminated_markers.append({
                        "marker": marker,
                        "main_profile": "/".join(p.allele for p in contamination.main_profile_peaks),
                        "contamination_peaks": ", ".join(f"{p.allele}({p.relative_height:.1f}%)" 
                                                       for p in contamination.contamination_peaks),
                        "relative_distance": contamination.relative_distance
                    })
                else:
                    variant_info["contamination"] = None
            else:
                variant_info["contamination"] = None
            
            # Create marker results
            marker_results = {
                "allele_count": stats["allele_count"],
                "median_height": stats["median_height"],
                "dye": dye,
                "std_height": stats["std_height"],
                "height_limits": self._get_height_limits(dye),
                "variants": {
                    variant_key: variant_info
                }
            }
            
            results["LocusResults"][marker] = marker_results
        
        # Add sample contamination summary
        total_valid_markers = len(valid_markers)
        total_contaminated = len(contaminated_markers)
        results["SampleContamination"] = {
            "contamination_rate": round(total_contaminated / total_valid_markers * 100, 1) if total_valid_markers > 0 else 0.0,
            "contaminated_markers": contaminated_markers,
            "total_valid_markers": total_valid_markers,
            "total_contaminated_markers": total_contaminated
        }
        
        return results
    
    def save_results(self,
                    results: Dict[str, Any],
                    output_dir: str | Path) -> None:
        """
        Save results to a JSON file.
        
        Args:
            results: Results dictionary to save
            output_dir: Directory to save the results in
        """
        output_dir = Path(output_dir)
        output_path = output_dir / f"{results['SampleParameters']['SampleId']}.STR_analysis.json"
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2) 