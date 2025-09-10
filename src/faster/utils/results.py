import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
import numpy as np
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

class ResultGenerator:
    """Generates standardized JSON results for STR analysis."""
    
    def __init__(self, config_path: Optional[str | Path] = None):
        """
        Initialize the result generator.
        
        Args:
            config_path: Path to marker configuration file (JSON). If None, uses default config.
        """
        # Load marker configuration
        if config_path is None:
            config_path = Path(__file__).parent.parent / 'config' / 'marker_info.json'
        
        with open(config_path) as f:
            config = json.load(f)
            
        self.marker_info = {}
        for marker, info in config['markers'].items():
            pos = info['position'].split(':')
            chr_pos = pos[0]
            start_end = pos[1].split('-')
            self.marker_info[marker] = {
                "chr": chr_pos,
                "start": int(start_end[0]),
                "end": int(start_end[1]),
                "motif": info['motif']
            }
            
        self.dye_cutoffs = config['dye_cutoffs']
        self.sex_markers_config = config.get('sex_markers', {}) # Load sex_markers configuration

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

    def _determine_sex(self, peaks_by_marker: Dict[str, pd.DataFrame]) -> str:
        """
        Determine sex based on the presence and genotype of sex-specific markers.

        Args:
            peaks_by_marker: Dictionary mapping marker names to peak DataFrames.

        Returns:
            str: "Male", "Female", or "Uncertain".
        """
        sex_calls = []

        for marker_name, marker_props in self.sex_markers_config.items():
            marker_data = peaks_by_marker.get(marker_name)
            
            if marker_props['type'] == 'Y':
                if marker_data is not None and not marker_data.empty:
                    sex_calls.append("Male") # Presence of Y marker indicates Male
                else:
                    sex_calls.append("Female") # Absence of Y marker indicates Female
            elif marker_props['type'] == 'XY' and marker_name == "AMEL": # Specifically for AMEL
                if marker_data is not None and not marker_data.empty:
                    # Simplified genotype check for AMEL based on allele names
                    # Assumes AMEL alleles are 'X' and 'Y' or similar.
                    # A more robust check might involve looking at the actual genotype call if available.
                    alleles = set(marker_data['allele'].astype(str).str.upper())
                    if 'Y' in alleles and 'X' in alleles:
                        sex_calls.append("Male")
                    elif 'X' in alleles and 'Y' not in alleles : #Only X is present
                        sex_calls.append("Female")
                    else:
                        sex_calls.append("Uncertain")
                    # If only Y is present (unlikely for AMEL) or other combinations, it could be uncertain or an anomaly.
                    # For now, if not clearly Male or Female, we don't add a call, letting the consensus logic handle it.

        if not sex_calls or peaks_by_marker.get("AMEL") is None:
            return "Uncertain" # No sex marker data found

        # Check for consensus
        if all(call == "Male" for call in sex_calls):
            return "Male"
        elif all(call == "Female" for call in sex_calls):
            return "Female"
        else:
            return "Uncertain" # Conflicting or insufficient data among markers

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
        
        # Determine sex
        determined_sex = self._determine_sex(peaks_by_marker)
        
        results = {
            "LocusResults": {},
            "SampleParameters": {
                "SampleId": clean_sample_name,
                "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "sex": determined_sex # Add determined sex
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
                "peaks": {
                    "is_contaminated": False,
                    "main_profile_peaks": [],
                    "contamination_peaks": []
                }
            }
            
            # Add contamination information if present
            if marker in contamination_by_marker:
                contamination_result = contamination_by_marker[marker]
                # Unpack tuple return value (contamination_info, join_points)
                contamination = contamination_result[0] if isinstance(contamination_result, tuple) else contamination_result
                
                # Always add main profile information if available
                if contamination and contamination.main_profile_peaks:
                    variant_info["peaks"] = {
                        "is_contaminated": contamination.is_contaminated if contamination else False,
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
                        ] if contamination and contamination.is_contaminated else []
                    }
                    
                    # Add to contaminated markers list only if actually contaminated
                    if contamination and contamination.is_contaminated:
                        contaminated_markers.append({
                            "marker": marker,
                            "main_profile": "/".join(p.allele for p in contamination.main_profile_peaks),
                            "contamination_peaks": ", ".join(f"{p.allele}({p.relative_height:.1f}%)" 
                                                           for p in contamination.contamination_peaks),
                            "relative_distance": contamination.relative_distance
                        })
                else:
                    # If no contamination info but we have peaks, use them as main profile
                    variant_info["peaks"] = {
                        "is_contaminated": False,
                        "main_profile_peaks": [
                            {
                                "allele": p["allele"],
                                "height": float(p["height"]),
                                "size": float(p["size"]),
                                "relative_height": float(p["height"] / sorted_peaks[0]["height"]) if sorted_peaks else 0
                            }
                            for p in sorted_peaks
                        ],
                        "contamination_peaks": []
                    }
            else:
                # If no contamination info but we have peaks, use them as main profile
                variant_info["peaks"] = {
                    "is_contaminated": False,
                    "main_profile_peaks": [
                        {
                            "allele": p["allele"],
                            "height": float(p["height"]),
                            "size": float(p["size"]),
                            "relative_height": float(p["height"] / sorted_peaks[0]["height"]) if sorted_peaks else 0
                        }
                        for p in sorted_peaks
                    ],
                    "contamination_peaks": []
                }
            
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
    
    def save_results(self, results: Dict, output_dir: str):
        """Save analysis results to output directory.
        
        Args:
            results: Analysis results dictionary
            output_dir: Output directory path
        """
        # Create output directory if it doesn't exist
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save JSON results
        sample_id = results["SampleParameters"]["SampleId"]
        json_path = output_dir / f"{sample_id}.STR_analysis.json"
        
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Results saved to: {output_dir}")
        logger.info("---") 