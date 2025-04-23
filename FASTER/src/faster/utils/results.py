import json
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np
from datetime import datetime

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