import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import logging
import math
from ..models.data_models import (
    MarkerGenotype, 
    VectorProperties, 
    VectorizationResult,
    VectorComparisonResult,
    CartesianCoordinates
)

logger = logging.getLogger(__name__)

class GenotypeVectorizer:
    def __init__(self, config_path: Optional[str] = None):
        """Initialize the GenotypeVectorizer class.
        Loads the list of markers and normalization/weighting info from config.
        
        Args:
            config_path: Optional path to marker configuration file
        """
        if config_path is None:
            config_path = str(Path(__file__).parent.parent / 'config' / 'marker_info.json')
            
        try:
            with open(config_path) as f:
                config = json.load(f)
        except FileNotFoundError:
             logger.error(f"Configuration file not found at: {config_path}")
             raise
        except json.JSONDecodeError:
            logger.error(f"Error decoding JSON from configuration file: {config_path}")
            raise
            
        self.vectorize_markers = config.get('vectorize_markers', []) 
        if not self.vectorize_markers:
            logger.warning("No markers specified for vectorization ('vectorize_markers') in config file.")
        
        # Load normalization ranges from the config file
        self.marker_ranges = config.get('normalization_ranges', {}) 
        if not self.marker_ranges:
             logger.warning("Normalization ranges ('normalization_ranges') not found in config file. Using default range [0, 50] for normalization.")

        # Load marker weights from config file, ensuring all vectorize_markers have a weight
        default_weight = 1.0
        loaded_weights = config.get('marker_weights', {}) # Load weights if they exist
        self.marker_weights = {}
        for marker in self.vectorize_markers:
            self.marker_weights[marker] = loaded_weights.get(marker, default_weight)
            
        self.vector_size = 2  # Each marker has 2 alleles
        
    def _parse_genotype(self, genotype: str) -> Optional[Tuple[float, float]]:
        """Parse genotype string into two allele numbers, always sorted in ascending order.
        Handles potential decimal values by taking the integer part.
        Returns None if parsing fails or alleles are non-numeric.
        
        Args:
            genotype: Genotype string (e.g., "10/8", "15.3", "14/14.1")
            
        Returns:
            Tuple of (allele1, allele2) as floats, sorted, or None.
        """
        try:
            if '/' in genotype:
                alleles_str = genotype.split('/')
                # Take integer part before converting to float
                allele1 = float(int(float(alleles_str[0]))) 
                allele2 = float(int(float(alleles_str[1])))
            else:
                # Single allele case - homozygous
                allele1 = allele2 = float(int(float(genotype)))
                
            return tuple(sorted([allele1, allele2]))
        except (ValueError, IndexError):
            # Handle parsing errors or non-numeric alleles
            logger.warning(f"Could not parse genotype: {genotype}")
            return None

    def _combine_exhunter_genotypes(self, variants: Dict) -> Optional[Tuple[float, float]]:
        """Combine multiple ExpansionHunter variant genotypes into a single average genotype.
        
        Args:
            variants: Dictionary of variants from ExpansionHunter
            
        Returns:
            Combined genotype as (allele1, allele2), sorted, or None if no valid genotypes found.
        """
        valid_alleles = []
        for variant in variants.values():
            genotype_str = variant.get('Genotype')
            if genotype_str:
                parsed_alleles = self._parse_genotype(genotype_str)
                if parsed_alleles:
                    valid_alleles.extend(parsed_alleles)
        
        if not valid_alleles:
            return None
            
        # Simple average for demonstration; could be more sophisticated
        # This logic might need refinement depending on how EH genotypes should be combined
        # For now, just take the min and max observed allele values across variants
        if len(valid_alleles) >= 2:
             return tuple(sorted([min(valid_alleles), max(valid_alleles)]))
        elif len(valid_alleles) == 1:
             return tuple(sorted([valid_alleles[0], valid_alleles[0]])) # Homozygous
        else: 
             return None

    def _polar_to_cartesian(self, magnitude: float, angle_radians: float) -> CartesianCoordinates:
        """Convert polar coordinates to Cartesian coordinates.
        
        Args:
            magnitude: Distance from origin (r)
            angle_radians: Angle in radians (θ)
            
        Returns:
            CartesianCoordinates object with x and y values
        """
        x = magnitude * np.cos(angle_radians)
        y = magnitude * np.sin(angle_radians)
        return CartesianCoordinates(x=round(float(x), 6), y=round(float(y), 6))

    def _cartesian_to_polar(self, x: float, y: float) -> Tuple[float, float]:
        """Convert Cartesian coordinates to polar coordinates.
        
        Args:
            x: X coordinate
            y: Y coordinate
            
        Returns:
            Tuple of (magnitude, angle_radians)
        """
        magnitude = np.sqrt(x**2 + y**2)
        angle = np.arctan2(y, x)
        # Ensure angle is positive (0 to 2π)
        if angle < 0:
            angle += 2 * np.pi
        return round(float(magnitude), 6), round(float(angle), 6)

    def _normalize_allele(self, allele: float, marker: str) -> float:
        """Normalize allele value to range [0,1] based on marker-specific ranges.
        Clips values outside the defined range to 0 or 1.
        Uses default range [0, 50] if marker range is missing or invalid in config.
        
        Args:
            allele: Raw allele value (integer part)
            marker: Marker name
            
        Returns:
            Normalized allele value between 0 and 1
        """
        marker_range = self.marker_ranges.get(marker)
        
        # Validate the loaded range
        if isinstance(marker_range, (list, tuple)) and len(marker_range) == 2 and all(isinstance(x, (int, float)) for x in marker_range):
            min_val, max_val = marker_range
        else:
            # Use default if marker not found or range is invalid
            if marker_range is not None: # Log only if invalid format, not just missing
                 logger.warning(f"Invalid or missing range for marker {marker} in config: {marker_range}. Using default [0, 50].")
            min_val, max_val = 0, 50
        
        if max_val == min_val: # Avoid division by zero
            return 0.5 # Midpoint if range is zero
            
        normalized = (allele - min_val) / (max_val - min_val)
        
        # Clip values to be strictly within [0, 1]
        return max(0.0, min(1.0, normalized))

    def _process_marker_genotype(self, marker: str, allele1: float, allele2: float) -> List[float]:
        """Process a marker's genotype with normalization and weighting.
        
        Args:
            marker: Marker name
            allele1: First allele value (integer part)
            allele2: Second allele value (integer part)
            
        Returns:
            List of processed values [normalized_weighted_allele1, normalized_weighted_allele2]
        """
        # Normalize alleles
        norm_allele1 = self._normalize_allele(allele1, marker)
        norm_allele2 = self._normalize_allele(allele2, marker)
        
        # Apply marker-specific weight
        weight = self.marker_weights.get(marker, 1.0)
        
        # Return weighted and normalized values
        return [norm_allele1 * weight, norm_allele2 * weight]

    def _calculate_vector_properties(self, vector: np.ndarray) -> VectorProperties:
        """Calculate magnitude and angle of the processed vector.
        Magnitude is the L2 norm of the entire flattened vector.
        Angle is the simple mean angle of allele pairs treated as 2D points (x=processed_allele1, y=processed_allele2).
        """
        if vector.size == 0:
            return VectorProperties(magnitude=0.0, angle_radians=0.0, angle_degrees=0.0)
            
        # Calculate L2 norm (magnitude) of the entire flattened vector
        magnitude = float(np.linalg.norm(vector))
        
        # Reshape into nx2 matrix for angle calculation
        points = vector.reshape(-1, 2)
        
        # Avoid calculating angles if all points are at the origin (e.g., after normalization/weighting)
        if np.all(points == 0):
             mean_angle = 0.0
        else:
            # Calculate angles for each point (processed allele pair)
            # arctan2(y, x) where y=processed_allele2, x=processed_allele1
            angles = np.arctan2(points[:, 1], points[:, 0])
            # Ensure angles are positive (0 to 2π)
            angles = np.where(angles < 0, angles + 2*np.pi, angles)
            # Calculate the simple mean angle
            mean_angle = float(np.mean(angles))
            
        return VectorProperties(
            magnitude=round(magnitude, 6),
            angle_radians=round(mean_angle, 6),
            angle_degrees=round(math.degrees(mean_angle), 6)
        )

    def vectorize_str(self, data: Dict, sample_id: str) -> VectorizationResult:
        """Vectorize STR analysis results using normalized and weighted allele values.
        
        Args:
            data: Dictionary containing STR analysis results
            sample_id: Sample identifier
            
        Returns:
            VectorizationResult object
        """
        marker_genotypes = []
        processed_vectors = [] # Store processed (normalized, weighted) values
        
        for marker in self.vectorize_markers:
            locus_result = data.get("LocusResults", {}).get(marker)
            if not locus_result:
                logger.warning(f"Marker {marker} not found for sample {sample_id}")
                continue
                
            variants = locus_result.get("variants", {})
            if not variants:
                logger.warning(f"No variants for marker {marker} in sample {sample_id}")
                continue
                
            variant = list(variants.values())[0]
            genotype = variant.get("genotype")
            if not genotype:
                 logger.warning(f"No genotype for marker {marker} in sample {sample_id}")
                 continue
            
            parsed_alleles = self._parse_genotype(genotype)
            
            if parsed_alleles is not None:
                allele1, allele2 = parsed_alleles
                marker_genotypes.append(
                    MarkerGenotype(
                        marker=marker,
                        allele1=allele1, # Store raw parsed alleles
                        allele2=allele2
                    )
                )
                # Process alleles (normalize + weight) and add to vector list
                processed_values = self._process_marker_genotype(marker, allele1, allele2)
                processed_vectors.extend(processed_values)
            else:
                logger.info(f"Skipping marker {marker} due to unparsable genotype: {genotype}")
        
        # Calculate properties based on the processed vector
        vector_array = np.array(processed_vectors, dtype=float)
        vector_props = self._calculate_vector_properties(vector_array)
        
        # Calculate Cartesian coordinates from the final polar properties
        cartesian = self._polar_to_cartesian(
            vector_props.magnitude,
            vector_props.angle_radians
        )
        
        return VectorizationResult(
            sample_id=sample_id,
            source_type="str",
            markers=marker_genotypes, # Stores raw alleles
            vector_properties=vector_props, # Based on processed vector
            markers_used=len(marker_genotypes),
            compact_vector=f"{vector_props.magnitude:.6f},{vector_props.angle_radians:.6f}",
            cartesian_coordinates=cartesian
        )

    def vectorize_eh(self, data: Dict, sample_id: str) -> VectorizationResult:
        """Vectorize ExpansionHunter results using normalized and weighted allele values.
        
        Args:
            data: Dictionary containing ExpansionHunter results
            sample_id: Sample identifier
            
        Returns:
            VectorizationResult object
        """
        marker_genotypes = []
        processed_vectors = [] # Store processed values
        
        for marker in self.vectorize_markers:
            locus_result = data.get("LocusResults", {}).get(marker)
            if not locus_result:
                logger.warning(f"Marker {marker} not found in EH results for {sample_id}")
                continue
                
            variants = locus_result.get("Variants")
            if not variants:
                logger.warning(f"No variants for marker {marker} in EH results for {sample_id}")
                continue
                
            combined_alleles = self._combine_exhunter_genotypes(variants)
            
            if combined_alleles is not None:
                allele1, allele2 = combined_alleles
                marker_genotypes.append(
                    MarkerGenotype(
                        marker=marker,
                        allele1=allele1, # Store raw combined alleles
                        allele2=allele2
                    )
                )
                # Process alleles (normalize + weight) and add to vector list
                processed_values = self._process_marker_genotype(marker, allele1, allele2)
                processed_vectors.extend(processed_values)
            else:
                logger.info(f"Skipping marker {marker} due to no valid combined genotype from EH variants")
        
        # Calculate properties based on the processed vector
        vector_array = np.array(processed_vectors, dtype=float)
        vector_props = self._calculate_vector_properties(vector_array)
        
        cartesian = self._polar_to_cartesian(
            vector_props.magnitude,
            vector_props.angle_radians
        )
        
        return VectorizationResult(
            sample_id=sample_id,
            source_type="eh",
            markers=marker_genotypes, # Stores raw alleles
            vector_properties=vector_props, # Based on processed vector
            markers_used=len(marker_genotypes),
            compact_vector=f"{vector_props.magnitude:.6f},{vector_props.angle_radians:.6f}",
            cartesian_coordinates=cartesian
        )

    def save_vector(self, result: VectorizationResult, output_path: str):
        """Save vectorization results to a JSON file.
        
        Args:
            result: VectorizationResult object
            output_path: Path to save the output file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            # Use pydantic's json method for better serialization
            f.write(result.model_dump_json(indent=2))

    def load_vector(self, file_path: str) -> VectorizationResult:
        """Load vectorization results from a JSON file.
        
        Args:
            file_path: Path to the vector file
            
        Returns:
            VectorizationResult object
        """
        with open(file_path) as f:
            data = json.load(f)
        # Use pydantic's parse_obj for robust loading
        return VectorizationResult.model_validate(data)

    def compare_vectors(self, vector1_path: str, vector2_path: str) -> VectorComparisonResult:
        """Compare two vectors using loaded vector data.
        
        Args:
            vector1_path: Path to first vector file
            vector2_path: Path to second vector file
            
        Returns:
            VectorComparisonResult object
        """
        vec1 = self.load_vector(vector1_path)
        vec2 = self.load_vector(vector2_path)
        
        # Compare Cartesian coordinates for Euclidean distance
        x1 = vec1.cartesian_coordinates.x
        y1 = vec1.cartesian_coordinates.y
        x2 = vec2.cartesian_coordinates.x
        y2 = vec2.cartesian_coordinates.y
        euclidean_dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        # Compare polar coordinates for magnitude and angle differences
        mag1 = vec1.vector_properties.magnitude
        angle1 = vec1.vector_properties.angle_radians
        mag2 = vec2.vector_properties.magnitude
        angle2 = vec2.vector_properties.angle_radians
        
        mag_diff = abs(mag1 - mag2)
        angle_diff = abs(angle1 - angle2)
        # Ensure angle difference is the smaller arc (0 to pi)
        if angle_diff > np.pi:
            angle_diff = 2 * np.pi - angle_diff
        
        # Calculate normalized similarity score (0-1)
        # Avoid division by zero if magnitudes are 0
        max_mag = max(mag1, mag2, 1e-9) # Add small epsilon
        mag_similarity = 1 - (mag_diff / max_mag)
        angle_similarity = 1 - (angle_diff / np.pi) # Normalize angle diff to [0,1]
        similarity_score = (mag_similarity + angle_similarity) / 2
        
        return VectorComparisonResult(
            vector1=vec1,
            vector2=vec2,
            euclidean_distance=round(float(euclidean_dist), 6),
            magnitude_difference=round(mag_diff, 6),
            angle_difference_radians=round(angle_diff, 6),
            angle_difference_degrees=round(math.degrees(angle_diff), 6),
            similarity_score=round(similarity_score, 6)
        )

    def print_comparison(self, vector1_path: str, vector2_path: str):
        """Print a formatted comparison of two vectors.
        
        Args:
            vector1_path: Path to first vector file
            vector2_path: Path to second vector file
        """
        result = self.compare_vectors(vector1_path, vector2_path)
        
        print("\nVector Comparison Results:")
        print("=" * 50)
        print(f"Vector 1: {Path(vector1_path).name} ({result.vector1.sample_id}, {result.vector1.source_type})")
        print(f"  Magnitude: {result.vector1.vector_properties.magnitude:.6f}")
        print(f"  Angle: {result.vector1.vector_properties.angle_degrees:.6f}°")
        print(f"  Cartesian: (x={result.vector1.cartesian_coordinates.x:.6f}, y={result.vector1.cartesian_coordinates.y:.6f})")
        print(f"  Markers used: {result.vector1.markers_used}")
        
        print(f"\nVector 2: {Path(vector2_path).name} ({result.vector2.sample_id}, {result.vector2.source_type})")
        print(f"  Magnitude: {result.vector2.vector_properties.magnitude:.6f}")
        print(f"  Angle: {result.vector2.vector_properties.angle_degrees:.6f}°")
        print(f"  Cartesian: (x={result.vector2.cartesian_coordinates.x:.6f}, y={result.vector2.cartesian_coordinates.y:.6f})")
        print(f"  Markers used: {result.vector2.markers_used}")
        
        print("\nDifferences:")
        print(f"  Euclidean Distance: {result.euclidean_distance:.6f}")
        print(f"  Magnitude Difference: {result.magnitude_difference:.6f}")
        print(f"  Angular Difference: {result.angle_difference_degrees:.6f}°")
        print(f"\nSimilarity Score (0-1): {result.similarity_score:.6f}")
        print("=" * 50) 