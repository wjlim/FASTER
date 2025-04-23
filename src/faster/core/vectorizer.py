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
    def __init__(self):
        """Initialize the GenotypeVectorizer class."""
        self.markers = []
        self.vector_size = 2  # Each marker has 2 alleles

    def _parse_genotype(self, genotype: str) -> Tuple[float, float]:
        """Parse genotype string into two allele numbers, always sorted in ascending order.
        Skip if allele contains non-numeric characters (e.g. 'X').
        
        Args:
            genotype: Genotype string in format "allele1/allele2" (e.g. "10/8") or single "allele"
            
        Returns:
            Tuple of (allele1, allele2) as floats, sorted in ascending order.
            For example, "10/8" will return (8.0, 10.0)
            Returns None if any allele contains non-numeric characters
        """
        try:
            if '/' in genotype:
                alleles = genotype.split('/')
                allele1 = float(alleles[0].split('.')[0])  # Remove decimal part
                allele2 = float(alleles[1].split('.')[0])  # Remove decimal part
            else:
                # Single allele case - homozygous
                allele1 = allele2 = float(genotype.split('.')[0])
                
            # Always return alleles sorted in ascending order (smaller number first)
            return tuple(sorted([allele1, allele2]))
        except ValueError:
            # Return None if conversion to float fails (non-numeric allele)
            return None

    def _combine_exhunter_genotypes(self, variants: Dict) -> Tuple[float, float]:
        """Combine multiple ExpansionHunter variant genotypes into a single genotype.
        
        Args:
            variants: Dictionary of variants from ExpansionHunter
            
        Returns:
            Combined genotype as (allele1, allele2), sorted in ascending order
        """
        total_allele1 = 0.0
        total_allele2 = 0.0
        
        for variant in variants.values():
            if 'Genotype' not in variant:
                continue
                
            # Get sorted alleles (smaller first)
            allele1, allele2 = self._parse_genotype(variant['Genotype'])
            total_allele1 += allele1
            total_allele2 += allele2
            
        # Return sorted alleles
        return tuple(sorted([total_allele1, total_allele2]))

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

    def _calculate_vector_properties(self, vector: np.ndarray) -> VectorProperties:
        """Calculate magnitude and angle of the 44-dimensional vector.
        
        Args:
            vector: 44-dimensional numpy array (22 markers x 2 alleles)
            
        Returns:
            VectorProperties object containing magnitude and angles
        """
        # Reshape into 22x2 matrix
        points = vector.reshape(-1, 2)
        
        # Calculate magnitude (Euclidean distance from origin)
        magnitude = float(np.sqrt(np.sum(points ** 2)))
        
        # Calculate angles for each point (in radians)
        angles = np.arctan2(points[:, 1], points[:, 0])
        # Ensure angles are positive (0 to 2π)
        angles = np.where(angles < 0, angles + 2*np.pi, angles)
        # Take average angle
        mean_angle = float(np.mean(angles))
        
        return VectorProperties(
            magnitude=round(magnitude, 6),
            angle_radians=round(mean_angle, 6),
            angle_degrees=round(math.degrees(mean_angle), 6)
        )

    def vectorize_str(self, data: Dict, sample_id: str) -> VectorizationResult:
        """Vectorize STR analysis results.
        
        Args:
            data: Dictionary containing STR analysis results
            sample_id: Sample identifier
            
        Returns:
            VectorizationResult object
        """
        marker_genotypes = []
        vectors = []
        
        for marker, info in data["LocusResults"].items():
            variant = list(info["variants"].values())[0]
            genotype = variant["genotype"]
            result = self._parse_genotype(genotype)
            
            # Skip markers with non-numeric alleles
            if result is not None:
                allele1, allele2 = result
                marker_genotypes.append(
                    MarkerGenotype(
                        marker=marker,
                        allele1=allele1,
                        allele2=allele2
                    )
                )
                vectors.extend([allele1, allele2])
            else:
                logger.info(f"Skipping marker {marker} due to non-numeric allele in genotype: {genotype}")
        
        vector = np.array(vectors)
        vector_props = self._calculate_vector_properties(vector)
        
        # Calculate Cartesian coordinates
        cartesian = self._polar_to_cartesian(
            vector_props.magnitude,
            vector_props.angle_radians
        )
        
        return VectorizationResult(
            sample_id=sample_id,
            source_type="str",
            markers=marker_genotypes,
            vector_properties=vector_props,
            markers_used=len(marker_genotypes),
            compact_vector=f"{vector_props.magnitude:.6f},{vector_props.angle_radians:.6f}",
            cartesian_coordinates=cartesian
        )

    def vectorize_eh(self, data: Dict, sample_id: str) -> VectorizationResult:
        """Vectorize ExpansionHunter results.
        
        Args:
            data: Dictionary containing ExpansionHunter results
            sample_id: Sample identifier
            
        Returns:
            VectorizationResult object
        """
        marker_genotypes = []
        vectors = []
        
        for marker, info in data["LocusResults"].items():
            if "Variants" not in info:
                continue
                
            try:
                allele1, allele2 = self._combine_exhunter_genotypes(info["Variants"])
                marker_genotypes.append(
                    MarkerGenotype(
                        marker=marker,
                        allele1=allele1,
                        allele2=allele2
                    )
                )
                vectors.extend([allele1, allele2])
            except (ValueError, TypeError):
                logger.info(f"Skipping marker {marker} due to non-numeric allele")
                continue
        
        vector = np.array(vectors)
        vector_props = self._calculate_vector_properties(vector)
        
        # Calculate Cartesian coordinates
        cartesian = self._polar_to_cartesian(
            vector_props.magnitude,
            vector_props.angle_radians
        )
        
        return VectorizationResult(
            sample_id=sample_id,
            source_type="eh",
            markers=marker_genotypes,
            vector_properties=vector_props,
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
        # Convert to JSON and save
        with open(output_path, 'w') as f:
            json.dump(result.dict(), f, indent=2)

    def load_vector(self, file_path: str) -> VectorizationResult:
        """Load vectorization results from a JSON file.
        
        Args:
            file_path: Path to the vector file
            
        Returns:
            VectorizationResult object
        """
        with open(file_path) as f:
            data = json.load(f)
        return VectorizationResult.parse_obj(data)

    def compare_vectors(self, vector1_path: str, vector2_path: str) -> VectorComparisonResult:
        """Compare two vectors in polar coordinates.
        
        Args:
            vector1_path: Path to first vector file
            vector2_path: Path to second vector file
            
        Returns:
            VectorComparisonResult object
        """
        # Load vectors
        vec1 = self.load_vector(vector1_path)
        vec2 = self.load_vector(vector2_path)
        
        # Get properties
        mag1 = vec1.vector_properties.magnitude
        angle1 = vec1.vector_properties.angle_radians
        mag2 = vec2.vector_properties.magnitude
        angle2 = vec2.vector_properties.angle_radians
        
        # Calculate differences
        mag_diff = abs(mag1 - mag2)
        angle_diff = abs(angle1 - angle2)
        # Ensure angle difference is the smaller arc
        if angle_diff > np.pi:
            angle_diff = 2 * np.pi - angle_diff
            
        # Calculate Euclidean distance in polar coordinates
        x1 = mag1 * np.cos(angle1)
        y1 = mag1 * np.sin(angle1)
        x2 = mag2 * np.cos(angle2)
        y2 = mag2 * np.sin(angle2)
        euclidean_dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        # Calculate normalized similarity score (0-1)
        max_mag = max(mag1, mag2)
        mag_similarity = 1 - (mag_diff / max_mag)
        angle_similarity = 1 - (angle_diff / np.pi)
        similarity_score = (mag_similarity + angle_similarity) / 2
        
        return VectorComparisonResult(
            vector1=vec1,
            vector2=vec2,
            euclidean_distance=round(euclidean_dist, 6),
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