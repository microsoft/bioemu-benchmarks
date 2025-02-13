from dataclasses import dataclass

import mdtraj
import numpy as np


@dataclass
class ProjectionParameters:
    """
    Class for storing projection parameters.

    Attributes:
        sqrt_inv_cov: Projection matrix (square root of inverse covariance matrix). Shape is
          [num_features x num_features], where `num_features` is the number of features in the final
          coordinate transformation.
        mean: Mean of features used to generate projection (shape [num_features]).
    """

    sqrt_inv_cov: np.ndarray
    mean: np.ndarray


class FeatureSettings:
    """
    Setting for generating features from molecular dynamics trajectories.

    Attributes:
        n_trim: Trim N residues from start and end of protein chain.
        exclude_neighbors: Exclude all neighbors within N residues along the sequence from distance
          computation.
        effective_distance: Parameter used for computing contact maps.
    """

    n_trim: int = 2
    exclude_neighbors: int = 2
    effective_distance: float = 0.8


def project_samples(
    samples: dict[str, mdtraj.Trajectory], projection_params: dict[str, ProjectionParameters]
) -> dict[str, np.ndarray]:
    """
    Compute projections for a set of sample trajectories.

    Args:
        samples: Dictionary of trajectories with test case IDs as keys.
        projection_params: Dictionary of parameters used for projections with test case IDs as keys.

    Returns:
        Dictionary of projections with test case IDs as keys.
    """

    sample_projections: dict[str, np.ndarray] = {}
    for test_case in samples:
        # Compute features.
        features = compute_features(samples[test_case], feature_settings=FeatureSettings())
        # Get projection.
        projection = project_features(features, projection_params[test_case])
        sample_projections[test_case] = projection

    return sample_projections


def project_features(
    features: np.ndarray, projection_parameters: ProjectionParameters
) -> np.ndarray:
    """
    Compute feature projection.

    Args:
        features: Contact map features of samples.
        projection_parameters: Parameters for projecting.

    Returns:
        Projected features.
    """
    projection = features - projection_parameters.mean
    projection = projection @ projection_parameters.sqrt_inv_cov
    return projection


def compute_features(
    trajectory: mdtraj.Trajectory, feature_settings: FeatureSettings
) -> np.ndarray:
    """
    Compute contact map features from a trajectory.

    Args:
        trajectory: Trajectory of samples.
        feature_settings: Settings for features.

    Returns:
        Array with computed features.
    """

    # Extract C alpha coordinates and apply trimming.
    ca_coordinates = get_ca_coordinates(trajectory, n_trim=feature_settings.n_trim)
    # Compute C alpha distance matrices and exclude neighbours.
    distance_matrices = compute_distance_matrices(
        ca_coordinates, exclude_neighbors=feature_settings.exclude_neighbors
    )

    # Compute contact map features.
    features = distance_matrices / feature_settings.effective_distance
    features = np.minimum(np.exp(-features), 1.0)

    # Extract upper triangular part (unique contact maps).
    dim_features = features.shape[-1]
    idx_i, idx_j = np.triu_indices(dim_features)
    features = features[:, idx_i, idx_j]

    return features


def compute_distance_matrices(coordinates: np.ndarray, exclude_neighbors: int = 2) -> np.ndarray:
    """
    Compute distance matrices along a set of coordinates with the option to exclude close neighbors.
    Uses units of nm.

    Args:
        coordinates: Coordinates extracted from trajectory. Should be of shape
          [n_samples x n_ca x 3].
        exclude_neighbors: Minimum neighbor distance to consider assuming atoms are ordered
          according to sequence.

    Returns:
        Distance matrix array of the shape [n_samples x n_ca x n_ca].
    """
    distance_matrices = np.linalg.norm(
        np.expand_dims(coordinates, 2) - np.expand_dims(coordinates, 1), axis=-1
    )
    entry_idx = np.arange(coordinates.shape[1])
    neighbor_mask = np.abs(entry_idx[:, None] - entry_idx[None, :]) <= exclude_neighbors
    distance_matrices[:, neighbor_mask] = 0.0
    return distance_matrices


def get_ca_coordinates(trajectory: mdtraj.Trajectory, n_trim: int = 0) -> np.ndarray:
    """
    Extract C alpha coordinates from trajectory.

    Args:
        trajectory: Samples trajectory.
        n_trim: Trim N residues at start and end of protein chain.

    Returns:
        C alpha coordinates in nm with the shape [n_samples x (n_ca - 2*n_trim) x 3].
    """
    atom_indices = trajectory.top.select(
        f"name CA and resid {n_trim} to {trajectory.n_residues-1-n_trim}"
    )
    return trajectory.xyz[:, atom_indices]
