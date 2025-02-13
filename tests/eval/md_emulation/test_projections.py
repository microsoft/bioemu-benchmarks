import numpy as np

from bioemu_benchmarks.eval.md_emulation.projections import (
    FeatureSettings,
    compute_distance_matrices,
    compute_features,
    get_ca_coordinates,
    project_features,
    project_samples,
)


def test_get_ca_coordinates(samples_traj, samples_ca_coordinates):
    """Test C alpha coordinate extraction and trimming."""
    ca_coordinates = get_ca_coordinates(samples_traj, n_trim=2)

    # Check shape is consistent with selection.
    num_samples = len(samples_traj)
    num_ca = len(samples_traj.top.select("name CA")) - 4  # remove trim on both ends.
    assert ca_coordinates.shape == (num_samples, num_ca, 3)

    # Check if correct coordinate values were extracted.
    np.testing.assert_allclose(ca_coordinates, samples_ca_coordinates)


def test_compute_distance_matrices(samples_ca_coordinates, samples_distance_matrices):
    """Test distance matrix computation."""
    distance_matrices = compute_distance_matrices(samples_ca_coordinates, exclude_neighbors=2)

    # Check shape.
    num_samples, num_atoms, _ = samples_ca_coordinates.shape
    assert distance_matrices.shape == (num_samples, num_atoms, num_atoms)

    # Check if first and second off diagonal (exclude neighbors) are 0. Assumes sequence is ordered.
    for k in [1, 2]:
        off_diag = np.diagonal(distance_matrices, offset=k, axis1=1, axis2=2)
        assert not np.any(off_diag)

    # Check if values are correct.
    np.testing.assert_allclose(distance_matrices, samples_distance_matrices)


def test_compute_features(samples_traj, samples_features):
    """Test feature computation."""
    features = compute_features(samples_traj, FeatureSettings())

    # Check if shape matches (features do not exclude diagonal).
    num_samples, num_atoms, _ = get_ca_coordinates(
        samples_traj, n_trim=FeatureSettings().n_trim
    ).shape
    assert features.shape == (num_samples, int(num_atoms * (num_atoms + 1) / 2))

    # Check against reference values.
    np.testing.assert_allclose(features, samples_features)


def test_project_features(samples_features, samples_projections, projection_params):
    """Test basic feature projection function."""
    projections = project_features(samples_features, projection_parameters=projection_params)

    # Check shape.
    assert projections.shape == (samples_features.shape[0], projection_params.sqrt_inv_cov.shape[1])

    # Check values.
    np.testing.assert_allclose(projections, samples_projections)


def test_project_samples(samples_traj, samples_projections, projection_params):
    """Test trajectory transformation."""
    # Build dummy dictionaries for testing.
    projections = project_samples(
        samples=dict(x=samples_traj),
        projection_params=dict(x=projection_params),
    )

    # Check output and shapes.
    assert len(projections) == 1 and "x" in projections
    target_projections = projections["x"]
    assert target_projections.shape == (
        samples_projections.shape[0],
        projection_params.sqrt_inv_cov.shape[1],
    )

    # Check values.
    np.testing.assert_allclose(target_projections, samples_projections)
