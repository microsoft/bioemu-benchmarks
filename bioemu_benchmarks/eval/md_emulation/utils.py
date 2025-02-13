from pathlib import Path

import numpy as np

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.md_emulation.projections import ProjectionParameters
from bioemu_benchmarks.paths import MD_EMULATION_ASSET_DIR


def load_reference_projections() -> dict[str, np.ndarray]:
    """
    Load reference projections from file.

    Args:

    Returns:
        Dictionary with test case identifiers as keys and the target projections arrays as entries.
    """
    benchmark = Benchmark.MD_EMULATION
    projection_path = (
        Path(MD_EMULATION_ASSET_DIR) / benchmark.value.lower() / "reference_projections.npz"
    )
    projections = np.load(projection_path)
    return dict(projections)


def load_projection_parameters() -> dict[str, ProjectionParameters]:
    """
    Load parameters used for projecting samples.

    Args:
        benchmark: MD emulation benchmark.

    Returns:
        Dictionary of projection parameters using test case identifiers as keys.
    """
    benchmark = Benchmark.MD_EMULATION
    parameter_dir = Path(MD_EMULATION_ASSET_DIR) / benchmark.value.lower()
    projection_sqrt_inv_cov = dict(np.load(parameter_dir / "projections_sqrt_inv_cov.npz"))
    projection_mean = dict(np.load(parameter_dir / "projections_mean.npz"))

    assert set(projection_mean) == set(
        projection_sqrt_inv_cov
    ), "Mismatch in systems found for projection parameters."

    projection_params: dict[str, ProjectionParameters] = {}
    for test_case in projection_sqrt_inv_cov.keys():
        projection_params[test_case] = ProjectionParameters(
            sqrt_inv_cov=projection_sqrt_inv_cov[test_case],
            mean=projection_mean[test_case],
        )
    return projection_params
