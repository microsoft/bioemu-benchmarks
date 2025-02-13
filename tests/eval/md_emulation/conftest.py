from dataclasses import dataclass
from pathlib import Path

import mdtraj
import numpy as np
import pytest

from bioemu_benchmarks.eval.folding_free_energies.free_energies import K_BOLTZMANN
from bioemu_benchmarks.eval.md_emulation.projections import ProjectionParameters
from bioemu_benchmarks.eval.md_emulation.state_metric import compute_density_2D, histogram_bin_edges
from bioemu_benchmarks.eval.md_emulation.utils import (
    load_projection_parameters,
    load_reference_projections,
)

from ... import TEST_DATA_DIR


@dataclass
class FreeEnergyErrors:
    energies_predicted: np.ndarray
    energies_target: np.ndarray
    mae_base: float
    rmse_base: float
    mae_minimized: float
    rmse_minimized: float


@dataclass
class DensityMetricTargets:
    score_mae: float
    score_nonzero_mae: float
    score_rmse: float
    score_nonzero_rmse: float
    score_nonzero_coverage: float


@pytest.fixture
def samples_path() -> Path:
    test_data_path = Path(TEST_DATA_DIR) / "samples_example" / "md_emulation"
    return test_data_path


@pytest.fixture
def samples_traj(samples_path) -> mdtraj.Trajectory:
    top_path = samples_path / "cath1_1bl0A02.pdb"
    xtc_path = samples_path / "cath1_1bl0A02.xtc"

    traj = mdtraj.load_xtc(xtc_path, top=top_path)
    return traj


@pytest.fixture
def samples_ca_coordinates() -> np.ndarray:
    test_data_path = Path(TEST_DATA_DIR) / "md_emulation" / "test_cath1_1bl0A02_ca_coordinates.npy"
    return np.load(test_data_path)


@pytest.fixture
def samples_distance_matrices() -> np.ndarray:
    test_data_path = Path(TEST_DATA_DIR) / "md_emulation" / "test_cath1_1bl0A02_distmat.npy"
    return np.load(test_data_path)


@pytest.fixture
def samples_features() -> np.ndarray:
    test_data_path = Path(TEST_DATA_DIR) / "md_emulation" / "test_cath1_1bl0A02_features.npy"
    return np.load(test_data_path)


@pytest.fixture
def samples_projections() -> np.ndarray:
    test_data_path = Path(TEST_DATA_DIR) / "md_emulation" / "test_cath1_1bl0A02_projections.npy"
    return np.load(test_data_path)


@pytest.fixture
def projection_params() -> ProjectionParameters:
    all_projections = load_projection_parameters()
    return all_projections["cath1_1bl0A02"]


@pytest.fixture
def reference_projections() -> np.ndarray:
    all_reference = load_reference_projections()
    return all_reference["cath1_1bl0A02"][::1000]


@pytest.fixture
def free_energy_targets(reference_projections, samples_projections) -> FreeEnergyErrors:
    bins_x = histogram_bin_edges(reference_projections[:, 0], 5, 0.1)
    bins_y = histogram_bin_edges(reference_projections[:, 1], 5, 0.1)

    density_target = compute_density_2D(reference_projections, bins_x, bins_y)
    density_samples = compute_density_2D(samples_projections, bins_x, bins_y)

    kBT = 300.0 * K_BOLTZMANN

    def _clamp_density(density, cutoff=4.0):
        energy_min = -kBT * np.log(np.max(density))
        p_cut = np.exp(-(energy_min + cutoff) / kBT)
        return np.maximum(density, p_cut)

    density_target_c = _clamp_density(density_target)
    density_samples_c = _clamp_density(density_samples)

    energy_target = -kBT * np.log(density_target_c).flatten()
    energy_samples = -kBT * np.log(density_samples_c).flatten()

    return FreeEnergyErrors(
        energies_predicted=energy_samples,
        energies_target=energy_target,
        mae_base=0.3860775974382061,
        rmse_base=0.8687197895923561,
        mae_minimized=0.3121663464279002,
        rmse_minimized=0.8594285807027948,
    )


@pytest.fixture
def density_metric_targets() -> DensityMetricTargets:
    return DensityMetricTargets(
        score_mae=1.0811314147747684,
        score_rmse=1.3776751970232277,
        score_nonzero_mae=0.6341615432056811,
        score_nonzero_rmse=0.711393908412758,
        score_nonzero_coverage=0.8,
    )
