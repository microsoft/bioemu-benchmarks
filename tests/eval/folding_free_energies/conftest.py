from dataclasses import dataclass
from pathlib import Path

import mdtraj
import numpy as np
import pandas as pd
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.folding_free_energies.utils import (
    load_reference,
)

from ... import TEST_DATA_DIR


@dataclass
class FNCTestData:
    fnc: np.ndarray
    sequence: str
    test_case: str
    target_dg: float
    target_ddg: float | None
    p_fold_thr: float
    steepness: float
    temperature: float


@pytest.fixture
def fnc_test_data_wt() -> FNCTestData:
    """FNC data for 1TG0."""
    return FNCTestData(
        fnc=np.array([0.5556425, 0.9077428, 0.07344922, 0.85962605, 0.9426347]),
        test_case="1TG0",
        sequence="EVPFKVVAQFPYKSDYEDDLNFEKDQEIIVTSVEDAEWYFGEYQDSNGDVIEGIFPKSFVAVQG",
        target_dg=-0.645082176012767,
        target_ddg=None,
        p_fold_thr=0.5,
        steepness=10.0,
        temperature=295.0,
    )


@pytest.fixture
def fnc_test_data_mutant() -> FNCTestData:
    """
    FNC data for 1TG0__D45C_N47P. NOTE: the FNC values here are adapted for testing. Original
    values result in dGs almost identical to wild type.
    """
    return FNCTestData(
        fnc=np.array([0.5872849, 0.08238664, 0.6484105, 0.88388467, 0.9725501]),
        sequence="EVPFKVVAQFPYKSDYEDDLNFEKDQEIIVTSVEDAEWYFGEYQCSPGDVIEGIFPKSFVAVQG",
        test_case="1TG0__D45C_N47P",
        target_dg=-0.6771580745606026,
        target_ddg=-0.03207589854783566,
        p_fold_thr=0.5,
        steepness=10.0,
        temperature=295.0,
    )


@pytest.fixture(params=["wt", "mutant"])
def fnc_test_data_parameterized(fnc_test_data_wt, fnc_test_data_mutant, request) -> FNCTestData:
    if request.param == "wt":
        return fnc_test_data_wt
    else:
        return fnc_test_data_mutant


@pytest.fixture
def system_info() -> pd.DataFrame:
    """Convenience fixture for getting system info data frame."""
    return pd.read_csv(Path(Benchmark.FOLDING_FREE_ENERGIES.asset_dir) / "system_info.csv")


@pytest.fixture
def samples_path() -> Path:
    test_data_path = Path(TEST_DATA_DIR) / "samples_example" / "folding_free_energies"
    return test_data_path


@pytest.fixture
def samples_traj_wt(samples_path) -> mdtraj.Trajectory:
    top_path = samples_path / "test_1TG0.pdb"
    xtc_path = samples_path / "test_1TG0.xtc"

    traj = mdtraj.load_xtc(xtc_path, top=top_path)
    return traj


@pytest.fixture
def reference_conformation_wt(fnc_test_data_wt) -> mdtraj.Trajectory:
    """Reference structure of 1TG0 wild type."""
    return load_reference(fnc_test_data_wt.test_case)


@pytest.fixture
def reference_conformation_mutant(fnc_test_data_mutant) -> mdtraj.Trajectory:
    """Reference structure of 1TG0__D45C_N47P mutant."""
    return load_reference(
        fnc_test_data_mutant.test_case,
    )


@pytest.fixture
def reference_contacts_mutant() -> tuple[np.ndarray, np.ndarray]:
    """Reference contacts and distances of 1TG0 wild type."""
    test_contacts_path = (
        Path(TEST_DATA_DIR) / "folding_free_energies" / "test_1TG0__D45C_N47P_contacts.npz"
    )
    test_contact_data = np.load(test_contacts_path)
    pair_idx = test_contact_data["pair_indices"]
    pair_dist = test_contact_data["pair_distances"]
    return pair_idx, pair_dist


@pytest.fixture
def test_free_energy_results() -> pd.DataFrame:
    """Dataframe with reference results for testing."""
    test_free_energy_results_path = (
        Path(TEST_DATA_DIR) / "folding_free_energies" / "test_results.csv"
    )
    return pd.read_csv(test_free_energy_results_path)


@pytest.fixture
def test_metrics_dg() -> dict[str, float]:
    """Target metrics for dG."""
    metrics = {
        "mae": 1.1010791628428436,
        "pearson_corrcoef": -0.9999999999999998,
        "spearman_corrcoef": -0.9999999999999999,
    }
    return metrics


@pytest.fixture
def test_metrics_ddg() -> dict[str, float]:
    """Target metrics for ddG. Correlations are NaN since single test datapoint is present."""
    metrics = {"mae": 2.062185936444976, "pearson_corrcoef": np.nan, "spearman_corrcoef": np.nan}
    return metrics
