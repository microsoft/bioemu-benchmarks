from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.folding_free_energies.evaluate import (
    FoldingFreeEnergiesResults,
    evaluate_folding_free_energies,
)
from bioemu_benchmarks.samples import IndexedSamples, find_samples_in_dir
from tests.eval.folding_free_energies.conftest import TARGET_MAE

from . import TEST_DATA_DIR

_expected_output_files = [
    "results.pkl",
    "contact_scores.npz",
    "results_metrics.csv",
    "results_systems.csv",
    "scatter_ddG.png",
    "scatter_dG.png",
]


@pytest.fixture
def samples_path() -> Path:
    test_data_path = Path(TEST_DATA_DIR) / "samples_example" / "folding_free_energies"
    return test_data_path


def test_evaluate_folding_free_energies(samples_path: Path, tmp_path: Path):
    """Check if evaluator produces expected output files."""
    sequence_samples = find_samples_in_dir(samples_path)
    indexed_samples = IndexedSamples.from_benchmark(
        Benchmark.FOLDING_FREE_ENERGIES, sequence_samples=sequence_samples
    )

    results = evaluate_folding_free_energies(indexed_samples=indexed_samples, temperature_K=295.0)

    results.save_results(output_dir=tmp_path)
    results.plot(output_dir=tmp_path)
    results.to_pickle(tmp_path / "results.pkl")

    for expected_file in _expected_output_files:
        assert (tmp_path / expected_file).exists(), f"Missing expected file {expected_file}."

    # Perform check if metrics were computed and stored properly.
    metric_results = pd.read_csv(tmp_path / "results_metrics.csv")
    np.testing.assert_allclose(metric_results.mae, TARGET_MAE, rtol=1e-6)

    # Test if pickle serialization works.
    loaded_results = FoldingFreeEnergiesResults.from_pickle(tmp_path / "results.pkl")
    pd.testing.assert_frame_equal(loaded_results.metrics, results.metrics)
    pd.testing.assert_frame_equal(
        loaded_results.free_energies_per_system, results.free_energies_per_system
    )
    for test_case in results.fnc_per_system:
        np.testing.assert_allclose(
            loaded_results.fnc_per_system[test_case],
            results.fnc_per_system[test_case],
        )
