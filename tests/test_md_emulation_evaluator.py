from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.md_emulation.evaluate import MDEmulationResults, evaluate_md_emulation
from bioemu_benchmarks.samples import IndexedSamples, find_samples_in_dir
from tests.eval.md_emulation.test_evaluate import TARGET_METRICS

from . import TEST_DATA_DIR

_expected_output_files = [
    "results.pkl",
    "results_projections.npz",
    "results_metrics.csv",
    "metrics.png",
    "projections.png",
]


@pytest.fixture
def samples_path() -> Path:
    test_data_path = Path(TEST_DATA_DIR) / "samples_example" / "md_emulation"
    return test_data_path


def test_evaluate_folding_free_energies(samples_path: Path, tmp_path: Path):
    """Check if evaluator produces expected output files."""
    sequence_samples = find_samples_in_dir(samples_path)
    indexed_samples = IndexedSamples.from_benchmark(
        Benchmark.MD_EMULATION, sequence_samples=sequence_samples
    )

    results = evaluate_md_emulation(
        indexed_samples=indexed_samples, temperature_K=300, random_seed=42
    )

    results.save_results(output_dir=tmp_path)
    results.plot(output_dir=tmp_path)
    results.to_pickle(tmp_path / "results.pkl")

    for expected_file in _expected_output_files:
        assert (tmp_path / expected_file).exists(), f"Missing expected file {expected_file}."

    # Perform check if metrics were computed and stored properly.
    metric_results = pd.read_csv(tmp_path / "results_metrics.csv")
    metric_results = metric_results.set_index("test_case")
    metrics_system = metric_results.loc["cath1_1bl0A02"]
    for metric in TARGET_METRICS:
        np.testing.assert_allclose(metrics_system[metric], TARGET_METRICS[metric])

    # Test if pickle serialization works.
    loaded_results = MDEmulationResults.from_pickle(tmp_path / "results.pkl")
    pd.testing.assert_frame_equal(loaded_results.metrics, results.metrics)
    for test_case in results.sample_projections:
        np.testing.assert_allclose(
            loaded_results.sample_projections[test_case],
            results.sample_projections[test_case],
        )
