import numpy as np

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.md_emulation.evaluate import evaluate_md_emulation
from bioemu_benchmarks.samples import IndexedSamples, find_samples_in_dir

TARGET_METRICS = {
    "mae": 0.941111545976,
    "rmse": 1.209342521795,
    "coverage": 0.438759689922,
}


def test_evaluate_md_emulation(samples_path, samples_projections):
    """Check if evaluate function works as intended."""
    sequence_samples = find_samples_in_dir(samples_path)
    indexed_samples = IndexedSamples.from_benchmark(Benchmark.MD_EMULATION, sequence_samples)

    md_emulation_results = evaluate_md_emulation(indexed_samples, temperature_K=300, random_seed=42)

    # Check basic shapes. Metrics also compute average, so there is one more.
    assert len(md_emulation_results.metrics) == 2
    assert len(md_emulation_results.sample_projections) == 1

    # Check if projections for test system exist and compare against target value.
    assert "cath1_1bl0A02" in md_emulation_results.sample_projections
    np.testing.assert_allclose(
        md_emulation_results.sample_projections["cath1_1bl0A02"], samples_projections
    )

    # Compare metric values against targets.
    metrics_system = md_emulation_results.metrics.loc["cath1_1bl0A02"]
    for metric in TARGET_METRICS:
        np.testing.assert_allclose(metrics_system[metric], TARGET_METRICS[metric])
