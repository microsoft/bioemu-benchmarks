import os

import numpy as np
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.multiconf.align import setup_tm_align
from bioemu_benchmarks.eval.multiconf.evaluate import (
    MULTICONF_METRIC_TYPES,
    MetricType,
    evaluate_multiconf,
    get_closest_sample_per_reference,
)
from bioemu_benchmarks.eval.multiconf.unfolding_evaluate import evaluate_singleconf_unfolding
from bioemu_benchmarks.paths import MULTICONF_ASSET_DIR
from bioemu_benchmarks.samples import IndexedSamples, SequenceSample, find_samples_in_dir
from tests.test_samples import DEFAULT_MULTICONF_BENCHMARK, DEFAULT_MULTICONF_BENCHMARK_SAMPLE_DIR

from ... import TEST_DATA_DIR

# Multiconf metrics
_expected_value_global: dict[MetricType, float] = {
    MetricType.RMSD: 7.442272,
    MetricType.TMSCORE: 0.62649,
    MetricType.LDDT: 0.68465437,
    MetricType.DSSP_ACC: 0.83333333,
    MetricType.CONTACT_DISTANCE: 0.78787879,
}

_expected_value_local: dict[MetricType, float] = {
    MetricType.RMSD: 4.686166,
    MetricType.TMSCORE: 0.30437168,
    MetricType.LDDT: 0.76430888,
    MetricType.DSSP_ACC: 0.85714286,
    MetricType.CONTACT_DISTANCE: 0.6122449,
}


# Singleconf unfolding mean sample metrics
_expected_value_local_unfold_mean: dict[MetricType, float] = {MetricType.FNC_UNFOLD_F: 0.2972972}


@pytest.mark.parametrize(
    "is_local_eval, expected_vals",
    [(False, _expected_value_global), (True, _expected_value_local)],
)
def test_evaluate_multiconf(is_local_eval, expected_vals, test_test_case: str = "Q699R5") -> None:
    """Tests that we can compute multiconf metrics using a test sample and reference directory on blob"""
    setup_tm_align()
    sequence_samples = find_samples_in_dir(DEFAULT_MULTICONF_BENCHMARK_SAMPLE_DIR)
    indexed_samples = IndexedSamples.from_benchmark(
        benchmark=DEFAULT_MULTICONF_BENCHMARK, sequence_samples=sequence_samples
    )
    benchmark_dir = os.path.join(MULTICONF_ASSET_DIR, "ood60")
    references_dir = os.path.join(benchmark_dir, "reference")
    references_localresidinfo_dir = (
        os.path.join(benchmark_dir, "local_residinfo") if is_local_eval else None
    )

    evals_per_test_case = evaluate_multiconf(
        indexed_samples=indexed_samples,
        references_dir=references_dir,
        references_localresidinfo_dir=references_localresidinfo_dir,
        metric_types=MULTICONF_METRIC_TYPES,
    )

    assert test_test_case in evals_per_test_case

    reference_metrics_test_test_case = evals_per_test_case[
        test_test_case
    ].metrics_between_references
    sample_metrics_test_test_case = evals_per_test_case[test_test_case].metrics_against_references

    # Check nothing has changed in how metrics are computed

    for metric_type in expected_vals.keys():
        assert reference_metrics_test_test_case is not None  # shut up mypy
        assert np.isclose(
            reference_metrics_test_test_case[metric_type][(0, 1)],
            expected_vals[metric_type],
            atol=0.01,
        )

        closest_samples = evals_per_test_case[test_test_case].closest_samples
        assert closest_samples is not None  # shut up mypy
        closest_samples_mt = closest_samples[metric_type]
        assert len(closest_samples_mt) == sample_metrics_test_test_case[metric_type].shape[1]

    # Check that collation has happened correctly for samples
    for metric_type in expected_vals.keys():
        rf_mt = evals_per_test_case[test_test_case].references_names
        assert rf_mt is not None  # shut up mypy
        assert sample_metrics_test_test_case[metric_type].shape == (
            5,
            len(rf_mt),
        )
        assert ~np.any(np.isnan(sample_metrics_test_test_case[metric_type]))


def test_evaluate_singleconfunfolding(test_test_case: str = "O88273"):
    """Tests that we can compute singleconf unfolding metrics using test sample and reference directories"""
    benchmark = Benchmark.SINGLECONF_LOCALUNFOLDING
    sequence_samples = find_samples_in_dir(
        os.path.join(TEST_DATA_DIR, "samples_example", "singleconf_localunfolding")
    )

    indexed_samples = IndexedSamples.from_benchmark(
        benchmark=benchmark,
        sequence_samples=sequence_samples,
    )
    benchmark_dir = os.path.join(MULTICONF_ASSET_DIR, "localunfolding")

    references_dir = os.path.join(benchmark_dir, "reference")
    references_localresidinfo_dir = os.path.join(benchmark_dir, "local_residinfo")

    evals_per_test_case = evaluate_singleconf_unfolding(
        indexed_samples=indexed_samples,
        references_dir=references_dir,
        references_localresidinfo_dir=references_localresidinfo_dir,
    )
    test_case_metrics = evals_per_test_case[test_test_case].metrics_against_references

    for metric_type in _expected_value_local_unfold_mean:
        np.testing.assert_allclose(
            test_case_metrics[metric_type].mean(),
            _expected_value_local_unfold_mean[metric_type],
            atol=0.01,
        )


def test_closest_sample(metric_type: MetricType = MetricType.RMSD):
    """
    Tests that we can correctly log the best sample after computing metrics
    against arbitrary references
    """
    sequence_samples = [
        SequenceSample(trajectory_file="1K0M_1.63788d73.xtc", topology_file="1K0M_1.63788d73.pdb"),
        SequenceSample(trajectory_file="1RK4_1.9bf5961b.xtc", topology_file="1RK4_1.9bf5961b.pdb"),
    ]

    references = ["foo/ref1.pdb", "bar/ref2.pdb"]

    metrics_against_references = {
        metric_type: np.array([[1.0, 2.0], [2.0, 2.0], [3.0, 1.5], [2.0, 1.5], [2.0, 0.75]])
    }

    topology_ids = np.array([0, 0, 0, 1, 1], dtype=int)

    closest_samples = get_closest_sample_per_reference(
        sequence_samples,
        metrics_against_references=metrics_against_references,
        references=references,
        topology_ids=topology_ids,
    )
    best_ss, best_indices, best_vals, closest_refs = (
        [cs.sequence_sample for cs in closest_samples[metric_type]],
        [cs.frame_idx for cs in closest_samples[metric_type]],
        [cs.metric_value for cs in closest_samples[metric_type]],
        [cs.reference_path for cs in closest_samples[metric_type]],
    )
    assert best_ss == sequence_samples
    assert best_indices == [0, 1]
    assert best_vals == [1.0, 0.75]
    assert closest_refs == references
