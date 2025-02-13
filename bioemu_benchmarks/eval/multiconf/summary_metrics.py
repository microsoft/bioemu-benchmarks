import operator

import numpy as np

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.multiconf.evaluate import (
    METRIC_OPERATOR_BETTER,
    MetricType,
    TestCaseResult,
    split_holo_apo,
)
from bioemu_benchmarks.eval.multiconf.results import MulticonfResults

METRICS_MIN_VAL: dict[MetricType, float] = {
    MetricType.RMSD: 0.0,
    MetricType.TMSCORE: 0.0,
    MetricType.LDDT: 0.0,
    MetricType.CONTACT_DISTANCE: 0.0,
    MetricType.DSSP_ACC: 0.0,
    MetricType.FNC_UNFOLD_U: 0.0,
    MetricType.FNC_UNFOLD_F: 0.0,
}
METRICS_MAX_VAL: dict[MetricType, float] = {
    MetricType.RMSD: 10.0,  # arbitrary, unbounded
    MetricType.TMSCORE: 1.0,
    MetricType.LDDT: 1.0,
    MetricType.CONTACT_DISTANCE: 10.0,  # arbitrary, unbounded
    MetricType.DSSP_ACC: 1.0,
    MetricType.FNC_UNFOLD_U: 1.0,
    MetricType.FNC_UNFOLD_F: 1.0,
}


def coverage_bootstrap(
    results: dict[str, np.ndarray],
    metric_type: MetricType,
    nsuccess: int = 1,
    nbootstrap: int = 20,
    nsample: int = 1000,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Takes a dictionary of test cases mapped to n_samples x n_references arrays
    of metrics and computes coverage at different thresholds using bootstrapping.

    Args:
        results: Dictionary of test cases mapped to n_samples x n_referencs arrays of metrics
        metric_type: The metric reported
        nsuccess: Number of successes required to count a sample as covered
        nbootstrap: Number of bootstrap resamples
        nsample: Number of samples to bootstrap

    Returns:
        A tuple of the thresholds and the coverages at each threshold, respectively.
    """
    n_samples: list[np.ndarray] = []
    for _ in range(nbootstrap):
        resampled_results: dict[str, np.ndarray] = {}
        for s in results.keys():
            x = results[s]
            indices = np.random.randint(0, x.shape[0], size=nsample)
            resampled_results[s] = x[indices]
        x_range, coverages = coverage(resampled_results, metric_type, nsuccess)
        n_samples.append(coverages)

    return x_range, np.array(n_samples)


def coverage(
    results: dict[str, np.ndarray], metric_type: MetricType, nsuccess: int = 1
) -> tuple[np.ndarray, np.ndarray]:
    """
    Takes a dictionary of test cases mapped to n_samples x n_referencs arrays of metrics and computed coverage at different thresholds

    Args:
        results: Dictionary of test cases mapped to n_samples x n_referencs arrays of metrics
        metric_type: The metric reported
        nsuccess: Number of successes required to count a sample as covered
        bootstrap: Whether to bootstrap resample
        nsample: Number of samples to bootstrap, cannot be None if bootstrap is True

    Returns:
        A tuple of the thresholds and the coverage at each threshold, respectively.

    """
    xrange = (METRICS_MIN_VAL[metric_type], METRICS_MAX_VAL[metric_type])
    xrange = np.linspace(xrange[0], xrange[1], num=100)
    n_below = np.zeros_like(xrange)

    for s in results.keys():
        x = results[s]

        for j, r in enumerate(xrange):
            C = METRIC_OPERATOR_BETTER[metric_type](x, r).astype(int).sum(axis=0)

            n_below[j] += np.mean(C >= nsuccess)

    # normalize
    n_below /= len(results.keys())
    return xrange, n_below


def k_recall(
    results: dict[str, np.ndarray], metric_type: MetricType, k: int = 1
) -> dict[str, float]:
    """
    Computes the k-recall value for a given metric, which we define as the mean of the
    best k values.  'Best' can be smallest or biggest depending on the metric_type.

    Args:
        results: Dictionary of test cases mapped to n_samples x n_references arrays of metrics
        metric_type: The metric reported
        k: Number of samples to use to compute the recall

    Returns:
        dict[str, float]: Dictionary of test cases mapped to the k-recall
    """
    recalls = {}
    for s in results.keys():
        x = np.sort(results[s], axis=0)
        operator_metric = METRIC_OPERATOR_BETTER[metric_type]

        if operator_metric == operator.gt:
            x = x[::-1]

        recall = np.mean(x[:k])
        recalls[s] = recall
    return recalls


def k_recall_bootstrap(
    results: dict[str, np.ndarray],
    metric_type: MetricType,
    k: int = 1,
    nbootstrap: int = 20,
    nsample: int = 1000,
) -> dict[str, tuple[float, float]]:
    """
    Computes the k-recall metric for a given metric with bootstrapping.

    Args:
        results: Dictionary of test cases mapped to n_samples x n_referencs arrays of metrics
        metric_type: The metric reported
        k: Number of samples to use to compute the recall
        sigma: Number of standard deviations to report
        nbootstrap: Number of bootstrap resamples
        nsample: Number of samples to bootstrap

    Returns:
        Dictionary of test cases mapped to a tuple of the mean and standard deviation of the k-recall
    """
    all_recalls: dict[str, list[float]] = {s: [] for s in results.keys()}
    for _ in range(nbootstrap):
        resampled_results = {}
        for s in results.keys():
            x = results[s]
            indices = np.random.randint(0, x.shape[0], size=nsample)
            resampled_results[s] = x[indices]

        recalls = k_recall(resampled_results, metric_type, k)
        for s in recalls.keys():
            all_recalls[s].append(recalls[s])

    return {s: (np.mean(all_recalls[s]), np.std(all_recalls[s])) for s in all_recalls.keys()}


def _list_metrics(results: dict[str, TestCaseResult]) -> list[MetricType]:
    """List the types of metrics in the results.

    Args:
        results: Dictionary of test cases mapped to TestCaseResult objects.
    """
    metrics_per_test_case = [x.metrics_against_references.keys() for x in results.values()]
    assert all(
        set(x) == set(metrics_per_test_case[0]) for x in metrics_per_test_case
    ), "All test cases should have the same metrics"
    return list(metrics_per_test_case[0])


def compute_coverage_and_k_recall(
    evals_per_test_case: dict[str, TestCaseResult],
    benchmark: Benchmark,
) -> MulticonfResults:
    """
    Compute coverage and k-recall for a collection of test cases.

    Args:
        evals_per_test_case: Dictionary of test cases mapped to TestCaseResult objects
        benchmark: The benchmark being evaluated

    Returns:
        MulticonfResults: The compiled results

    """
    metric_types = _list_metrics(evals_per_test_case)

    multiconf_results = MulticonfResults(
        per_system=evals_per_test_case,
        coverage={
            benchmark.value: {
                metric_type: coverage_bootstrap(
                    {
                        test_case: _eval.metrics_against_references[metric_type]
                        for test_case, _eval in evals_per_test_case.items()
                    },
                    metric_type=metric_type,
                )
                for metric_type in metric_types
            }
        },
        krecall={
            benchmark.value: {
                metric_type: k_recall_bootstrap(
                    {
                        test_case: _eval.metrics_against_references[metric_type]
                        for test_case, _eval in evals_per_test_case.items()
                    },
                    metric_type=metric_type,
                )
                for metric_type in metric_types
            }
        },
        benchmark=benchmark,
    )
    return multiconf_results


def compute_coverage_and_k_recall_crypticpocket(
    evals_per_test_case: dict[str, TestCaseResult],
    benchmark: Benchmark,
) -> MulticonfResults:
    """Compute coverage and k-recall for a collection of test cases for the cryptic pocket benchmark."""
    metric_types = _list_metrics(evals_per_test_case)
    holo_evals_per_test_case: dict[str, TestCaseResult] = {}
    apo_evals_per_test_case: dict[str, TestCaseResult] = {}

    for test_case in evals_per_test_case.keys():
        holo_evals_per_test_case[test_case], apo_evals_per_test_case[test_case] = split_holo_apo(
            evals_per_test_case[test_case]
        )

    multiconf_results = MulticonfResults(
        per_system=evals_per_test_case,
        coverage={
            f"{benchmark.value}_holo": {
                metric_type: coverage_bootstrap(
                    {
                        test_case: _eval.metrics_against_references[metric_type]
                        for test_case, _eval in holo_evals_per_test_case.items()
                    },
                    metric_type=metric_type,
                )
                for metric_type in metric_types
            },
            f"{benchmark.value}_apo": {
                metric_type: coverage_bootstrap(
                    {
                        test_case: _eval.metrics_against_references[metric_type]
                        for test_case, _eval in apo_evals_per_test_case.items()
                    },
                    metric_type=metric_type,
                )
                for metric_type in metric_types
            },
        },
        krecall={
            f"{benchmark.value}_holo": {
                metric_type: k_recall_bootstrap(
                    {
                        test_case: _eval.metrics_against_references[metric_type]
                        for test_case, _eval in holo_evals_per_test_case.items()
                    },
                    metric_type=metric_type,
                )
                for metric_type in metric_types
            },
            f"{benchmark.value}_apo": {
                metric_type: k_recall_bootstrap(
                    {
                        test_case: _eval.metrics_against_references[metric_type]
                        for test_case, _eval in apo_evals_per_test_case.items()
                    },
                    metric_type=metric_type,
                )
                for metric_type in metric_types
            },
        },
        benchmark=benchmark,
    )

    return multiconf_results
