from pathlib import Path

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.multiconf.evaluate import (
    MetricType,
    evaluate_multiconf,
)
from bioemu_benchmarks.eval.multiconf.results import MulticonfResults
from bioemu_benchmarks.eval.multiconf.summary_metrics import (
    compute_coverage_and_k_recall,
    compute_coverage_and_k_recall_crypticpocket,
)
from bioemu_benchmarks.eval.multiconf.unfolding_evaluate import (
    evaluate_singleconf_unfolding,
)
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.samples import IndexedSamples
from bioemu_benchmarks.utils import StrPath

LOGGER = get_logger(__name__)


def evaluate_local_unfolding(
    indexed_samples: IndexedSamples, asset_dir: StrPath
) -> MulticonfResults:
    """Evaluate the local unfolding benchmark."""
    evals_per_test_case = evaluate_singleconf_unfolding(
        indexed_samples=indexed_samples,
        references_dir=Path(asset_dir) / "reference",
        references_localresidinfo_dir=Path(asset_dir) / "local_residinfo",
    )
    multiconf_results = compute_coverage_and_k_recall(
        evals_per_test_case, Benchmark.SINGLECONF_LOCALUNFOLDING
    )
    return multiconf_results


def evaluate_multiconf_local(
    indexed_samples: IndexedSamples,
    asset_dir: StrPath,
    benchmark: Benchmark,
    metric_types: list[MetricType] | None = None,
) -> MulticonfResults:
    """Evaluate a local multiconf benchmark. 'local' means that structures are aligned and compared using only subsets of their residues."""
    evals_per_test_case = evaluate_multiconf(
        indexed_samples=indexed_samples,
        references_dir=Path(asset_dir) / "reference",
        metric_types=metric_types,
        references_localresidinfo_dir=Path(asset_dir) / "local_residinfo",
    )
    multiconf_results = compute_coverage_and_k_recall(evals_per_test_case, benchmark=benchmark)
    return multiconf_results


def evaluate_crypticpocket(
    indexed_samples: IndexedSamples,
    asset_dir: StrPath,
    metric_types: list[MetricType] | None = None,
) -> MulticonfResults:
    """Evaluate the cryptic pocket benchmark. It's like a local multiconf benchmark, but the results are reported differently,
    distinguishing between coverage of holo and apo states.
    """
    evals_per_test_case = evaluate_multiconf(
        indexed_samples=indexed_samples,
        references_dir=Path(asset_dir) / "reference",
        metric_types=metric_types,
        references_localresidinfo_dir=Path(asset_dir) / "local_residinfo",
    )
    multiconf_results = compute_coverage_and_k_recall_crypticpocket(
        evals_per_test_case, benchmark=Benchmark.MULTICONF_CRYPTICPOCKET
    )
    return multiconf_results


def evaluate_multiconf_global(
    indexed_samples: IndexedSamples,
    asset_dir: StrPath,
    benchmark: Benchmark,
    metric_types: list[MetricType] | None = None,
) -> MulticonfResults:
    """Evaluate a global multiconf benchmark."""
    evals_per_test_case = evaluate_multiconf(
        indexed_samples=indexed_samples,
        references_dir=Path(asset_dir) / "reference",
        metric_types=metric_types,
    )
    multiconf_results = compute_coverage_and_k_recall(evals_per_test_case, benchmark=benchmark)
    return multiconf_results
