from collections.abc import Callable
from functools import partial

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.folding_free_energies.evaluate import evaluate_folding_free_energies
from bioemu_benchmarks.eval.md_emulation.evaluate import evaluate_md_emulation
from bioemu_benchmarks.multiconf_evaluator import (
    evaluate_crypticpocket,
    evaluate_local_unfolding,
    evaluate_multiconf_global,
    evaluate_multiconf_local,
)
from bioemu_benchmarks.results import BenchmarkResults
from bioemu_benchmarks.samples import IndexedSamples

Evaluator = Callable[[IndexedSamples], BenchmarkResults]


def evaluator_from_benchmark(benchmark: Benchmark, **kwargs) -> Evaluator:
    """
    Utility function for getting evaluators depending on the benchmark type.

    Args:
        benchmark: Benchmark type.

    Returns:
        Evaluator for benchmark.
    """

    if benchmark == Benchmark.SINGLECONF_LOCALUNFOLDING:
        return partial(evaluate_local_unfolding, asset_dir=benchmark.asset_dir, **kwargs)
    elif benchmark == Benchmark.MULTICONF_OOD60:
        return partial(
            evaluate_multiconf_local, asset_dir=benchmark.asset_dir, benchmark=benchmark, **kwargs
        )
    elif benchmark == Benchmark.MULTICONF_CRYPTICPOCKET:
        return partial(evaluate_crypticpocket, asset_dir=benchmark.asset_dir, **kwargs)
    elif benchmark in [Benchmark.MULTICONF_OODVAL, Benchmark.MULTICONF_DOMAINMOTION]:
        return partial(
            evaluate_multiconf_global, asset_dir=benchmark.asset_dir, benchmark=benchmark, **kwargs
        )
    elif benchmark == Benchmark.FOLDING_FREE_ENERGIES:
        return partial(evaluate_folding_free_energies, **kwargs)
    elif benchmark == Benchmark.MD_EMULATION:
        return partial(evaluate_md_emulation, **kwargs)
    else:
        raise ValueError(f"Unrecognized benchmark {benchmark}")
