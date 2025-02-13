import argparse
import json
from pathlib import Path
from typing import Literal

import pandas as pd

from bioemu_benchmarks.benchmarks import (
    Benchmark,
)
from bioemu_benchmarks.evaluator_utils import evaluator_from_benchmark
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.results import BenchmarkResults
from bioemu_benchmarks.samples import (
    IndexedSamples,
    NoSamples,
    filter_unphysical_samples,
    find_samples_in_dir,
)
from bioemu_benchmarks.utils import StrPath

LOGGER = get_logger(__name__)

BENCHMARK_CHOICES: list[str] = [b.value for b in Benchmark] + ["all"]


def benchmarks_from_choices(
    choices: list[
        Literal[
            "multiconf_ood60",
            "multiconf_oodval",
            "multiconf_domainmotion",
            "multiconf_crypticpocket",
            "singleconf_localunfolding",
            "folding_free_energies",
            "md_emulation",
            "all",
        ]
    ],
) -> list[Benchmark]:
    """
    Map benchmark selection to glorious Enums.

    Args:
        choices: List of benchmark string identifiers. This can contain `all`, which will return
          all available benchmarks without further selections.

    Returns:
        List of selected benchmark types.
    """
    if "all" in choices:
        return [b for b in Benchmark]
    else:
        return [Benchmark(bench) for bench in set(choices)]


def get_sampling_specs(benchmark: Benchmark) -> pd.DataFrame:
    """
    Load test case IDs, sequences, set up number of recommended samples and store in data frame.

    Args:
        benchmark: Benchmark type.

    Returns:
        Dataframe with benchmark test cases, sequences and number of samples for each system.
    """
    benchmark.metadata["num_samples"] = benchmark.default_samplesize
    return benchmark.metadata


def save_sampling_specs(benchmarks: list[Benchmark], csv_file: str) -> None:
    """
    Collect sampling specifications for all systems, merge and write to CSV file.

    Args:
        benchmarks: List of benchmark types.
        csv_file: Output CSV file.
    """
    sampling_specs: list[pd.DataFrame] = []

    for benchmark in benchmarks:
        sampling_specs.append(get_sampling_specs(benchmark))

    merged_specs = pd.concat(sampling_specs)
    merged_specs.to_csv(csv_file, index=False)

    LOGGER.info(f"Sampling specifications written to {csv_file}")


def run_benchmarks(
    benchmarks: list[Benchmark],
    sample_dirs: list[StrPath],
    output_dir: Path,
    overwrite: bool = False,
    filter_samples: bool = True,
) -> None:
    """
    Perform evaluation for all requested benchmarks and store results.

    Args:
        benchmarks: Benchmarks to be evaluated.
        sample_dirs: Directories which should be searched for samples.
        output_dir: Directory to which results are stored. Each benchmark will have its own
          subdirectory.
        overwrite: Overwrite previous results.
        filter_samples: Filter samples by basic criteria (structure clashes, etc).
    """
    from bioemu_benchmarks.utils import BIOEMU_HEADER

    print(BIOEMU_HEADER + "\n")

    LOGGER.info("Running benchmarks:")
    for bench in benchmarks:
        LOGGER.info(f"  - {bench.value}")
    LOGGER.info("")

    aggregate_results: dict[str, dict[str, float]] = {}
    sequence_samples = [x for d in sample_dirs for x in find_samples_in_dir(d)]

    for idx, benchmark in enumerate(benchmarks):
        LOGGER.info(f"Running benchmark {idx+1} of {len(benchmarks)}: {benchmark.value} ...")

        results_dir = output_dir / benchmark.value
        results_dir.mkdir(exist_ok=overwrite, parents=True)

        LOGGER.info("  Loading data...")
        try:
            indexed_samples = IndexedSamples.from_benchmark(
                benchmark, sequence_samples=sequence_samples
            )

            # Filter samples and save statistics.
            if filter_samples:
                indexed_samples, filter_stats = filter_unphysical_samples(indexed_samples)
                filter_stats_mean = {k: v.mean() for k, v in filter_stats.items()}
                with open(results_dir / "filter_statistics.json", "w") as out_file:
                    json.dump(filter_stats_mean, out_file, indent=2, sort_keys=True)
        except NoSamples:
            LOGGER.warning(f"No samples found for {benchmark.value}, skipping benchmark.")
            continue

        LOGGER.info("  Evaluating benchmark ...")
        evaluator = evaluator_from_benchmark(benchmark)
        results: BenchmarkResults = evaluator(indexed_samples)

        LOGGER.info("  Extracting aggregate metrics ...")
        aggregate_results[benchmark.value] = results.get_aggregate_metrics()

        LOGGER.info("  Storing results ...")
        results.save_results(results_dir)
        results.plot(results_dir)
        results.to_pickle(results_dir / "results.pkl")

    # Collect aggregate metrics.
    with open(output_dir / "benchmark_metrics.json", "w") as benchmark_out:
        json.dump(aggregate_results, benchmark_out, indent=2, sort_keys=True)

    LOGGER.info(f"Aggregate results stored to {output_dir / 'benchmark_metrics.json'}")
    LOGGER.info("Finished.")


def parse_arguments() -> argparse.Namespace:
    """Set up parser and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perform BioEmu benchmarks on a set of samples.")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    # Subparser for sample sequence printing.
    parser_print = subparsers.add_parser(
        "specs",
        help="Store sample specifications to CSV.",
        description="Store sample specifications to CSV.",
    )
    parser_print.add_argument("output_csv", type=str, help="Output CSV file.")
    parser_print.add_argument(
        "--benchmarks",
        "-b",
        nargs="+",
        type=str,
        choices=BENCHMARK_CHOICES,
        required=True,
        help="List of benchmarks to evaluate. If set to `all`, "
        "all available benchmarks will be run.",
    )

    # Subparser for running evaluations.
    parser_eval = subparsers.add_parser(
        "eval",
        help="Evaluate benchmarks on a set of samples.",
        description="Evaluate benchmarks on a set of samples.",
    )
    parser_eval.add_argument(
        "output_dir",
        type=Path,
        help="Output directory to which results will be written. "
        "Each benchmark will have its own subdirectory.",
    )
    parser_eval.add_argument(
        "--benchmarks",
        "-b",
        nargs="+",
        type=str,
        choices=BENCHMARK_CHOICES,
        required=True,
        help="List of benchmarks to evaluate. If set to `all`, "
        "all available benchmarks will be run.",
    )
    parser_eval.add_argument(
        "--sample_dirs",
        "-s",
        nargs="+",
        type=Path,
        required=True,
        help="Directory or list of directories containing samples. Each sample "
        "should have a *.pdb and *.xtc file associated with it.",
    )
    parser_eval.add_argument("--overwrite", action="store_true", help="Overwrite previous results.")
    parser_eval.add_argument(
        "--skip_filtering",
        action="store_true",
        help="Do not filter structures for unphysical features like clashes.",
    )

    args = parser.parse_args()

    return args


def cli():
    args = parse_arguments()

    # Set up benchmark enums.
    target_benchmarks = benchmarks_from_choices(choices=args.benchmarks)

    if args.cmd == "specs":
        save_sampling_specs(benchmarks=target_benchmarks, csv_file=args.output_csv)
    elif args.cmd == "eval":
        run_benchmarks(
            benchmarks=target_benchmarks,
            sample_dirs=args.sample_dirs,
            output_dir=args.output_dir,
            overwrite=args.overwrite,
            filter_samples=not args.skip_filtering,
        )
    else:
        raise NotImplementedError(f"Unrecognized command {args.cmd}")


if __name__ == "__main__":
    cli()
