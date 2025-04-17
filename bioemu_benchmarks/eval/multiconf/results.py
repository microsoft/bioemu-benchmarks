import os
from dataclasses import dataclass

import h5py
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from bioemu_benchmarks.benchmarks import MULTICONF_BENCHMARKS, Benchmark
from bioemu_benchmarks.eval.multiconf.evaluate import MetricType, TestCaseResult
from bioemu_benchmarks.eval.multiconf.plot import (
    METRICS_SUCCESS_THRESHOLD,
    METRICS_UNIT,
    plot_2D_free_energy_landscapes_in_grid,
    plot_coverage_bootstrap,
    plot_free_energy_landscapes_by_fnc_in_grid,
)
from bioemu_benchmarks.utils import StrPath

BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS: dict[Benchmark, dict[MetricType, float]] = {
    Benchmark.MULTICONF_CRYPTICPOCKET: {MetricType.RMSD: 1.5}
}
from bioemu_benchmarks.results import BenchmarkResults


@dataclass
class MulticonfResults(BenchmarkResults):
    """Class that contain results of a multiconf benchmark evaluation

    per_system: multiconf evaluation results by test case
    coverage: dictionary of provided labels to overall coverage at different thresholds
    krecall: dictionary of provided labels to krecall results indexed by test case and metricy type
    """

    per_system: dict[str, TestCaseResult]
    coverage: dict[
        str, dict[MetricType, tuple[np.ndarray, np.ndarray]]
    ]  # map from a label to coverage results
    krecall: dict[
        str, dict[MetricType, dict[str, tuple[float, float]]]
    ]  # map from a label to krecall results

    def save_to_h5(self, outfile: StrPath) -> None:
        """Saves main metrics as a hierarchical H5 file.
        Note: H5 does not currently support deserialization
        to `MulticonfResults`.

        Args:
            outfile: Path to target H5 file
        """
        h5 = h5py.File(outfile, mode="w", libver="latest")

        # Dump coverage
        for label, coverage_results in self.coverage.items():
            for metric_type, (thresholds, values) in coverage_results.items():
                mgr = h5.create_group(name=f"coverage_{label}_{metric_type.value}")
                mgr.create_dataset(name="thresholds", data=thresholds)
                mgr.create_dataset(name="values", data=values)

        # Dump krecall
        for label, krecall_results in self.krecall.items():
            for metric_type, test_case_dict in krecall_results.items():
                mgr = h5.create_group(name=f"krecall_{label}_{metric_type.value}")
                for test_case, (km, kstd) in test_case_dict.items():
                    mgr.create_dataset(name=test_case, data=np.array([km, kstd]))

        # Dump system metrics
        sgr = h5.create_group(name="sample_metrics")
        if self.benchmark in MULTICONF_BENCHMARKS:
            rgr = h5.create_group(name="reference_metrics")

        for test_case, meval in self.per_system.items():
            sgr_testcase = sgr.create_group(name=test_case)
            if self.benchmark in MULTICONF_BENCHMARKS:
                rgr_testcase = rgr.create_group(name=test_case)
                sgr_testcase.create_dataset(name="topology_ids", data=meval.topology_ids)

            # Sample metrics
            for metric_type, sample_arr in meval.metrics_against_references.items():
                sgr_testcase.create_dataset(name=metric_type.value, data=sample_arr)

            if meval.references_names is not None:
                # Reference metrics
                rgr_testcase.create_dataset(name="references_names", data=meval.references_names)
                if meval.metrics_between_references is not None:
                    for metric_type, mbr in meval.metrics_between_references.items():
                        if mbr is not None:
                            for (ref_i, ref_j), value in mbr.items():
                                rgr_testcase.create_dataset(
                                    name=f"{metric_type.value}_{ref_i}_{ref_j}",
                                    data=value,
                                )

    def plot(self, output_dir: StrPath) -> None:
        """
        Main entrypoint to generate plots

        Args:
            output_dir: Directory to save plot results

        """

        per_system_evals = self.per_system

        benchmark_values = (
            [f"{self.benchmark.value}_holo", f"{self.benchmark.value}_apo"]
            if self.benchmark == Benchmark.MULTICONF_CRYPTICPOCKET
            else [self.benchmark.value]
        )

        for benchmark_value in benchmark_values:
            for metric_type in tqdm(
                self.coverage[benchmark_value].keys(), desc="Plotting results..."
            ):
                # if a benchmark specific success threshold is defined, use it

                success_threshold = None
                if self.benchmark in BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS:
                    if metric_type in BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS[self.benchmark]:
                        success_threshold = BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS[self.benchmark][
                            metric_type
                        ]

                # plot coverage

                if not os.path.isdir(output_dir):
                    os.makedirs(output_dir)

                for label, coverage_results in self.coverage.items():
                    thresholds, coverages = coverage_results[metric_type]
                    fig = plt.figure(figsize=(3, 3))
                    ax = plt.gca()
                    plot_coverage_bootstrap(
                        thresholds=thresholds,
                        coverages=coverages,
                        ax=ax,
                        metric_type=metric_type,
                        success_threshold=success_threshold,
                    )

                    plt.xlabel(metric_type.value.upper() + " " + METRICS_UNIT[metric_type])
                    plt.ylabel("0.1% Coverage")

                    fig.savefig(
                        os.path.join(output_dir, f"{label}_{metric_type.value}_coverage.png"),
                        dpi=300,
                        bbox_inches="tight",
                    )

                # extract metrics per sample

                sample_metrics = {
                    test_case: _eval.metrics_against_references[metric_type]
                    for test_case, _eval in per_system_evals.items()
                }

                # plot free energies

                if metric_type not in [MetricType.FNC_UNFOLD_F, MetricType.FNC_UNFOLD_U]:
                    fig = plot_2D_free_energy_landscapes_in_grid(
                        sample_metrics, metric_type, success_threshold=success_threshold
                    )
                    fig.savefig(
                        os.path.join(
                            output_dir,
                            f"{benchmark_value}_{metric_type.value}_free_energy.png",
                        ),
                        dpi=300,
                        bbox_inches="tight",
                    )
                else:
                    fig = plot_free_energy_landscapes_by_fnc_in_grid(sample_metrics)
                    fig.savefig(
                        os.path.join(
                            output_dir,
                            f"{benchmark_value}_{metric_type.value}_free_energy.png",
                        ),
                        dpi=300,
                        bbox_inches="tight",
                    )

    def save_closest_samples(self, output_dir: StrPath) -> None:
        """
        Saves closest samples to references in PDB format under `results_dir`

        Args:
            output_dir: Directory to save closest sample results.
                         PDB files will be stored under the `closest` subdirectory
                         of `results_dir`
        """
        assert (
            not self.benchmark == Benchmark.SINGLECONF_LOCALUNFOLDING
        ), "Save close samples not supported for SingleconfBenchmark"
        closest_dir = os.path.join(output_dir, "closest")
        os.makedirs(closest_dir, exist_ok=True)

        for test_case, meval in self.per_system.items():
            if meval.closest_samples is not None:
                for _metric_type, closest_samples_l in meval.closest_samples.items():
                    for close_sample in closest_samples_l:
                        close_sample.save_to_pdb(test_case=test_case, closest_dir=closest_dir)

    def save_results(self, output_dir: StrPath) -> None:
        """
        Utility function for saving metrics as HDF5 file and closest samples in PDB format at the
        same time.

        Args:
            output_dir: Directory to which result outputs should be saved.
        """
        self.save_to_h5(os.path.join(output_dir, "results.h5"))

        if self.benchmark != Benchmark.SINGLECONF_LOCALUNFOLDING:
            self.save_closest_samples(output_dir)

    def get_aggregate_metrics(self) -> dict[str, float]:
        """
        Compute aggregate metrics for multiconf benchmarks. These are the averages of k-recall
        values and their standard deviation, as well as the area under the coverage curve and the
        coverage at the benchmark / metric specific threshold.

        Returns:
            Dictionary of aggregate metrics.
        """
        aggregate_metrics: dict[str, float] = {}

        benchmark_values = (
            [f"{self.benchmark.value}_holo", f"{self.benchmark.value}_apo"]
            if self.benchmark == Benchmark.MULTICONF_CRYPTICPOCKET
            else [self.benchmark.value]
        )

        for benchmark_value in benchmark_values:
            # De-nest nested dicts.
            krecall_metrics = self.krecall[benchmark_value]
            coverage_metrics = self.coverage[benchmark_value]

            # Select suffix for plotting in crypticpocket benchmark.
            if self.benchmark == Benchmark.MULTICONF_CRYPTICPOCKET:
                suffix = f"_{benchmark_value.split("_")[-1]}"
            else:
                suffix = ""

            for metric in krecall_metrics:
                assert metric in coverage_metrics

                # De-nest at metric level.
                krecall = krecall_metrics[metric]
                coverage = coverage_metrics[metric]

                # Get coverage for metric.
                if (
                    self.benchmark in BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS
                    and metric in BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS[self.benchmark]
                ):
                    threshold = BENCHMARK_SPECIFIC_SUCCESS_THRESHOLDS[self.benchmark][metric]
                else:
                    threshold = METRICS_SUCCESS_THRESHOLD[metric]

                # Extract coordinates and compute mean over system axis.
                coverage_coordinates, coverage_per_system = coverage
                coverage_mean = np.mean(coverage_per_system, axis=0)

                # Get area under the curve.
                coverage_auc = np.trapz(coverage_mean, coverage_coordinates)
                # Get value at metric / benchmark threshold.
                coverage_threshold = np.interp(threshold, coverage_coordinates, coverage_mean)

                aggregate_metrics.update(
                    {
                        f"krecall_mean_{metric.value}{suffix}": float(
                            np.mean([x[0] for x in krecall.values()])
                        ),
                        f"krecall_std_{metric.value}{suffix}": float(
                            np.mean([x[1] for x in krecall.values()])
                        ),
                        f"coverage_auc_{metric.value}{suffix}": float(coverage_auc),
                        f"coverage_threshold_{metric.value}{suffix}": float(coverage_threshold),
                    }
                )

        return aggregate_metrics
