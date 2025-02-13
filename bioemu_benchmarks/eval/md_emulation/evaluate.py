from dataclasses import dataclass
from pathlib import Path

import mdtraj
import numpy as np
import pandas as pd

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.md_emulation.plot import plot_metrics, plot_projections
from bioemu_benchmarks.eval.md_emulation.projections import project_samples
from bioemu_benchmarks.eval.md_emulation.state_metric import (
    DistributionMetricSettings,
    compute_state_metrics,
)
from bioemu_benchmarks.eval.md_emulation.utils import (
    load_projection_parameters,
    load_reference_projections,
)
from bioemu_benchmarks.results import BenchmarkResults
from bioemu_benchmarks.samples import IndexedSamples
from bioemu_benchmarks.utils import StrPath


@dataclass
class MDEmulationResults(BenchmarkResults):
    """
    Data class for collecting MD emulation benchmark results.

    Attributes:
        sample_projections: Dictionary containing projections computed for samples. Keys are the
          test case IDs.
        metrics: Pandas data frame collecting aggregate metrics.
        temperature_K: Temperature used for analysis in units of Kelvin.
        random_seed: Random seed used for analysis (mainly for resampling in computing free energy
          metric).
    """

    sample_projections: dict[str, np.ndarray]
    metrics: pd.DataFrame
    temperature_K: float
    random_seed: int

    def save_results(self, output_dir: StrPath) -> None:
        """
        Save individual evaluator results in accessible files (txt, csv, npz).

        Args:
            output_dir: Directory to which result outputs should be saved.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.metrics.to_csv(output_dir / "results_metrics.csv")
        np.savez(output_dir / "results_projections.npz", **self.sample_projections)

    def plot(self, output_dir: StrPath, max_energy: float = 7.0) -> None:
        """
        Generate plots associated with benchmark and write to output directory.

        Args:
            output_dir: Directory where plots will be written.
            max_energy: Upper energy cutoff for 2D plots in kcal/mol.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        reference_projections = load_reference_projections()

        fig_projections = plot_projections(
            self.sample_projections,
            reference_projections,
            temperature_K=self.temperature_K,
            max_energy=max_energy,
        )
        fig_metrics = plot_metrics(
            self.metrics,
            label_map={
                "mae": r"MAE (kcal/mol) $\downarrow$",
                "rmse": r"RMSE (kcal/mol) $\downarrow$",
                "coverage": r"coverage$\uparrow$",
            },
        )

        fig_projections.savefig(output_dir / "projections.png")
        fig_metrics.savefig(output_dir / "metrics.png")

    def get_aggregate_metrics(self) -> dict[str, float]:
        """
        Collect aggregate mean absolute error, root mean squared error and coverage metrics.

        Returns:
            Dictionary of aggregate metrics.
        """
        return dict(self.metrics.loc["mean"])


def evaluate_md_emulation(
    indexed_samples: IndexedSamples,
    temperature_K: float = 300.0,
    random_seed: int = 42,
) -> MDEmulationResults:
    """
    Load samples, compute projections and compare free energy surfaces spanned by projections.

    Args:
        indexed_samples: `IndexedSamples` containing samples, preferably filtered.
        temperature_K: Temperature used for computing free energies in Kelvin.
        random_seed: Random seed used for analysis (mainly for resampling in computing free energy
          metric).

    Returns:
        MD emulation benchmark results data class.
    """
    # Load samples and convert lists to single trajectories.
    samples = indexed_samples.get_all_trajs()
    samples = {k: mdtraj.join(v) for k, v in samples.items()}

    # Load reference projections.
    reference_projections = load_reference_projections()

    # Load projection parameters.
    projection_params = load_projection_parameters()
    # Compute projections on sample structures.
    sample_projections = project_samples(samples, projection_params)

    # Sort outputs by test case ID for consistent outputting.
    sample_projections = dict(sorted(sample_projections.items()))

    # Compute state metrics.
    metrics = compute_state_metrics(
        sample_projections,
        reference_projections,
        temperature_K=temperature_K,
        random_seed=random_seed,
        n_resample=DistributionMetricSettings.n_resample,
        sigma_resample=DistributionMetricSettings.sigma_resample,
        num_bins=DistributionMetricSettings.num_bins,
        energy_cutoff=DistributionMetricSettings.energy_cutoff,
        padding=DistributionMetricSettings.padding,
    )

    results = MDEmulationResults(
        sample_projections=sample_projections,
        benchmark=Benchmark.MD_EMULATION,
        metrics=metrics,
        temperature_K=temperature_K,
        random_seed=random_seed,
    )

    return results
