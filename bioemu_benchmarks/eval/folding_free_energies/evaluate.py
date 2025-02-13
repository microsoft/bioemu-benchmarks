from dataclasses import dataclass
from pathlib import Path

import mdtraj
import numpy as np
import pandas as pd

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.folding_free_energies.analysis import (
    analyze_ddg,
    analyze_dg,
    compute_confidence_intervals_ddG,
    compute_confidence_intervals_dG,
)
from bioemu_benchmarks.eval.folding_free_energies.fraction_native_contacts import (
    FNCSettings,
    get_fnc_from_samples_trajectory,
)
from bioemu_benchmarks.eval.folding_free_energies.free_energies import compute_dg_ddg_from_fnc
from bioemu_benchmarks.eval.folding_free_energies.utils import load_reference
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.results import BenchmarkResults
from bioemu_benchmarks.samples import IndexedSamples
from bioemu_benchmarks.utils import StrPath

LOGGER = get_logger(__name__)


@dataclass
class FoldingFreeEnergiesResults(BenchmarkResults):
    """
    Data class for free energy evaluation results.

    Attributes:
        fnc_per_system: Dictionary containing fraction of native contacts for each sample. Keys are
          the test case identifiers.
        metrics: Pandas dataframe collecting aggregate metrics.
        free_energies_per_system: Pandas dataframe collecting per system results (free energies,
          confidence intervals).
        temperature_K: Temperature used for evaluation (in Kelvin).
    """

    fnc_per_system: dict[str, np.ndarray]
    metrics: pd.DataFrame
    free_energies_per_system: pd.DataFrame
    temperature_K: float

    def save_results(self, output_dir: StrPath) -> None:
        """
        Save individual evaluator results in accessible files (txt, csv, npz).

        Args:
            output_dir: Directory to which result outputs should be saved.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.free_energies_per_system.to_csv(output_dir / "results_systems.csv", index=False)
        self.metrics.to_csv(output_dir / "results_metrics.csv")
        np.savez(output_dir / "contact_scores.npz", **self.fnc_per_system)

    def plot(self, output_dir: StrPath) -> None:
        """
        Generate plots associated with benchmark and write to output directory.

        Args:
            output_dir: Directory where plots will be written.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Compute overall metrics and generate plots.
        _, fig_dg = analyze_dg(self.free_energies_per_system)
        _, fig_ddg = analyze_ddg(self.free_energies_per_system)

        fig_dg.savefig(output_dir / "scatter_dG.png")
        fig_ddg.savefig(output_dir / "scatter_ddG.png")

    def get_aggregate_metrics(self) -> dict[str, float]:
        """
        Collect aggregate mean absolute errors and correlation coefficients for dG and ddG.

        Returns:
            Dictionary of all aggregate metrics.
        """
        aggregate_metrics: dict[str, float] = {}
        for idx, row in self.metrics.iterrows():
            aggregate_metrics.update({f"{idx}_{k}": v for k, v in dict(row).items()})

        return aggregate_metrics


def evaluate_folding_free_energies(
    indexed_samples: IndexedSamples,
    temperature_K: float = 295,
) -> FoldingFreeEnergiesResults:
    """
    Load samples, compute free energies, measure errors compared to experiment.

    Args:
        indexed_samples: `IndexedSamples` containing samples, preferably filtered.
        temperature_K: Temperature used for computing free energies in Kelvin.

    Returns:
        Folding free energies results data class.
    """
    benchmark = Benchmark.FOLDING_FREE_ENERGIES
    # Load samples and convert lists to single trajectories.
    samples = indexed_samples.get_all_trajs()
    samples = {k: mdtraj.join(v) for k, v in samples.items()}

    # 1) Compute fraction of native contacts.
    dict_fnc: dict[str, np.ndarray] = {}
    for test_case in samples:
        # Load reference structure.
        reference_structure = load_reference(test_case)
        sample_traj = samples[test_case]

        # Compute FNC contact scores.
        contact_scores = get_fnc_from_samples_trajectory(
            sample_traj,
            reference_structure,
            sequence_separation=FNCSettings.sequence_separation,
            contact_cutoff=FNCSettings.contact_cutoff,
            contact_beta=FNCSettings.contact_beta,
            contact_lambda=FNCSettings.contact_lambda,
            contact_delta=FNCSettings.contact_delta,
        )
        dict_fnc[test_case] = contact_scores

    # 3) Compute dGs and ddGs and get experimental targets.
    system_info = pd.read_csv(Path(benchmark.asset_dir) / "system_info.csv")
    free_energies_per_system = compute_dg_ddg_from_fnc(
        dict_fnc=dict_fnc, system_info=system_info, temperature=temperature_K
    )

    # Add confidence intervals.
    free_energies_per_system = compute_confidence_intervals_dG(free_energies_per_system)
    free_energies_per_system = compute_confidence_intervals_ddG(free_energies_per_system)

    # Compute overall metrics and generate plots.
    metrics_dg, _ = analyze_dg(free_energies_per_system)
    metrics_ddg, _ = analyze_ddg(free_energies_per_system)

    metrics = pd.DataFrame(
        [metrics_dg, metrics_ddg],
        index=["dG", "ddG"],
    ).rename_axis("benchmark")

    folding_energies_results = FoldingFreeEnergiesResults(
        benchmark=benchmark,
        temperature_K=temperature_K,
        fnc_per_system=dict_fnc,
        free_energies_per_system=free_energies_per_system,
        metrics=metrics,
    )

    return folding_energies_results
