from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import bisect

from bioemu_benchmarks.eval.folding_free_energies.free_energies import K_BOLTZMANN


@dataclass
class DistributionMetricSettings:
    """
    Data class for collecting distribution metric settings.

    Attributes:
        n_resample: Resample projections using this many points.
        sigma_resample: Standard deviation of noise added to resamples points (mainly used to
            improve binning stability. Formally corresponds to Gaussian convolution of discretized
            density.
        num_bins: Number of bins used for discretizing density.
        energy_cutoff: Energy cutoff used for computing metric (units are kcal/mol).
        padding: Padding used for discretization grid.
    """

    n_resample: int = 1000000
    sigma_resample: float = 0.25
    num_bins: int = 50
    energy_cutoff: float = 4.0
    padding: float = 0.5


def histogram_bin_edges(x: np.ndarray, num_bins: int, padding: float | None = 0.5) -> np.ndarray:
    """
    Generate histogram bin edges for provided 1D array. Creates `num_bins` + 1 edges between minimum
    and maximum of array using optional padding.

    Args:
        x: Array used for bin edge computation.
        num_bins: Number of bins.
        padding: If provided, upper and lower limits will be extended by this value times the grid
          spacing.

    Returns:
        Array of bin edges.
    """
    x_min = np.min(x)
    x_max = np.max(x)

    if padding is not None:
        delta_x = (x_max - x_min) / (num_bins + 1)
        x_min = x_min - padding * delta_x
        x_max = x_max + padding * delta_x

    return np.linspace(x_min, x_max, num_bins + 1)


def compute_density_2D(x: np.ndarray, edges_x: np.ndarray, edges_y: np.ndarray) -> np.ndarray:
    """
    Auxiliary function for computing normalized density from a two dimensional array of the shape
    given bin edges.

    Args:
        x: Array of the shape [num_samples x 2].
        edges_x: Bin edges along first dimension.
        edges_y: Bin edges along second dimension.

    Returns:
        Discretized density of the shape [num_bins_x x num_bins_y].
    """
    density, _, _ = np.histogram2d(x[:, 0], x[:, 1], bins=(edges_x, edges_y), density=True)
    return density


def resample_with_noise(
    x: np.ndarray,
    num_samples: int,
    sigma: float,
    rng: np.random.Generator | int | None = None,
) -> np.ndarray:
    """
    Resample `num_samples` points from `x` adding Gaussian noise with standard deviation `sigma`.
    This is e.g. used to make binning for statistics more robust.

    Args:
        x: Array of data to resample. Each row corresponds to a data point (shape [N x d]).
        num_samples: Number of samples to redraw.
        sigma: Standard deviation of Gaussian noise added to resamples points.
        rng: Optional random generator or integer seed for reproducibility.

    Returns:
        Array of resamples, noised points with shape [num_samples x d].
    """
    rng = np.random.default_rng(rng)
    indices = np.arange(x.shape[0])
    sel = rng.choice(indices, size=num_samples)
    return x[sel] + sigma * rng.standard_normal((num_samples, x.shape[1]))


def compute_rmse(
    energies_pred: np.ndarray, energies_target: np.ndarray, minimize: bool = True
) -> float:
    """
    Compute root mean squared error between two arrays containing energies. Optionally, energy
    difference can be minimized with respect to an integration constant between both arrays.

    Args:
        energies_pred: Array of predicted energies (shape [num_grid]).
        energies_target: Array of reference energies (shape [num_grid]).
        minimize: Minimize error by finding an optimal scalar shift between values.

    Returns:
        Computed error.
    """

    if minimize:
        # Optimal shift for RMSE / MSE is difference in means.
        energy_shift: float = np.mean(energies_target) - np.mean(energies_pred)
    else:
        energy_shift = 0.0

    energies_difference = energies_pred - energies_target + energy_shift

    return np.sqrt(np.mean(energies_difference**2))


def compute_mae(
    energies_pred: np.ndarray,
    energies_target: np.ndarray,
    minimize: bool = True,
) -> float:
    """
    Compute mean absolute error between two arrays containing energies. Optionally, energy
    difference can be minimized with respect to an integration constant between both arrays.

    Args:
        energies_pred: Array of predicted energies (shape [num_samples]).
        energies_target: Array of reference energies (shape [num_samples]).
        minimize: Minimize error by finding an optimal scalar shift between values.

    Returns:
        Computed error.
    """
    if minimize:
        # Optimal shift needs to be determined with a short numerical optimization.
        def mae_derivative(delta_energies: float) -> float:
            return np.sum(np.sign(energies_pred - energies_target + delta_energies))

        limit_lower = np.min(energies_pred) - np.max(energies_target)
        limit_upper = np.max(energies_pred) - np.min(energies_target)
        energy_shift = bisect(mae_derivative, limit_lower, limit_upper, disp=False)
    else:
        energy_shift = 0.0

    energies_difference = energies_pred - energies_target + energy_shift

    return np.mean(np.abs(energies_difference))


class DistributionMetrics2D:
    def __init__(
        self,
        reference_projections: np.ndarray,
        n_resample: int = 1000000,
        sigma_resample: float = 0.25,
        num_bins: int = 50,
        energy_cutoff: float = 4.0,
        temperature_K: float = 300.0,
        padding: float = 0.5,
        random_seed: int | None = None,
    ):
        """
        Class for computing free energy mean absolute and root mean squared errors between reference
        and sample densities computed from low dimensional projections based on protein structures.

        This routine follows the overall procedure:

            1) Resample data and add Gaussian noise.
            2) Discretize resampled data and normalize to density.
            3) Either select bin coordinates where reference probabilities are greater than a cutoff
               and clamp low sample probabilities to this cutoff (`score`) or
               select bin coordinates where there are reference and sample densities
               (`score_nonzero`).
            4) Compute free energies on those bins.
            5) Compute minimized mean absolute error and root mean squared error (optimizing global
               energy offset).

        Args:
            reference_projections: Array of reference projections (shape [N x 2]).
            n_resample: Resample projections using this many points.
            sigma_resample: Standard deviation of noise added to resamples points (mainly used to
                improve binning stability. Formally corresponds to Gaussian convolution of discretized
                density.
            num_bins: Number of bins used for discretizing density.
            energy_cutoff: Energy cutoff used for computing metric (units are kcal/mol).
            temperature_K: Temperature used for analysis in Kelvin.
            padding: Padding used for discretization grid.
            random_seed: Random seed for resampling.
        """
        self.n_resample = n_resample
        self.sigma_resample = sigma_resample
        self.kBT = temperature_K * K_BOLTZMANN
        self.energy_cutoff = energy_cutoff
        self.random_seed = random_seed

        # Resample reference projections.
        reference_projections_noised = resample_with_noise(
            reference_projections, n_resample, sigma_resample, rng=self.random_seed
        )
        # Get bin edges for discretization.
        self.edges_x = histogram_bin_edges(
            reference_projections_noised[:, 0], num_bins, padding=padding
        )
        self.edges_y = histogram_bin_edges(
            reference_projections_noised[:, 1], num_bins, padding=padding
        )
        # Compute discretized density.
        self.density_ref = compute_density_2D(
            reference_projections_noised, self.edges_x, self.edges_y
        )

        # Compute reference density mask based on energy cutoff,
        p_cutoff = self._compute_density_cutoff(self.density_ref)
        self.low_energy_mask = self.density_ref > p_cutoff

    def _compute_density_cutoff(self, density: np.ndarray) -> float:
        """
        Auxiliary function for computing probability cutoff based on energy threshold.

        Args:
            density: Density for which cutoff should be computed.

        Returns:
            Density cutoff corresponding to energy cutoff.
        """
        energy_min = -self.kBT * np.log(np.max(density))
        return np.exp(-(energy_min + self.energy_cutoff) / self.kBT)

    def score(self, sample_projections: np.ndarray) -> tuple[float, float]:
        """
        Compute free energy errors between sample and reference densities. Region where sample free
        energies are below the probability cutoff are set to the cutoff.

        Args:
            sample_projections: Sample projections with shape [M x 2].

        Returns:
            Optimized mean absolute and root mean squared errors.
        """

        # Resample sample projections and compute discretized density.
        sample_projections_noised = resample_with_noise(
            sample_projections, self.n_resample, self.sigma_resample, rng=self.random_seed
        )
        sample_density = compute_density_2D(sample_projections_noised, self.edges_x, self.edges_y)

        # Pad low probability regions.
        p_cutoff = self._compute_density_cutoff(sample_density)
        density_padded = np.maximum(sample_density, p_cutoff)

        # Compute free energy surfaces based on low probability mask.
        energy_samples = -self.kBT * np.log(density_padded[self.low_energy_mask])
        energy_ref = -self.kBT * np.log(self.density_ref[self.low_energy_mask])

        # Compute minimal errors between surfaces.
        mae_min = compute_mae(energy_samples, energy_ref, minimize=True)
        rmse_min = compute_rmse(energy_samples, energy_ref, minimize=True)

        return mae_min, rmse_min

    def score_nonzero(self, sample_projections: np.ndarray) -> tuple[float, float, float]:
        """
        Compute free energy errors between sample and reference densities scores in regions where
        there are reference data points.

        Args:
            sample_projections: Sample projections with shape [M x 2].

        Returns:
            Optimized mean absolute, root mean squared errors and coverage (fraction of grid states
            where there are samples).
        """
        # Resample sample projections and compute discretized density.
        sample_projections_noised = resample_with_noise(
            sample_projections, self.n_resample, self.sigma_resample, rng=self.random_seed
        )
        sample_density = compute_density_2D(sample_projections_noised, self.edges_x, self.edges_y)

        # Select subset of states where the reference energy is low and where there are samples in
        # the model density.
        common_mask = np.logical_and(self.low_energy_mask, sample_density > 0)

        # Compute free energy surfaces based on common mask.
        energy_ref = -self.kBT * np.log(self.density_ref[common_mask])
        energy_samples = -self.kBT * np.log(sample_density[common_mask])

        # Compute minimal errors between surfaces and coverage of samples counted based on above
        # threshold criterion.
        mae_min = compute_mae(energy_samples, energy_ref, minimize=True)
        rmse_min = compute_rmse(energy_samples, energy_ref, minimize=True)
        coverage = np.count_nonzero(common_mask) / np.count_nonzero(self.low_energy_mask)

        return mae_min, rmse_min, coverage


def compute_state_metrics(
    sample_projections: dict[str, np.ndarray],
    reference_projections: dict[str, np.ndarray],
    temperature_K: float = 300.0,
    random_seed: int = 42,
    n_resample: int = 1000000,
    sigma_resample: float = 0.25,
    num_bins: int = 50,
    energy_cutoff: float = 4.0,
    padding: float = 0.5,
) -> pd.DataFrame:
    """
    Compute state metrics (free energy mean absolute and root mean squared errors in kcal/mol and
    state coverages) for individual systems and on average.

    Args:
        sample_projections: Sample projections.
        reference_projections: Target projections.
        temperature_K: Temperature used for evaluation in Kelvin (used to scale free energy
          surfaces).
        random_seed: Random seed for reproducibility (used for resampling and noising).
        n_resample: Resample projections using this many points.
        sigma_resample: Standard deviation of noise added to resamples points (mainly used to
            improve binning stability. Formally corresponds to Gaussian convolution of discretized
            density.
        num_bins: Number of bins used for discretizing density.
        energy_cutoff: Energy cutoff used for computing metric (units are kcal/mol).
        padding: Padding used for discretization grid.

    Returns:
        Dataframe with per system and aggregate metrics.
    """

    results_mae: dict[str, float] = {}
    results_rmse: dict[str, float] = {}
    results_coverage: dict[str, float] = {}
    for test_case in sample_projections:
        # Compute non zero distribution score with standard settings.
        metric = DistributionMetrics2D(
            reference_projections[test_case],
            random_seed=random_seed,
            temperature_K=temperature_K,
            n_resample=n_resample,
            sigma_resample=sigma_resample,
            num_bins=num_bins,
            energy_cutoff=energy_cutoff,
            padding=padding,
        )
        mae, rmse, coverage = metric.score_nonzero(sample_projections[test_case])

        # Collect state metrics.
        results_mae[test_case] = mae
        results_rmse[test_case] = rmse
        results_coverage[test_case] = coverage

    # Compute means and add entries.
    def _dict_mean(result_dict: dict[str, float]) -> float:
        return np.mean(list(result_dict.values()))

    results_mae["mean"] = _dict_mean(results_mae)
    results_rmse["mean"] = _dict_mean(results_rmse)
    results_coverage["mean"] = _dict_mean(results_coverage)

    # Create dataframe and rename index column.
    results_df = pd.DataFrame(
        [results_mae, results_rmse, results_coverage], index=["mae", "rmse", "coverage"]
    ).T
    results_df.index.name = "test_case"

    return results_df
