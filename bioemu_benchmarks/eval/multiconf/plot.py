import copy
from collections import defaultdict

import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binned_statistic_2d

from bioemu_benchmarks.eval.multiconf.evaluate import MetricType
from bioemu_benchmarks.logger import get_logger

DEFAULT_FONTSIZE = 12
matplotlib.rcParams.update({"font.size": DEFAULT_FONTSIZE})
matplotlib.rc("xtick", labelsize=DEFAULT_FONTSIZE)
matplotlib.rc("ytick", labelsize=DEFAULT_FONTSIZE)

LOGGER = get_logger(__name__)


METRICS_SUCCESS_THRESHOLD: dict[MetricType, float] = {
    MetricType.RMSD: 3.0,
    MetricType.TMSCORE: 0.75,
    MetricType.LDDT: 0.75,
    MetricType.CONTACT_DISTANCE: 1.0,
    MetricType.DSSP_ACC: 0.75,
    MetricType.FNC_UNFOLD_U: 0.3,
    MetricType.FNC_UNFOLD_F: 0.7,
}

METRICS_DENSITY_VMIN: dict[MetricType, float] = {
    MetricType.RMSD: -6.0,
    MetricType.TMSCORE: -6.0,
    MetricType.LDDT: -6.0,
    MetricType.CONTACT_DISTANCE: -6.0,
    MetricType.DSSP_ACC: -6.0,
}

METRICS_DENSITY_VMAX: dict[MetricType, float] = {
    MetricType.RMSD: 2.0,
    MetricType.TMSCORE: 8.0,
    MetricType.LDDT: 8.0,
    MetricType.CONTACT_DISTANCE: 2.0,
    MetricType.DSSP_ACC: 8.0,
}

METRICS_UNIT: dict[MetricType, str] = defaultdict(lambda: "")
METRICS_UNIT[MetricType.RMSD] = "(Å)"


def plot_free_energy(
    x: np.ndarray,
    numbins: int = 20,
    figsize: tuple[int, int] = (5, 5),
    max_energy: float = 10,
    levels: int = 20,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    cbar: bool = True,
    kT: float = 1,
    ax: matplotlib.axes.Axes | None = None,
):
    """
    Plot the free energy landscape defined by a 2D dataset

    Args:
        x: array of shape (n_samples, 2)
        numbins: number of bins for the 2D histogram
        figsize: figure size
        max_energy: maximum energy to plot
        levels: number of levels for the contour plot
        x_range: x range (min, max) for the plot
        y_range: y range (min, max) for the plot
        cbar: whether to plot the colorbar
        kT: temperature factor for the energy
        ax: Axes to plot on
    """
    if x_range is None:
        x_range = (min(x[:, 0]), max(x[:, 0]))
    if y_range is None:
        y_range = (min(x[:, 1]), max(x[:, 1]))

    # shape x,y for plot
    grid_x = np.linspace(x_range[0], x_range[1], numbins)
    grid_x = np.stack([grid_x for _ in range(numbins)])
    grid_y = np.transpose(grid_x)

    # Bin the data
    binned = binned_statistic_2d(
        x[:, 0],
        x[:, 1],
        None,
        "count",
        bins=[numbins, numbins],
        range=[x_range, y_range],
    )
    grid_z = binned.statistic  # .T

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.subplots(1, 1)

    # Generate the contour plot
    if ax is None:
        plt.figure(figsize=figsize)
    energy = -np.log(grid_z + 1e-6)
    energy -= energy.min()

    # scale energy
    energy = kT * energy
    energy = np.minimum(energy, max_energy + 1)
    cmap = copy.copy(plt.cm.turbo)
    cmap.set_over(color="w")

    if ax is None:
        plt.contourf(grid_x, grid_y, energy, cmap=cmap, levels=levels, vmin=0, vmax=max_energy)
        plt.xlim(x_range[0], x_range[1])
        plt.ylim(y_range[0], y_range[1])
        plt.clim(0, max_energy)
        if cbar:
            cbar_ = plt.colorbar()
            cbar_.ax.set_ylim(0, max_energy)
            cbar_.set_label("Energy (kcal/mol)")
    else:
        ax.contourf(grid_x, grid_y, energy, cmap=cmap, levels=levels, vmin=0, vmax=max_energy)
        ax.set_xlim(x_range[0], x_range[1])
        ax.set_ylim(y_range[0], y_range[1])


def plot_free_energy_with_threshold(
    x: np.ndarray,
    success_threshold: float,
    max_range: float,
    tick_spacing: int | None = 2,
    numbins: int = 20,
    figsize: tuple[int, int] = (5, 5),
    max_energy: float = 10,
    levels: int = 20,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    cbar: bool = False,
    kT: float = 1,
    line_color: str = "black",
    ax: matplotlib.axes.Axes | None = None,
):
    """
    Generates a square free energy plot with a success threshold plotted to the specified maximum values

    Args:
        x: array of shape (n_samples, 2)
        success_threshold: success threshold to plot
        max_range: maximum range to plot
        tick_spacing: spacing for the ticks
        numbins: number of bins for the 2D histogram
        figsize: figure size
        max_energy: maximum energy to plot
        levels: number of levels for the contour plot
        x_range: x range (min, max) for the plot
        y_range: y range (min, max) for the plot
        cbar: whether to plot the colorbar
        kT: temperature factor for the energy
        line_color: color of the success threshold line
        ax: Axes to plot on
    """

    plot_free_energy(
        x,
        numbins=numbins,
        figsize=figsize,
        max_energy=max_energy,
        levels=levels,
        x_range=x_range,
        y_range=y_range,
        cbar=cbar,
        kT=kT,
        ax=ax,
    )
    if ax is None:
        plt.axhline(y=success_threshold, color=line_color, linestyle="--")
        plt.axvline(x=success_threshold, color=line_color, linestyle="--")
        if tick_spacing is not None:
            plt.xticks(np.arange(0, max_range, tick_spacing))
            plt.yticks(np.arange(0, max_range, tick_spacing))

    else:
        ax.axhline(y=success_threshold, color=line_color, linestyle="--")
        ax.axvline(x=success_threshold, color=line_color, linestyle="--")
        if tick_spacing is not None:
            ax.set_xticks(np.arange(0, max_range, tick_spacing))
            ax.set_yticks(np.arange(0, max_range, tick_spacing))


def plot_2D_free_energy_landscapes_in_grid(
    results: dict[str, np.ndarray],
    metric: MetricType,
    numbins: int = 50,
    max_range_multiplier: float = 2.0,
    success_threshold: float | None = None,
) -> matplotlib.figure.Figure:
    """
    Plots the reference vs reference energy landscape for each system in a provided dictionary

    Args:
        results: Dictionary of test cases mapped to 2 x n_samples arrays of a particular metric
        metric: The metric reported
        numbins: number of bins for the 2D histogram
        max_range_multiplier: multiplier for the maximum range determined by sample RMSD bounds
        success_threshold: threshold for success
    """

    if success_threshold is None:
        success_threshold = METRICS_SUCCESS_THRESHOLD[metric]

    systems = [s for s in results.keys()]
    systems.sort()
    n_results = len(systems)
    fig, axes = plt.subplots(int(np.ceil(n_results / 3.0)), 3, figsize=(10, 1.5 * n_results))
    axes = axes.T.flatten()
    for i, s in enumerate(systems):
        x = results[s]
        d = np.maximum(x[np.argmin(x[:, 0]), 1], x[np.argmin(x[:, 1]), 0])
        max_range = max_range_multiplier * d

        # plot PMF + scatter
        plot_free_energy_with_threshold(
            x,
            success_threshold,
            numbins=numbins,
            max_range=max_range,
            levels=20,
            max_energy=10,
            x_range=(0, max_range),
            y_range=(0, max_range),
            ax=axes[i],
            line_color="red",
            tick_spacing=None,
        )

        axes[i].set_title(s)
    return fig


def smoothed_1d_free_energy(
    x: np.ndarray,
    noise: float = 0.025,
    range: tuple[float, float] = (0, 1),
    bins: int = 100,
    noise_amplification: int = 1000,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Smooths a 1D free energy landscape by adding noise to the data prior to computing binned free energies

    Args:
        x: 1D array of values
        noise: standard deviation of the noise to add
        range: range of the free energy landscape
        bins: number of bins for the histogram
        noise_amplification: number of times to apply noise to the data

    Returns:
        A tuple of the x values and the smoothed free energies
    """
    x_noised = np.concatenate([f + noise * np.random.randn(noise_amplification) for f in x], axis=0)
    hist, bin_edges = np.histogram(x_noised, bins=bins, range=range, density=True)
    x = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return x, -np.log(hist)


def plot_smoothed_1d_free_energy(
    x: np.ndarray,
    noise: float = 0.025,
    range: tuple[float, float] = (0, 1),
    bins: int = 100,
    noise_amplification: int = 1000,
    kT: float = 0.6,
    color: str = "black",
    ax: matplotlib.axes.Axes | None = None,
):
    """
    Plots a smoothed 1D free energy landscape

    Args:
        x: 1D array of values
        noise: standard deviation of the noise to add
        range: range of the free energy landscape
        bins: number of bins for the histogram
        noise_amplification: number of times to apply noise to the data
        kT: temperature factor for the energy
        color: color of the line
        ax: Axes to plot on optionally
    """
    x, y = smoothed_1d_free_energy(
        x, noise=noise, range=range, bins=bins, noise_amplification=noise_amplification
    )
    y = kT * y
    y[np.logical_not(np.isfinite(y))] = np.max(y[np.isfinite(y)]) + 0.1
    ymin = np.floor(np.min(y) - 1.0)
    if ax is None:
        plt.figure(figsize=(5, 4))
        plt.plot(x, y, color=color, linewidth=3)
        plt.fill_between(x, ymin, y, color=color, alpha=0.2)
        plt.xlim(range[0], range[1])
        plt.ylim(ymin, np.max(y) - 0.1)
    else:
        ax.plot(x, y, color=color, linewidth=3)
        ax.fill_between(x, ymin, y, color=color, alpha=0.2)
        ax.set_xlim(range[0], range[1])
        ax.set_ylim(ymin, np.max(y) - 0.1)


def plot_free_energy_landscapes_by_fnc_in_grid(results: dict[str, np.ndarray]):
    """
    Plots free energy landscape based on compute fraction of native contacts for different systems

    Args:
        results: Dictionary of test case ids mapped to FNC arrays of length n_samples
    """
    systems = [s for s in results.keys()]
    systems.sort()
    n_results = len(systems)
    fig, axes = plt.subplots(int(np.ceil(n_results / 3.0)), 3, figsize=(10, 1.5 * n_results))
    axes = axes.T.flatten()
    for i, s in enumerate(systems):
        x = results[s]

        plot_smoothed_1d_free_energy(x, ax=axes[i])

        axes[i].set_title(s)

        if i <= n_results // 3:
            axes[i].set_ylabel("free energy (kcal/mol)")
        if (i - 1) == (n_results // 3) * 2:
            axes[i].set_xlabel("fraction of native contacts")
    return fig


def plot_coverage_bootstrap(
    thresholds: np.ndarray,
    coverages: np.ndarray,
    metric_type: MetricType,
    nsigma: int = 1,
    color: str = "black",
    label: str | None = None,
    ax: matplotlib.axes.Axes | None = None,
    success_threshold: float | None = None,
) -> float:
    """
    Takes a dictionary of test cases mapped to n_samples x n_references arrays of metrics and plots
    coverage at different thresholds using bootstrapping

    Args:
        results: Dictionary of test cases mapped to n_samples x n_references arrays of metrics
        metric_type: The metric reported
        nsuccess: Number of successes required to count a sample as covered
        nbootstrap: Number of bootstrap resamples
        nsample: Number of samples to bootstrap
        nsigma: Number of standard deviations to plot
        color: Color of the line
        label: Label of the line
        ax: Axes to plot on
        success_threshold: Threshold to plot as a dashed line

    Returns:
        Coverage at metric success threshold
    """

    if success_threshold is None:
        success_threshold = METRICS_SUCCESS_THRESHOLD[metric_type]

    ns_mean = coverages.mean(axis=0)

    if ax is None:
        ax = plt.gca()
    ax.plot(thresholds, ns_mean, color=color, linewidth=1.5, label=label)
    lower = np.maximum(0, ns_mean - nsigma * coverages.std(axis=0))
    upper = np.minimum(1, ns_mean + nsigma * coverages.std(axis=0))
    ax.fill_between(thresholds, lower, upper, color=color, alpha=0.2)

    i = np.argmin(np.abs(thresholds - success_threshold))
    ax.plot(
        [thresholds[i], thresholds[i]],
        [0, ns_mean[i]],
        color=color,
        linestyle="--",
        linewidth=1,
    )
    ax.plot(
        [0, thresholds[i]],
        [ns_mean[i], ns_mean[i]],
        color=color,
        linestyle="--",
        linewidth=1,
    )
    ax.set_xlim(0, thresholds[-1])
    ax.set_ylim(0, 1.05)

    return ns_mean[i]
