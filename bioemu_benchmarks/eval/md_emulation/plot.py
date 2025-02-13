import copy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import binned_statistic_2d

from bioemu_benchmarks.eval.folding_free_energies.free_energies import K_BOLTZMANN

ENERGY_METRICS: list[str] = ["mae", "rmse"]


def plot_free_energy_on_axes(
    axes: Axes,
    projections: np.ndarray,
    num_bins: int = 20,
    max_energy: float = 10.0,
    levels: int = 10,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    add_colorbar: bool = False,
    colorbar_axis: Axes | None = None,
    kBT: float = 1.0,
) -> None:
    """
    Convert sample projections to free energy surfaces and plot a contour plot on a given axes
    object.

    Args:
        axes: matplotlib axes object.
        projections: Projections to plot. Should be array of shape [n_samples x 2].
        num_bins: Number of bins used to discretize the free energy surface.
        max_energy: Clamp the free energy surface to this value. Units are the same as `kBT`.
        levels: Number of levels used in the contour plot.
        x_range: Optional range of x axis values.
        y_range: Optional range of y axis values.
        add_colorbar: Add color bar to plot (requires `colorbar_axis`).
        colorbar_axis: Axes object to which color bar should be added.
        kBT: Product of Boltzmann constant and temperature. Defines units in the plot.
    """

    # Extract x and y coordinates.
    projections_x = projections[:, 0]
    projections_y = projections[:, 1]

    # Set up plot ranges if none are provided.
    if x_range is None:
        x_range, _ = (min(projections_x), max(projections_x))
    if y_range is None:
        y_range = (min(projections_y), max(projections_y))

    # Set up x and y grids for plot.
    grid_x, grid_y = np.mgrid[
        x_range[0] : x_range[1] : num_bins * 1j, y_range[0] : y_range[1] : num_bins * 1j  # type: ignore
    ]

    # Discretize data distribution.
    binned = binned_statistic_2d(
        projections_x,
        projections_y,
        None,
        "count",
        bins=[num_bins, num_bins],
        range=[x_range, y_range],
    )

    # Generate the contour plot. The offset here is for numerical stability, the energy limit used
    # for plotting is handled purely by `max_energy`.
    # p = exp(-e/kT) => e = -kT ln(p)
    energy = -kBT * np.log(binned.statistic + 1e-12)

    # Remove minimum energy and clamp to maximum for consistent plotting.
    energy -= energy.min()
    energy = np.minimum(energy, max_energy + 1)

    # Change cutoff color to white.
    cmap = copy.copy(plt.cm.turbo)
    cmap.set_over(color="w")

    # Add contours.
    contours = axes.contourf(
        grid_x, grid_y, energy, cmap=cmap, levels=levels, vmin=0, vmax=max_energy
    )
    contours.set_clim(0, max_energy)

    # Add color bar if requested.
    if add_colorbar:
        assert colorbar_axis is not None, "Colorbar axis required for color bar."
        cbar_ = axes.figure.colorbar(contours, cax=colorbar_axis)
        cbar_.ax.set_ylim(0, max_energy)


def _get_axis_ranges(data: list[np.ndarray], buffer: float = 0.1) -> list[tuple[float, float]]:
    """
    Utility function for finding minium and maximum values for each dimension across a list of
    arrays. Used to determine common plotting ranges for samples and reference projections.

    Args:
        data: List of arrays for which ranges should be computed. They need not match in the first
          dimension, but the feature dimension should be the same, e.g [[N x F], [M x F], ...].
        buffer: Buffer used to extend the ranges for nicer plotting. It specifies the fraction of
          the total range added to each side of the range.

    Returns:
        List of tuples containing minimum and maximum values. The length of the list corresponds to
        the number of feature dimensions (F in the above example).
    """
    data_mins = np.min(np.array([np.min(x, axis=0) for x in data]), axis=0)
    data_maxs = np.max(np.array([np.max(x, axis=0) for x in data]), axis=0)

    # Apply buffer to ranges.
    buffer_zone = buffer * (data_maxs - data_mins)
    data_mins -= buffer_zone
    data_maxs += buffer_zone

    return list(zip(data_mins, data_maxs))


def plot_projections(
    sample_projections: dict[str, np.ndarray],
    reference_projections: dict[str, np.ndarray],
    sort_alphabetically: bool = True,
    num_bins: int = 40,
    temperature_K: float = 300.0,
    max_energy: float = 10,
    levels: int = 10,
) -> Figure:
    """
    Plot the free energy surfaces of reference and model sample projections side by side for each
    system.

    Args:
        sample_projections: Dictionary of projections computed from model samples. Keys are the
          test case IDs, the general shape is [num_samples_model x 2].
        reference_projections: Dictionary of reference projections. Keys are the test case IDs, the
          general shape is [num_samples_ref x 2].
        sort_alphabetically: Use alphabetic order (by test case ID) for plotting.
        num_bins: Number of bins used to discretize free energy surfaces.
        temperature_K: Temperature in Kelvin used for computing free energy surfaces.
        max_energy: Maximum energy considered for plotting. Everything higher will be clamped to
          this value. Units are kcal/mol, range should be adjusted when changing temperature or the
          number of samples.
        levels: Number of levels used in free energy contour plot.

    Returns:
        Figure object containing the plots.
    """
    # Get effective Boltzmann factor and number of systems.
    kBT = K_BOLTZMANN * temperature_K
    n_systems = len(sample_projections)

    # Set up basic subplot structure (last, narrow column is used for color bars).
    fig, axs = plt.subplots(
        n_systems,
        3,
        figsize=(6, 2.55 * n_systems),
        gridspec_kw={"width_ratios": [1, 1, 0.05]},
    )

    if n_systems == 1:
        axs = [axs]

    # Sort plots alphabetically by test_case.
    if sort_alphabetically:
        sample_projections = dict(sorted(sample_projections.items()))

    for n, test_case in enumerate(sample_projections):
        # Get sample and reference projections for system.
        reference_projection = reference_projections[test_case]
        sample_projection = sample_projections[test_case]

        # Get plotting axes and share axes of first two subplots.
        ax = axs[n]
        ax[0].sharex(ax[1])
        ax[0].sharey(ax[1])
        ax[1].tick_params(labelleft=False)

        # Get axis ranges.
        x_range, y_range = _get_axis_ranges([reference_projection, sample_projection], buffer=0.5)

        # Plot reference surface in first column.
        plot_free_energy_on_axes(
            ax[0],
            reference_projection,
            add_colorbar=False,
            kBT=kBT,
            max_energy=max_energy,
            levels=levels,
            num_bins=num_bins,
            x_range=x_range,
            y_range=y_range,
        )
        # Plot sample surface in second and color bar in third column.
        plot_free_energy_on_axes(
            ax[1],
            sample_projection,
            add_colorbar=True,
            colorbar_axis=ax[2],
            kBT=kBT,
            max_energy=max_energy,
            levels=levels,
            num_bins=num_bins,
            x_range=x_range,
            y_range=y_range,
        )

        # Set up column title, use test case ID as y axis label and specify tick format.
        ax[0].set_title("MD reference")
        ax[0].set_ylabel(f"{test_case}")
        ax[0].xaxis.set_major_formatter(FormatStrFormatter("%4.1f"))
        ax[0].yaxis.set_major_formatter(FormatStrFormatter("%4.1f"))

        # Set up column title for samples and specify tick format.
        ax[1].set_title(f"model ({sample_projection.shape[0]} samples)")
        ax[1].xaxis.set_major_formatter(FormatStrFormatter("%4.1f"))
        ax[1].yaxis.set_major_formatter(FormatStrFormatter("%4.1f"))

    fig.tight_layout()

    return fig


def plot_metrics(
    metrics: pd.DataFrame,
    label_map: dict[str, str],
    energy_hline: float | None = 1.0,
) -> sns.axisgrid.PairGrid:
    """
    Plot different metrics for each system.

    Args:
        metrics: Data frame with metric information (index should be test case ID, metric values in
          columns).
        label_map: Dictionary mapping column name in data frame to label used in plot.
        energy_hline: If provided, plot energy limit in energy related metrics (currently MAE and
          RMSE). Uses the same units as evaluation energies (default is kcal/mol).

    Returns:
        Seaborn axis grid with plots.
    """
    # Convert index to column for easier handling.
    metrics = metrics.reset_index()

    # Set up initial axis grid.
    width = 5
    height = 2.2 * len(label_map)
    axis_grid = sns.PairGrid(
        metrics,
        y_vars=label_map.keys(),
        x_vars=["test_case"],
        height=height,
        aspect=width / (height / len(label_map)),
    )

    # Add plot using stripplot.
    axis_grid.map(
        sns.stripplot,
        size=7,
        orient="v",
        jitter=False,
        linewidth=1,
        edgecolor="none",
        color="k",
        legend="brief",
    )

    # Remove borders and labels.
    sns.despine(top=True, right=False, left=True, bottom=False)
    axis_grid.set(xlabel="", ylabel="")

    # Use consistent x limit for all plots.
    xlim = axis_grid.axes[0, 0].get_xlim()

    for x, ax in zip(label_map, axis_grid.axes.flat):
        ax.set(title=None, ylim=(0, metrics[x].max() * 1.2), xlim=xlim)
        ax.xaxis.grid(True)
        ax.yaxis.grid(False)
        ax.tick_params(axis="x", labelrotation=90)
        ax.set_ylabel(label_map[x])
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.yaxis.set_major_formatter(lambda x, pos: "{:3.2f}".format(x))

        # Add energy threshold where relevalt.
        if energy_hline is not None and x in ENERGY_METRICS:
            ax.hlines(energy_hline, *xlim, linewidth=0.5, color="grey")

    # Finalize figure.
    axis_grid.figure.set_size_inches(9, 1.85 * len(label_map) + 1.0)
    axis_grid.figure.tight_layout()

    return axis_grid
