import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta, spearmanr

from bioemu_benchmarks.eval.folding_free_energies.free_energies import K_BOLTZMANN
from bioemu_benchmarks.logger import get_logger

LOGGER = get_logger(__name__)


def _clopper_pearson_confidence_interval(
    predicted_dG: np.ndarray | float,
    num_samples: np.ndarray | float,
    temperature: np.ndarray | float,
    confidence: float = 0.95,
    epsilon: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Estimate confidence interval of dG values using Clopper--Pearson method.

    Args:
        predicted_dG: Array of predicted dG values of a series of systems of the shape
          [num_samples], or a single value.
        num_samples: Array of number of samples used for computing the predicted dGs (shape
          [num_samples]), or a single value.
        temperature: Array of temperature for predicting the dG values in Kelvin with the shape
          [num_samples], or a single value.
        confidence: Percentage of confidence interval.
        epsilon: Small offset for numerical stability.

    Return:
        The upper bound and lower bound of dG values representing the requested confidence interval.
    """
    ratio = np.exp(-predicted_dG / (K_BOLTZMANN * temperature))  # ratio = n_folded / n_unfolded
    p_folded = ratio / (
        1 + ratio
    )  # p_folded * N is subject to a binomial distribution, asymptotically normal
    k = np.round(p_folded * num_samples).astype(int)  # Number of folded samples

    # Calculate the confidence interval using the Clopper-Pearson method
    alpha = 1 - confidence

    p_folded_lower = np.zeros(len(k)) + epsilon
    p_folded_lower[k > 0] = beta.ppf(alpha / 2, k, num_samples - k + 1)[k > 0]  # if k > 0 else 0
    p_folded_upper = np.ones(len(k)) - epsilon
    p_folded_upper[k < num_samples] = beta.ppf(1 - alpha / 2, k + 1, num_samples - k)[
        k < num_samples
    ]  # if k < N else 1

    # Convert back to dGs.
    dG_upper = -K_BOLTZMANN * temperature * np.log(p_folded_lower / (1 - p_folded_lower))
    dG_lower = -K_BOLTZMANN * temperature * np.log(p_folded_upper / (1 - p_folded_upper))
    return dG_upper, dG_lower


def compute_confidence_intervals_dG(
    free_energy_df: pd.DataFrame, confidence: float = 0.95
) -> pd.DataFrame:
    """
    Compute confidence intervals for experimental and predicted dGs and store in dataframe.

    Args:
        free_energy_df: Input dataframe containing experimental info and model predictions.
        confidence: Confidence used for computing interval.

    Returns:
        Input dataframe with confidence intervals added.
    """
    # Extract free energy data.
    p_dGs = free_energy_df.dg_pred.values.astype(float)
    t_dGs = free_energy_df.dg_exp.values.astype(float)
    ci95_low = free_energy_df.dg_ci95_low.values.astype(float)
    ci95_high = free_energy_df.dg_ci95_high.values.astype(float)
    n_samples = free_energy_df.num_samples.values.astype(int)
    temperature = free_energy_df.temperature.values.astype(float)

    # Calculate the experimental error bars
    free_energy_df["exp_errors_dg_lower"] = t_dGs - ci95_high
    free_energy_df["exp_errors_dg_upper"] = ci95_low - t_dGs

    # Calculate the 95% confidence interval of the predicted dG
    p_dGs_upper, p_dGs_lower = _clopper_pearson_confidence_interval(
        p_dGs, n_samples, temperature, confidence=confidence
    )
    err_lower = p_dGs - p_dGs_lower
    err_upper = p_dGs_upper - p_dGs

    free_energy_df["model_errors_dg_lower"] = err_lower * (err_lower > 0)
    free_energy_df["model_errors_dg_upper"] = err_upper * (err_upper > 0)

    return free_energy_df


def compute_confidence_intervals_ddG(
    free_energy_df: pd.DataFrame, confidence: float = 0.975
) -> pd.DataFrame:
    """
    Compute confidence intervals for predicted ddGs and store in dataframe.

    Args:
        free_energy_df: Input dataframe containing experimental info and model predictions.
        confidence: Confidence used for computing interval.

    Returns:
        Input dataframe with ddG confidence intervals added.
    """

    # Filter for missing ddG values.
    missing_ddg = free_energy_df.ddg_pred.isna()
    dG_df = free_energy_df[~missing_ddg]

    p_ddGs = dG_df.ddg_pred.values.astype(float)
    p_dGs_mut = dG_df.dg_pred.values.astype(float)
    temperature = dG_df.temperature.values.astype(float)
    n_samples = dG_df.num_samples.values.astype(int)

    def get_wt_dg(x):
        try:
            dg = free_energy_df[free_energy_df.name == x.name_wt].dg_pred.values[0]
            return dg
        except IndexError:
            return None

    # Extract wild type dGs.
    p_dGs_wt = dG_df.apply(get_wt_dg, axis=1)

    # Calculate the 97.5% confidence interval of the predicted dG
    p_dGs_wt_upper, p_dGs_wt_lower = _clopper_pearson_confidence_interval(
        p_dGs_wt, n_samples, temperature, confidence=confidence
    )
    p_dGs_mut_upper, p_dGs_mut_lower = _clopper_pearson_confidence_interval(
        p_dGs_mut, n_samples, temperature, confidence=confidence
    )
    p_ddGs_lower = p_dGs_mut_lower - p_dGs_wt_upper
    p_ddGs_upper = p_dGs_mut_upper - p_dGs_wt_lower
    err_lower = p_ddGs - p_ddGs_lower
    err_upper = p_ddGs_upper - p_ddGs
    model_errors_ddG_masked_lower = err_lower * (err_lower > 0)
    model_errors_ddG_masked_upper = err_upper * (err_upper > 0)

    def unmask_errors(errors: np.ndarray) -> np.ndarray:
        unmasked = np.zeros(missing_ddg.shape)
        unmasked[~missing_ddg] = errors
        unmasked[missing_ddg] = np.NaN
        return unmasked

    free_energy_df["model_errors_ddg_lower"] = unmask_errors(model_errors_ddG_masked_lower)
    free_energy_df["model_errors_ddg_upper"] = unmask_errors(model_errors_ddG_masked_upper)

    return free_energy_df


def scatter_plot_with_errorbar(
    model_pred: np.ndarray,
    exp_targets: np.ndarray,
    quantity: str,
    axis_range: np.ndarray,
    model_errors: np.ndarray | None = None,
    exp_errors: np.ndarray | None = None,
    markersize: float = 3.0,
    label_mae: float | None = None,
    label_correlation: float | None = None,
) -> matplotlib.figure.Figure:
    """
    Create a scatter plot for predicted and target free energy values.

    Args:
        model_pred: Array of predicted values.
        exp_targets: Array of target values.
        quantity: Quantity analyzed, used for axis labeling.
        axis_range: Array of shape [2] used to determined x and y axis upper and lower bounds.
        model_errors: Optional errors for model predictions (used for error bars if provided).
        exp_errors: Optional experimental errors (used for error bars if provided).
        markersize: Marker size used in scatter plot.
        label_mae: Optional MAE value for free energy error used in labeling.
        label_correlation: Optional correlation used in labeling.

    Returns:
        Figure containing scatter plot.
    """
    fig, ax = plt.subplots(figsize=(5, 5))

    # Make the scatter plot with error bars.
    ax.errorbar(
        model_pred,
        exp_targets,
        xerr=model_errors,
        yerr=exp_errors,
        fmt="o",
        capsize=2,
        markersize=markersize,
        markeredgewidth=0.5,
        markeredgecolor="cornflowerblue",
        color="cornflowerblue",
        ecolor="cornflowerblue",
        elinewidth=0.75,
        alpha=0.7,
    )

    # Create label annotations.
    if label_mae is not None or label_correlation is not None:
        label = []
        if label_mae is not None:
            label.append(f"Error {label_mae:.2f} kcal/mol")
        if label_correlation is not None:
            label.append(f"Correlation {label_correlation:.2f}")
        ax.text(axis_range[0] + 0.5, axis_range[1] - 0.75, "\n".join(label))

    ax.plot(axis_range, axis_range, color="grey", ls=":")
    # Plot the shaded area between y=x-1 and y=x+1 with no border line (denoting area within 1 kcal
    # /mol accuracy).
    ax.fill_between(
        axis_range, axis_range - 1, axis_range + 1, color="silver", alpha=0.2, linewidth=0
    )

    # Make the plot square.
    ax.set_xlim(axis_range)
    ax.set_ylim(axis_range)
    ax.set_xlabel(f"Predicted {quantity} (kcal/mol)")
    ax.set_ylabel(f"Experimental {quantity} (kcal/mol)")
    ax.set_aspect("equal", adjustable="box")  # set equal aspect ratio

    # Set the tick labels to be the same for x and y axis.
    ax.set_xticks(np.arange(np.ceil(axis_range[0]), np.round(axis_range[1]) + 1, 1))
    ax.set_yticks(np.arange(np.ceil(axis_range[0]), np.round(axis_range[1]) + 1, 1))
    fig.tight_layout()

    return fig


def compute_error_metrics(model_pred: np.ndarray, exp_targets: np.ndarray) -> dict[str, float]:
    # Compute mean absolute error and correlation coefficients
    mae = np.mean(np.abs(model_pred - exp_targets))
    pearson_corrcoef = np.corrcoef(model_pred, exp_targets)[0, 1]
    spearman_corrcoef = spearmanr(model_pred, exp_targets)[0]

    metrics: dict[str, float] = {
        "mae": mae,
        "pearson_corrcoef": pearson_corrcoef,
        "spearman_corrcoef": spearman_corrcoef,
    }

    return metrics


def analyze_dg(
    free_energy_df: pd.DataFrame,
) -> tuple[dict[str, float], matplotlib.figure.Figure]:
    """
    Analyze dG predictions computing accuracy metrics and generating scatter plot.

    Args:
        free_energy_df: Data frame containing free energy results, experimental targets and
          uncertainties.

    Returns:
        Dictionary collecting different error metrics and analysis figure.
    """
    # Remove wild type systems which only serve as reference for ddGs.
    free_energy_df = free_energy_df.infer_objects()
    dG_df = free_energy_df[~free_energy_df.wt_only_reference.astype("bool")]

    # Extract free energy data.
    p_dGs = dG_df.dg_pred.values.astype(float)
    t_dGs = dG_df.dg_exp.values.astype(float)

    assert len(p_dGs) > 0, "No data found for delta G benchmark."

    # Get errors for error bars.
    exp_errors = [dG_df["exp_errors_dg_lower"], dG_df["exp_errors_dg_upper"]]

    # Get the 95% confidence interval of the predicted dG.
    model_errors = [dG_df["model_errors_dg_lower"], dG_df["model_errors_dg_upper"]]

    # Check for outliers. These are with thresholds 0 or 1, where folded and unfolded states are not
    # well represented.
    outliers = free_energy_df[free_energy_df.threshold.isin([0, 1])]
    if len(outliers) > 0:
        LOGGER.info(f"Found {len(outliers)} outliers for delta G benchmark:")
        LOGGER.info(outliers.name)

    # Compute correlation coefficients and MAE.
    error_metrics = compute_error_metrics(model_pred=p_dGs, exp_targets=t_dGs)

    error_plot = scatter_plot_with_errorbar(
        model_pred=p_dGs,
        exp_targets=t_dGs,
        model_errors=model_errors,
        exp_errors=exp_errors,
        quantity="ΔG",
        axis_range=np.array([-5.5, 2]),
        label_mae=error_metrics["mae"],
        label_correlation=error_metrics["pearson_corrcoef"],
    )

    return error_metrics, error_plot


def analyze_ddg(
    free_energy_df: pd.DataFrame,
) -> tuple[dict[str, float], matplotlib.figure.Figure]:
    """
    Analyze ddG predictions computing accuracy metrics and generating scatter plot.

    Args:
        free_energy_df: Data frame containing free energy results and experimental targets.

    Returns:
        Dictionary collecting different error metrics and analysis figure.
    """
    # Filter out wild type ddGs which will be trivially zero and should be excluded from analysis.
    dG_df = free_energy_df[~free_energy_df.ddg_pred.isna()]

    p_ddGs = dG_df.ddg_pred.values.astype(float)
    t_ddGs = dG_df.ddg_exp.values.astype(float)
    assert len(p_ddGs) > 0, "No data found for delta delta G benchmark."

    model_errors = [dG_df["model_errors_ddg_lower"], dG_df["model_errors_ddg_upper"]]
    error_metrics = compute_error_metrics(model_pred=p_ddGs, exp_targets=t_ddGs)

    error_plot = scatter_plot_with_errorbar(
        model_pred=p_ddGs,
        exp_targets=t_ddGs,
        model_errors=model_errors,
        exp_errors=None,
        quantity="ΔΔG",
        axis_range=np.array([-2, 5]),
        markersize=3,
        label_mae=error_metrics["mae"],
        label_correlation=error_metrics["pearson_corrcoef"],
    )

    return error_metrics, error_plot
