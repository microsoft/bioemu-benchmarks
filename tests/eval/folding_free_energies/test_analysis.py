import numpy as np
from matplotlib.figure import Figure
from scipy.stats import spearmanr

from bioemu_benchmarks.eval.folding_free_energies.analysis import (
    analyze_ddg,
    analyze_dg,
    compute_confidence_intervals_ddG,
    compute_confidence_intervals_dG,
    compute_error_metrics,
)

REQUIRED_COLS_DG_CONFIDENCE = [
    "exp_errors_dg_lower",
    "exp_errors_dg_upper",
    "model_errors_dg_lower",
    "model_errors_dg_upper",
]

REQUIRED_COLS_DDG_CONFIDENCE = [
    "model_errors_ddg_lower",
    "model_errors_ddg_upper",
]


def test_compute_confidence_intervals_dG(test_free_energy_results):
    """Test dG confidence interval computation."""
    # Remove confidence intervals for testing.
    dummy_df = test_free_energy_results.drop(columns=REQUIRED_COLS_DG_CONFIDENCE)

    # Compute intervals.
    results_df = compute_confidence_intervals_dG(dummy_df)

    for col in REQUIRED_COLS_DG_CONFIDENCE:
        # Check if columns were added.
        assert col in results_df.columns
        # Check if computed values match target.
        np.testing.assert_allclose(test_free_energy_results[col], results_df[col])


def test_compute_confidence_intervals_ddG(test_free_energy_results):
    """Test ddG confidence interval computation."""
    # Remove confidence intervals for testing.
    dummy_df = test_free_energy_results.drop(columns=REQUIRED_COLS_DDG_CONFIDENCE)

    # Compute intervals.
    results_df = compute_confidence_intervals_ddG(dummy_df)

    for col in REQUIRED_COLS_DDG_CONFIDENCE:
        # Check if columns were added.
        assert col in results_df.columns

        # Get NaN masks and check consistency.
        nan_test = test_free_energy_results[col].isna()
        nan_target = results_df[col].isna()
        np.testing.assert_allclose(nan_test, nan_target)

        # Check if computed values match target.
        np.testing.assert_allclose(
            test_free_energy_results[col][~nan_test], results_df[col][~nan_target]
        )


def test_compute_error_metrics():
    """Test error metric computation."""
    model_pred = np.random.randn(10)
    exp_targets = np.random.randn(10)

    # Compute target mean absolute error and correlation coefficients
    mae = np.mean(np.abs(model_pred - exp_targets))
    pearson_corrcoef = np.corrcoef(model_pred, exp_targets)[0, 1]
    spearman_corrcoef = spearmanr(model_pred, exp_targets)[0]

    metrics = compute_error_metrics(model_pred, exp_targets)

    assert "mae" in metrics
    assert "pearson_corrcoef" in metrics
    assert "spearman_corrcoef" in metrics

    np.testing.assert_allclose(metrics["mae"], mae)
    np.testing.assert_allclose(metrics["pearson_corrcoef"], pearson_corrcoef)
    np.testing.assert_allclose(metrics["spearman_corrcoef"], spearman_corrcoef)


def test_analyze_dg(test_free_energy_results, test_metrics_dg):
    """Test analysis routine for dG."""
    results, fig = analyze_dg(test_free_energy_results)

    # Check for presence of results and compare to reference.
    for key in test_metrics_dg:
        assert key in results
        np.testing.assert_allclose(results[key], test_metrics_dg[key], rtol=1e-6)

    # Check if figure was generated.
    assert isinstance(fig, Figure)


def test_analyze_ddg(test_free_energy_results, test_metrics_ddg):
    """Test analysis routine for ddG."""
    results, fig = analyze_ddg(test_free_energy_results)

    for key in test_metrics_ddg:
        # Check if results were generated.
        assert key in results

        # Check for matching NaNs (correlation) and numerical values.
        if np.isnan(test_metrics_ddg[key]):
            assert np.isnan(results[key])
        else:
            np.testing.assert_allclose(results[key], test_metrics_ddg[key])

    # Check if figure was generated.
    assert isinstance(fig, Figure)
