import numpy as np

from bioemu_benchmarks.eval.md_emulation.state_metric import (
    DistributionMetrics2D,
    compute_density_2D,
    compute_mae,
    compute_rmse,
    compute_state_metrics,
    histogram_bin_edges,
    resample_with_noise,
)

NUM_BINS = 5
PADDING = 0.1
NUM_SAMPLES = 100
SIGMA_RESAMPLE = 0.1
ENERGY_CUTOFF = 4.0


def test_histogram_bin_edges():
    """Test histogram bin generation."""
    test_data = np.random.randn(50)
    test_min = np.min(test_data)
    test_max = np.max(test_data)
    test_delta = (test_max - test_min) / (NUM_BINS + 1)

    # Test edges with no padding.
    edges_nopad = histogram_bin_edges(test_data, num_bins=NUM_BINS, padding=None)
    edges_target = np.linspace(test_min, test_max, NUM_BINS + 1)
    np.testing.assert_allclose(edges_nopad, edges_target)

    # Test edges with padding.
    edges_pad = histogram_bin_edges(test_data, num_bins=NUM_BINS, padding=PADDING)
    edges_target = np.linspace(
        test_min - test_delta * PADDING, test_max + test_delta * PADDING, NUM_BINS + 1
    )
    np.testing.assert_allclose(edges_pad, edges_target)


def test_compute_density_2D(samples_projections):
    """Test computation of densities from samples."""
    bins_x = histogram_bin_edges(samples_projections[:, 0], NUM_BINS, PADDING)
    bins_y = histogram_bin_edges(samples_projections[:, 1], NUM_BINS, PADDING)

    density = compute_density_2D(samples_projections, bins_x, bins_y)

    # Check for proper normalization.
    bin_area = (bins_x[1] - bins_x[0]) * (bins_y[1] - bins_y[0])
    np.testing.assert_allclose(density.sum() * bin_area, 1.0)

    # Compare directly against target output.
    target_density, _, _ = np.histogram2d(
        samples_projections[:, 0], samples_projections[:, 1], bins=(bins_x, bins_y), density=True
    )
    np.testing.assert_allclose(density, target_density)


def test_resample_with_noise(samples_projections):
    """Test resampling."""
    # Check basic shapes.
    resamples_noseed_nonoise = resample_with_noise(
        samples_projections, NUM_SAMPLES, sigma=0.0, rng=1
    )
    assert resamples_noseed_nonoise.shape == (NUM_SAMPLES, 2)

    # Using different seed should yield different samples.
    resamples_seed_nonoise = resample_with_noise(
        samples_projections, NUM_SAMPLES, sigma=0.0, rng=42
    )
    assert resamples_seed_nonoise.shape == (NUM_SAMPLES, 2)
    assert not np.allclose(resamples_noseed_nonoise, resamples_seed_nonoise)

    # Adding noise but keeping seed the same should still yield different samples.
    resamples_seed_noise = resample_with_noise(
        samples_projections, NUM_SAMPLES, sigma=SIGMA_RESAMPLE, rng=42
    )
    assert not np.allclose(resamples_seed_nonoise, resamples_seed_noise)

    # Using the same setting with same seed should yield the same samples.
    resamples_seed_noise_2 = resample_with_noise(
        samples_projections, NUM_SAMPLES, sigma=SIGMA_RESAMPLE, rng=42
    )
    np.testing.assert_allclose(resamples_seed_noise, resamples_seed_noise_2)

    # Using same seed but different noise settings should yield the same samples but noised in one
    # case. This is a generous check if the assumption holds and the noise is in a reasonable range.
    assert np.std(resamples_seed_noise - resamples_seed_nonoise) < 2 * SIGMA_RESAMPLE


def test_compute_rmse(free_energy_targets):
    """Test basic and minimum RMSE computation."""
    # Get energies for comparison.
    energies_predicted = free_energy_targets.energies_predicted
    energies_target = free_energy_targets.energies_target

    # Compute basic and optimized RMSE.
    rmse_base = compute_rmse(energies_predicted, energies_target, minimize=False)
    rmse_min = compute_rmse(energies_predicted, energies_target, minimize=True)

    def rmse(x, y):
        return np.sqrt(np.mean((x - y) ** 2))

    # Compare basic RMSE to online computation and target.
    np.testing.assert_allclose(rmse(energies_predicted, energies_target), rmse_base)
    np.testing.assert_allclose(rmse_base, free_energy_targets.rmse_base)

    # Compare minimized RMSE to online computation and target.
    energies_pred_centered = energies_predicted - np.mean(energies_predicted)
    energies_target_centered = energies_target - np.mean(energies_target)
    np.testing.assert_allclose(
        rmse(energies_pred_centered, energies_target_centered), free_energy_targets.rmse_minimized
    )
    np.testing.assert_allclose(rmse_min, free_energy_targets.rmse_minimized)


def test_compute_mae(free_energy_targets):
    """Test basic and minimum MAE computation."""
    # Get energies for comparison.
    energies_predicted = free_energy_targets.energies_predicted
    energies_target = free_energy_targets.energies_target

    # Compute basic and optimized MAE.
    mae_base = compute_mae(energies_predicted, energies_target, minimize=False)
    mae_min = compute_mae(energies_predicted, energies_target, minimize=True)

    # Compare basic MAE to online computation and target.
    np.testing.assert_allclose(np.mean(np.abs(energies_predicted - energies_target)), mae_base)
    np.testing.assert_allclose(mae_base, free_energy_targets.mae_base)

    # Compare minimized MAE to target value.
    np.testing.assert_allclose(mae_min, free_energy_targets.mae_minimized)


def test_distribution_metrics_2D(
    reference_projections, samples_projections, density_metric_targets
):
    """Test distribution metric class."""
    # Set up metric.
    metric = DistributionMetrics2D(
        reference_projections,
        random_seed=42,
        temperature_K=300,
        n_resample=NUM_SAMPLES,
        sigma_resample=SIGMA_RESAMPLE,
        num_bins=NUM_BINS,
        energy_cutoff=4.0,
        padding=PADDING,
    )
    # Compute base score and compare against target values.
    score_mae, score_rmse = metric.score(samples_projections)
    np.testing.assert_allclose(score_mae, density_metric_targets.score_mae)
    np.testing.assert_allclose(score_rmse, density_metric_targets.score_rmse)

    # Compute nonzero scores and compare against target values.
    score_nonzero_mae, score_nonzero_rmse, score_nonzero_coverage = metric.score_nonzero(
        samples_projections
    )
    np.testing.assert_allclose(score_nonzero_mae, density_metric_targets.score_nonzero_mae)
    np.testing.assert_allclose(score_nonzero_rmse, density_metric_targets.score_nonzero_rmse)
    np.testing.assert_allclose(
        score_nonzero_coverage, density_metric_targets.score_nonzero_coverage
    )


def test_compute_state_metrics(reference_projections, samples_projections, density_metric_targets):
    """Test computation of state metrics."""

    # TODO: relying on seeding for these tests might be wishful thinking. Maybe just remove
    #   resampling randomness.

    results = compute_state_metrics(
        dict(x=samples_projections),
        dict(x=reference_projections),
        n_resample=NUM_SAMPLES,
        sigma_resample=SIGMA_RESAMPLE,
        num_bins=NUM_BINS,
        padding=PADDING,
        energy_cutoff=ENERGY_CUTOFF,
        temperature_K=300,
        random_seed=42,
    )

    # Should have mean and results.
    assert len(results) == 2

    results_system = results.loc["x"]
    np.testing.assert_allclose(results_system.mae, density_metric_targets.score_nonzero_mae)
    np.testing.assert_allclose(results_system.rmse, density_metric_targets.score_nonzero_rmse)
    np.testing.assert_allclose(
        results_system.coverage, density_metric_targets.score_nonzero_coverage
    )
