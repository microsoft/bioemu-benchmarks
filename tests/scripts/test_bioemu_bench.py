import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.scripts.bioemu_bench import (
    BENCHMARK_CHOICES,
    benchmarks_from_choices,
    get_sampling_specs,
    save_sampling_specs,
)

ROOT_DIR = Path(__file__).parents[2]
TEST_SAMPLES_DIR = ROOT_DIR / "tests" / "test_data" / "samples_example"
SCRIPT_PATH = ROOT_DIR / "bioemu_benchmarks" / "scripts" / "bioemu_bench.py"
EXPECTED_SCRIPT_EVAL_OUTPUTS: dict[str, list[str]] = {
    "folding_free_energies": [
        "contact_scores.npz",
        "results_metrics.csv",
        "results.pkl",
        "results_systems.csv",
        "scatter_ddG.png",
        "scatter_dG.png",
    ],
    "md_emulation": [
        "metrics.png",
        "projections.png",
        "results_metrics.csv",
        "results.pkl",
        "results_projections.npz",
    ],
    "multiconf_ood60": [
        "closest",
        "multiconf_ood60_rmsd_coverage.png",
        "multiconf_ood60_rmsd_free_energy.png",
        "results.h5",
        "results.pkl",
    ],
    "singleconf_localunfolding": [
        "results.h5",
        "results.pkl",
        "singleconf_localunfolding_fnc_unfold_f_coverage.png",
        "singleconf_localunfolding_fnc_unfold_f_free_energy.png",
        "singleconf_localunfolding_fnc_unfold_u_coverage.png",
        "singleconf_localunfolding_fnc_unfold_u_free_energy.png",
    ],
    "multiconf_crypticpocket": [
        "results.h5",
        "results.pkl",
        "closest",
        "multiconf_crypticpocket_apo_rmsd_free_energy.png",
        "multiconf_crypticpocket_apo_rmsd_coverage.png",
        "multiconf_crypticpocket_holo_rmsd_free_energy.png",
        "multiconf_crypticpocket_holo_rmsd_coverage.png",
    ],
    "multiconf_domainmotion": [
        "results.h5",
        "results.pkl",
        "closest",
        "multiconf_domainmotion_rmsd_coverage.png",
        "multiconf_domainmotion_rmsd_free_energy.png",
    ],
}


@pytest.mark.parametrize("benchmark_choice", BENCHMARK_CHOICES)
def test_benchmarks_from_choices(benchmark_choice):
    """Test mapping of script options to benchmarks."""
    benchmarks = benchmarks_from_choices([benchmark_choice])

    if benchmark_choice == "all":
        assert set(benchmarks) == set([b for b in Benchmark])
    else:
        assert set(benchmarks) == set([Benchmark(benchmark_choice)])


@pytest.mark.parametrize("benchmark", [b for b in Benchmark])
def test_get_sampling_specs(benchmark):
    """Test set up of sample specifications based on benchmarks."""
    sampling_specs = get_sampling_specs(benchmark)

    # Specs should at least have these columns.
    assert set(["test_case", "sequence", "num_samples"]).issubset(sampling_specs.columns)

    # Check if required columns match targets.
    assert sampling_specs.test_case.equals(benchmark.metadata.test_case)
    assert sampling_specs.sequence.equals(benchmark.metadata.sequence)
    assert np.all(sampling_specs.num_samples == benchmark.default_samplesize)


@pytest.mark.parametrize("benchmark_choice", BENCHMARK_CHOICES)
def test_save_sampling_specs(benchmark_choice, tmpdir):
    """Test saving of sampling specifications to file."""
    tmp_csv = tmpdir / "test_specs.csv"

    benchmarks = benchmarks_from_choices([benchmark_choice])
    save_sampling_specs(benchmarks, tmp_csv)

    # Check if output file was created.
    assert tmp_csv.exists()

    # Load file and check against target values.
    test_specs = pd.read_csv(tmp_csv)

    target_specs = []
    for benchmark in benchmarks:
        target_specs.append(get_sampling_specs(benchmark))

    target_specs = pd.concat(target_specs)
    target_specs.reset_index(drop=True, inplace=True)

    pd.testing.assert_frame_equal(test_specs, target_specs)


def test_run_bioemu_bench_eval(tmpdir):
    """Check if CLI script runs and generates expected outputs."""

    # Target benchmarks are the ones for which test samples exist.
    target_benchmarks = list(EXPECTED_SCRIPT_EVAL_OUTPUTS.keys())
    target_sample_dirs = [TEST_SAMPLES_DIR / benchmark for benchmark in target_benchmarks]

    # Generate command running on all benchmarks with test samples. Filtering needs to be disabled
    # since multiconf evaluators would complain otherwise.
    command = (
        ["python", SCRIPT_PATH, "eval", tmpdir]
        + ["--benchmarks"]
        + target_benchmarks
        + ["--sample_dirs"]
        + target_sample_dirs
        + ["--skip_filtering"]
    )

    # Set PYTHONPATH for script.
    env = os.environ.copy()
    env["PYTHONPATH"] = ROOT_DIR

    # Run script and check exit status.
    results = subprocess.run(command, env=env)
    assert results.returncode == 0

    # Check if collected metric file was generated.
    assert (tmpdir / "benchmark_metrics.json").exists()

    # Check if all individual benchmark files were written.
    for benchmark in target_benchmarks:
        target_output_files = EXPECTED_SCRIPT_EVAL_OUTPUTS[benchmark]
        for target_file in target_output_files:
            target_path = tmpdir / benchmark / target_file
            assert target_path.exists()


def test_run_bioemu_bench_specs(tmpdir):
    """Test specs mode of CLI script."""
    target_benchmarks = list(EXPECTED_SCRIPT_EVAL_OUTPUTS.keys())
    output_csv = tmpdir / "test_specs.csv"

    # Mode should be specs, target benchmarks are the ones for which test samples exists.
    command = ["python", SCRIPT_PATH, "specs", output_csv] + ["--benchmarks"] + target_benchmarks

    # Set PYTHONPATH for script.
    env = os.environ.copy()
    env["PYTHONPATH"] = ROOT_DIR

    # Run script and check exit status.
    results = subprocess.run(command, env=env)
    assert results.returncode == 0

    # Check if output file was generated.
    assert output_csv.exists()

    # Check if output file contains target columns and if number of systems is as expected.
    test_specs = pd.read_csv(output_csv)
    assert set(["test_case", "sequence", "num_samples"]).issubset(test_specs.columns)
    num_systems = sum([len(get_sampling_specs(Benchmark(b))) for b in target_benchmarks])
    assert len(test_specs) == num_systems
