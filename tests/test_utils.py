from pathlib import Path

import mdtraj
import numpy as np
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.samples import IndexedSamples, find_samples_in_dir
from bioemu_benchmarks.utils import _compute_filtered_contacts

from . import TEST_DATA_DIR

CLASH_DISTANCE = 1.0
CA_FILTER_CUTOFF = 2 * 2.5 + CLASH_DISTANCE  # 2 x Ca-O in frame + clash


@pytest.fixture
def samples_path() -> Path:
    test_data_path = Path(TEST_DATA_DIR) / "samples_example" / "md_emulation"
    return test_data_path


def test_compute_filtered_contacts(samples_path):
    """
    Check in mdtraj native `compute_contacts` and internal routine with filtering based on Ca
    distances lead to the same masks.
    """
    sequence_samples = find_samples_in_dir(samples_path)
    indexed_samples = IndexedSamples.from_benchmark(
        Benchmark.MD_EMULATION, sequence_samples=sequence_samples
    )

    samples_traj = indexed_samples.get_trajs_for_test_case("cath1_1bl0A02")[0]

    # Make atoms in first frame clash.
    samples_traj.xyz[0, -5:] -= 0.4
    has_clashes_expected = np.array([False] + [True] * (len(samples_traj) - 1))

    def _get_frames_non_clash(distances: np.ndarray) -> np.ndarray:
        return np.all(
            mdtraj.utils.in_units_of(distances, "nanometers", "angstrom") > CLASH_DISTANCE, axis=1
        )

    # Routine with prefiltering.
    distances_test, _ = _compute_filtered_contacts(
        samples_traj, ca_clash_cutoff=CA_FILTER_CUTOFF, sequence_separation=2
    )
    frames_non_clash_test = _get_frames_non_clash(distances_test)

    # Basic mdtraj routine.
    distances_target, _ = mdtraj.compute_contacts(samples_traj, periodic=False)
    frames_non_clash_target = _get_frames_non_clash(distances_target)

    np.testing.assert_equal(frames_non_clash_test, frames_non_clash_target)
    np.testing.assert_equal(frames_non_clash_test, has_clashes_expected)
