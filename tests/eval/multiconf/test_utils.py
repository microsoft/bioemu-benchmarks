import os

import mdtraj
import pytest

from bioemu_benchmarks.utils import filter_unphysical_traj

from ... import TEST_DATA_DIR


def test_filter_unphysical_traj() -> None:
    test_pdb = os.path.join(
        TEST_DATA_DIR, "samples_example/multiconf_ood60/K0JNC6_1d8ee045_seed0.pdb"
    )

    traj = mdtraj.load(test_pdb)
    traj = filter_unphysical_traj(traj, strict=True)
    assert traj.n_frames == 1

    # Modify one coord so that it's rejected. Generate a fake clash
    traj.xyz[:, 0, :] = traj.xyz[:, 20, :]

    with pytest.raises(AssertionError):
        filter_unphysical_traj(traj, strict=True)

    traj = mdtraj.load(test_pdb)

    # CA-CA distance test (CA is always written in the second position)
    traj = mdtraj.load(test_pdb)
    traj.xyz[:, 1, :] += 1000

    with pytest.raises(AssertionError):
        filter_unphysical_traj(traj, strict=True)

    # C-N distance test (C is always written in the third position)
    traj = mdtraj.load(test_pdb)
    traj.xyz[:, 2, :] += 1000

    with pytest.raises(AssertionError):
        filter_unphysical_traj(traj, strict=True)
