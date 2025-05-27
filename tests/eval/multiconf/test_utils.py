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

    # Modify one coord so that it's rejected. Generate a fake clash.
    # NOTE: our clash filter assumes rigid frames and prefilters based on Ca for sake of efficiency.
    # Because of this, we make the first and last full non-glycine residue clash (residues have
    # 5 atoms due to C beta).
    idx_no_glycine = traj.top.select("not resname GLY")
    traj.xyz[:, idx_no_glycine[:5], :] = traj.xyz[:, idx_no_glycine[-5:], :]

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
