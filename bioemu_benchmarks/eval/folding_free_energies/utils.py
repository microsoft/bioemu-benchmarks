from pathlib import Path

import mdtraj

from bioemu_benchmarks.benchmarks import Benchmark


def load_reference(test_case: str) -> mdtraj.Trajectory:
    """
    Load reference structure based on test case ID and benchmark type.

    Args:
        test_case: which test case to load the reference structure for.

    Returns:
        Loaded reference structure in mdtraj Trajectory format.
    """
    path = list(Path(Benchmark.FOLDING_FREE_ENERGIES.asset_dir).glob(f"**/{test_case}.pdb"))
    assert len(path) == 1, f"Expected 1 reference structure for {test_case}, found {len(path)}."
    return mdtraj.load_pdb(path[0])
