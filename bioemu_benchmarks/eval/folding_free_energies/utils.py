from pathlib import Path

import mdtraj

from bioemu_benchmarks.benchmarks import Benchmark


def load_reference(test_case: str) -> mdtraj.Trajectory:
    """
    Load reference structure based on test case ID and benchmark type.
    Note that we load the wild type reference structure for mutants as well.

    Args:
        test_case: which test case to load the reference structure for.

    Returns:
        Loaded reference structure in mdtraj Trajectory format.
    """
    name_wt = test_case.split("__")[0]  # Get the wild type name from the test case ID.
    path = list(Path(Benchmark.FOLDING_FREE_ENERGIES.asset_dir).glob(f"**/{name_wt}.pdb"))
    assert len(path) == 1, f"Expected 1 reference structure for {test_case}, found {len(path)}."
    return mdtraj.load_pdb(path[0])
