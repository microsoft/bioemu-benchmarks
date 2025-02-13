import os
import stat
import subprocess
import tarfile
from subprocess import CompletedProcess

import mdtraj
import numpy as np
from Bio import Align

from bioemu_benchmarks.aa_utils import get_aa1code_from_aa3code
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.paths import UTILS_DIR
from bioemu_benchmarks.utils import StrPath, download

LOGGER = get_logger(__name__)


def setup_tm_align() -> StrPath:
    """
    Sets up TM-align under a default directory, returns location of binary.
    """
    os.makedirs(UTILS_DIR, exist_ok=True)
    us_align_exec = os.path.join(UTILS_DIR, "bin", "USalign")
    if not os.path.exists(us_align_exec):
        # Download from source
        LOGGER.info("Downloading tm-score...")
        us_align_src = "https://anaconda.org/bioconda/usalign/2024.07.30/download/linux-64/usalign-2024.07.30-h503566f_1.tar.bz2"
        us_align_dst = os.path.join(UTILS_DIR, os.path.basename(us_align_src))
        download(us_align_src, us_align_dst)

        with tarfile.open(us_align_dst, "r") as usalign_tar:
            usalign_tar.extractall(UTILS_DIR)

        os.chmod(
            us_align_exec,
            stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH,
        )
    return us_align_exec


def tm_align(
    mobile_pdbfile: StrPath,
    target_pdbfile: StrPath,
    outfile: StrPath,
    us_align_exec: StrPath | None = None,
) -> CompletedProcess:
    """Wrapper over TM Score executable. Performs sequence alignment between proteins
    in `pdbfile_i` and `pdbfile_j` before structural alignment + scoring

    Args:
        mobile_pdbfile: PDB file of mobile protein
        target_pdbfile: PDB file of reference protein
        outfile: Aligned mobile PDB file
        us_align_exec (optional): Location to `USalign` executable
                      If not provided it is downloaded from `USALIGN_REMOTE`. Defaults to None.
    """

    if us_align_exec is None:
        us_align_exec = setup_tm_align()

    outfile_name = os.path.splitext(os.path.basename(outfile))[0]

    res = subprocess.run(
        args=[
            us_align_exec,
            mobile_pdbfile,
            target_pdbfile,
            "-mol",
            "prot",
            "-TMscore",
            "5",
            "-o",
            os.path.join(os.path.dirname(outfile), outfile_name),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if res.returncode != 0:
        raise ValueError(f"Something went wrong during alignment: {res.stderr.decode()}")
    return res


def range_to_all_matching_resids(matching_ranges: np.ndarray) -> list[tuple[int, int]]:
    """
    Converts a list of aligned residue index ranges between two sequences into a list of tuples of explicit aligned
    residue indices.

    Args:
        matching_ranges : List of aligned residue ranges, as present in the `Bio.Alignment.aligned` attribute
                            (e.g., (e.g. [[[0, 154], [0, 154]], [[156, 168], [156, 168]], ...]))

    Returns:
        A list of tuples with explicitly aligned indices between two sequences(e.g. [(0, 0), (1, 1), ...])
    """
    ranges_i, ranges_j = matching_ranges

    matching_resids: list[tuple[int, int]] = []
    for range_row_i, range_row_j in zip(ranges_i, ranges_j):
        start_i, end_i = range_row_i
        start_j, end_j = range_row_j

        start_to_end_i = range(start_i, end_i)
        start_to_end_j = range(start_j, end_j)
        matching_resids.extend(list(zip(start_to_end_i, start_to_end_j)))
    return matching_resids


def seq_pairwise_align_trajs(
    traj_i: mdtraj.Trajectory, traj_j: mdtraj.Trajectory
) -> list[tuple[int, int]]:
    """Gets matching pairs of residues between `traj_i` and `traj_j` via
    global pairwise sequence alignment.

    Args:
        traj_i: First trajectory
        traj_j: Second trajectory

    Returns:
        List of tuples with matching residue indices between `traj_i` and `traj_j`.
    """

    def _get_seq_resid_traj(traj: mdtraj.Trajectory) -> tuple[list[int], str]:
        residues = list(traj.topology.residues)
        resids = [r.resSeq for r in residues]
        resnames = [r.name for r in residues]
        return resids, "".join([get_aa1code_from_aa3code(r) for r in resnames])

    resid_i, seq_i = _get_seq_resid_traj(traj_i)
    resid_j, seq_j = _get_seq_resid_traj(traj_j)
    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-0.5)
    alignment = aligner.align(seq_i, seq_j)[
        0
    ]  # Get first alignment (max score) if more than one available
    matching_ranges = alignment.aligned
    matching_resid_zero_idx = range_to_all_matching_resids(matching_ranges)
    return [(resid_i[i], resid_j[j]) for i, j in matching_resid_zero_idx]


def fmt_resseq_selection(resids: list[int]) -> str:
    """mdtraj resid selection string generator"""
    base_str = "resSeq {} or "
    current_str = ""
    for resid in resids:
        current_str += base_str.format(resid)
    return current_str.rstrip("or ")


def align_traj_to_ref_using_range(
    ref_traj: mdtraj.Trajectory, traj: mdtraj.Trajectory, matching_resids: list[tuple[int, int]]
) -> mdtraj.Trajectory:
    """
    Structurally aligns `traj` to `ref_traj` using the matching residue indices provided
    under `matching_resids`

    Args:
        ref_traj: Reference trajectory instance
        traj: Mobile trajectory instance
        matching_resids: List of matching sequence-aligned resids (idx_i: idx_j) between
                         `ref_traj` and `traj`.
    """
    # Align sample trajectory to reference trajectory using user-defined alignment residues
    sel_alignment_ref = fmt_resseq_selection([m[0] for m in matching_resids])
    sel_alignment_traj = fmt_resseq_selection([m[1] for m in matching_resids])

    atomsel_alignment_ref = ref_traj.topology.select(sel_alignment_ref)
    atomsel_alignment_traj = traj.topology.select(sel_alignment_traj)

    assert len(atomsel_alignment_ref) == len(
        atomsel_alignment_traj
    ), "Uneven alignment atom selections"

    return traj.superpose(
        ref_traj,
        atom_indices=atomsel_alignment_traj,
        ref_atom_indices=atomsel_alignment_ref,
    )
