import os
import tempfile
from collections import defaultdict
from itertools import combinations

import mdtraj
import numpy as np

from bioemu_benchmarks.eval.multiconf.align import (
    fmt_resseq_selection,
    seq_pairwise_align_trajs,
    tm_align,
)
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.utils import StrPath

LOGGER = get_logger(__name__)


def tm_score(pdbfile_i: StrPath, pdbfile_j: StrPath, us_align_exec: StrPath | None = None) -> float:
    """TM score between `pdbfile_i` and `pdbfile_j`. Caution: this performs
    structural alignment of `pdbfile_i` to `pdbfile_j`.

    Args:
        pdbfile_i: Mobile PDB file
        pdbfile_j: Static (target) PDB file
        us_align_exec (optional): location of `USalign` binary.
    """
    _seq_aln_error_msg = "There is no alignment between the two structures"
    with tempfile.TemporaryDirectory() as tempdir:
        res = tm_align(
            pdbfile_i,
            pdbfile_j,
            outfile=os.path.join(tempdir, "out.pdb"),
            us_align_exec=us_align_exec,
        )

    stdout = res.stdout.decode().split("\n")
    if stdout[0].startswith(_seq_aln_error_msg):
        LOGGER.error(_seq_aln_error_msg)
        return 0.0
    # 'User-specified initial alignment: TM/Lali/rmsd = 0.75292,  295,  4.028'
    tm_score_i = float(stdout[16].split("=")[1].split("(")[0])
    tm_score_j = float(stdout[17].split("=")[1].split("(")[0])
    return max(tm_score_i, tm_score_j)


def tm_score_local(pairwise_distances: np.ndarray) -> np.ndarray:
    """
    TM-score implementation over a local region.
    Formula taken from https://en.wikipedia.org/wiki/Template_modeling_score#The_TM-score_equation

    Args:

        pairwise_distances: Pairwise distances between backbone atoms in the region of interest. (n_frames = n_samples, n_paired_residues * 4)

    """
    assert pairwise_distances.ndim == 2
    # TM-score usually only uses C-alpha atoms
    # Assumes references and samples have backbone atoms in the following order (N, CA, C, O)
    ca_pairwise_distances = pairwise_distances[:, 1::4]  # [n_frames, n_residues]
    ltarget = ca_pairwise_distances.shape[1]
    d0 = 1.24 * np.cbrt(ltarget - 15) - 1.8
    summand = 1 / (1 + (ca_pairwise_distances / d0) ** 2)
    return 1 / ltarget * summand.sum(axis=-1)  # [n_frames]


def dssp_match_accuracy(
    traj_i: mdtraj.Trajectory,
    traj_j: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]] | None = None,
) -> np.ndarray:
    """Secondary structure agreement between structures in `traj_i` and `traj_j`

    Args:
        traj_i: First trajectory.
        traj_j: Second trajectory.
        matching_resids (optional): Matching resSeq indices between the trajectories
    """
    assert (
        traj_i.n_frames == traj_j.n_frames
        or (traj_i.n_frames == 1 and traj_j.n_frames > 1)
        or (traj_i.n_frames > 1 and traj_j.n_frames == 1)
    ), "at least one trajectory must have either only one frame (i.e a single reference conformation) or both need to have the same number of frames"

    def _get_zero_resid_map(traj: mdtraj.Trajectory) -> dict[int, int]:
        """mdtraj returns dssp in zero-indexed arrays"""
        zero_indexed_resid: dict[int, int] = {}
        for idx, residue in enumerate(traj.top.residues):
            zero_indexed_resid[residue.resSeq] = idx
        return zero_indexed_resid

    dssp_i, dssp_j = mdtraj.compute_dssp(traj_i), mdtraj.compute_dssp(traj_j)

    if matching_resids is None:
        matching_resids = seq_pairwise_align_trajs(traj_i[0], traj_j[0])

    zero_resid_map_i, zero_resid_map_j = (
        _get_zero_resid_map(traj_i[0]),
        _get_zero_resid_map(traj_j[0]),
    )
    zero_resid_i = [zero_resid_map_i[resid] for resid in [r[0] for r in matching_resids]]
    zero_resid_j = [zero_resid_map_j[resid] for resid in [r[1] for r in matching_resids]]

    dssp_i, dssp_j = dssp_i[:, zero_resid_i], dssp_j[:, zero_resid_j]
    return np.mean(dssp_i == dssp_j, axis=1)


def lddt(
    traj_i: mdtraj.Trajectory,
    traj_j: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]] | None = None,
    inclusion_radius: float = 15.0,
    thresholds: list[float] | None = None,
) -> np.ndarray:
    """
    Local distance difference test between two trajectories, as reported in Mariani et al. (2013)

    Args:
        traj_i: First trajectory. Note: this trajectory is used as reference to choose what gets
                atoms get included via `inclusion_radius`.
        traj_j: Second trajectory.
        matching_resids: Matching resSeq indices between the trajectories
        inclusion_radius: (in angstroms) chooses the local cutoff radius in which distances
                          are considered
        thresholds: (in angstroms) determines the difference thresholds over which the proportion of matched
                    distances are averaged. Defaults to [0.5, 1.0, 2.0, 4.0]
    """

    assert (
        traj_i.n_frames == 1 and traj_j.n_frames >= 1
    ) or (), "traj_i must have one frame only, and traj_j must have at least one frame"

    def _get_atom_indices_same_residue(traj: mdtraj.Trajectory) -> set[tuple[int, int]]:
        """Gets pairs of atom indices belonging to the same residue in `traj`"""
        atom_indices_same_residue: list[tuple[int, int]] = []

        for residue in traj.topology.residues:
            atom_indices = [atom.index for atom in residue.atoms]
            atom_indices_same_residue.extend(list(combinations(atom_indices, 2)))
        return set(atom_indices_same_residue)

    if matching_resids is None:
        matching_resids = seq_pairwise_align_trajs(traj_i, traj_j)

    if thresholds is None:
        thresholds = [0.5, 1.0, 2.0, 4.0]

    thresholds = np.array(thresholds)
    sel_i = fmt_resseq_selection([m[0] for m in matching_resids])
    sel_j = fmt_resseq_selection([m[1] for m in matching_resids])

    atomsel_i = set(traj_i.topology.select(sel_i))
    atomsel_j = set(traj_j.topology.select(sel_j))

    # Make pairs of atoms belonging to same residue
    atom_pairs_i = list(combinations(range(traj_i.n_atoms), 2))
    atom_pairs_j = list(combinations(range(traj_j.n_atoms), 2))

    same_residue_pairs_i = _get_atom_indices_same_residue(traj_i)
    same_residue_pairs_j = _get_atom_indices_same_residue(traj_j)

    # Exclude atom pairs from pair distances if they belong to same residue
    valid_pairs_i = [
        pair
        for pair in atom_pairs_i
        if (pair not in same_residue_pairs_i) and all([p in atomsel_i for p in pair])
    ]
    valid_pairs_j = [
        pair
        for pair in atom_pairs_j
        if (pair not in same_residue_pairs_j) and all([p in atomsel_j for p in pair])
    ]

    # Compute distances, *10 for conversion between nm and angstroms
    dist_i, dist_j = (
        mdtraj.compute_distances(traj_i, atom_pairs=valid_pairs_i) * 10,
        mdtraj.compute_distances(traj_j, atom_pairs=valid_pairs_j) * 10,
    )
    # Exclude distances above the inclusion radius
    valid_dist_indices = (dist_i <= inclusion_radius)[0, :]
    dist_i, dist_j = dist_i[:, valid_dist_indices], dist_j[:, valid_dist_indices]

    # Compute LDDT at different thresholds
    assert thresholds is not None  # shup up mypy
    dist_diff = np.repeat(np.abs(dist_i - dist_j)[np.newaxis, :], repeats=len(thresholds), axis=0)
    return (dist_diff < thresholds[:, np.newaxis, np.newaxis]).mean(axis=-1).mean(axis=0)  # type: ignore


def compute_contacts(
    traj_i: mdtraj.Trajectory,
    traj_j: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]] | None = None,
    reference_resid_pairs: list[tuple[int, int]] | None = None,
    threshold: float = 8.0,
    exclude_n_neighbours: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes contact maps for two trajectories over a set of common residues specified by `matching_resids`.

    Args:
        traj_i: First trajectory. Note: this trajectory is used as reference to determine which residues
                are sequence neighbours.
        traj_j: Second trajectory.
        matching_resids: Matching resSeq indices between the trajectories
        reference_resid_pairs: If passed, contacts will only be computed between pairs of specified
                            residue resSeqs.
        threshold: Distance threshold (in angstroms) to define contacts.
        exclude_n_neighbours: If set to > 0, contact maps will be excluded between residues whose residue index
                            differences is fewer or equal to this number.

    Returns
        * `contacts_i`, `contacts_j`: Arrays with contact maps
        * A list of pairs of residue indices corresponding to `traj_i`
          over which the contact map was computed
        * A list of pairs of zero-based residue indices corresponding to `traj_i`
          over which the contact map was computed
    """
    if matching_resids is None:
        matching_resids = seq_pairwise_align_trajs(traj_i, traj_j)

    # Get zero-indexed residue indices from matching_resids
    def _get_zero_index_resid_map(traj: mdtraj.Trajectory, resids: list[int]) -> dict[int, int]:
        resids_set = set(resids)
        resids_to_zeroidx_resid: dict[int, int] = {}
        for residue in traj.topology.residues:
            if residue.resSeq in resids_set:
                resids_to_zeroidx_resid[residue.resSeq] = residue.index
        return resids_to_zeroidx_resid

    resids_i, resids_j = [m[0] for m in matching_resids], [m[1] for m in matching_resids]
    resid_i_to_j = {k: v for k, v in matching_resids}

    # Exclude neighbouring residues using `traj_i` as a reference.
    valid_resid_combs_i = [
        (r, l) for r, l in list(combinations(resids_i, 2)) if abs(r - l) >= exclude_n_neighbours
    ]

    if reference_resid_pairs is not None:
        reference_resid_pairs_set = set(reference_resid_pairs)
        valid_resid_combs_i = [
            (r, l)
            for r, l in valid_resid_combs_i
            if (r, l) in reference_resid_pairs_set or (l, r) in reference_resid_pairs_set
        ]

    valid_resid_combs_i = np.array(valid_resid_combs_i)
    valid_resid_combs_j = np.array(
        [(resid_i_to_j[r], resid_i_to_j[l]) for r, l in valid_resid_combs_i]
    )

    # Translate resids to zero-based resids to use `compute_contacts`
    resid_to_zeroidx_resid_i = _get_zero_index_resid_map(traj_i, resids_i)
    resid_to_zeroidx_resid_j = _get_zero_index_resid_map(traj_j, resids_j)

    valid_zeroidx_resid_combs_i = np.array(
        [
            (resid_to_zeroidx_resid_i[resid_r], resid_to_zeroidx_resid_i[resid_l])
            for resid_r, resid_l in valid_resid_combs_i
        ]
    )
    valid_zeroidx_resid_combs_j = np.array(
        [
            (resid_to_zeroidx_resid_j[resid_r], resid_to_zeroidx_resid_j[resid_l])
            for resid_r, resid_l in valid_resid_combs_j
        ]
    )

    distances_i, _ = mdtraj.compute_contacts(
        traj_i, scheme="ca", contacts=valid_zeroidx_resid_combs_i
    )
    distances_j, _ = mdtraj.compute_contacts(
        traj_j, scheme="ca", contacts=valid_zeroidx_resid_combs_j
    )

    # Convert to angstrom
    distances_i *= 10.0
    distances_j *= 10.0

    contacts_i, contacts_j = (
        (distances_i < threshold).astype(int),
        (distances_j < threshold).astype(int),
    )
    return contacts_i, contacts_j, valid_resid_combs_i, valid_zeroidx_resid_combs_i


def fraction_native_contacts(
    traj_i: mdtraj.Trajectory,
    traj_j: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]] | None = None,
    reference_resid_pairs: list[tuple[int, int]] | None = None,
    threshold: float = 8.0,
    exclude_n_neighbours: int = 0,
) -> np.ndarray:
    """Computes fraction of native contacts between a reference trajectory `traj_i` and another one `traj_j`,
    over a region of matching residues `matching_resids`. Only positive contacts in `traj_i` are used for
    comparison.

    Args:
        traj_i: Reference trajectory
        traj_j: Trajectory to compare against
        matching_resids: Matching resSeq indices between the trajectories.
        reference_resid_pairs: If passed, contacts will only be computed between pairs of specified
                            residue resSeqs
        threshold: Distance threshold (in Angstrom) to consider two residues as contacting.
        exclude_n_neighbours: If set to > 0, contact maps will be excluded between residues whose residue index
                            differences is fewer or equal to this number.
    """
    contacts_i, contacts_j, _, _ = compute_contacts(
        traj_i,
        traj_j,
        matching_resids=matching_resids,
        reference_resid_pairs=reference_resid_pairs,
        threshold=threshold,
        exclude_n_neighbours=exclude_n_neighbours,
    )
    native_contact_indices = np.where(contacts_i[0, :] == 1)[0]
    return np.mean(
        contacts_i[:, native_contact_indices] == contacts_j[:, native_contact_indices], axis=1
    )


def contact_distances(
    traj_i: mdtraj.Trajectory,
    traj_j: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]] | None = None,
    threshold: float = 8.0,
) -> np.ndarray:
    """
    Measures differences over 'hard' contact maps, as defined by a distance threshold.
    Range: [0, +inf), where 0 represents a perfect match between `traj_i` and `traj_j` over `matching_resids`.
    Neighbouring residue indices (e.g., `idx + 1, idx + 2`) are discarded from the computation.

    Args:
        traj_i: First trajectory. Note: this trajectory is used as reference to determine which residues
                are sequence neighbours.
        traj_j: Second trajectory.
        matching_resids: Matching resSeq indices between the trajectories.
        threshold: Distance threshold (in angstroms) to define contacts.
    """
    contacts_i, contacts_j, _, valid_zeroidx_resid_combs_i = compute_contacts(
        traj_i, traj_j, matching_resids=matching_resids, threshold=threshold
    )
    unmatched_contacts = np.abs(contacts_i - contacts_j)

    # Get the total number of mismatched contacts, per residue.
    sum_per_resid: dict[int, np.ndarray] = defaultdict(
        lambda: np.zeros(unmatched_contacts.shape[0])
    )
    for (resid_r, resid_l), mc in zip(valid_zeroidx_resid_combs_i, unmatched_contacts.T):
        sum_per_resid[resid_r] += mc
        sum_per_resid[resid_l] += mc

    sum_per_resid_arr: np.ndarray = np.vstack(list(sum_per_resid.values()))

    return sum_per_resid_arr.mean(axis=0)
