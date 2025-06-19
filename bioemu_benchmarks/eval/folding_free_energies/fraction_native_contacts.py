from dataclasses import dataclass
from itertools import combinations

import mdtraj
import numpy as np
from Bio import pairwise2
from scipy.special import expit

from bioemu_benchmarks.logger import get_logger

LOGGER = get_logger(__name__)


@dataclass
class FNCSettings:
    """
    Data class for collecting fraction of native contact settings.

    Attributes:
        sequence_separation: Minimum separation of residues for which contacts are computed.
        contact_cutoff: Maximum cutoff distance (in Angstrom) used for contact computation.
        contact_beta: Scaling factor for contact score computations (see `_compute_contact_score`).
        contact_delta: Offset for contact score computations (see `_compute_contact_score`).
        contact_lambda: Scaling off offset in contact score computation (see
          `_compute_contact_score`).
    """

    sequence_separation: int = 3
    contact_cutoff: float = 10.0
    contact_beta: float = 5.0
    contact_delta: float = 0.0
    contact_lambda: float = 1.2


def _compute_reference_contacts(
    *, reference_conformation: mdtraj.Trajectory, sequence_separation: int, contact_cutoff: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute native contacts in reference conformation based on given parameters. Uses `mdtraj`
    for atom selection and distance computations. Contacts are computed between heavy atoms only.

    Args:
        reference_conformation: Reference structure in mdtraj Trajectory format.
        sequence_separation: Minimum separation between residues in sequence to be considered for
          contact computation.
        contact_cutoff: Distance cutoff. Only contacts between atoms below this cutoff are
          considered.

    Returns:
        Tuple of arrays containing the contact indices and contact distances.
    """
    # Select only heavy atoms for contact computation
    heavy_atoms = reference_conformation.topology.select_atom_indices("heavy")

    heavy_pair_indices: list[tuple[int, int]] = []
    for i, j in combinations(heavy_atoms, 2):
        res_position_i = reference_conformation.topology.atom(i).residue.index
        res_position_j = reference_conformation.topology.atom(j).residue.index

        # Check if residues are sufficiently far apart in the sequence:
        if abs(res_position_i - res_position_j) > sequence_separation:
            heavy_pair_indices.append((i, j))

            # Make symmetric for per-residue resolution of contacts:
            heavy_pair_indices.append((j, i))

    # Compute the distances between the valid heavy pairs.
    heavy_pair_indices = np.array(heavy_pair_indices)
    heavy_pairs_distances = mdtraj.compute_distances(reference_conformation, heavy_pair_indices)[0]

    # Convert distances from nanometers to Angstrom
    heavy_pairs_distances = mdtraj.utils.in_units_of(
        heavy_pairs_distances, "nanometers", "angstrom"
    )

    # Filter according to cutoff
    heavy_pair_indices = heavy_pair_indices[heavy_pairs_distances <= contact_cutoff]
    heavy_pairs_distances = heavy_pairs_distances[heavy_pairs_distances <= contact_cutoff]

    return heavy_pair_indices, heavy_pairs_distances


def _get_aligned_indices(seq_alignment_1: str, seq_alignment_2: str) -> list[int]:
    """
    Compute the indices of the aligned residues in sequence 1 without gaps in the alignment.

    E.g. seq1=ABCDE, seq2=GABDF, then seq_alignment_1='-ABCDE-', seq_alignment_2='GAB-D-F'.
    The dashes are gaps and not counted when incrementing the index.

    Args:
        seq_alignment_1: First of the aligned sequences.
        seq_alignment_2: Second sequence.

    Returns:
        List of indices.
    """

    aligned_indices = []
    n = 0
    for i, s in enumerate(seq_alignment_1):
        if s != "-":
            if seq_alignment_2[i] != "-":
                aligned_indices.append(n)
            n += 1
    return aligned_indices


def _get_sequence_index_map(samples_sequence: str, reference_sequence: str) -> np.ndarray:
    """
    Perform pairwise alignment of reference and sample sequence to get the indices mapping reference
    atoms to sample atoms.

    Args:
        samples_sequence: Sequence of sample structures.
        reference_sequence: Reference sequence.

    Returns:
        Array mapping reference atom indices to sample atom indices. Atoms with no mapping are
        assigned -1.
    """
    # TODO: potential duplicate with `get_pairwise_align_traj `.
    alignments = pairwise2.align.globalxx(samples_sequence, reference_sequence)
    aligned_indices_sample = _get_aligned_indices(alignments[0].seqA, alignments[0].seqB)
    aligned_indices_ref = _get_aligned_indices(alignments[0].seqB, alignments[0].seqA)
    assert len(aligned_indices_sample) == len(aligned_indices_ref)

    ref_to_samples_map = np.full(len(reference_sequence), -1, dtype=np.int64)
    for ref, smp in zip(aligned_indices_ref, aligned_indices_sample):
        ref_to_samples_map[ref] = smp

    return ref_to_samples_map


def _compute_contact_score(
    *,
    samples_contact_distances: np.ndarray,
    reference_contact_distances: np.ndarray,
    contact_delta: float,
    contact_beta: float,
    contact_lambda: float,
) -> np.ndarray:
    """
    Compute contact scores for all pairs of contacts using the equation:

    .. math::

        q = \frac{1}{N_\mathrm{contacts}} \sum_{c}^{N_\mathrm{contacts}} \frac{1}{1 + \exp(-\beta(d_c - \lambda (d^\mathrm{ref}_c + \delta)))}

    Args:
        samples_contact_distances: Array of contact distances in samples with shape
          [num_samples x num_contacts].
        reference_contact_distances: Array of reference contact distances with shape [num_contacts].
        contact_delta: Offset for contact score determination.
        contact_beta: Scaling of exponential term.
        contact_lambda: Scaling for reference contact to account for fluctuations.

    Returns:
        Contact scores for all pairwise interactions.
    """

    q_ij = expit(
        -contact_beta
        * (
            samples_contact_distances
            - contact_lambda * (reference_contact_distances[None, :] + contact_delta)
        )
    )
    return np.mean(q_ij, axis=-1)


def get_fnc_from_samples_trajectory(
    samples: mdtraj.Trajectory,
    reference_conformation: mdtraj.Trajectory,
    sequence_separation: int = 3,
    contact_cutoff: float = 10.0,
    contact_beta: float = 5.0,
    contact_lambda: float = 1.2,
    contact_delta: float = 0.0,
) -> np.ndarray:
    """
    Compute fraction of native contact scores for samples.

    Args:
        samples: Samples in mdtraj Trajectory format.
        reference_conformation: Reference structure in mdtraj Trajectory format.
        sequence_separation: Minimum separation of residues for which contacts are computed.
        contact_cutoff: Maximum cutoff distance (in Angstrom) used for contact computation.
        contact_beta: Scaling factor for contact score computations (see `_compute_contact_score`).
        contact_delta: Offset for contact score computations (see `_compute_contact_score`).
        contact_lambda: Scaling off offset in contact score computation (see
          `_compute_contact_score`).

    Returns:
        Array containing fraction of native contacts score for each sample.
    """

    # Select CA atoms in reference.
    ca_indices = reference_conformation.top.select("name CA")
    reference_conformation_ca = reference_conformation.atom_slice(ca_indices)

    # Select CA atoms in samples.
    ca_indices_samples = samples.top.select("name CA")
    samples_ca = samples.atom_slice(ca_indices_samples)

    # Compute reference contacts.
    reference_contact_indices, reference_contact_distances = _compute_reference_contacts(
        reference_conformation=reference_conformation_ca,
        sequence_separation=sequence_separation,
        contact_cutoff=contact_cutoff,
    )

    # Perform sequence alignment for residue mapping.
    reference_sequence = reference_conformation.top.to_fasta()[0]
    samples_sequence = samples_ca.top.to_fasta()[0]
    ref_to_samples_map = _get_sequence_index_map(samples_sequence, reference_sequence)

    # Remove unmatched pairs with -1 in it
    reference_contact_distances = reference_contact_distances[
        np.all(ref_to_samples_map[reference_contact_indices] != -1, axis=1)
    ]
    reference_contact_indices = reference_contact_indices[
        np.all(ref_to_samples_map[reference_contact_indices] != -1, axis=1)
    ]

    # Map contact indices to current samples.
    aligned_contact_indices = ref_to_samples_map[reference_contact_indices]

    # Compute sample distances and convert to Angstrom.
    samples_contact_distances = mdtraj.geometry.compute_distances(
        samples_ca, aligned_contact_indices, periodic=False
    )
    samples_contact_distances = mdtraj.utils.in_units_of(
        samples_contact_distances, "nanometers", "angstrom"
    )

    # Compute contact score.
    contact_score = _compute_contact_score(
        samples_contact_distances=samples_contact_distances,
        reference_contact_distances=reference_contact_distances,
        contact_beta=contact_beta,
        contact_lambda=contact_lambda,
        contact_delta=contact_delta,
    )

    return contact_score
