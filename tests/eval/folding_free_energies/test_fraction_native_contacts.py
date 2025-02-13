import mdtraj
import numpy as np

from bioemu_benchmarks.eval.folding_free_energies.fraction_native_contacts import (
    FNCSettings,
    _compute_contact_score,
    _compute_reference_contacts,
    _get_sequence_index_map,
    get_fnc_from_samples_trajectory,
)


def test_compute_reference_contacts(reference_conformation_mutant, reference_contacts_mutant):
    """Test determination of contact indices and distances."""
    pair_idx, pair_dist = _compute_reference_contacts(
        reference_conformation=reference_conformation_mutant,
        sequence_separation=FNCSettings.sequence_separation,
        contact_cutoff=FNCSettings.contact_cutoff,
    )

    # Compare to precomputed reference contacts.
    np.testing.assert_allclose(pair_idx, reference_contacts_mutant[0])
    np.testing.assert_allclose(pair_dist, reference_contacts_mutant[1])

    # Make sure heavy element selection works properly.
    def check_indices(index_slice):
        unique_idx = np.unique(index_slice)
        sliced_top = reference_conformation_mutant.atom_slice(unique_idx).top
        atomic_numbers = np.array([a.element.number for a in sliced_top.atoms])
        assert np.all(atomic_numbers > 1)

    check_indices(pair_idx[:, 0])
    check_indices(pair_idx[:, 1])


def test_get_sequence_index_map(fnc_test_data_wt, fnc_test_data_mutant):
    """Test generation of index sequence map."""
    sequence_index_map = _get_sequence_index_map(
        fnc_test_data_wt.sequence, fnc_test_data_mutant.sequence
    )
    # Convert to array for easier comparison.
    sequence_index_map = np.array(sequence_index_map)

    # Target map for this wild type and mutant should be sequential indices with -1 at the points
    # where sequence mismatches.
    target_map = np.arange(len(sequence_index_map))
    mask = np.array(
        [a != b for a, b in zip(fnc_test_data_wt.sequence, fnc_test_data_mutant.sequence)]
    )
    target_map[mask] = -1

    np.testing.assert_allclose(target_map, sequence_index_map)


def test_compute_contact_score(samples_traj_wt, reference_conformation_wt, fnc_test_data_wt):
    """Test computation of contact score from distances."""
    # Filter for CA atoms in reference contacts.
    reference_conformation_wt_ca = reference_conformation_wt.atom_slice(
        reference_conformation_wt.top.select("name CA")
    )
    reference_contact_indices, reference_contact_distances = _compute_reference_contacts(
        reference_conformation=reference_conformation_wt_ca,
        sequence_separation=FNCSettings.sequence_separation,
        contact_cutoff=FNCSettings.contact_cutoff,
    )

    # Select CA atoms in samples.
    ca_indices_samples = samples_traj_wt.top.select("name CA")
    samples_ca = samples_traj_wt.atom_slice(ca_indices_samples)

    # Compute sample distances and convert to Angstrom.
    samples_contact_distances = mdtraj.geometry.compute_distances(
        samples_ca, reference_contact_indices, periodic=False
    )
    samples_contact_distances = mdtraj.utils.in_units_of(
        samples_contact_distances, "nanometers", "angstrom"
    )

    # Compute contact score.
    contact_score = _compute_contact_score(
        samples_contact_distances=samples_contact_distances,
        reference_contact_distances=reference_contact_distances,
        contact_beta=FNCSettings.contact_beta,
        contact_lambda=FNCSettings.contact_lambda,
        contact_delta=FNCSettings.contact_delta,
    )

    # Check consistency.
    np.testing.assert_allclose(contact_score, fnc_test_data_wt.fnc)


def test_get_fnc_from_samples_trajectory(
    samples_traj_wt, reference_conformation_wt, fnc_test_data_wt
):
    """Test computation of contact score from trajecties."""
    contact_score = get_fnc_from_samples_trajectory(
        samples_traj_wt,
        reference_conformation_wt,
        sequence_separation=FNCSettings.sequence_separation,
        contact_cutoff=FNCSettings.contact_cutoff,
        contact_beta=FNCSettings.contact_beta,
        contact_lambda=FNCSettings.contact_lambda,
        contact_delta=FNCSettings.contact_delta,
    )
    # Check consistency.
    np.testing.assert_allclose(contact_score, fnc_test_data_wt.fnc)
