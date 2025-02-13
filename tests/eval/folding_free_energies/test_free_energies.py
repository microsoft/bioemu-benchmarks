import numpy as np

from bioemu_benchmarks.eval.folding_free_energies.free_energies import (
    _compute_dG,
    _compute_threshold,
    compute_dg_ddg_from_fnc,
)


def test_compute_dG(fnc_test_data_parameterized):
    """Test dG computation."""
    dG = _compute_dG(
        fnc_test_data_parameterized.fnc,
        threshold=fnc_test_data_parameterized.threshold,
        temperature=fnc_test_data_parameterized.temperature,
    )
    np.testing.assert_allclose(dG, fnc_test_data_parameterized.target_dg)


def test_compute_threshold(fnc_test_data_parameterized):
    """Test threshold computation."""
    threshold = _compute_threshold(fnc_test_data_parameterized.fnc)
    np.testing.assert_allclose(threshold, fnc_test_data_parameterized.threshold)


def test_compute_dg_ddg_from_fnc(fnc_test_data_wt, fnc_test_data_mutant, system_info):
    """Test computation of multiple dGs and ddGs."""
    dict_fnc = {
        fnc_test_data_wt.test_case: fnc_test_data_wt.fnc,
        fnc_test_data_mutant.test_case: fnc_test_data_mutant.fnc,
    }

    results_df = compute_dg_ddg_from_fnc(
        dict_fnc=dict_fnc,
        system_info=system_info,
        fixed_threshold=None,
        temperature=fnc_test_data_wt.temperature,
    )

    assert len(results_df) == 2

    # Check results for wild type (ddG should be NaN).
    results_wt = results_df[results_df.sequence == fnc_test_data_wt.sequence]
    np.testing.assert_allclose(results_wt.dg_pred.item(), fnc_test_data_wt.target_dg)
    assert np.isnan(results_wt.ddg_pred.item())
    np.testing.assert_allclose(results_wt.threshold.item(), fnc_test_data_wt.threshold)

    # Check results for mutant.
    results_mutant = results_df[results_df.sequence == fnc_test_data_mutant.sequence]
    np.testing.assert_allclose(results_mutant.dg_pred.item(), fnc_test_data_mutant.target_dg)
    np.testing.assert_allclose(results_mutant.ddg_pred.item(), fnc_test_data_mutant.target_ddg)
    np.testing.assert_allclose(results_mutant.threshold.item(), fnc_test_data_mutant.threshold)
