from typing import Any

import numpy as np
import pandas as pd

from bioemu_benchmarks.logger import get_logger

LOGGER = get_logger(__name__)

K_BOLTZMANN = 0.001987203599772605  # Boltzmann constant in kcal / mol / K


def _foldedness_from_fnc(fnc: np.ndarray, p_fold_thr: float, steepness: float) -> np.ndarray:
    """
    Compute foldedness from fraction of native contacts (FNC).

    Args:
        fnc: Fraction of native contacts.
        p_fold_thr: FNC that has foldedness 0.5.
        steepness: Steepness of the sigmoid function.

    Returns:
        Foldedness values.
    """
    return 1 / (1 + np.exp(-2 * steepness * (fnc - p_fold_thr)))


def _compute_dG(fnc: np.ndarray, temperature: float, p_fold_thr: float, steepness: float) -> float:
    """Compute dG from sigmoid of fraction of native contacts"""
    p_fold_from_stat = _foldedness_from_fnc(fnc, p_fold_thr=p_fold_thr, steepness=steepness).mean()
    p_fold_from_stat = np.clip(p_fold_from_stat, 1e-10, 1 - 1e-10)

    ratio = p_fold_from_stat / (1 - p_fold_from_stat)
    ratio = np.clip(ratio, 1e-10, 1e10)

    dG = -np.log(ratio) * K_BOLTZMANN * temperature  # default temperature 295 K
    return dG


def compute_dg_ddg_from_fnc(
    *,
    dict_fnc: dict[str, np.ndarray],
    system_info: pd.DataFrame,
    temperature: float = 295.0,
    p_fold_thr: float = 0.5,
    steepness: float = 10.0,
) -> pd.DataFrame:
    """
    Compute dG and ddG for a collection of systems based on their native contact scores.

    Args:
        dict_fnc: Dictionary with arrays containing fraction of native contact values for the
          different systems as entries and test case IDs as keys.
        system_info: Data frame containing benchmark information.
        temperature: Temperature used for free energy computation in Kelvin. Default: 295 K.
        p_fold_thr: Threshold for foldedness, i.e., FNC value that corresponds to foldedness 0.5,
            used as an offset in the sigmoid function. Default: 0.5.
        steepness: Steepness of the sigmoid function used to compute foldedness from FNC values.
            Default: 10.0.
    Returns:
        Data frame containing computed dG and ddGs, as well as experimental reference values and
        uncertainties.
    """
    free_energy_results: dict[str, dict[str, Any]] = {}

    # Compute corresponding dG from generated samples.
    for test_case in dict_fnc:
        # Load fraction of native contacts.
        fnc = dict_fnc[test_case]

        # Extract basic system info from benchmark definitions.
        sequence_dataframe = system_info[system_info.name == test_case]
        free_energy_results[test_case] = sequence_dataframe.to_dict(orient="records")[0]

        # Store temperature.
        free_energy_results[test_case]["temperature"] = temperature

        # Store number of samples and warn if below target value.
        num_samples_target = free_energy_results[test_case]["num_samples"]
        num_samples = len(fnc)
        if num_samples < 0.7 * num_samples_target:
            LOGGER.warning(
                f"Number of samples for {test_case} below recommendation "
                f"({num_samples}/{num_samples_target})."
            )
        free_energy_results[test_case]["num_samples"] = num_samples

        # Compute dG and store.
        dg_pred = _compute_dG(
            fnc, temperature=temperature, p_fold_thr=p_fold_thr, steepness=steepness
        )
        free_energy_results[test_case]["dg_pred"] = dg_pred

    # Compute ddGs from dGs.
    for test_case in free_energy_results:
        test_case_wt = free_energy_results[test_case]["name_wt"]

        # Do not compute ddGs within the same system.
        if test_case == test_case_wt:
            continue

        if test_case_wt not in free_energy_results:
            LOGGER.warning(f"Could not find wild type results for {test_case_wt} for ddG")
            continue

        free_energy_results[test_case]["ddg_pred"] = (
            free_energy_results[test_case]["dg_pred"] - free_energy_results[test_case_wt]["dg_pred"]
        )

    return pd.DataFrame(free_energy_results).T
