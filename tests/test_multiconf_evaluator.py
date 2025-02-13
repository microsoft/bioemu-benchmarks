import os
import tempfile

import h5py
import numpy as np
import pytest

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.eval.multiconf.evaluate import MULTICONF_METRIC_TYPES, MetricType
from bioemu_benchmarks.eval.multiconf.results import MulticonfResults
from bioemu_benchmarks.evaluator_utils import evaluator_from_benchmark
from bioemu_benchmarks.samples import IndexedSamples, find_samples_in_dir

from . import TEST_DATA_DIR

_EXPECTED_test_case_KRECALL = ("E1C7U0", MetricType.RMSD, [6.0333076, 4.7683716e-07])
_EXPECTED_EVALUATOR_RESULTS: dict[Benchmark, list[str]] = {
    Benchmark.MULTICONF_OOD60: [
        "multiconf_ood60_rmsd_free_energy.png",
        "multiconf_ood60_contact-distance_coverage.png",
        "multiconf_ood60_dssp_acc_free_energy.png",
        "multiconf_ood60_dssp_acc_coverage.png",
        "multiconf_ood60_rmsd_coverage.png",
        "multiconf_ood60_tm-score_free_energy.png",
        "multiconf_ood60_lddt_free_energy.png",
        "multiconf_ood60_lddt_coverage.png",
        "closest",
        "multiconf_results.h5",
        "multiconf_ood60_tm-score_coverage.png",
        "multiconf_ood60_contact-distance_free_energy.png",
    ],
    Benchmark.SINGLECONF_LOCALUNFOLDING: [
        "singleconf_localunfolding_fnc_unfold_f_free_energy.png",
        "multiconf_results.h5",
        "singleconf_localunfolding_fnc_unfold_f_coverage.png",
        "singleconf_localunfolding_fnc_unfold_u_coverage.png",
        "singleconf_localunfolding_fnc_unfold_u_free_energy.png",
    ],
}
_EXPECTED_H5_KEYS: dict[Benchmark, list[str]] = {
    Benchmark.MULTICONF_OOD60: [
        "coverage_multiconf_ood60_contact-distance",
        "coverage_multiconf_ood60_dssp_acc",
        "coverage_multiconf_ood60_lddt",
        "coverage_multiconf_ood60_rmsd",
        "coverage_multiconf_ood60_tm-score",
        "krecall_multiconf_ood60_contact-distance",
        "krecall_multiconf_ood60_dssp_acc",
        "krecall_multiconf_ood60_lddt",
        "krecall_multiconf_ood60_rmsd",
        "krecall_multiconf_ood60_tm-score",
        "reference_metrics",
        "sample_metrics",
    ],
    Benchmark.SINGLECONF_LOCALUNFOLDING: [
        "coverage_singleconf_localunfolding_fnc_unfold_f",
        "coverage_singleconf_localunfolding_fnc_unfold_u",
        "krecall_singleconf_localunfolding_fnc_unfold_f",
        "krecall_singleconf_localunfolding_fnc_unfold_u",
        "sample_metrics",
    ],
}


@pytest.mark.parametrize(
    "benchmark", [Benchmark.MULTICONF_OOD60, Benchmark.SINGLECONF_LOCALUNFOLDING]
)
def test_multiconf_evaluator(benchmark: Benchmark):
    """
    Tests that the multiconf evaluator dumps what it's supposed to in a
    specified results directory
    """
    _test_samples_basedir = {
        Benchmark.MULTICONF_OOD60: "multiconf_ood60",
        Benchmark.SINGLECONF_LOCALUNFOLDING: "singleconf_localunfolding",
    }
    kwargs = (
        dict(metric_types=MULTICONF_METRIC_TYPES) if benchmark == Benchmark.MULTICONF_OOD60 else {}
    )

    evaluator = evaluator_from_benchmark(benchmark, **kwargs)

    sequence_samples = find_samples_in_dir(
        os.path.join(TEST_DATA_DIR, "samples_example", _test_samples_basedir[benchmark])
    )

    indexed_samples = IndexedSamples.from_benchmark(
        benchmark=benchmark, sequence_samples=sequence_samples
    )

    with tempfile.TemporaryDirectory() as tempdir:
        evalresults = evaluator(indexed_samples)
        assert isinstance(evalresults, MulticonfResults)
        test_case, metric_type, expected_krecall = _EXPECTED_test_case_KRECALL
        if benchmark == Benchmark.MULTICONF_OOD60:
            np.testing.assert_allclose(
                list(evalresults.krecall[benchmark.value][metric_type][test_case]),
                expected_krecall,
                err_msg="Recall calculation changed!",
            )
            np.testing.assert_allclose(
                [evalresults.coverage[benchmark.value][metric_type][1][0][-1]],
                [0.8157894736842105],
                err_msg="Coverage calculation changed!",
            )

        # Test serialisation of evaluator / results class
        h5_file = os.path.join(tempdir, "multiconf_results.h5")
        evalresults.save_to_h5(h5_file)
        evalresults.plot(tempdir)
        if benchmark == Benchmark.MULTICONF_OOD60:
            evalresults.save_closest_samples(tempdir)

        for eo in _EXPECTED_EVALUATOR_RESULTS[benchmark]:
            assert os.path.exists(
                os.path.join(tempdir, eo)
            ), f"Missing expected output file/dir: {eo}"

        # Test H5 keys consistency
        h5file = h5py.File(h5_file, mode="r")
        for expected_key in _EXPECTED_H5_KEYS[benchmark]:
            assert (
                expected_key in h5file
            ), f"Expected result key {expected_key} not found in results H5"
        # Test serialization of evaluator results
        pkl_file = os.path.join(tempdir, "evalresults.pkl")
        evalresults.to_pickle(pkl_file)
        _new_evalresults = MulticonfResults.from_pickle(pkl_file)
