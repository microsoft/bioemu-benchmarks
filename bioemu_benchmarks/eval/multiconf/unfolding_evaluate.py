import os
from glob import glob
from itertools import combinations, product

import mdtraj
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

from bioemu_benchmarks.eval.multiconf.align import seq_pairwise_align_trajs
from bioemu_benchmarks.eval.multiconf.evaluate import LocalResidInfo, MetricType, TestCaseResult
from bioemu_benchmarks.eval.multiconf.metrics import fraction_native_contacts
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.samples import IndexedSamples, SequenceSample
from bioemu_benchmarks.utils import StrPath, tqdm_joblib

LOGGER = get_logger(__name__)


def evaluate_singleconf_unfolding_test_case(
    test_case: str,
    references_dir: StrPath,
    sequence_samples: list[SequenceSample],
    references_localresidinfo_dir: StrPath,
    filter_backbone: bool = True,
) -> TestCaseResult:
    """Evaluates a single test case. Args are the same as `evaluate_singleconf_unfolding`."""

    ref_pdbs = sorted(glob(os.path.join(references_dir, test_case, "*.pdb")))

    assert len(sequence_samples) == len(ref_pdbs) == 1

    sequence_sample = sequence_samples[0]
    ref_traj = mdtraj.load_pdb(ref_pdbs[0])

    sample_traj = sequence_sample.get_traj()
    residinfo_json = os.path.join(references_localresidinfo_dir, f"{test_case}.json")
    localresidinfo = LocalResidInfo.from_json(
        residinfo_json, n_residues=max([r.resSeq for r in ref_traj.topology.residues])
    )

    if filter_backbone:
        sample_traj = sample_traj.atom_slice(sample_traj.topology.select("backbone"))

    matching_resids = seq_pairwise_align_trajs(ref_traj, sample_traj[0])

    metrics: dict[MetricType, np.ndarray] = {}

    # Get folding resids
    folding_resids = set(localresidinfo.metric_resids)

    # We want to compute the FNC within the folding region and between the folding region and anything else
    # in the protein
    matching_resids_folding = [(i, j) for i, j in matching_resids if i in folding_resids]
    matching_resids_rest = [(i, j) for i, j in matching_resids if i not in folding_resids]

    reference_resid_pairs_folding = list(combinations([m[0] for m in matching_resids_folding], 2))
    reference_resid_pairs_folding_v_rest = list(
        product(
            [m[0] for m in matching_resids_folding],
            [m[0] for m in matching_resids_rest],
        )
    )

    reference_resid_pairs: list[tuple[int, int]] = (
        reference_resid_pairs_folding + reference_resid_pairs_folding_v_rest
    )
    assert reference_resid_pairs

    metrics[MetricType.FNC_UNFOLD_F] = fraction_native_contacts(
        ref_traj,
        sample_traj,
        matching_resids=matching_resids,
        reference_resid_pairs=reference_resid_pairs,
        exclude_n_neighbours=3,
    )
    metrics[MetricType.FNC_UNFOLD_U] = metrics[MetricType.FNC_UNFOLD_F]
    return TestCaseResult(
        test_case=test_case,
        metrics_against_references=metrics,
        references_names=None,
        metrics_between_references=None,
        closest_samples=None,
        topology_ids=None,
    )


def evaluate_singleconf_unfolding(
    indexed_samples: IndexedSamples,
    references_dir: StrPath,
    references_localresidinfo_dir: StrPath,
    n_jobs: int = -1,
) -> dict[str, TestCaseResult]:
    """
    Computes local unfolding metrics for a benchmark

    Args:
        indexed_samples: An `IndexedSamples` instance containing samples to evaluate.
        references_dir: Parent reference directory containing a bunch of PDB files representing
                        static references.
        references_localresidinfo_dir: Parent reference directory with local unfolding and contact region information
                        in JSON format.
        n_jobs: Number of processes used for evaluation

    Returns:
        A dictionary, indexed by test case, contaning `TestCaseResult` instances for
        each of them
    """

    with tqdm_joblib(
        tqdm(
            desc="Now computing metrics...", total=len(indexed_samples.test_case_to_sequencesamples)
        )
    ) as _pbar:
        evaluations: list[TestCaseResult] = Parallel(n_jobs=n_jobs)(
            delayed(evaluate_singleconf_unfolding_test_case)(
                test_case,
                references_dir=references_dir,
                sequence_samples=sequence_samples,
                references_localresidinfo_dir=references_localresidinfo_dir,
            )
            for test_case, sequence_samples in indexed_samples.test_case_to_sequencesamples.items()
        )
    return {ev.test_case: ev for ev in evaluations}
