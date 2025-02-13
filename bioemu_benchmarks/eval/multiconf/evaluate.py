import dataclasses
import json
import operator
import os
import tempfile
from dataclasses import dataclass
from enum import Enum
from glob import glob
from itertools import combinations
from typing import Callable

import mdtraj
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

from bioemu_benchmarks.eval.multiconf.align import (
    align_traj_to_ref_using_range,
    fmt_resseq_selection,
    seq_pairwise_align_trajs,
    setup_tm_align,
)
from bioemu_benchmarks.eval.multiconf.metrics import (
    contact_distances,
    dssp_match_accuracy,
    lddt,
    tm_score,
    tm_score_local,
)
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.paths import MULTICONF_ASSET_DIR
from bioemu_benchmarks.samples import IndexedSamples, SequenceSample
from bioemu_benchmarks.utils import StrPath, tqdm_joblib

LOGGER = get_logger(__name__)

NM2ANGS = 10


class MetricType(str, Enum):
    """Metrics supported by multiconf evaluation"""

    RMSD = "rmsd"
    TMSCORE = "tm-score"  # in the case of local metrics, alignment computed with Kabsch alg.
    CONTACT_DISTANCE = "contact-distance"  # per-residue average contact-map distance
    LDDT = "lddt"
    DSSP_ACC = "dssp_acc"  # secondary structure accuracy
    FNC_UNFOLD_U = "fnc_unfold_u"  # fraction of native contacts for evaluating unfoldedness
    FNC_UNFOLD_F = "fnc_unfold_f"  # fraction of native contacts for evaluating foldedness


METRIC_OPERATOR_BETTER: dict[MetricType, Callable] = {
    MetricType.RMSD: operator.lt,
    MetricType.TMSCORE: operator.gt,
    MetricType.LDDT: operator.gt,
    MetricType.CONTACT_DISTANCE: operator.lt,
    MetricType.DSSP_ACC: operator.gt,
    MetricType.FNC_UNFOLD_U: operator.lt,
    MetricType.FNC_UNFOLD_F: operator.gt,
}


@dataclass
class LocalResidInfo:
    alignment_resid_ranges: list[list[int | None]] | None  # ranges of residues to align upon
    metric_resid_ranges: list[list[int | None]]  # ranges of residues to compute metrics on
    n_residues: int  # maximum number of residue indices between all conformations of a single
    # test case id (or number of columns in their MSA)

    @classmethod
    def from_json(cls, json_file: StrPath, n_residues: int):
        """Parses alignment and metric residue information from JSON file"""
        with open(json_file) as json_handle:
            localresid_dict = json.load(json_handle)
        return cls(**localresid_dict, n_residues=n_residues)

    def parse_range(self, r: list[int | None]) -> tuple[int, int]:
        """
        Parses either an alignment or a metric range, in the format `[idx_i, idx_j]`, where
        `idx_x` can be `:`, representing _up to_, or _from_ in the Python slice sense, and
        returns explicit residue begin/end indices in the format `(begin_resid, end_resid)`.
        """
        assert len(r) == 2
        if r[0] is None:
            begin_resid = 1
        else:
            assert isinstance(r[0], int)  # shut up mypy
            begin_resid = r[0]

        if r[1] is None:
            end_resid = self.n_residues
        else:
            assert isinstance(r[1], int)  # shut up mypy
            end_resid = r[1]

        assert isinstance(begin_resid, int) and isinstance(end_resid, int)  # shut up mypy
        assert begin_resid <= end_resid
        return begin_resid, end_resid

    def to_explicit_resids(self, resid_ranges: list[list[int | None]]) -> list[int]:
        """Converts a list of resid ranges into an explicit residue index list"""
        explicit_resids: list[int] = []
        last_e_idx = -1
        for r in resid_ranges:
            b_idx, e_idx = self.parse_range(r)
            if e_idx > self.n_residues:
                explicit_resids.extend(list(range(b_idx, self.n_residues + 1)))
                return explicit_resids
            assert b_idx > last_e_idx  # ensure ascending order
            explicit_resids.extend(list(range(b_idx, e_idx + 1)))
            last_e_idx = e_idx
        return explicit_resids

    def __post_init__(self) -> None:
        self.alignment_resids: list[int] = []  # explicit alignment residue indices
        self.metric_resids: list[int] = []  # explicit metric residue indices

        if self.alignment_resid_ranges is not None:
            self.alignment_resids = self.to_explicit_resids(self.alignment_resid_ranges)

        self.metric_resids = self.to_explicit_resids(self.metric_resid_ranges)


@dataclass
class ClosestSample:
    """
    Keeps information about the closest sample to a given reference.

    Args:
        sequence_sample: XTC/PDB paths to sample file
        frame_idx: Frame index where closest sample in the XTC file is located
        metric_type: The type of metric that was computed
        metric_value: The associated numerical value of that metric against the reference
        reference_path: Path to the reference PDB file
    """

    sequence_sample: SequenceSample
    frame_idx: int
    metric_type: MetricType
    metric_value: float
    reference_path: StrPath

    def save_to_pdb(self, test_case: str, closest_dir: StrPath) -> None:
        dirpath = os.path.join(closest_dir, test_case, self.metric_type.value)

        os.makedirs(dirpath, exist_ok=True)
        filename = os.path.splitext(os.path.basename(self.reference_path))[0]
        outfile = os.path.join(dirpath, f"{filename}.pdb")

        traj = mdtraj.load(
            self.sequence_sample.trajectory_file, top=self.sequence_sample.topology_file
        )
        traj[self.frame_idx].save_pdb(outfile)


@dataclass
class TestCaseResult:
    """Results for a single multiconf test case.

    Args:
        test_case: Benchmark test case
        reference_names (optional): Names of the references the samples were evaluated
                        against
        metrics_between_references (optional): Metrics computed between references, indexed
                        by `MetricType`, then by the indices in `reference_names`
        metrics_against_references: Sample metrics computed against the references, indexed
                        by `MetricType`. Array shapes are [n_samples, n_references]
        closest_samples (optional): Closest sampled structure to each considered reference,
                        indexed by `MetricType`.
        topology_ids (optional): If different sampled topologies were used to evaluate the
                        same test case, this keeps track of which one was used for
                        which evaluation under the arrays used in
                        `metrics_against_references`.
    """

    test_case: str
    references_names: list[str] | None
    metrics_between_references: dict[MetricType, dict[tuple[int, int], float]] | None
    metrics_against_references: dict[MetricType, np.ndarray]
    closest_samples: (
        dict[
            MetricType,
            list[ClosestSample],
        ]
        | None
    )
    topology_ids: np.ndarray | None


def calc_metrics_global(
    ref_traj: mdtraj.Trajectory,
    sample_traj: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]],
    metric_types: list[MetricType],
    n_jobs: int = -1,
) -> dict[MetricType, np.ndarray]:
    """Computes additional metrics between reference and sample structures.
    Currently this returns TM-Score, LDDT and secondary structure agreement.
    All metrics, where appropriate, are in Angstrom units.

    Args:
        ref_traj: Reference trajectory instance
        sample_pdb: Sample trajectory instance
        matching_resids: List of matching residue ids between
                        `ref_traj` and `sample_traj`.
        metric_types: List of metrics to be computed.
        n_jobs: Number of workers used. Defaults to all available cores.


    Returns:
        A dictionary of metrics indexed by metric type.
    """

    sel_ref_traj = fmt_resseq_selection([m[0] for m in matching_resids])
    sel_sample_traj = fmt_resseq_selection([m[1] for m in matching_resids])

    metrics = {}

    ref_atomsel = ref_traj.topology.select(sel_ref_traj)
    sample_atomsel = sample_traj.topology.select(sel_sample_traj)

    assert len(ref_atomsel) == len(
        sample_atomsel
    ), f"Unequal atom selections, got {str(ref_atomsel)} and {str(sample_atomsel)}"

    # TODO: this code could be refactored so that we could iterate over `metric_types`
    # instead of checking whether the metrics are requested one by one
    # rmsd
    if MetricType.RMSD in metric_types:
        # Note: mdtraj uses this implementation to compute the minimum RMSD between two structures
        # https://theobald.brandeis.edu/qcp/
        # The RMSD calc does differ from a manual one if the structures are first aligned with another
        # software (e.g. pymol).
        # Note: While `matching_resids` takes care of aligning residues between reference and sample trajectories,
        # the calculation below does _not_ check that the reference and sample indices actually correspond
        # to the same atoms within residue. This needs to be guaranteed by cleaning `ref_traj` and `sample_traj`
        # accordingly. The current atom name ordering convention we're using within each residue is: [N, CA, C, O]
        metrics[MetricType.RMSD] = (
            mdtraj.rmsd(
                target=sample_traj,
                reference=ref_traj,
                atom_indices=sample_atomsel,
                ref_atom_indices=ref_atomsel,
                parallel=True,
            )
            * NM2ANGS
        )

    # tm score
    if MetricType.TMSCORE in metric_types:

        def _tm_score_wrapper_subprocess(
            sample_traj: mdtraj.Trajectory, ref_pdb: StrPath, frame_idx: int
        ):
            sample_pdb = os.path.join(tempdir, f"sample_{frame_idx}.pdb")
            sample_traj_one = sample_traj[frame_idx]
            sample_traj_one.save_pdb(sample_pdb)
            return tm_score(ref_pdb, sample_pdb)

        with tempfile.TemporaryDirectory() as tempdir:
            ref_pdb = os.path.join(tempdir, "ref.pdb")
            ref_traj.save_pdb(ref_pdb)
            metrics[MetricType.TMSCORE] = np.array(
                Parallel(n_jobs=n_jobs)(
                    delayed(_tm_score_wrapper_subprocess)(sample_traj, ref_pdb, frame_idx)
                    for frame_idx in range(sample_traj.n_frames)
                )
            )

    # dssp match accuracy
    if MetricType.DSSP_ACC in metric_types:
        metrics[MetricType.DSSP_ACC] = dssp_match_accuracy(
            ref_traj, sample_traj, matching_resids=matching_resids
        )

    # lddt
    if MetricType.LDDT in metric_types:
        metrics[MetricType.LDDT] = lddt(ref_traj, sample_traj, matching_resids=matching_resids)

    # contact-match
    if MetricType.CONTACT_DISTANCE in metric_types:
        metrics[MetricType.CONTACT_DISTANCE] = contact_distances(
            ref_traj, sample_traj, matching_resids=matching_resids
        )

    return metrics


def calc_metrics_local(
    ref_traj: mdtraj.Trajectory,
    sample_traj: mdtraj.Trajectory,
    matching_resids: list[tuple[int, int]],
    localresidinfo: LocalResidInfo,
    metric_types: list[MetricType],
) -> dict[MetricType, float]:
    """
    Computes metrics using a user-specified sequence range for superposition, and
    a potentially different user-defined sequence range to compute the metrics on.

    Args:
        ref_traj: Reference trajectory instance
        sample_pdb: Sample trajectory instance
        matching_resids: List of matching residue ids between
                                                 `ref_traj` and `sample_traj`.
        localresidinfo: instance containing sequence index ranges for
                                         alignment and metric calculation.
        metric_types: List of metrics to be computed


    Returns:
        A dictionary of metrics indexed by metric type
    """

    metrics = {}

    # Filter matching resids to user-defined ones for alignment
    # If no alignment residues have been defined, by default align on
    # everything
    matching_resids_alignment = (
        [m for m in matching_resids if m[0] in set(localresidinfo.alignment_resids)]
        if localresidinfo.alignment_resid_ranges is not None
        else matching_resids
    )

    sample_traj = align_traj_to_ref_using_range(
        ref_traj=ref_traj, traj=sample_traj, matching_resids=matching_resids_alignment
    )

    # Define matched residue ids for metric computation
    matching_resids_metric = [
        m for m in matching_resids if m[0] in set(localresidinfo.metric_resids)
    ]

    sel_metric_ref = fmt_resseq_selection([m[0] for m in matching_resids_metric])
    sel_metric_sample = fmt_resseq_selection([m[1] for m in matching_resids_metric])

    atomsel_metric_ref = ref_traj.topology.select(sel_metric_ref)
    atomsel_metric_sample = sample_traj.topology.select(sel_metric_sample)

    assert len(atomsel_metric_ref) == len(atomsel_metric_sample), "Uneven metric atom selections"

    ref_traj_filtered = ref_traj.atom_slice(atomsel_metric_ref)
    sample_traj_filtered = sample_traj.atom_slice(atomsel_metric_sample)

    # RMSD (manual)
    if MetricType.RMSD in metric_types or MetricType.TMSCORE in metric_types:
        pairwise_distances = np.sqrt(
            np.sum(
                ((sample_traj_filtered.xyz - ref_traj_filtered.xyz) * NM2ANGS) ** 2,
                axis=-1,
            )
        )  # [n_frames, n_residues * 4]

    if MetricType.RMSD in metric_types:
        metrics[MetricType.RMSD] = np.mean(pairwise_distances, axis=-1)  # [n_frames],

    # TM-score (manual, from https://en.wikipedia.org/wiki/Template_modeling_score)
    # TODO: Compute TMscore over the alignment region using tm-align
    if MetricType.TMSCORE in metric_types:
        metrics[MetricType.TMSCORE] = tm_score_local(pairwise_distances)

    # dssp accuracy
    if MetricType.DSSP_ACC in metric_types:
        metrics[MetricType.DSSP_ACC] = dssp_match_accuracy(
            ref_traj, sample_traj, matching_resids=matching_resids_metric
        )

    # lddt
    if MetricType.LDDT in metric_types:
        metrics[MetricType.LDDT] = lddt(
            ref_traj, sample_traj, matching_resids=matching_resids_metric
        )

    # contact-match
    if MetricType.CONTACT_DISTANCE in metric_types:
        metrics[MetricType.CONTACT_DISTANCE] = contact_distances(
            ref_traj, sample_traj, matching_resids=matching_resids_metric
        )
    return metrics


def get_metrics_against_references(
    sample_traj: mdtraj.Trajectory,
    ref_trajs: list[mdtraj.Trajectory],
    metric_types: list[MetricType],
    localresidinfo: LocalResidInfo | None = None,
    filter_backbone: bool = True,
) -> dict[MetricType, np.ndarray]:
    """Gets all-atom RMSDs of between the structures provided in `sample_traj`
    and reference PDB structures in `ref_pdbs`

    Args:
        sample_traj: Sample trajectory instance
        ref_trajs: Location of reference pdb files.
        metric_types: List of metrics to be computed
        localresidinfo: Object containing information about which residues to use for alignment
                        and which ones for metric computation in `ref_trajs`.
        filter_samples_backbone: Whether to filter `sample_traj` with mdtraj
                                selection `backbone`. Useful since samples are generated with
                                CB atoms, which can cause mismatch errors with reference samples.
    Returns:
        A dictionary indexed by type of metric containing a list the same length of `sample_traj` per
        each reference trajectory in `ref_trajs`
    """

    metrics: dict[MetricType, list[np.ndarray]] = {}

    # Sample traj needs to be filtered to exclude side-chain knob (CB) from backbone.
    if filter_backbone:
        sample_traj = sample_traj.atom_slice(sample_traj.topology.select("backbone"))
        ref_trajs = [
            ref_traj.atom_slice(ref_traj.topology.select("backbone")) for ref_traj in ref_trajs
        ]

    for ref_traj in ref_trajs:
        matching_resids = seq_pairwise_align_trajs(
            ref_traj, sample_traj[0]
        )  # all samples have same sequence

        if localresidinfo is None:
            metrics_ref = calc_metrics_global(
                ref_traj=ref_traj,
                sample_traj=sample_traj,
                matching_resids=matching_resids,
                metric_types=metric_types,
            )
        else:
            metrics_ref = calc_metrics_local(
                ref_traj=ref_traj,
                sample_traj=sample_traj,
                matching_resids=matching_resids,
                localresidinfo=localresidinfo,
                metric_types=metric_types,
            )

        for metric_type, metric_value in metrics_ref.items():
            metrics.setdefault(metric_type, []).append(metric_value)

    # Stack metrics across sample axis
    stacked_metrics: dict[MetricType, np.ndarray] = {}
    for metric_type, metric_list in metrics.items():
        stacked_metrics[metric_type] = np.vstack(metric_list).T

    return stacked_metrics


def get_metrics_between_references(
    ref_trajs: list[mdtraj.Trajectory],
    metric_types: list[MetricType],
    localresidinfo: LocalResidInfo | None = None,
) -> dict[MetricType, dict[tuple[int, int], float]]:
    """
    Computes metrics between reference structures

    Args:
        ref_trajs: Reference trajectories
    """
    metrics_between_references: dict[MetricType, dict[tuple[int, int], float]] = {}

    comparisons_trajs = list(combinations(ref_trajs, 2))
    comparisons_idxs = list(combinations(list(range(len(ref_trajs))), 2))
    metrics_comparisons = Parallel(n_jobs=-1)(
        delayed(get_metrics_against_references)(
            ref_traj_i,
            [ref_traj_j],
            metric_types=metric_types,
            filter_backbone=True,
            localresidinfo=localresidinfo,
        )
        for ref_traj_i, ref_traj_j in comparisons_trajs
    )

    for metric_comparison, comparison_idxs in zip(
        metrics_comparisons, comparisons_idxs, strict=True
    ):
        for metric_type, value_arr in metric_comparison.items():
            if metric_type not in metrics_between_references:
                metrics_between_references[metric_type] = {}
            assert len(value_arr) == 1
            metrics_between_references[metric_type][comparison_idxs] = value_arr[0]

    return metrics_between_references


def evaluate_test_case(
    test_case: str,
    references_dir: StrPath,
    sequence_samples: list[SequenceSample],
    metric_types: list[MetricType],
    references_localresidinfo_dir: StrPath | None = None,
) -> TestCaseResult:
    """Evaluates a single test case. Args are the same as `evaluate_multiconf`."""

    ref_pdbs = sorted(glob(os.path.join(references_dir, test_case, "*.pdb")))
    ref_trajs = [mdtraj.load_pdb(pdb) for pdb in ref_pdbs]
    assert (
        len(ref_pdbs) > 1
    ), f"Not enough references in specified directory for test case {test_case}: {ref_pdbs}!"

    # Get local information if requested and available, otherwise default to global eval
    if references_localresidinfo_dir is not None:
        max_resseq = max(
            [max([res.resSeq for res in traj.topology.residues]) for traj in ref_trajs]
        )
        residinfo_json = os.path.join(references_localresidinfo_dir, f"{test_case}.json")
        if not os.path.exists(residinfo_json):
            LOGGER.warning(
                f"Requested local evaluation for {test_case} but local info JSON not present. Defaulting to global eval."
            )
            localresidinfo = None
        else:
            localresidinfo = LocalResidInfo.from_json(residinfo_json, n_residues=max_resseq)
    else:
        localresidinfo = None

    assert all([os.path.exists(p) for p in ref_pdbs])

    metrics_between_references = get_metrics_between_references(
        ref_trajs, localresidinfo=localresidinfo, metric_types=metric_types
    )

    references_names = [str(os.path.basename(p)) for p in ref_pdbs]

    # Compute metrics per available traj (potentially different sequences per test_case)
    all_metrics_against_references: list[dict[MetricType, np.ndarray]] = []
    topology_ids = []
    for top_index, sequence_sample in enumerate(sequence_samples):
        sample_traj = sequence_sample.get_traj()
        all_metrics_against_references.append(
            get_metrics_against_references(
                sample_traj,
                ref_trajs,
                localresidinfo=localresidinfo,
                metric_types=metric_types,
            )
        )
        topology_ids.append(np.ones(sample_traj.n_frames, dtype=int) * top_index)

    # Flatten metrics against references into a single array
    flattened_metrics_against_references: dict[MetricType, np.ndarray] = {}

    for metric_type in all_metrics_against_references[0].keys():
        flattened_metrics_against_references[metric_type] = np.vstack(
            [mar[metric_type] for mar in all_metrics_against_references]
        )

    topology_ids = np.hstack(topology_ids)

    # Get info of closest samples
    closest_samples = get_closest_sample_per_reference(
        sequence_samples,
        metrics_against_references=flattened_metrics_against_references,
        references=ref_pdbs,
        topology_ids=topology_ids,
    )

    return TestCaseResult(
        test_case=test_case,
        references_names=references_names,
        metrics_between_references=metrics_between_references,
        metrics_against_references=flattened_metrics_against_references,
        closest_samples=closest_samples,
        topology_ids=topology_ids,
    )


def get_closest_sample_per_reference(
    sequence_samples: list[SequenceSample],
    metrics_against_references: dict[MetricType, np.ndarray],
    references: list[str],
    topology_ids: np.ndarray,
) -> dict[MetricType, list[ClosestSample]]:
    """
    Keeps track of the best samples, per metric type.

    Args:
        sequence_samples: A list of evaluated `SequenceSample` instances
        metrics_against_references: Sample metrics against references
        references: String (or path) identifier of the reference
        topology_ids: 1D array that keeps track which topology was used on
                     each row of the values available in `metrics_against_references`.

    Returns:
        A dictionary indexed by metric type, containing information about the closest samples
        to each reference.
    """

    closest_per_metric_type: dict[
        MetricType,
        list[ClosestSample],
    ] = {}
    uq_topology_ids = set(topology_ids)
    assert len(sequence_samples) == len(uq_topology_ids)

    for metric_type, mar in metrics_against_references.items():
        assert mar.shape[1] == len(references)
        operator_type = METRIC_OPERATOR_BETTER[metric_type]
        minmax_fun = np.max if operator_type == operator.gt else np.min
        argminmax_fun = np.argmax if operator_type == operator.gt else np.argmin

        frame_idxs = np.hstack(
            [np.arange(np.sum(topology_ids == uq_t)) for uq_t in uq_topology_ids]
        )

        best_idxs = argminmax_fun(mar, axis=0)
        best_topology_idxs = topology_ids[best_idxs]
        best_vals = list(mar[best_idxs, :])
        best_argvals = [argminmax_fun(v) for v in best_vals]

        best_frame_idxs = frame_idxs[best_idxs]

        best_ss = [sequence_samples[idx] for idx in best_topology_idxs]

        closest_per_metric_type[metric_type] = [
            ClosestSample(
                ss,
                frame_idx=idx,
                metric_type=metric_type,
                metric_value=minmax_fun(val),
                reference_path=references[argval],
            )
            for ss, idx, val, argval in zip(
                best_ss, best_frame_idxs, best_vals, best_argvals, strict=True
            )
        ]

    return closest_per_metric_type


MULTICONF_METRIC_TYPES = [
    MetricType.RMSD,
    MetricType.TMSCORE,
    MetricType.LDDT,
    MetricType.CONTACT_DISTANCE,
    MetricType.DSSP_ACC,
]


def evaluate_multiconf(
    indexed_samples: IndexedSamples,
    references_dir: StrPath,
    metric_types: list[MetricType] | None = None,
    references_localresidinfo_dir: StrPath | None = None,
    n_jobs: int = -1,
) -> dict[str, TestCaseResult]:
    """
    Args:
        indexed_samples: An `IndexedSamples` instance containing samples to evaluate for a given
                         benchmark.
        references_dir: Directory where reference PDB files are stored
        metric_types: Metrics to be computed
        references_localresidinfo_dir: Directory where local reference residue information files (JSON)
                                       are stored.
        n_jobs: Number of processes used for evaluation


    Returns:
        * A dictionary, indexed by test case identifier, containing `TestCaseResult` instances for
          each of them.
    """
    metric_types = metric_types or [MetricType.RMSD]
    assert set(metric_types).issubset(MULTICONF_METRIC_TYPES)
    setup_tm_align()

    with tqdm_joblib(
        tqdm(
            desc="Now computing metrics...", total=len(indexed_samples.test_case_to_sequencesamples)
        )
    ) as _pbar:
        multiconf_evals = Parallel(n_jobs=n_jobs)(
            delayed(evaluate_test_case)(
                test_case,
                references_dir=references_dir,
                sequence_samples=sequence_samples,
                references_localresidinfo_dir=references_localresidinfo_dir,
                metric_types=metric_types,
            )
            for test_case, sequence_samples in indexed_samples.test_case_to_sequencesamples.items()
        )

    return {mconfeval.test_case: mconfeval for mconfeval in multiconf_evals}


def split_holo_apo(multiconf_eval: TestCaseResult) -> tuple[TestCaseResult, TestCaseResult]:
    """
    Takes an eval on the cryptic pocket benchmark and splits the results into holo and apo, respectively
    """

    metadata = pd.read_csv(os.path.join(MULTICONF_ASSET_DIR, "crypticpocket/references.csv"))

    row = metadata[metadata["test_case"] == multiconf_eval.test_case].iloc[0]

    assert multiconf_eval.references_names is not None

    is_holo = [
        reference_name.split(".pdb")[0] == row["holo_pdbidchain"].upper()
        for reference_name in multiconf_eval.references_names
    ]

    holo_index = np.argmax(is_holo)

    multiconf_eval_holo = dataclasses.replace(multiconf_eval)

    multiconf_eval_holo.metrics_against_references = {
        k: v[:, holo_index : holo_index + 1]
        for k, v in multiconf_eval_holo.metrics_against_references.items()
    }

    apo_index = 1 - holo_index

    multiconf_eval_apo = dataclasses.replace(multiconf_eval)

    multiconf_eval_apo.metrics_against_references = {
        k: v[:, apo_index : apo_index + 1]
        for k, v in multiconf_eval_apo.metrics_against_references.items()
    }

    return multiconf_eval_holo, multiconf_eval_apo
