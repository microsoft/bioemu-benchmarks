from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import mdtraj
import numpy as np
import pandas as pd
from tqdm import tqdm

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.logger import get_logger
from bioemu_benchmarks.utils import StrPath, get_physical_traj_indices

LOGGER = get_logger(__name__)


@dataclass(frozen=True, eq=True)
class SequenceSample:
    """Paths to files containing samples for a single sequence"""

    topology_file: StrPath  # PDB
    trajectory_file: StrPath  # XTC
    frames_to_include: np.ndarray | None = (
        None  # Indices to keep in trajectory (e.g. [0, 4, 7]). None means keep all frames.
    )

    def get_traj(self) -> mdtraj.Trajectory:
        """
        Returns associated `Trajectory` instance
        """
        traj = mdtraj.load(self.trajectory_file, top=self.topology_file)

        if self.frames_to_include is not None:
            traj = traj[self.frames_to_include]
        return traj

    def __post_init__(self):
        if self.frames_to_include is not None:
            assert self.frames_to_include.size > 0


def assert_topology_has_backbone_atoms(top: mdtraj.Topology) -> None:
    """
    Assert that topology is not missing any backbone atoms.
    """
    for residue in top.residues:
        atom_names = set([atom.name for atom in residue.atoms])
        if not atom_names.issuperset(["N", "CA", "C", "O"]):
            raise MissingBackbone


class MissingBackbone(Exception):
    pass


class NoSamples(Exception):
    pass


class MissingTopology(Exception):
    pass


def find_samples_in_dir(samples_dir: StrPath) -> list[SequenceSample]:
    """
    Find pairs of .xtc and .pdb files in a directory.

    Samples for each sequence should be stored as a pair of files: a .pdb file containing the topology and a .xtc file containing a 'trajectory'
    where each 'frame' is a sampled conformation of the protein.

    We first look for .xtc files in the specified directory. For each .xtc file we first look for a
    .pdb file with the same name, or failing that, a 'topology.pdb' in the same directory.
    """
    sequence_samples = []
    xtc_files = [x for x in Path(samples_dir).glob("**/*.xtc")]
    for f in xtc_files:
        if f.with_suffix(".pdb").exists():
            # First look for a PDB file with the same name.
            sequence_samples.append(
                SequenceSample(topology_file=f.with_suffix(".pdb"), trajectory_file=f)
            )
        elif (f.parent / "topology.pdb").exists():
            # Failing that, look for 'topology.pdb' in the same directory.
            sequence_samples.append(
                SequenceSample(topology_file=f.parent / "topology.pdb", trajectory_file=f)
            )
        else:
            raise MissingTopology(f"Could not find .pdb file to use as topology for {f}.")
    return sequence_samples


def select_relevant_samples(
    sequence_samples: list[SequenceSample], relevant_sequences: set[str] | set[None]
) -> list[SequenceSample]:
    seqs = [mdtraj.load_topology(ss.topology_file).to_fasta()[0] for ss in sequence_samples]

    if relevant_sequences == set([None]):
        return [files for files in sequence_samples]
    else:
        irrelevant_sequences_sampled = set(seqs).difference(relevant_sequences)
        if irrelevant_sequences_sampled:
            LOGGER.info(
                f"Ignoring samples for {len(irrelevant_sequences_sampled)} irrelevant sequences."
            )

        return [files for files, seq in zip(sequence_samples, seqs) if seq in relevant_sequences]


class IndexedSamples: 
    def __init__(self, test_case_to_sequencesamples: dict[str, list[SequenceSample]]):
        assert test_case_to_sequencesamples, "Empty input"
        self.test_case_to_sequencesamples = test_case_to_sequencesamples

    @classmethod
    def from_benchmark(
        cls, benchmark: Benchmark, sequence_samples: list[SequenceSample], test_case_list: list[str], include_relevant_sequences: bool =False
    ) -> "IndexedSamples":
        
        if include_relevant_sequences:
            benchmark_sequences = benchmark.metadata["sequence"]
        else:
            benchmark_sequences = [None]

        # Ignore irrelevant samples, i.e., samples of sequences that are not in the benchmark.
        sequence_samples = select_relevant_samples(sequence_samples, set(benchmark_sequences))

        sampled_sequences: set[str] = set()
        sequence_sample_to_test_cases: dict[SequenceSample, list[str]] = defaultdict(list)

        for sequence_sample in sequence_samples:
            top = mdtraj.load_topology(sequence_sample.topology_file)
            assert top.n_chains == 1, "Multi-chain sample PDB file not supported"
            sequence: str = top.to_fasta()[0]

            assert_topology_has_backbone_atoms(top)

            if include_relevant_sequences:
                assoc_test_cases = benchmark.metadata.loc[
                    benchmark.metadata["sequence"] == sequence
                ].test_case
                if isinstance(assoc_test_cases, str):
                    sequence_sample_to_test_cases[sequence_sample].append(assoc_test_cases)
                else:
                    assert isinstance(assoc_test_cases, pd.Series)
                    sequence_sample_to_test_cases[sequence_sample].extend(assoc_test_cases)

               
            else:
                assoc_test_cases = test_case_list
                sequence_sample_to_test_cases[sequence_sample]= assoc_test_cases

            sampled_sequences.add(sequence)



        # Check if any relevant sequences have been sampled.
        if len(sampled_sequences) == 0:
            raise NoSamples("No samples found for benchmark.")

        # Report if missing sampled sequences
        if include_relevant_sequences:
            missing_sequences = set(benchmark_sequences).difference(sampled_sequences)
            if len(missing_sequences) > 0:
                LOGGER.warning(
                    f"Missing samples for {len(missing_sequences)} sequence(s) "
                    f"for this benchmark: {missing_sequences}"
                )

        # Build map from test case to samples
        test_case_to_sequencesamples: dict[str, list[SequenceSample]] = defaultdict(list)
        for sequence_sample, test_cases in sequence_sample_to_test_cases.items():
            for test_case in test_cases:
                test_case_to_sequencesamples[test_case].append(sequence_sample)

        return cls(test_case_to_sequencesamples)

    def get_trajs_for_test_case(self, test_case: str) -> list[mdtraj.Trajectory]:
        """Returns trajectory instances associated with a particular identifier"""
        system_samples = self.test_case_to_sequencesamples[test_case]
        return [ss.get_traj() for ss in system_samples]

    def get_all_trajs(self) -> dict[str, list[mdtraj.Trajectory]]:
        """Returns all available (sampled) trajectory instances for this benchmark."""
        test_case_to_trajs = {}
        for test_case in self.test_case_to_sequencesamples.keys():
            test_case_to_trajs[test_case] = self.get_trajs_for_test_case(test_case)
        return test_case_to_trajs

    def __repr__(self):
        return f"IndexedSamples(test_cases={list(self.test_case_to_sequencesamples)})"


def filter_unphysical_sequencesample(sequence_sample: SequenceSample, **kwargs) -> SequenceSample:
    """
    Adds physically valid frame information to `sequence_sample` instance.
    """

    assert sequence_sample.frames_to_include is None
    valid_indices = get_physical_traj_indices(sequence_sample.get_traj(), **kwargs)
    return SequenceSample(
        topology_file=sequence_sample.topology_file,
        trajectory_file=sequence_sample.trajectory_file,
        frames_to_include=valid_indices,
    )


def filter_unphysical_samples(
    indexed_samples: IndexedSamples,
) -> tuple[IndexedSamples, dict[str, np.ndarray]]:
    """
    Filters 'unphysical' samples froma given definition of samples. Filtered
    trajectories are saved on a temporary directory.

    Samples are excluded according to the default criteria described on:
    `eval.multiconf.utils.filter_unphysical_traj`

    Args:
        indexed_samples: A `IndexedSamples` instance for a given benchmark.

    Returns:
        * Another `IndexedSamples` instance with filtered samples
        * A dictionary indexed by test identifier detailing the percentage of
        samples that have remained after filtering, per each of the associated
        `SequenceSample` instances.
    """

    test_case_to_trajs = indexed_samples.get_all_trajs()

    # Get original sample sizes
    test_case_to_sizes: dict[str, list[int]] = {
        test_case: [traj.n_frames for traj in trajs]
        for test_case, trajs in test_case_to_trajs.items()
    }

    # Filter unphysical samples
    test_case_to_physical_sequencesamples: dict[str, list[SequenceSample]] = {
        test_case: [filter_unphysical_sequencesample(ss) for ss in sequence_samples]
        for test_case, sequence_samples in tqdm(
            indexed_samples.test_case_to_sequencesamples.items(),
            desc="Filtering unphysical samples...",
        )
    }


    # Get filtered sample sizes
    test_case_to_filtered_sizes: dict[str, list[int]] = defaultdict(list)
    for test_case, sequence_samples in test_case_to_physical_sequencesamples.items():
        for ss in sequence_samples:
            assert ss.frames_to_include is not None  # shut up mypy
            test_case_to_filtered_sizes[test_case].append(len(ss.frames_to_include))

    # Warn user if all samples have been filtered for a given test_case
    all_filtered_test_cases: list[str] = []

    for test_case, sizes in test_case_to_filtered_sizes.items():
        for size in sizes:
            if size == 0:
                all_filtered_test_cases.append(test_case)
    if all_filtered_test_cases:
        LOGGER.warning(
            f"Filtered all samples for the following test cases: {all_filtered_test_cases}. Will be ignored from analyses"
        )

    filtered_indexed_samples = IndexedSamples(test_case_to_physical_sequencesamples)

    # Get stats and return them
    test_case_to_perc_kept_samples: dict[str, np.ndarray] = {}  # {test_case: [0.8, 0.6, ...]}
    for test_case in test_case_to_sizes.keys():
        o_sizes = np.array(test_case_to_sizes[test_case])
        f_sizes = np.array(test_case_to_filtered_sizes[test_case])
        assert len(o_sizes) == len(f_sizes)
        test_case_to_perc_kept_samples[test_case] = f_sizes / o_sizes
    return filtered_indexed_samples, test_case_to_perc_kept_samples
