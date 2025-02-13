import pickle
from abc import ABC, abstractmethod
from dataclasses import dataclass

from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.utils import StrPath


@dataclass
class BenchmarkResults(ABC):
    """
    Base data class for collecting benchmark results.
    """

    benchmark: Benchmark

    def to_pickle(self, outfile: StrPath) -> None:
        """
        Save results class as pickle file.

        Args:
            outfile: Path to output file.
        """
        with open(outfile, "wb") as handle:
            pickle.dump(self, handle)

    @classmethod
    def from_pickle(cls, infile: StrPath):
        """
        Load results class from pickle file.

        Args:
            infile: Path to pickle file.

        Returns:
            Loaded results class.
        """
        with open(infile, "rb") as handle:
            read_pkl = pickle.load(handle)

        return cls(**read_pkl.__dict__)

    @abstractmethod
    def save_results(self, output_dir: StrPath) -> None:
        """
        Save individual evaluator results in accessible files (txt, csv, npz).

        Args:
            output_dir: Directory to which result outputs should be saved.
        """
        ...

    @abstractmethod
    def plot(self, output_dir: StrPath) -> None:
        """
        Generate plots associated with benchmark data and write to output directory.

        Args:
            output_dir: Directory where plots will be written.
        """
        ...

    @abstractmethod
    def get_aggregate_metrics(self) -> dict[str, float]:
        """Collect and return aggregate metrics in flat dictionary."""
        ...
