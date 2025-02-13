import os
from enum import Enum
from functools import cached_property
from io import StringIO
from typing import Literal

import fire
import pandas as pd

from bioemu_benchmarks.paths import (
    FOLDING_FREE_ENERGY_ASSET_DIR,
    MD_EMULATION_ASSET_DIR,
    MULTICONF_ASSET_DIR,
)


class Benchmark(str, Enum):
    MULTICONF_OOD60 = "multiconf_ood60"
    MULTICONF_OODVAL = "multiconf_oodval"
    MULTICONF_DOMAINMOTION = "multiconf_domainmotion"
    MULTICONF_CRYPTICPOCKET = "multiconf_crypticpocket"
    SINGLECONF_LOCALUNFOLDING = "singleconf_localunfolding"
    FOLDING_FREE_ENERGIES = "folding_free_energies"
    MD_EMULATION = "md_emulation"

    @cached_property
    def asset_dir(self) -> str:
        """
        Directory with assets for each benchmark
        """
        if self in MULTICONF_BENCHMARKS or self == Benchmark.SINGLECONF_LOCALUNFOLDING:
            return os.path.join(MULTICONF_ASSET_DIR, self.value.split("_")[1])
        elif self == Benchmark.FOLDING_FREE_ENERGIES:
            return os.path.join(FOLDING_FREE_ENERGY_ASSET_DIR, "folding_free_energies")
        elif self == Benchmark.MD_EMULATION:
            return os.path.join(MD_EMULATION_ASSET_DIR, "md_emulation")
        else:
            raise ValueError("Benchmark not recognised")

    @cached_property
    def metadata(self) -> pd.DataFrame:
        """
        Returns sequence info for the benchmark
        """
        metadata = pd.read_csv(os.path.join(self.asset_dir, "testcases.csv"))
        return metadata

    @cached_property
    def default_samplesize(self) -> list[int]:
        """
        Gets recommended sample sizes for each of the available
        benchmarks
        """
        if self in MULTICONF_BENCHMARKS or self == Benchmark.SINGLECONF_LOCALUNFOLDING:
            return [4000] * len(self.metadata)
        elif self == Benchmark.MD_EMULATION:
            return [10000] * len(self.metadata)
        elif self == Benchmark.FOLDING_FREE_ENERGIES:
            ss_metadata = pd.read_csv(
                os.path.join(
                    self.asset_dir,
                    "system_info.csv",
                )
            )
            ss_metadata.set_index("sequence", inplace=True)

            sample_sizes = []
            for seq in self.metadata["sequence"]:
                sample_sizes.append(ss_metadata.loc[ss_metadata.index == seq].num_samples.iloc[0])
            return sample_sizes
        else:
            raise ValueError("Unrecognised benchmark name")


MULTICONF_BENCHMARKS: list[Benchmark] = [
    Benchmark.MULTICONF_OOD60,
    Benchmark.MULTICONF_OODVAL,
    Benchmark.MULTICONF_DOMAINMOTION,
    Benchmark.MULTICONF_CRYPTICPOCKET,
]


def main(
    benchmark_name: Literal[
        "multiconf_ood60",
        "multiconf_oodval",
        "multiconf_domainmotion",
        "multiconf_crypticpocket",
        "singleconf_localunfolding",
        "folding_free_energies",
        "md_emulation",
    ],
) -> None:
    """
    Main entrypoint to obtain benchmark sequences and recommended sample
    sizes. Prints out a CSV with sequence and sample size information.

    Args:
        benchmark_name: The name of the benchmark to be evaluated
        num_samples (optional): If set, the generated CSV will have the `num_samples`
                    field set to this value. Otherwise a benchmark-specific default
                    is set.
    """
    benchmark = Benchmark(benchmark_name)
    output = StringIO()
    benchmark.metadata["default_sample_size"] = benchmark.default_samplesize
    benchmark.metadata.to_csv(output, index=False)
    print(output.getvalue())  # print to stdout


if __name__ == "__main__":
    fire.Fire(main)
