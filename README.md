# Biomolecular Emulator - Benchmarks (BioEmu-Benchmarks)

Accompanying benchmark code for the BioEmu [paper](https://www.science.org/doi/10.1126/science.adv9817). For the BioEmu sampling code please check [here](https://www.github.com/microsoft/bioemu).

## Table of Contents
- [Installation](#installation)
- [Benchmarks](#available-benchmarks)
- [Usage](#usage)
- [Citation](#citation)
- [Get in touch](#get-in-touch)

## Installation

`bioemu-benchmarks` is provided as a pip-installable package, requiring Python >= 3.10:

```bash
pip install bioemu-benchmarks
```

Alternatively, an `environment.yml` file is also provided in case you prefer to clone the repository and create a conda environment from scratch.

## Available benchmarks

`bioemu-benchmarks` implements the following benchmarks for evaluating the emulation performance of models:
* `multiconf_ood60`: Measures local protein conformational changes on a set that is different from the training set in terms of sequence similarity.
* `multiconf_domainmotion`: Measures global protein domain motions. Main metric measured is global RMSD.
* `singleconf_localunfolding`: Measures local protein unfolding. Main metric measured is a custom fraction of native contacts calculation on a pre-defined set of residues.
* `multiconf_crypticpocket`: Measures pocket backbone changes upon ligand binding. Default metric is local RMSD on a predefined set of residues in, or close to, a given binding pocket. 
* `md_emulation`: Measure the match between model samples and target molecular dynamics distribution on low dimensional free energies surfaces. Includes free energy MAE and RMSE, as well as coverage of the samples as metrics.
* `folding_free_energies`: Measures the ability to predict folding free energies ($\Delta G$ and $\Delta\Delta G$) based on the provided samples. Main metrics are the MAEs, as well as correlation coefficients on both metrics. 

For details of the different benchmarks, please refer to the SI of the accompanying [publication](https://www.science.org/doi/10.1126/science.adv9817).

## Usage

`bioemu-benchmarks` provides both a Python interface as well as a regular CLI. The CLI is intended to provide quick access to the different benchmark, while the python API provides more flexibility.

### Sample format
To run the benchmarks in this repository you will need to prepare samples in `.xtc` format so that they can be loaded as `mdtraj.Trajectory` objects.
Each `.xtc` file needs an accompanying `.pdb` file which defines the topology. When loading samples, for each `.xtc` file, the code will look for either a `.pdb` file of the same name or failing that, a `topology.pdb` file in the same directory.  For example, you can store your samples like this
```bash
my_samples/
├── foo.pdb
├── foo.xtc
├── bar.pdb
├── bar.xtc
...
```
or like this (which [bioemu](https://github.com/microsoft/bioemu) will do)
```bash
my_samples/
├── foo/
├──── samples.xtc
├──── topology.pdb
├── bar/
├──── samples.xtc
├──── topology.pdb
```
In order to know which sequences to sample for each benchmark, users can check the `testcases.csv` file under the `assets/<benchmark_type>/<benchmark>` folder on this repo. Alternatively, either the python API or CLI can be used to get sample specifications (see below for examples).

### Bash CLI

Upon installation, the `bioemu-bench` benchmark script is added to the PATH. This script provides a simple CLI for running benchmarks, collecting results and getting benchmark specifications (sequences to samples and recommended number of samples).

#### Loading samples and running benchmarks

To run a single or multiple benchmarks, the `eval` mode of the script is used:
```bash
bioemu-bench eval <output_dir> --benchmarks / -b [...] --sample_dirs / -s [...]
```
`<output_dir>` is the directory to which results will be written and `--benchmarks` specifies the benchmarks to evaluate on the given samples. It can take a single or list of benchmarks as input (for available benchmarks see `bioemu-bench eval --help`). There also is the option to specify `--benchmarks all`, in which case all benchmarks are run. The `--sample_dirs` option takes the path to a directory (or multiple directories) from which samples will be loaded (expecting the format described above).

`bioemu-bench eval` will collect results in the <output_dir>, with each requested benchmark getting its own subdirectory:
```
<output_dir>
├── benchmark_metrics.json
├── domainmotion
│   ├── ...
│   └── results.pkl
...   ...
├── folding_free_energies
│   ├── ...
│   └── scatter_dG.png
└── md_emulation
    ├── ...
    └── results_projections.npz
```
The file `benchmark_metrics.json` collects aggregate metrics for all benchmarks in `json` format.

#### Getting sample specifications

The `specs` mode of `bioemu-bench` can be used to generate a CSV file collecting sequence information and recommended number of samples for the requested benchmarks:
```bash
bioemu-bench specs <output_csv> --benchmarks/-b [...] 
```
`<output_csv>` is the output CSV generated by the script. `--benchmarks` once again specifies the benchmarks for which this information should be generated.

### Python API

The python API provides a set of tools for evaluating samples according to a benchmark, the central ones being:
* `Benchmark`: an Enum that defines the available benchmarks in the repository.
* `IndexedSamples`: used for loading, validating and optionally filtering samples via the `filter_unphysical_samples` function available under that module.
* `Evaluator` functions: These define the evaluations to be done and are called on a instance of `IndexedSamples`. `evaluator_utils.py` provides a function for retrieving the evaluator function for each `Benchmark`.
* `BenchmarkResults` classes: Each `Evaluator` returns a `BenchmarksResults` instance as an output. This classes collect the benchmark results and offer utilities for storing results (`save_results`), plotting (`plot`) and assessing aggregate metrics (`get_aggregate_metrics`).

#### Loading samples and running benchmarks

An example for performing the `ood60` benchmark on a set of samples would look like the following:

```python
from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.samples import IndexedSamples, filter_unphysical_samples, find_samples_in_dir
from bioemu_benchmarks.evaluator_utils import evaluator_from_benchmark

# Specify the benchmark you want to run (e.g., OOD60)
benchmark = Benchmark.MULTICONF_OOD60

# This validates samples
sequence_samples = find_samples_in_dir("/path/to/your/sample_dir")
samples = IndexedSamples.from_benchmark(benchmark=benchmark, sequence_samples=sequence_samples)

# Filter unphysical-looking samples from getting evaluated
samples, _sample_stats = filter_unphysical_samples(samples)

# Instanstiate an evaluator for a given benchmark
evaluator = evaluator_from_benchmark(benchmark=benchmark)

# `results` has methods for plotting / computing summary metrics
results = evaluator(samples)
results.plot('/path/to/result/plots') 
results.save_results('/path/to/result/metrics/')
```

#### Getting sample specifications

The python API also offers a way to get information on the sequences to sample for a benchmark from its `metadata` attribute:
```python
from bioemu_benchmarks.benchmarks import Benchmark

# Returns a pandas df with info about the Multiconf OOD60 benchmark
metadata = Benchmark.MULTICONF_OOD60.metadata
```

## Citation
If you are using our code or model, please consider citing our work:
```bibtex
@article{bioemu2025,
  title={Scalable emulation of protein equilibrium ensembles with generative deep learning},
  author={Lewis, Sarah and Hempel, Tim and Jim{\'e}nez-Luna, Jos{\'e} and Gastegger, Michael and Xie, Yu and Foong, Andrew YK and Satorras, Victor Garc{\'\i}a and Abdin, Osama and Veeling, Bastiaan S and Zaporozhets, Iryna and others},
  journal={Science},
  pages={eadv9817},
  year={2025},
  publisher={American Association for the Advancement of Science}
}
```

## Licensing
The code of this project is licensed under the MIT License. See the LICENSE file for details.
The accompanying dataset is licensed under the Community Data License Agreement – Permissive, Version 2.0 (CDLA-Permissive-2.0). See the DATASET_LICENSE for details.

## Trademarks

This project may contain trademarks or logos for projects, products, or services. Authorized use of Microsoft 
trademarks or logos is subject to and must follow 
[Microsoft's Trademark & Brand Guidelines](https://www.microsoft.com/en-us/legal/intellectualproperty/trademarks/usage/general).
Use of Microsoft trademarks or logos in modified versions of this project must not cause confusion or imply Microsoft sponsorship.
Any use of third-party trademarks or logos are subject to those third-party's policies.

## Get in touch
If you have any questions not covered here, please create an issue or contact the BioEmu team by writing to the corresponding author on our [paper](https://www.science.org/doi/10.1126/science.adv9817).