---
# For reference on dataset card metadata, see the spec: https://github.com/huggingface/hub-docs/blob/main/datasetcard.md?plain=1
# Doc / guide: https://huggingface.co/docs/hub/datasets-cards
{{ card_data }}
---

# Dataset Card for {{ pretty_name | default("BioEmu 1.0 Benchmark", true) }}

The BioEmu 1.0 Benchmark dataset contains protein structure data and experimental folding free energies for evaluating models that predict protein structure ensembles.

{{ dataset_summary | default("", true) }}

## Dataset Details

### Dataset Description

The BioEmu 1.0 Benchmark dataset delivers data for three different benchmarks that are used in the [BioEmu manuscript](https://www.biorxiv.org/content/10.1101/2024.12.05.626885). Please refer to the manuscript for full details.
The dataset contains the following parts:

- Protein structures for multiconf benchmark: Contains selected and cleaned structures from the [PDB](https://www.rcsb.org) for proteins that have multiple conformations.
- Folding free energies and relative folding free energies from mutations for dG benchmark, extracted from [wet-lab experiments](https://doi.org/10.1038/s41586-023-06328-6), and reference structures from the [PDB](https://www.rcsb.org).
- Projections and ground truth for MAE benchmark: Projection matrices and projected MD trajectories from long MD simulations.

{{ dataset_description | default("", true) }}

- **Curated by:** Sarah Lewis, Tim Hempel, José Jiménez-Luna, Michael Gastegger, Yu Xie, Andrew Y. K. Foong, Victor García Satorras, Osama Abdin, Bastiaan S. Veeling, Iryna Zaporozhets, Yaoyi Chen, Soojung Yang, Arne Schneuing, Jigyasa Nigam, Federico Barbero, Vincent Stimper, Andrew Campbell, Jason Yim, Marten Lienen, Yu Shi, Shuxin Zheng, Hannes Schulz, Usman Munir, Cecilia Clementi, Frank Noé {{ curators | default("[More Information Needed]", true)}}
- **Funded by:** Microsoft Research AI for Science {{ funded_by | default("[More Information Needed]", true)}}
- **Shared by:** Microsoft Research AI for Science {{ shared_by | default("[More Information Needed]", true)}}
- **License:** CDLA-2.0 (https://cdla.dev/permissive-2-0/) {{ license | default("[More Information Needed]", true)}}


### Dataset Sources

The dataset is shared as part of the BioEmu benchmark repository.

- **Repository:** https://github.com/microsoft/bioemu-benchmarks/ {{ repo | default("[More Information Needed]", true)}}
- **Paper:** https://www.biorxiv.org/content/10.1101/2024.12.05.626885 {{ paper | default("[More Information Needed]", true)}}

## Uses

The dataset can be used to benchmark models that emulate the thermodynamic ensemble of protein structures (i.e., Boltzmann distribution) given an amino acid sequence.
Definitions of the benchmarks and related metrics are given in this repository.

### Direct Use

Examples of direct usages include but are not limited to evaluating a model's capabilities
- to capture multiple protein conformations in the multi-conf benchmark
- to predict protein stabilities with the folding free energy benchmark
- as an MD emulator with the MAE benchmark


{{ direct_use | default("[More Information Needed]", true)}}

### Out-of-Scope Use
Out-of-scope use cases include but are not limited to benchmarking models with higher than atomistic precision or optimizing models for cases that are not in this benchmark (e.g., multi-chain conformational dynamics).
The dataset is intended for research and experimental purposes. Further testing/development are needed before considering its application in real-world scenarios.
<!-- This section addresses misuse, malicious use, and uses that the dataset will not work well for. -->

{{ out_of_scope_use | default("[More Information Needed]", true)}}

## Dataset Structure
Raw data for the datasets is located in the repo under the `assets` directory.  More specifically:
- The static multiconformational benchmarks are located under the `multiconf_benchmark_x.y` subdirectory, where `x.y` denotes the version of the benchmark that is shipped with a particular version of the repository. A `README.md` with more detailed information about each benchmark and its structure is also included.
- The folding free energy benchmark is located under the `folding_free_energies_benchmark_x.y` subdirectory, with the same `x.y` version encoding convention as for the multiconf benchmark. The directory includes a `README.md` with a detailed description of further directory structure and file contents.
- The MD emulation MAE benchmark can be found in the `md_emulation_benchmark_x.y` subdirectory, with `x.y` used as a version identifier. Details on the directory substructure as well as file contents cre provided in the `README.md` file included in the directory.

{{ dataset_structure | default("[More Information Needed]", true)}}

## Dataset Creation

### Curation Rationale
The dataset was created to define a benchmark for deep learning models that predict conformational plasticity in proteins.
<!-- Motivation for the creation of this dataset. -->

{{ curation_rationale_section | default("[More Information Needed]", true)}}

### Source Data

- PDB structures for multiconf & free energy benchmarks downloaded from the PDB ([https://www.rcsb.org]) and processed internally
- Free energy differences from mega-scale dataset from [https://zenodo.org/records/7992926] for free energy benchmark, with corresponding native state structures from the PDB ([https://www.rcsb.org]).
- Projection matrices and projected MD trajectories derived from in-house MD data

#### Data Collection and Processing

All data in this benchmark has been collected to represent structural heterogeneity in proteins. For details, check the [BioEmu manuscript](https://www.biorxiv.org/content/10.1101/2024.12.05.626885).
<!-- This section describes the data collection and processing process such as data selection criteria, filtering and normalization methods, tools and libraries used, etc. -->

{{ data_collection_and_processing_section | default("[More Information Needed]", true)}}

## Bias, Risks, and Limitations
This benchmark is designed for single-chain protein systems and does not reflect conformational dynamics of multi-chain systems.
Multiconf and dG benchmark data are biased towards the set of proteins that can be resolved by currently available experimental methods. 
MAE benchmark data provides structural ensembles from molecular dynamics, which might be biased towards starting structures and by the used empirical force fields.
The benchmark is limited by the low number of systems that are either structurally resolved in experimental studies or extensively sampled with MD.
Given the comparably small number of systems, relying exclusively on this data may harm model performance.

We do not recommend using this benchmark in commercial or real-world applications without further testing and development.
It is being released for (non-commercial) research purposes.

<!-- This section is meant to convey both technical and sociotechnical limitations. -->

{{ bias_risks_limitations | default("[More Information Needed]", true)}}

### Recommendations
We recommend combining the benchmark data presented here with other data, e.g., with single-structured data for protein structure prediction benchmarks.

Users are responsible for sourcing their datasets responsibly. This could include securing appropriate copyrights, appropriate data preprocessing and validation, or ensuring the anonymization of data prior to use in research.

<!-- This section is meant to convey recommendations with respect to the bias, risk, and technical limitations. -->

{{ bias_recommendations | default("Users should be made aware of the risks, biases and limitations of the dataset. More information needed for further recommendations.", true)}}

## Citation
Please cite the [BioEmu manuscript](https://www.biorxiv.org/content/10.1101/2024.12.05.626885) when using the benchmark or dataset.
<!-- If there is a paper or blog post introducing the dataset, the APA and Bibtex information for that should go in this section. -->

**BibTeX:**
```
@misc{lewis_scalable_2024,
  title = {Scalable Emulation of Protein Equilibrium Ensembles with Generative Deep Learning},
  author = {Lewis, Sarah and Hempel, Tim and Jim{\'e}nez Luna, Jos{\'e} and Gastegger, Michael and Xie, Yu and Foong, Andrew Y. K. and Garc{\'i}a Satorras, Victor and Abdin, Osama and Veeling, Bastiaan S. and Zaporozhets, Iryna and Chen, Yaoyi and Yang, Soojung and Schneuing, Arne and Nigam, Jigyasa and Barbero, Federico and Stimper, Vincent and Campbell, Andrew and Yim, Jason and Lienen, Marten and Shi, Yu and Zheng, Shuxin and Schulz, Hannes and Munir, Usman and Clementi, Cecilia and No{\'e}, Frank},
  year = {2024},
  doi = {10.1101/2024.12.05.626885},
  archiveprefix = {BioRXiv},
  url = {https://www.biorxiv.org/content/10.1101/2024.12.05.626885}
}
```

{{ citation_bibtex | default("[More Information Needed]", true)}}

**APA:**
Lewis, S., Hempel, T., Jiménez Luna, J., Gastegger, M., Xie, Y., Foong, A. Y. K., García Satorras, V., Abdin, O., Veeling, B. S., Zaporozhets, I., Chen, Y., Yang, S., Schneuing, A., Nigam, J., Barbero, F., Stimper, V., Campbell, A., Yim, J., Lienen, M., … Noé, F. (2024). Scalable emulation of protein equilibrium ensembles with generative deep learning. Molecular Biology. https://doi.org/10.1101/2024.12.05.626885

{{ citation_apa | default("[More Information Needed]", true)}}


## Dataset Card Authors
- Tim Hempel (timhempel@microsoft.com)
- Jose Salvador Jimenez Luna (jjimenezluna@microsoft.com)
- Michael Gastegger (mgastegger@microsoft.com)

{{ dataset_card_authors | default("[More Information Needed]", true)}}

## Dataset Card Contact
- Franke Noe (franknoe@microsoft.com)
- Sarah Lewis (sarahlewis@microsoft.com)
- Jose Salvador Jimenez Luna (jjimenezluna@microsoft.com)
- Michael Gastegger (mgastegger@microsoft.com)

{{ dataset_card_contact | default("[More Information Needed]", true)}}
