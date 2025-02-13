# Definitions for the multi-conformation benchmarks v0.1

Each benchmark comes with its own subdirectory, indexed by test case identifier, containing the following information:

* A `indexed_samples.csv` file containing sequence information for each evaluated test case identifier (which in most cases maps to a Uniprot Id.)
* A `references.csv` file detailing the PDB Ids, `label_asym_ids`, and test case identifiers.
* `reference`: Contains processed backbone PDB files as used in our evaluation pipeline. Note that these PDB files may have a different residue numbering as the ones deposited in the PDB.
* `local_residinfo`: For those benchmarks featuring local evaluations, a set of JSON files are provided detailing ranges of residue numbers where (i) structural alignments are performed (`alignment_resid_ranges`), and (ii) the metric of interest is computed (`metric_resid_ranges`).
