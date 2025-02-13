# Definitions for the MD emulation benchmark v0.1

The MD emulation benchmark has the following substructure:
* `indexed_samples.csv`: contains sequence information for each evaluated test case identifier. 
* `projections_mean.npz` and `projections_sqrt_inv_cov.npz`: contain parameters for computing TICA projections based on sample feature coordinates for each system indexed by test case ID.
* `reference_projections.npz`: contains reference projections computed for the original MD data.