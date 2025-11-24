# BioEmu results

These are the results from BioEmu checkpoint version 1.2. NOTE: the model version associated with our [publication](https://www.science.org/doi/10.1126/science.adv9817) is BioEmu 1.1, whose results are shown in [BIOEMU_RESULTS_v1.1](BIOEMU_RESULTS_v1.1.md).

Results as generated from samples using the public release of [`bioemu`](https://github.com/microsoft/bioemu).
The released code also uses a faster sampler by default (DPM with 30 steps) instead of the sampler used for the paper (Heun sampler with 100 steps).
Here we deposit benchmark results.

## Multiconf benchmarks

### Domain motion

![bioemu-multiconf-domain-coverage](/repo_assets/v1.2/bioemu_domainmotion_rmsd_coverage.png)

![bioemu-multiconf-domain-free](/repo_assets/v1.2/bioemu_domainmotion_rmsd_free_energy.png)

### Cryptic pocket

![bioemu-multiconf-cryptic-coverage-holo](/repo_assets/v1.2/bioemu_crypticpocket_holo_rmsd_coverage.png)

![bioemu-multiconf-cryptic-coverage-apo](/repo_assets/v1.2/bioemu_crypticpocket_apo_rmsd_coverage.png)

![bioemu-multiconf-cryptic-free](/repo_assets/v1.2/bioemu_crypticpocket_rmsd_free_energy.png)

### Local unfolding

![bioemu-multiconf-unfolding-coverage-f](/repo_assets/v1.2/bioemu_localunfolding_fnc_unfold_f_coverage.png)

![bioemu-multiconf-unfolding-coverage-u](/repo_assets/v1.2/bioemu_localunfolding_fnc_unfold_u_coverage.png)

![bioemu-multiconf-unfolding-free-u](/repo_assets/v1.2/bioemu_localunfolding_fnc_unfold_u_free_energy.png)
 
### OOD60

![bioemu-ood60-domain-coverage](/repo_assets/v1.2/bioemu_ood60_rmsd_coverage.png)

![bioemu-ood60-domain-free](/repo_assets/v1.2/bioemu_ood60_rmsd_free_energy.png)

## Folding free energy benchmark

![bioemu-folding-free-energies](/repo_assets/v1.2/bioemu_folding_free_energies_scatter.png)
![bioemu-folding-free-energies-change](/repo_assets/v1.2/bioemu_folding_free_energy_changes_scatter.png)

## CATH MD benchmark

![bioemu-md-emulation](/repo_assets/v1.2/bioemu_md_emulation_mae_free_energy.png)
