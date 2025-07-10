# BioEmu results

Results as generated from samples using the public release of [`bioemu`](https://github.com/microsoft/bioemu).
The released code also uses a faster sampler by default (DPM with 30 steps) instead of the sampler used for the paper (Heun sampler with 100 steps).
Here we deposit benchmark results.

## Multiconf benchmarks

### Domain motion

![bioemu-multiconf-domain-coverage](/repo_assets/bioemu_domainmotion_rmsd_coverage.png)

![bioemu-multiconf-domain-free](/repo_assets/bioemu_domainmotion_rmsd_free_energy.png)

### Cryptic pocket

![bioemu-multiconf-cryptic-coverage-holo](/repo_assets/bioemu_crypticpocket_holo_rmsd_coverage.png)

![bioemu-multiconf-cryptic-coverage-apo](/repo_assets/bioemu_crypticpocket_apo_rmsd_coverage.png)

![bioemu-multiconf-cryptic-free](/repo_assets/bioemu_crypticpocket_rmsd_free_energy.png)

### Local unfolding

![bioemu-multiconf-unfolding-coverage-f](/repo_assets/bioemu_localunfolding_fnc_unfold_f_coverage.png)

![bioemu-multiconf-unfolding-coverage-u](/repo_assets/bioemu_localunfolding_fnc_unfold_u_coverage.png)

![bioemu-multiconf-unfolding-free-u](/repo_assets/bioemu_localunfolding_fnc_unfold_u_free_energy.png)
 
### OOD60

![bioemu-ood60-domain-coverage](/repo_assets/bioemu_ood60_rmsd_coverage.png)

![bioemu-ood60-domain-free](/repo_assets/bioemu_ood60_rmsd_free_energy.png)

## Folding free energy benchmark

![bioemu-folding-free-energies](/repo_assets/bioemu_folding_free_energies_scatter.png)
![bioemu-folding-free-energies-change](/repo_assets/bioemu_folding_free_energy_changes_scatter.png)

## CATH MD benchmark

![bioemu-md-emulation](/repo_assets/bioemu_md_emulation_mae_free_energy.png)
