# BioEmu-1 results

Results as generated from samples using the public release of [`bioemu`](https://github.com/microsoft/bioemu) slightly differ from those in the [preprint](https://www.biorxiv.org/content/10.1101/2024.12.05.626885v1.abstract). The main reason for this is that the public release of `bioemu` uses Evoformer embeddings as computed by [Colabfold](https://github.com/sokrypton/ColabFold) for ease of use, which involve both a different sequence search method (mmseqs2 vs. hhblits) and sequence database (BFD vs. Uniclust30) compared to the original AlphaFold2 code that we had used in the preprint.
The released code also uses a faster sampler by default (DPM with 30 steps) instead of the sampler used for the paper (Heun sampler with 100 steps).
In the interest of transparency, here we deposit benchmark results for our internal model (right) as well as those obtained using the public release of the code (left).

## Multiconf benchmarks

### Domain motion

| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-domain-coverage](/repo_assets/bioemu_domainmotion_rmsd_coverage.png) |  ![paper-multiconf-domain-coverage](/repo_assets/paper_domainmotion_rmsd_coverage.png) |


| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-domain-free](/repo_assets/bioemu_domainmotion_rmsd_free_energy.png) | ![paper-multiconf-domain-coverage](/repo_assets/paper_domainmotion_rmsd_free_energy.png) |



### Cryptic pocket

| BioEmu-1 (holo) | Paper (holo) |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-cryptic-coverage-holo](/repo_assets/bioemu_crypticpocket_holo_rmsd_coverage.png) | ![paper-multiconf-cryptic-coverage-holo](/repo_assets/paper_crypticpocket_holo_rmsd_coverage.png) |


| BioEmu-1 (apo) | Paper (apo) |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-cryptic-coverage-apo](/repo_assets/bioemu_crypticpocket_apo_rmsd_coverage.png) | ![paper-multiconf-cryptic-coverage-apo](/repo_assets/paper_crypticpocket_apo_rmsd_coverage.png) |

| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-cryptic-free](/repo_assets/bioemu_crypticpocket_rmsd_free_energy.png) | ![paper-multiconf-cryptic-free](/repo_assets/paper_crypticpocket_rmsd_free_energy.png) |



### Local unfolding

| BioEmu-1 (folding) | Paper (folding) |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-unfolding-coverage-f](/repo_assets/bioemu_localunfolding_fnc_unfold_f_coverage.png) |  ![paper-multiconf-unfolding-coverage-f](/repo_assets/paper_localunfolding_fnc_unfold_f_coverage.png) |

| BioEmu-1 (unfolding) | Paper (unfolding) |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-unfolding-coverage-u](/repo_assets/bioemu_localunfolding_fnc_unfold_u_coverage.png) |  ![paper-multiconf-unfolding-coverage-u](/repo_assets/paper_localunfolding_fnc_unfold_u_coverage.png) |


| BioEmu-1 (folding) | Paper (folding) |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-unfolding-free-f](/repo_assets/bioemu_localunfolding_fnc_unfold_f_free_energy.png) |   ![paper-multiconf-unfolding-free-f](/repo_assets/paper_localunfolding_fnc_unfold_f_free_energy.png) |


| BioEmu-1 (unfolding) | Paper (unfolding) |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-multiconf-unfolding-free-u](/repo_assets/bioemu_localunfolding_fnc_unfold_u_free_energy.png) |   ![paper-multiconf-unfolding-free-u](/repo_assets/bioemu_localunfolding_fnc_unfold_u_free_energy.png) |

 
### OOD60

| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-ood60-domain-coverage](/repo_assets/bioemu_ood60_rmsd_coverage.png)  |   ![paper-ood60-domain-coverage](/repo_assets/paper_ood60_rmsd_coverage.png) |


| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-ood60-domain-free](/repo_assets/bioemu_ood60_rmsd_free_energy.png)  |    ![paper-ood60-domain-free](/repo_assets/bioemu_ood60_rmsd_free_energy.png) |



## Folding free energy benchmark

| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-folding-free-energies](/repo_assets/bioemu_folding_free_energies_scatter.png)  |   ![paper-folding-free-energies](/repo_assets/paper_folding_free_energies_scatter.png) |

## CATH MD benchmark

| BioEmu-1 | Paper |
| ----------------------------------- | ----------------------------------- |
| ![bioemu-md-emulation](/repo_assets/bioemu_md_emulation_mae_free_energy.png)  |   ![paper-md-emulation](/repo_assets/paper_md_emulation_mae_free_energy.png) |
