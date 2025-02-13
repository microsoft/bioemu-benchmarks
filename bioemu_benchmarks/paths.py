import os

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

ASSETS_DIR = os.path.join(ROOT_DIR, "assets")


MD_EMULATION_ASSET_DIR = os.path.join(ASSETS_DIR, "md_emulation_benchmark_0.1")
FOLDING_FREE_ENERGY_ASSET_DIR = os.path.join(ASSETS_DIR, "folding_free_energies_benchmark_0.1")
MULTICONF_ASSET_DIR = os.path.join(ASSETS_DIR, "multiconf_benchmark_0.1")
UTILS_DIR = os.path.join(os.path.expanduser("~"), ".bioemubenchmarks")
