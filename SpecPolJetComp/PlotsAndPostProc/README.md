# PlotsAndPostProc

Post-processing and plotting code for the Swift J1727.8-1613 spectropolarimetric analysis. Fair warning: the script names in `code/` got a little... enthusiastic. I apologise, I temporarily LOST MY MIND.

---

## Directory Structure

### `code/`

The plotting scripts. Named with the energy of someone who had been staring at MeerKAT data for too long.

| Script / Notebook | Produces |
|---|---|
| `SUMMARY_PANEL_PLOTTER.py` | Figure 6 |
| `BIG_FIGURE_PLOTTER.py` | Figure 7 |
| `POL_PROPERTIES_PLOTTER.py` | Figures 8–9 and 15 |
| `Method_plots.ipynb` | Methods figures (various) |
| `Results_plots.ipynb` | Results figures (various) |

The notebooks cover the remaining methods and results plots and should be fairly self-explanatory — they're reasonably well-commented considering the circumstances.

`PULL_XF_TO_FILES.py` is a utility for pulling cross-hand phase solutions into the format expected by the plotting scripts.

---

### `files/`

Input data files, organised by type:

- `crosshand_phase_files/` — per-epoch cross-hand phase solutions from 3C286 and J1733-1304
- `LC_files/` — radio light curve data (`SW1727_Radio.csv`) and an extrapolation diagnostic
- `Xray_files/` — MAXI X-ray light curves in the 2–6 keV and 6–20 keV bands

---

### `fits_images/`

MeerKAT MFS Stokes I and linear polarisation FITS images for two epochs (November 2023 and February 2024), used in the imaging figures.

---

### `plots/`

Output directory for generated figures, split into `methods/` and `results/` subdirectories.

---

### `SSA_MC_analysis/`

The Swift J1727.8-1613 tuned synchrotron self-absorption (SSA) Monte Carlo analysis, adapted from [CITE PLACEHOLDER]. The script execution order is:

1. `ssa_tools.py` — defines the SSA physics and sampling utilities
2. `estimate_mass.py` — runs the MC analysis and produces the main output samples
3. `fit_beta_dist.py` *(optional)* — fits a beta distribution to the distance posterior and compares against the [Burridge et al. (2025)](https://arxiv.org/abs/2502.06448) distance prior


> **Note:** If you want to re-run the full pipeline from scratch, you will unfortunately need to re-process the polarisation data first, as the intermediate `.pkl` files are too large to host on GitHub. See [QUfitRMsynth](../QUfitRMsynth) for the relevant processing code.
