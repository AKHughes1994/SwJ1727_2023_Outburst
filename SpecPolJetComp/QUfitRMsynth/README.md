# QUfitRMsynth

Broadband QU-fitting and RM synthesis pipeline for the Swift J1727.8-1613 MeerKAT spectropolarimetric dataset.

> **Note on included files:** Only the results for the favoured models are included in `results/`. Given the number of model types fit per epoch and the size of the individual posterior sample `.pkl` files, including everything was not practical — only the JSON summaries and plots for the favoured models are tracked in the repository. Bash scripts are provided to reprocess everything from scratch (`run_all_model.sh`) or just the favoured models (`run_best_models.sh`). If doing so, **run `faraday_synthesis.py` first** — its RM-CLEAN outputs are used to set the Faraday depth prior range for each epoch before fitting.

The data were reduced and polarisation spectra extracted using [polkat](https://github.com/AKHughes1994/polkat), a MeerKAT-specific polarisation calibration and extraction pipeline. The output per-epoch QU text files live in `files/J1727/QU_text/` and are the starting point for everything here.

---

## Overview

The pipeline decomposes broadband complex polarisation spectra $P(\lambda^2) = Q + iU$ into a sum of phenomenological Faraday components (thin screens, thick screens, power-law screens) using Bayesian nested sampling via [dynesty](https://dynesty.readthedocs.io). Model comparison is then done via the Bayesian evidence to select the simplest model consistent with the data.

For a detailed description of the fitting formalism, prior choices, post-processing, and how to actually run the pipeline end-to-end, see [PLACEHOLDER — fitting pipeline guide] (to be written).

---

## The Five Core Modules (`code/`)

These five files form the library that everything else imports. They are not meant to be run directly. Brief descriptions are given below; the [PLACEHOLDER — fitting pipeline guide] goes into considerably more detail.

| Module | Role |
|---|---|
| `faraday_utils.py` | Physical model components and spectral utilities. Implements thin screens, power-law screens, and thick screens (super-Gaussian approach following Anderson et al. 2016), along with RM synthesis, RMTF computation, and polarization quantity helpers. |
| `faraday_data.py` | Data I/O and Stokes I fitting. Handles loading QU text files (both physical-unit and fractional-polarisation formats), channel flagging, power-law and polynomial total-intensity fits, and data validation. |
| `faraday_model.py` | Bayesian fitting engine. Wraps dynesty nested sampling, builds the likelihood from the composite model, handles prior transforms (including Cartesian ↔ polar conversion and RM ordering constraints), and stores/serialises fit results. |
| `faraday_plot.py` | Publication-quality visualisation. Fit result panels, posterior corner plots, convergence diagnostics, FDF diagnostic plots, and presentation/paper summary figures — all with consistent formatting. |
| `faraday_processing.py` | Posterior post-processing. Detects modes in the (potentially multi-modal) posterior via Gaussian KDE, handles circular/axial angle topology, computes 68% HDIs per mode, and converts Cartesian $(a, b)$ samples to derived polar coordinates $(p_0, \psi_0)$. |

---

## Scripts (`code/`)

### `faraday_obs_fit.py`

The main workhorse script. Invoked from the command line as:

```bash
python faraday_obs_fit.py <MODEL_TYPE> <path/to/QU_text_file.txt>
```

where `MODEL_TYPE` is a string like `S`, `ST`, `SST`, `PP`, etc. (S = thin screen, T = thick screen, P = power-law screen). It loads the epoch data, constructs the fit model and priors, runs the nested sampler, post-processes the posterior, and writes results and plots to the relevant subdirectories.

### `run_all_model.sh`

Calls `faraday_obs_fit.py` for every model type and every epoch — including models that are not ultimately favoured by the evidence. Running the full suite is necessary for model comparison via `compare_models.py`, but fair warning: it will take anywhere from several days to a week on a typical machine. If you have access to a cluster, parallelise across epochs rather than running this sequentially. The script is included as a complete record of what was actually run to produce the published results.

### `run_best_models.sh`

Runs only the model types that were favoured by the Bayesian evidence for each epoch. This is the script to use if you want to regenerate the posterior samples (`.pkl` files) and plots without re-running the full model comparison suite. Much faster than `run_all_model.sh` and sufficient for reproducing the final results.

The fit and corner plots for all favoured models, along with the JSON summary files containing the posterior parameter estimates, are already included in `plots/` and `results/` respectively — so you only need to run this script if you want to regenerate them from scratch or inspect the full posterior samples.

### `faraday_synthesis.py`

Runs RM synthesis and RM-CLEAN on all epochs using the [CIRADA RM-Tools](https://github.com/CIRADA-Tools/RM-Tools) `rmsynth1d` and `rmclean1d` command-line tools. The outputs (dirty and CLEANed FDFs, RMTF) are written to `rmsynth/`. This is independent of the QU-fitting and can be run before or after it.

### `faraday_example_simulate.py`

Generates a synthetic polarised dataset from a user-specified multi-component model, then fits it with the same pipeline used on the real data. Useful for validating the fitting machinery, testing prior sensitivity, and producing example figures. The true model, fit model, and all configuration parameters are defined at the top of the file. All outputs are written to `../simulated/`.

### `compare_models.py`

Reads the per-epoch nested sampling results and applies the model-selection criterion: the simplest model whose log-evidence is within $\Delta \ln Z \leq 10$ of the best-fitting model is selected. For each epoch, it writes a summary text file to `files/J1727/summary/` containing the selected model(s), their evidence, goodness-of-fit statistics, and posterior parameter estimates.

---

## Directory Structure

```
QUfitRMsynth/
├── code/               ← modules, scripts, and bash runner
├── files/
│   ├── J1727/
│   │   ├── QU_text/        ← per-epoch input spectra; fully included in the repo
│   │   ├── polkat_json/    ← polkat intermediate outputs (one compute environment)
│   │   ├── wapiti_output/  ← polkat intermediate outputs (second compute environment)
│   │   └── summary/        ← per-epoch model-selection summary files (from compare_models.py)
│   ├── 3C286/          ← calibrator QU spectra
│   └── J1733/          ← secondary calibrator QU spectra
├── results/
│   └── J1727/          ← per-epoch nested sampling outputs (JSON + pkl)
├── plots/
│   └── J1727/          ← per-epoch diagnostic plots
└── rmsynth/
    └── J1727/          ← RM synthesis and CLEAN outputs
```

`polkat_json/` and `wapiti_output/` are outputs from the same [polkat](https://github.com/AKHughes1994/polkat) pipeline run on two different compute environments (the names reflect the machines, not different software). The `QU_text/` files are the standardised spectra extracted from both and are the only inputs needed to re-run the fitting — they are included in full.

---

## Dependencies

A complete pinned environment is provided in `requirements.txt` (Python 3.12.3). To reproduce it:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Two packages sit outside PyPI and need separate installation:

**[polkat](https://github.com/AKHughes1994/polkat)** — the MeerKAT polarisation calibration and QU spectrum extraction pipeline that produced the input data. Ionospheric RM corrections were derived as part of the polkat pipeline using [Spinifex](https://github.com/lofar-astron/spinifex); the per-epoch outputs are stored in `files/spinifex_iono/`. The extracted `QU_text/` spectra are already included in the repo and are all you need to re-run the fitting. polkat is only relevant if you want to start completely from scratch and re-process the raw MeerKAT visibilities — in that case, refer to the polkat documentation or get in touch with me directly.

**[CIRADA RM-Tools](https://github.com/CIRADA-Tools/RM-Tools)** — provides the `rmsynth1d` and `rmclean1d` command-line tools used by `faraday_synthesis.py`. The package itself (`RM-Tools==1.4.9`) is included in `requirements.txt`, but make sure the CLI tools are on your `PATH` after installation. Only needed if you are re-running RM synthesis.

All other runtime dependencies (`dynesty`, `numpy`, `scipy`, `matplotlib`, `astropy`, `corner`, `numba`, etc.) are covered by `requirements.txt`.
