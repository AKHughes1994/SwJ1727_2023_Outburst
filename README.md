# Description

This repository contains the analysis I led (as part of a large collaboration) for the 2023-2024+ outburst of the black hole X-ray binary Swift J1727.8-1613. Each subdirectory contains data and analysis scripts; if the data are used, please cite the comprising papers listed within each subdirectory.

## Subdirectories:

---

### LightCurves

This directory contains data compiled and presented in Hughes et al. 2025a. The dataset includes over 200 epochs of integrated radio flux densities from seven facilities, with radio frequencies ranging from 0.3-300 GHz.

**Main components:**

1. **Radio data**: Available in CSV format at `files/SWJ1727_Radio.csv`.
2. **Plotting and analysis scripts**: Included in the provided Jupyter notebook.
3. **X-ray data**: The MAXI/GSC X-ray data required to create the Hardness-Intensity diagram are available in the `files` directory.

---

### RadioXrayCorrelation

This directory contains data and analysis for the radio-X-ray correlation ("LrLx"), as reported in Hughes et al. 2025b.

**Components:**

1. **X-ray analysis pipeline** (optional): 
   - All necessary files to reproduce the figures/analysis are located in the `SWJ1727_Files` directory. However, for those interested in re-running the pipeline, follow these steps:
        a. Run `get_xrt_from_pipeline.py` (requires the Swift-tools Python API and HEASoft).
        b. Execute `fit_xrt_spectra_rise.py` to fit the rise of the outburst.
        c. Execute `fit_xrt_spectra_decay.py` to fit the decay of the outburst.
        d. Run `output_xrt_txt.py` to convert JSON outputs from steps b/c into `.txt` files used in the notebook.
2. **X-ray and radio data**: Available in the `SWJ1727_Files` directory.
3. **Analysis and plotting scripts**: Included in the provided Jupyter notebook(s).

The notebooks are designed to be self-explanatory. However, if anything is unclear, please leave a note or start a discussion. The radio data is derived from the LightCurves directory.

---

### SpecPolJetComp

This directory contains data and analysis for the time-domain radio spectropolarimetry of Swift J1727.8-1613, as reported in Hughes et al. 2026a. The dataset is used to identify evolving plasma conditions in the jet through time-resolved measurements of the radio spectrum and linear polarization properties.

**Main components:**

1. **QUfitRMsynth**: Scripts and data for QU-fitting and RM synthesis analysis of the polarized radio emission.
2. **PlotsAndPostProc**: Post-processing scripts and plotting routines for generating the figures presented in the paper.
