# SpecPolJetComp

Spectropolarimetric analysis and jet comparison code for the Swift J1727.8-1613 2023 outburst.

## Citation

If you use this code, please cite:

```bibtex
PLACEHOLDER — citation will be added upon publication
```

---

### [`QUfitRMsynth/`](QUfitRMsynth/)

Broadband QU-fitting and RM synthesis pipeline. Takes the per-epoch MeerKAT polarisation spectra as input, fits phenomenological Faraday rotation models (thin screens, thick screens, power-law screens) via Bayesian nested sampling, and runs RM synthesis and CLEAN for comparison. Contains all fitting code, input data, results, and plots. See the [README](QUfitRMsynth/README.md) and [EXAMPLE](QUfitRMsynth/EXAMPLE.md) in that directory for full details.

---

### [`PlotsAndPostProc/`](PlotsAndPostProc/)

Post-processing and figure generation. Reads the QU-fitting outputs and other multi-wavelength data (radio light curves, X-ray, FITS images) to produce the publication figures and run the SSA Monte Carlo analysis. See the [README](PlotsAndPostProc/README.md) in that directory for a description of the scripts and what each one produces.
