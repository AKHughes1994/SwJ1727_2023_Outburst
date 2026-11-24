# SpecPolJetComp

Spectropolarimetric analysis and jet comparison code for the Swift J1727.8-1613 2023 outburst.

## Citation

If you use this code, please cite:

```bibtex
@ARTICLE{2026MNRAS.551g1283H,
       author = {{Hughes}, A.~K. and {Fender}, R.~P. and {Sivakoff}, G.~R. and {Cowie}, F.~J. and {Heywood}, I. and {Matthews}, J.~H. and {Savard}, K. and {Carotenuto}, F. and {Russell}, T.~D. and {Wood}, C.~M. and et al.},
        title = "{Spectropolarimetric detection of baryonic mass loading in a transient relativistic jet: application to the black hole X-ray binary Swift J1727.8{\ensuremath{-}}1613}",
      journal = {\mnras},
     keywords = {black hole physics, polarization, ISM: jets and outflows, radio continuum: stars, X-rays: binaries, High Energy Astrophysical Phenomena},
         year = 2026,
        month = sep,
       volume = {551},
       number = {1},
          eid = {stag1283},
        pages = {stag1283},
          doi = {10.1093/mnras/stag1283},
archivePrefix = {arXiv},
       eprint = {2607.05182},
 primaryClass = {astro-ph.HE},
       adsurl = {https://scixplorer.org/abs/2026MNRAS.551g1283H},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

---

### [`QUfitRMsynth/`](QUfitRMsynth/)

Broadband QU-fitting and RM synthesis pipeline. Takes the per-epoch MeerKAT polarisation spectra as input, fits phenomenological Faraday rotation models (thin screens, thick screens, power-law screens) via Bayesian nested sampling, and runs RM synthesis and CLEAN for comparison. Contains all fitting code, input data, results, and plots. See the [README](QUfitRMsynth/README.md) and [EXAMPLE](QUfitRMsynth/EXAMPLE.md) in that directory for full details.

> **Using this on your own sources?** `QUfitRMsynth/` is the frozen record of the Swift J1727.8-1613 analysis. An updated, blank version of the pipeline, without the J1727-specific data and easier to tailor to your own sources, is available at [specpol-template](https://github.com/AKHughes1994/specpol-template).

---

### [`PlotsAndPostProc/`](PlotsAndPostProc/)

Post-processing and figure generation. Reads the QU-fitting outputs and other multi-wavelength data (radio light curves, X-ray, FITS images) to produce the publication figures and run the SSA Monte Carlo analysis. See the [README](PlotsAndPostProc/README.md) in that directory for a description of the scripts and what each one produces.
