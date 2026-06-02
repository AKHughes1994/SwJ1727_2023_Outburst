# QU-Fitting: Worked Example

A step-by-step guide to loading data, building a model, specifying priors, and running a Bayesian nested sampling fit. This document is intended as a practical reference; a more detailed treatment of the formalism, flagging options, post-processing, and model comparison is in [PLACEHOLDER — fitting pipeline guide].

All examples assume you are running from `code/` and that the five core modules (`faraday_data`, `faraday_utils`, `faraday_model`, `faraday_plot`, `faraday_processing`) are importable from the same directory.

---

## 1. Loading Data

Data files are plain text with columns `freq[Hz]  I  Q  U  dI  dQ  dU` in physical flux units (e.g. mJy — the units do not matter as long as they are consistent, since the loader fits Stokes I and divides Q and U by it to produce a fractional polarisation spectrum).

```python
from faraday_data import DataLoader

data, I_fitter = DataLoader.load_physical_with_fit(
    '../files/J1727/QU_text/SwiftJ1727_WAPITI_20231006.txt',
    I_model='power-law',   # fit Stokes I with a power-law I ∝ λ^α
    poly_degree=4,         # only used if I_model='polynomial'
    central_freq=1.28,     # reference frequency in GHz (MeerKAT L-band centre)
)
```

`data` is a `PolarizationData` object containing fractional `Q`, `U` and their propagated errors. `I_fitter` holds the Stokes I fit parameters.

If you want to pass a channel flagging configuration, see `FlaggingConfig` in `faraday_data.py` and `flagging_config.toml` for the per-epoch overrides actually used in this paper — this is handled automatically by `faraday_obs_fit.py` and you normally do not need to configure it by hand.

---

## 2. Building a Model

A `CompositeModel` is a sum of one or more Faraday components. Three component types are available:

| Class | Symbol | Description |
|---|---|---|
| `ThinComponent` | `S` | Thin Faraday screen. Two free amplitudes `(a, b)` and a rotation measure `φ_rm`. |
| `ThickComponent` | `T` | Thick screen / internal Faraday dispersion (super-Gaussian FDF). Amplitudes `(a, b)`, peak depth `φ_peak`, dispersion `σ_φ`, and number of cells `N`. |
| `ThinPowerLawComponent` | `P` | Thin screen with wavelength-dependent depolarisation. Adds spectral index `β` and a reference wavelength `λ²_0`. |

The initial parameter values passed to each component are arbitrary when fitting real data — dynesty samples the full prior volume from scratch and ignores them entirely. The names you assign (`S1`, `T1`, `P1`, etc.) are used as keys in the prior dictionary.

```python
from faraday_utils import CompositeModel, ThinComponent, ThickComponent, ThinPowerLawComponent

model = CompositeModel()

# Add a thin Faraday screen (e.g., the ISM foreground)
model.add_component(ThinComponent(a=0.1, b=0.0, phi_rm=0.0, name='S1'))

# Add a thick screen (e.g., emission from within the jet)
model.add_component(ThickComponent(a=0.1, b=0.0, phi_peak=0.0,
                                   sigma_phi=1.0, N=2.0, name='T1'))
```

For a pure power-law model:

```python
model = CompositeModel()
model.add_component(ThinPowerLawComponent(a=0.1, b=0.0, phi_rm=0.0,
                                          beta=0.0, lambda_sq_0=0.09, name='P1'))
```

---

## 3. Defining Priors

Priors are passed as a plain Python dictionary. The top-level keys are the component names you assigned above. Each component entry is itself a dict specifying the prior on each free parameter.

### 3.1 Available prior types

Scalar parameters (`phi_rm`, `phi_peak`, `sigma_phi`, `N`, `beta`, `a`, `b`) accept a tuple:

| Syntax | Description |
|---|---|
| `('uniform', min, max)` | Uniform between `min` and `max` |
| `('log_uniform', min, max)` | Log-uniform (Jeffreys prior); both bounds must be > 0 |
| `('gaussian', mean, std)` | Unbounded Gaussian |
| `('truncated_gaussian', mean, std, min, max)` | Gaussian truncated to `[min, max]` |
| `('fixed', value)` | Fix parameter to `value` (removes it from the fit entirely) |

### 3.2 Amplitude priors via `polar_conversion`

Rather than placing priors directly on the Cartesian amplitudes `(a, b)`, it is far more natural to work in polar coordinates `(p0, ψ₀)` where `p0 = √(a²+b²)` is the fractional polarisation amplitude and `ψ₀` is the intrinsic EVPA. This is controlled by the `polar_conversion` sub-dict:

```python
'polar_conversion': {
    'p0_bounds':   (p0_min, p0_max, p0_distribution),
    'psi_0_bounds': (psi_min, psi_max, angle_prior),
}
```

**`p0_distribution`** options:

| Value | Description |
|---|---|
| `'radial_uniform'` | Uniform in `p0` — i.e., `p(p0) ∝ 1`. Recommended: physically motivated and avoids the ring bias of uniform Cartesian priors. |
| `'area_uniform'` | Uniform in area of the Poincaré disc — i.e., `p(p0) ∝ p0`. Equivalent to uniform Cartesian `(a, b)`. |
| `'log_uniform'` | Log-uniform in `p0`; useful if amplitudes span orders of magnitude. |

**`angle_prior`** options:

| Value | Description |
|---|---|
| `'uniform'` | Flat over `[-π/2, π/2]`. Use for sources where you have no prior knowledge of the EVPA. |
| `'peaked'` | Weakly informative prior that peaks at `peaked_angle_deg` (default 0°) and its orthogonal direction (90° away), with concentration `peaked_k` (default 1.0). Recommended for sources where the intrinsic EVPA is unlikely to be wildly large; reduces the sampler's tendency to explore unphysical EVPA solutions. |

### 3.3 Recommended priors (as used in this paper)

These are taken directly from `faraday_obs_fit.py` and reflect the priors actually used to produce the published results. The `phi_rm` / `phi_peak` range `[phi_min, phi_max]` is derived automatically from the RMclean output for each epoch (default fallback `±160 rad/m²`).

```python
import numpy as np

phi_min = -160.0   # rad/m² — set per-epoch from RMclean output in practice
phi_max =  160.0

THIN_PRIOR = {
    'type': 'ThinComponent',
    'polar_conversion': {
        'p0_bounds':    (0.0, 0.50, 'radial_uniform'),
        'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked'),
        'peaked_k':          1.0,   # concentration parameter
        'peaked_angle_deg':  0.0,   # peak direction (degrees)
    },
    'phi_rm': ('uniform', phi_min, phi_max),
}

THICK_PRIOR = {
    'type': 'ThickComponent',
    'polar_conversion': {
        'p0_bounds':    (0.0, 0.50, 'radial_uniform'),
        'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked'),
    },
    'phi_peak':  ('uniform',     phi_min, phi_max),
    'sigma_phi': ('log_uniform', 0.1,     50.0),    # rad/m²; upper bound ~ 2× RMTF FWHM
    'N':         ('log_uniform', 2.0,     50.0),    # number of Faraday cells
    'effective_width_mode': False,
}

POWERLAW_PRIOR = {
    'type': 'ThinPowerLawComponent',
    'polar_conversion': {
        'p0_bounds':    (0.0, 0.50, 'radial_uniform'),
        'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked'),
    },
    'phi_rm': ('uniform', phi_min, phi_max),
    'beta':   ('uniform', -5.0,    5.0),
}
```

### 3.4 ISM RM constraint

For models containing a thin screen, the first `S` component is constrained to the expected Galactic ISM foreground RM using a truncated Gaussian. The values below are those used for Swift J1727; adjust `PHI_ISM` and `PHI_SYS` for your source.

```python
PHI_ISM   = -0.5   # Galactic ISM Faraday depth (rad/m²)
PHI_SYS   = -0.3   # Systematic instrumental RM offset (rad/m²)
iono_rm   =  0.0   # Ionospheric contribution from Spinifex; set per-epoch

phi_center = PHI_ISM + PHI_SYS + iono_rm

ism_phi_prior = ('truncated_gaussian',
                 phi_center,          # mean
                 0.65,                # std (rad/m²)
                 phi_center - 2.3,    # lower truncation (≈ 3.5σ)
                 phi_center + 2.3)    # upper truncation

priors['S1']['phi_rm'] = ism_phi_prior
```

### 3.5 Tying `phi_rm` across components

For all-power-law models (`PP`, `PPP`, …), it is often appropriate to tie the rotation measure across components — i.e., all components share a single `φ_rm` with separate amplitude and spectral index. Enable this via:

```python
priors = {'tied_phi_rm': ('uniform', phi_min, phi_max)}
# ... then add each component prior WITHOUT a 'phi_rm' key
```

and pass `tie_phi_rm=True` to `FitSetup` (see below). `faraday_obs_fit.py` auto-enables this for any model consisting entirely of `P` components.

### 3.6 Assembling the prior dict

```python
priors = {}
priors['S1'] = THIN_PRIOR.copy()
priors['T1'] = THICK_PRIOR.copy()

# Optionally override one parameter after the fact:
priors['S1']['phi_rm'] = ism_phi_prior
```

---

## 4. FitSetup

`FitSetup` validates that the model, priors, and data are mutually consistent and pre-computes the prior transform used by dynesty.

```python
from faraday_model import FitSetup

# lambda_sq_ref: the λ² value at which p0 is defined (i.e., the reference band).
# Using the band centre is standard; here 0.05 m² ≈ 1.28 GHz.
lambda_sq_ref = (0.2998 / 1.28) ** 2   # ≈ 0.0549 m²

setup = FitSetup(
    model=model,
    priors=priors,
    data=data,
    enforce_ordering=True,   # enforce φ_rm ordering between same-type components
    tie_phi_rm=False,        # set True for all-P models
    lambda_sq_ref=lambda_sq_ref,
    p0_physical_max=0.7,     # hard upper bound on fractional polarisation
)
```

`enforce_ordering=True` (recommended) imposes `φ_rm(S1) < φ_rm(S2) < …` (and similarly for `T`, `P` components) to break the label-switching degeneracy when multiple components of the same type are present. This is implemented via order-statistics in the prior transform, so there is no likelihood penalty.

---

## 5. Running the Fit

`FaradayFitter` wraps dynesty. The `dynesty_kwargs` argument is **mandatory** — no defaults are provided to ensure you are always explicit about the sampler settings.

### 5.1 Parallelisation

By default (`use_pool=True`, `n_cores=None`), the fitter uses all available CPUs minus two. Override with `n_cores=N` if needed.

### 5.2 Static sampler (recommended for most single models)

```python
from faraday_model import FaradayFitter

dynesty_kwargs = {
    'sampler': {
        'nlive':  1000,         # number of live points
        'bound':  'multi',      # bounding method ('multi' = multi-ellipsoid)
        'sample': 'rwalk',      # sampling method
        'walks':  150,          # steps per random walk
    },
    'run': {
        'dlogz':   0.1,         # stopping criterion on log-evidence
        'maxiter': None,        # no iteration cap
        'maxcall': None,        # no likelihood call cap
    },
}

fitter = FaradayFitter(setup, sampler='static', dynesty_kwargs=dynesty_kwargs)
results = fitter.fit()
```

### 5.3 Dynamic sampler (used in this paper)

The dynamic sampler adaptively allocates live points to maximise posterior accuracy. It is slower per iteration but more efficient for complex posteriors. This is what `faraday_obs_fit.py` uses:

```python
dynesty_kwargs = {
    'sampler': {
        'bound':  'multi',
        'sample': 'rwalk',
        'walks':  150,
    },
    'run': {
        'dlogz_init':   0.01,
        'nlive_init':   450,      # = N_LIVE / 2
        'nlive_batch':  150,      # = N_LIVE / 6
        'n_effective':  10_000,
        'maxiter':      None,
        'maxcall':      None,
    },
}

fitter = FaradayFitter(setup, sampler='dynamic', dynesty_kwargs=dynesty_kwargs,
                       n_cores=15)
results = fitter.fit()
```

### 5.4 Saving results

```python
results.save('my_results.pkl')           # full pickle (samples, weights, etc.)
results.save_json('my_results.json')     # human-readable summary
```

---

## 6. Complete Worked Example

A self-contained script for a single epoch, fitting a thin screen + thick screen (`ST`) model:

```python
import numpy as np
from faraday_data import DataLoader
from faraday_utils import CompositeModel, ThinComponent, ThickComponent
from faraday_model import FitSetup, FaradayFitter
from faraday_processing import process_posterior_modes, ProcessingConfig
from faraday_plot import plot_fit_results, plot_corner
from pathlib import Path

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
data_file  = '../files/J1727/QU_text/SwiftJ1727_WAPITI_20231006.txt'
output_dir = Path('../results/J1727/example_ST')
output_dir.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------
data, I_fitter = DataLoader.load_physical_with_fit(
    data_file,
    I_model='power-law',
    central_freq=1.28,
)

# -------------------------------------------------------------------------
# Build model
# -------------------------------------------------------------------------
model = CompositeModel()
model.add_component(ThinComponent(a=0.1, b=0.0, phi_rm=0.0,   name='S1'))
model.add_component(ThickComponent(a=0.1, b=0.0, phi_peak=0.0,
                                   sigma_phi=1.0, N=2.0,       name='T1'))

# -------------------------------------------------------------------------
# Priors
# -------------------------------------------------------------------------
phi_min, phi_max = -160.0, 160.0
lambda_sq_ref    = (0.2998 / 1.28) ** 2

priors = {
    'S1': {
        'type': 'ThinComponent',
        'polar_conversion': {
            'p0_bounds':    (0.0, 0.50, 'radial_uniform'),
            'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked'),
        },
        # ISM RM constraint: truncated Gaussian around known foreground
        'phi_rm': ('truncated_gaussian', -0.8, 0.65, -3.1, 1.5),
    },
    'T1': {
        'type': 'ThickComponent',
        'polar_conversion': {
            'p0_bounds':    (0.0, 0.50, 'radial_uniform'),
            'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked'),
        },
        'phi_peak':  ('uniform',     phi_min, phi_max),
        'sigma_phi': ('log_uniform', 0.1,     50.0),
        'N':         ('log_uniform', 2.0,     50.0),
        'effective_width_mode': False,
    },
}

# -------------------------------------------------------------------------
# Setup + fit
# -------------------------------------------------------------------------
setup = FitSetup(model, priors, data,
                 enforce_ordering=True,
                 lambda_sq_ref=lambda_sq_ref)

dynesty_kwargs = {
    'sampler': {'nlive': 1000, 'bound': 'multi', 'sample': 'rwalk', 'walks': 150},
    'run':     {'dlogz': 0.1, 'maxiter': None, 'maxcall': None},
}
fitter  = FaradayFitter(setup, sampler='static', dynesty_kwargs=dynesty_kwargs)
results = fitter.fit()

# -------------------------------------------------------------------------
# Post-process and save
# -------------------------------------------------------------------------
pconfig = ProcessingConfig(
    filename_prefix='example_ST',
    detect_multiple_modes=True,
    plot_diagnostics=True,
    plot_dir=output_dir / 'modes',
)
results = process_posterior_modes(results, pconfig)
results.save(output_dir / 'example_ST.pkl')
results.save_json(output_dir / 'example_ST.json')
results.print_summary()

# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------
plot_fit_results(results, I_fitter=I_fitter,
                 filename='example_ST_fit.png', output_dir=output_dir)
plot_corner(samples=results.samples_processed,
            param_names=results.param_names_processed,
            param_summary=results.param_summary,
            processed_parameters=results.processed_parameters,
            weights=results.weights,
            setup=results.setup,
            filename='example_ST_corner.png', output_dir=output_dir)
```

For running every model type across every epoch (for full model comparison), see `run_all_model.sh`. To regenerate only the favoured models and their posterior samples, use `run_best_models.sh`.

---

## 7. Simulating Data

The pipeline includes a self-contained simulation script that generates a synthetic polarised dataset from a user-defined model, fits it, and produces the full suite of diagnostic plots. This is useful for validating the fitting machinery, testing prior sensitivity, and producing demonstration figures.

The script is `code/faraday_example_simulate.py`. All configuration lives at the top of the file — you do not need to touch anything below the configuration block to run a simulation.

### 7.1 Defining the true model

Components are specified as lists of parameter dictionaries. Add, remove, or comment out entries to change the model:

```python
# Thin screen components (S1, S2, ...)
THIN_COMPONENTS = [
    {'p0': 0.02, 'psi_0_deg': -2.6, 'phi_rm': -5.76},   # S1
]

# Power-law components (P1, P2, ...)
POWERLAW_COMPONENTS = []  # empty = no power-law components

# Thick screen components (T1, T2, ...)
THICK_COMPONENTS = [
    {'p0': 0.05, 'psi_0_deg': -35.0, 'phi_peak': 53.0, 'sigma_phi': 35.0, 'N': 40},  # T1
]
```

### 7.2 Choosing what model to fit

```python
MODEL_TO_FIT = None   # None = fit the exact same model as the truth
                      # or set e.g. 'ST', 'S', 'T' to fit a different model
```

Setting `MODEL_TO_FIT` to something other than `None` lets you test whether the sampler can recover the correct model even when fitted with a misspecified one.

### 7.3 Other key settings

```python
QUICK_MODE = False   # True = skip corner/convergence/FDF plots (faster for quick checks)
```

Data generation parameters (frequency range, number of channels, SNR, spectral index) are set inside `run_test_scenario()` at the clearly marked `SECTION 2` block.

### 7.4 Running the script

```bash
cd code/
python faraday_example_simulate.py
```

All outputs (synthetic data file, fit results, plots) are written to `../simulated/`.

---

## 8. Caveats and Implementation Notes

### 8.1 The reference wavelength `lambda_sq_ref`

The Cartesian sampler parameters `(a, b)` are the intrinsic complex polarisation amplitude — i.e., the amplitude before any Faraday rotation or depolarisation is applied. The full model is:

$$P(\lambda^2) = (a + ib) \times K(\lambda^2)$$

where $K$ encodes the Faraday rotation and depolarisation of each component type.

When `lambda_sq_ref` is set in `FitSetup`, the parameters are instead defined at the reference band: `(a, b)` become the *observed* complex amplitude at $\lambda^2_\text{ref}$, so that

$$a_\text{ref} + ib_\text{ref} = P(\lambda^2_\text{ref}) = (a_\text{int} + ib_\text{int}) \times K(\lambda^2_\text{ref})$$

and the intrinsic values are recovered internally via

$$a_\text{int} + ib_\text{int} = \frac{a_\text{ref} + ib_\text{ref}}{K(\lambda^2_\text{ref})}$$

The transfer factor $K$ differs by component type:

- **Thin screen:** $K(\lambda^2_\text{ref}) = e^{2i\phi_\text{rm}\lambda^2_\text{ref}}$ — a pure phase, so $|K| = 1$ and $p_{0,\text{ref}} = p_{0,\text{int}}$. The reference band has no effect on the amplitude prior for thin components.

- **Thick screen:** $K(\lambda^2_\text{ref}) = \int F_\text{norm}(\phi)\, e^{2i\phi\lambda^2_\text{ref}}\, d\phi$ — complex with $|K| < 1$ due to depolarisation at the reference band. Here $p_{0,\text{ref}} < p_{0,\text{int}}$, and the prior is placed on the physically observable quantity (the amplitude at the reference band) rather than the depolarisation-corrected intrinsic amplitude.

Using `lambda_sq_ref` is recommended whenever the dataset contains thick components, as it makes the prior on `p0` directly correspond to what is measured at the band centre. Set it to the band-centre wavelength squared — for MeerKAT L-band centred at 1.28 GHz this is `(0.2998 / 1.28)² ≈ 0.055 m²`.

### 8.2 Thick component numerical solver

For $N \neq 2$ (i.e., non-Gaussian FDFs), $P(\lambda^2)$ is computed by numerically Fourier-transforming the super-Gaussian FDF over a finite Faraday depth grid. The grid extent and point count are set by power-law fits calibrated against a sweep over $\sigma_\phi = 50\ \text{rad/m}^2$ in the MeerKAT L-band window (0.85–1.8 GHz); see `faraday_utils.py` lines 461–472:

```python
# Grid parameters from swept calibration on σ_φ=50 rad/m², band 0.85–1.8 GHz.
# Fit model: f(N) = a * N^b + c  (pinned through N=2, N_PHI inflated 20%).
# To recalibrate, rerun optimize_num_int.py and update the six constants below.
_nsig_a, _nsig_b, _nsig_c = 19.556, -2.190, 1.316
_nphi_a, _nphi_b, _nphi_c = 7252.291, -4.310, 55.407

n_sigma   = _nsig_a * self.N ** _nsig_b + _nsig_c   # grid half-width in σ_φ units
phi_range = n_sigma * self.sigma_phi
n_phi     = ...                                       # number of grid points
```

These constants are tuned for the Faraday thickness values and frequency coverage of this dataset. If you are working with significantly larger $\sigma_\phi$ values or a different frequency band, the grid may be too coarse to maintain the same numerical precision — the integrand becomes more oscillatory at large Faraday depths or at shorter wavelengths. In that case you will need to recalibrate by rerunning `optimize_num_int.py` and updating the six constants.

I plan to extend the pipeline to MeerKAT L+S-band data and will update the grid calibration accordingly when those data are in hand. If you are applying this code outside the L-band regime in the meantime, treat thick component results with some caution.
