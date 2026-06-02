# Set matplotlib to use non-GUI backend BEFORE any imports
import matplotlib
matplotlib.use('Agg')

import re
import sys
import numpy as np
from pathlib import Path
from astropy.time import Time
try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        raise ImportError("Python ≥3.11 required for tomllib, or install tomli: pip install tomli")
from faraday_utils import (ThinComponent, ThickComponent, ThinPowerLawComponent, CompositeModel, 
                           create_lambda_sq_grid)
from faraday_model import FitSetup, FaradayFitter, load_results
from faraday_plot import (plot_data_diagnostic, plot_fit_results, plot_corner,
                          plot_convergence, plot_fdf_diagnostic, plot_flagging_diagnostic, plot_thick_component_band_response)
from faraday_data import (DataLoader, IFitter, PolarizationData, FlaggingConfig,
                          save_flagged_data)
from faraday_processing import (process_posterior_modes, ProcessingConfig)


def load_iono_rm(iono_loc, source_name, data_file=None):
    """
    Parse a spinifex output file and return the weighted-mean ionospheric RM
    for the given source (case-insensitive substring match on the first column).
    Last two columns of each matching row must be (iono_rm, iono_err) in rad/m².

    iono_loc : path to a single spinifex file OR a directory of spinifex files.
               When a directory is given, the file whose leading Unix-timestamp
               stem is closest in MJD to the date encoded in data_file (YYYYMMDD)
               is selected; raises if no match falls within 1 day.
    data_file: required when iono_loc is a directory.
    """
    iono_path = Path(iono_loc)

    if iono_path.is_dir():
        if data_file is None:
            raise ValueError("load_iono_rm: data_file required when iono_loc is a directory")

        date_match = re.search(r'(\d{8})', Path(data_file).stem)
        if not date_match:
            raise ValueError(
                f"Cannot extract YYYYMMDD date from data filename: {data_file}"
            )
        ds = date_match.group(1)
        data_mjd = Time(f"{ds[:4]}-{ds[4:6]}-{ds[6:8]}T00:00:00", format='isot',
                        scale='utc').mjd
        print(f"CONSTRAIN_ISM_RM: data date {ds[:4]}-{ds[4:6]}-{ds[6:8]} (MJD {data_mjd:.3f})")

        best_file, best_diff, best_file_mjd = None, np.inf, None
        for f in sorted(iono_path.glob('*.txt')):
            ts_match = re.match(r'^(\d+)', f.stem)
            if not ts_match:
                continue
            file_mjd = Time(float(ts_match.group(1)), format='unix', scale='utc').mjd
            diff = abs(file_mjd - data_mjd)
            if diff < best_diff:
                best_diff = diff
                best_file = f
                best_file_mjd = file_mjd

        if best_file is None:
            raise ValueError(f"No spinifex files found in {iono_loc}")
        if best_diff > 1.0:
            raise ValueError(
                f"Closest spinifex file is {best_diff:.2f} days from data date "
                f"(MJD {data_mjd:.3f}): {best_file}"
            )
        matched_date = Time(best_file_mjd, format='mjd', scale='utc').to_value('isot')[:10]
        print(f"CONSTRAIN_ISM_RM: matched {best_file.name}  "
              f"[{matched_date}, MJD {best_file_mjd:.3f}, Δ{best_diff:.3f} d]")
        iono_path = best_file

    source_lower = source_name.lower()
    matches = []
    with open(iono_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if source_lower in parts[0].lower():
                try:
                    matches.append((float(parts[-2]), float(parts[-1])))
                except (ValueError, IndexError):
                    pass

    if not matches:
        raise ValueError(
            f"No ionospheric RM found for source matching '{source_name}' in {iono_path}"
        )
    rms  = np.array([m[0] for m in matches])
    errs = np.array([m[1] for m in matches])
    w    = 1.0 / errs**2
    return float(np.sum(w * rms) / np.sum(w)), float(1.0 / np.sqrt(np.sum(w)))


def phi_range_from_rmclean(data_file_arg, threshold_frac, padding, fallback_phi_center=None):
    """
    Derive phi_min/phi_max from non-negligible CLEAN components in the RMSynthesis
    FDFmodel.dat output.  Mirrors the files→rmsynth path substitution used elsewhere.

    Returns (phi_min, phi_max) with `padding` added on each side.
    Returns None if the model file cannot be found.
    If the file exists but all components are zero, falls back to
    (fallback_phi_center - padding, fallback_phi_center + padding) when provided,
    otherwise returns None.
    """
    data_path = Path(data_file_arg).resolve()
    data_id_local = data_path.stem

    if 'files' not in data_path.parts:
        return None

    files_idx = data_path.parts.index('files')
    base_dir = Path(*data_path.parts[:files_idx])
    subdir_parts = [p for p in data_path.parts[files_idx + 1:-1] if p != 'QU_text']
    subdir = Path(*subdir_parts) if subdir_parts else Path('.')

    rmclean_file = (base_dir / 'rmsynth' / subdir / data_id_local
                    / f"{data_id_local}_rmsynth_FDFmodel.dat")

    if not rmclean_file.exists():
        return None

    raw = np.loadtxt(rmclean_file)
    phi_arr = raw[:, 0]
    amp = np.sqrt(raw[:, 1]**2 + raw[:, 2]**2)

    max_amp = np.max(amp)
    has_components = max_amp > 0 and np.any(amp > threshold_frac * max_amp)

    if not has_components:
        if fallback_phi_center is not None:
            return fallback_phi_center - padding, fallback_phi_center + padding
        return None

    mask = amp > threshold_frac * max_amp
    return float(np.min(phi_arr[mask])) - padding, float(np.max(phi_arr[mask])) + padding


# =============================================================================
# CONFIGURATION
# =============================================================================
# Command line: python3 faraday_obs_fit.py [MODEL_TYPE] [DATA_FILE]
# MODEL_TYPE examples: 'S', 'SS', 'ST', 'PP', 'SPT', 'SSTT', etc.
# DATA_LOC_DEFAULT = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/J1727/QU_text/SwiftJ1727_WAPITI_20230906.txt'
DATA_LOC_DEFAULT = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/J1727/QU_text/SwiftJ1727_QU_202309016.txt'
MODEL_TYPE = sys.argv[1].upper() if len(sys.argv) > 1 else 'SS'
DATA_FILE_ARG = sys.argv[2] if len(sys.argv) > 2 else DATA_LOC_DEFAULT


MASTER_SEED = 25031995  # big day
_seed_seq = np.random.SeedSequence(MASTER_SEED)
rng = np.random.default_rng(MASTER_SEED)  # overridden per trial in the fit loop

N_CORES = 15
N_TRIALS = 1   # Repeat the fit N times with independent posterior draws; 1 = normal single run

TIE_RM = False
SKIP_IF_EXISTS = True  # If True, skip fitting and load existing .pkl if it exists

USE_RMCLEAN_PHI_RANGE = True   # Derive phi_min/phi_max from RMSynthesis CLEAN components
REQUIRE_RMCLEAN       = True  # If True and no RMclean file found, exit instead of falling back
RMCLEAN_THRESHOLD_FRAC = 1e-6  # Components > this fraction of peak amplitude are significant
RMCLEAN_PHI_PADDING    = 25.0  # rad/m² added to each side of the derived range

NLIVE_FACTOR = 150
SLICE_FACTOR = 3.0
SAMPLER = 'dynamic'

FLAGGING_CONFIG_FILE = Path(__file__).parent / 'flagging_config.toml'

# ISM RM constraint mode
CONSTRAIN_ISM_RM = True        # Enable to pin first S (or P) component's phi_rm to ISM+sys+iono
PHI_ISM          = -0.5          # Expected Galactic ISM Faraday depth (rad/m²)
PHI_SYS          = -0.3          # Systematic RM offset to add (rad/m²)
PHI_ISM_SIGMA    = +0.65          # Gaussian std for truncated prior (rad/m²)
PHI_ISM_RANGE    = 3.5 * PHI_ISM_SIGMA  # Truncation half-width; set 0 for a fixed prior
IONO_LOC         = Path(__file__).parent.parent / 'files' / 'spinifex_iono'  # spinifex dir or single file; None → iono_rm = 0
IONO_SOURCE      = 'J1727'      # Case-insensitive substring match against first column of spinifex file


# =============================================================================
# DERIVED SAMPLER CONFIGURATION
# =============================================================================
# Auto-enable TIE_RM for all-powerlaw models
if all(c == 'P' for c in MODEL_TYPE) and len(MODEL_TYPE) > 1:
    TIE_RM = True
    print(f"Auto-enabling TIE_RM for all-powerlaw model: {MODEL_TYPE} ({len(MODEL_TYPE)} components)")

# Parameter counts
# S (thin):     3 params (a, b, phi_rm)
# T (thick):    5 params (a, b, phi_peak, sigma_phi, N)
# P (powerlaw): 4 params (a, b, phi_rm, beta)
# TIE_RM: shared phi_rm saves (n_thin + n_powerlaw - 1) parameters
n_thin_components = MODEL_TYPE.count('S')
n_thick_components = MODEL_TYPE.count('T')
n_powerlaw_components = MODEL_TYPE.count('P')
n_rm_components = n_thin_components + n_powerlaw_components

if TIE_RM and n_rm_components > 1:
    n_parameters = (3 * n_thin_components + 5 * n_thick_components
                    + 4 * n_powerlaw_components - (n_rm_components - 1))
else:
    n_parameters = 3 * n_thin_components + 5 * n_thick_components + 4 * n_powerlaw_components

nlive = NLIVE_FACTOR * n_parameters
slices = 3 + int(SLICE_FACTOR * n_parameters)
walks = 25 + 15 * n_parameters

print(f"Model has {n_parameters} parameters ({n_thin_components} thin, {n_thick_components} thick, {n_powerlaw_components} powerlaw)")
print(f"Setting nlive = {NLIVE_FACTOR} * {n_parameters} = {nlive}")

if SAMPLER == 'static':
    FIT_KWARGS = {
        'sampler': {'nlive': nlive, 'bound': 'multi', 'sample': 'rslice', 'slices': slices, 'rstate': rng},
        'run': {'dlogz': 0.01, 'maxiter': None, 'maxcall': None}
    }
elif SAMPLER == 'dynamic':
    FIT_KWARGS = {
        'sampler': {'bound': 'multi', 'sample': 'rwalk', 'walks': walks, 'rstate': rng},
        # 'run': {'dlogz_init': 0.1, 'nlive_init': int(nlive / 3.0), 'nlive_batch': int(nlive / 6.0),
        #         'n_effective': 10_000, 'maxiter': None, 'maxcall': None}
        'run': {'dlogz_init': 0.01, 'nlive_init': int(nlive / 2.0), 'nlive_batch': int(nlive / 6.0),
                'n_effective': 10_000, 'maxiter': None, 'maxcall': None}
    }
else:
    raise ValueError(f"Invalid SAMPLER option: {SAMPLER}")

# =============================================================================
# PRIOR DEFINITIONS
# =============================================================================

phi_min = -160.0
phi_max = 160.0

if '1014' in DATA_FILE_ARG:
    phi_min = -250.0
    phi_max = 250.0

if '4U1630' in DATA_FILE_ARG:
    phi_min = 1000.0
    phi_max = 2000.0

if USE_RMCLEAN_PHI_RANGE:
    _rmclean_range = phi_range_from_rmclean(DATA_FILE_ARG, RMCLEAN_THRESHOLD_FRAC, RMCLEAN_PHI_PADDING,
                                            fallback_phi_center=PHI_ISM)
    if _rmclean_range is not None:
        phi_min, phi_max = _rmclean_range
        print(f"phi range from RMclean CLEAN components: [{phi_min:.1f}, {phi_max:.1f}] rad/m²")
    else:
        print("WARNING: RMclean model file not found.")
        if REQUIRE_RMCLEAN:
            print("REQUIRE_RMCLEAN=True — exiting.")
            sys.exit(1)
        print(f"Falling back to default phi range: [{phi_min:.1f}, {phi_max:.1f}] rad/m²")

sigma_phi_max = 50.0 # 2x the FWHM where you lose ~50% sensitivity to thick components, can adjust based on expected source properties and data quality

# Using polar conversion style for intuitive p0 specification
THIN_PRIOR = {
    'type': 'ThinComponent',
    'polar_conversion': {
        'p0_bounds': (0.0, 0.50, 'radial_uniform'),      # Uniform in p0 (no bias)
        'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked'),       # Unbounded angle
        'peaked_k':1.0, # Optional and equal to default,
        'peaked_angle_deg': 0.0 # Optional and equal to default
    },
    'phi_rm': ('uniform', phi_min, phi_max)
}

# For mixed S/T models, confine one S-component RM prior to a narrow window.
MIXED_ST_CONSTRAINED_PHI_RM = ('uniform', -20.0, 20.0)

# Thick component sigma_phi max adjustment
THICK_PRIOR = {
    'type': 'ThickComponent',
    'polar_conversion': {
        'p0_bounds': (0.0, 0.50, 'radial_uniform'),      # Uniform in p0 (no bias)
        'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked')       # Unbounded angle
        #'psi_0_bounds': (-np.pi/2, np.pi/2, 'uniform')       # For 4U1630 Unbounded angle
    },
    'phi_peak': ('uniform', phi_min, phi_max),
    'sigma_phi': ('log_uniform', 1e-1, sigma_phi_max),
    'N': ('log_uniform', 2.0, 50.0),
    'effective_width_mode': False
}

POWERLAW_PRIOR = {
    'type': 'ThinPowerLawComponent',
    'polar_conversion': {
        'p0_bounds': (0.0, 0.50, 'radial_uniform'),
        'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked')
    },
    'phi_rm': ('uniform', phi_min, phi_max),
    'beta': ('uniform', -5.0, 5.0),
}



# =============================================================================
# DATA & I/O PATHS
# =============================================================================
DATA_LOC = DATA_FILE_ARG
if not Path(DATA_LOC).is_file():
    raise FileNotFoundError(f"Data file not found: {DATA_LOC}")

data_id = Path(DATA_LOC).stem
data_path = Path(DATA_LOC).resolve()

# Smart directory structure: replace 'files' with 'results'/'plots' if it exists, else use CWD
if 'files' in data_path.parts:
    # Find where 'files' appears and replace it
    base_parts = []
    for part in data_path.parts:
        if part == 'files':
            break
        base_parts.append(part)
    base_dir = Path(*base_parts)
    
    # Get subdirectory structure after 'files' (e.g., 'J1727/QU_text')
    files_idx = data_path.parts.index('files')
    subdir_parts = data_path.parts[files_idx + 1:-1]  # Everything between 'files' and filename
    
    # Filter out 'QU_text' from subdirectory path
    subdir_parts = [part for part in subdir_parts if part != 'QU_text']
    subdir = Path(*subdir_parts) if subdir_parts else Path()
    
    RESULTS_BASE = base_dir / 'results' / subdir
    PLOTS_BASE = base_dir / 'plots' / subdir
else:
    # No 'files' directory found, use parent directory of data file
    RESULTS_BASE = data_path.parent
    PLOTS_BASE = data_path.parent

# Create data_id-specific subdirectories for better organization
RESULTS_DIR = RESULTS_BASE / data_id
PLOTS_DIR = PLOTS_BASE / data_id
if N_TRIALS == 1:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Results directory: {RESULTS_DIR}")
    print(f"Plots directory: {PLOTS_DIR}")


# =============================================================================
# FLAGGING CONFIGURATION
# =============================================================================
print("Loading data to calculate Stokes I median...")
data_for_jitter = np.loadtxt(DATA_LOC)
stokes_i_median = np.median(data_for_jitter[:, 1])

# Load flagging parameters from TOML; apply the first matching override for data_id
with open(FLAGGING_CONFIG_FILE, 'rb') as _f:
    _cfg = tomllib.load(_f)

_fparams = dict(_cfg['default'])
for _override in _cfg.get('overrides', []):
    if any(m in data_id for m in _override['match']):
        _fparams.update({k: v for k, v in _override.items() if k != 'match'})
        print(f"Flagging: applied override (match={_override['match']})")
        break
else:
    print("Flagging: using default configuration")

JITTER_FRAC         = _fparams['jitter_frac']
sigma_clip_extreme  = _fparams['sigma_clip_extreme']
sigma_clip_v        = _fparams['sigma_clip_v']
sigma_clip_errors   = _fparams['sigma_clip_errors']
sigma_clip_spectrum = _fparams['sigma_clip_spectrum']
spectrum_clip_sided = _fparams['spectrum_clip_sided']
spectrum_deg        = _fparams['spectrum_deg']
spectrum_stokes     = _fparams['spectrum_stokes']
freq_ranges         = _fparams['freq_ranges']

spectrum_jitter = JITTER_FRAC * stokes_i_median
print(f"Stokes I median: {stokes_i_median:.6f}")
print(f"Jitter fraction: {JITTER_FRAC} ({JITTER_FRAC*100:.1f}%)")
print(f"Calculated jitter: {spectrum_jitter:.6f}")

load_stokes_v = True

flagging = FlaggingConfig(
    sigma_clip_extreme=sigma_clip_extreme,
    freq_ranges=freq_ranges,
    clip_stokes_i=None,
    max_iter=10,
    sigma_clip_errors=sigma_clip_errors,
    sigma_clip_errors_space='log',
    sigma_clip_spectrum=sigma_clip_spectrum,
    spectrum_deg=spectrum_deg,
    spectrum_jitter=spectrum_jitter,
    spectrum_clip_sided=spectrum_clip_sided,
    spectrum_stokes=spectrum_stokes,
    sigma_clip_v=sigma_clip_v,
    clip_v_against='median',
    verbose=True,
)


# =============================================================================
# DATA LOADING
# =============================================================================
print("Loading data...")
load_result = DataLoader.load_physical_with_fit(DATA_LOC, I_model='power-law', poly_degree=4,
                                                 central_freq=1.28,
                                                 flagging=flagging, load_stokes_v=load_stokes_v)

if flagging is not None:
    data, I_fitter, flagging_info = load_result
    plot_flagging_diagnostic(flagging_info['raw_data'], flagging_info,
                             filename=f'{data_id}_flagging_diagnostic.png', output_dir=PLOTS_DIR)
    save_flagged_data(flagging_info, filename=f'{data_id}_flagged_data.txt', output_dir=RESULTS_DIR)
else:
    data, I_fitter = load_result

print(f"Power-law fit: I₀ = {I_fitter.fit_params[0]:.4f}, α = {I_fitter.fit_params[1]:.4f}")

plot_data_diagnostic(data, I_fitter=I_fitter, title=None,
                    filename=f"{data_id}_data.png", output_dir=PLOTS_DIR,)
                    
# =============================================================================
# MODEL SETUP - Automatic from MODEL_TYPE string
# =============================================================================
model = CompositeModel()

# Initialize priors with tied_phi_rm if TIE_RM is enabled
if TIE_RM:
    priors = {'tied_phi_rm': ('uniform', phi_min, phi_max)}
else:
    priors = {}

# Count each component type
s_count = 0
t_count = 0
p_count = 0
rm_component_count = 0  # Track total thin + powerlaw components for TIE_RM logic
has_mixed_st = ('S' in MODEL_TYPE) and ('T' in MODEL_TYPE)

# Parse MODEL_TYPE string and build components automatically
for char in MODEL_TYPE:
    if char == 'S':
        # Thin component
        s_count += 1
        rm_component_count += 1
        name = f'S{s_count}'
        model.add_component(ThinComponent(a=0.1, b=0.0, phi_rm=0.0, name=name))

        thin_prior = THIN_PRIOR.copy()
        if has_mixed_st and s_count == 1 and '4U1630' not in DATA_FILE_ARG and not TIE_RM:
            thin_prior['phi_rm'] = MIXED_ST_CONSTRAINED_PHI_RM
        priors[name] = thin_prior
        
    elif char == 'T':
        # Thick component
        t_count += 1
        name = f'T{t_count}'
        model.add_component(ThickComponent(a=0.1, b=0.0, phi_peak=0.0, sigma_phi=1.0, N=2.0, name=name))
        priors[name] = THICK_PRIOR.copy()
        
    elif char == 'P':
        # Powerlaw component
        p_count += 1
        rm_component_count += 1
        name = f'P{p_count}'
        model.add_component(ThinPowerLawComponent(a=0.1, b=0.0, phi_rm=0.0, beta=0.0, lambda_sq_0=0.09, name=name))
        priors[name] = POWERLAW_PRIOR.copy()
        
    else:
        raise ValueError(f"Invalid character '{char}' in MODEL_TYPE. Use only 'S', 'T', or 'P'.")

print(f"\nBuilt model: {MODEL_TYPE}")
print(f"  {s_count} thin components (S)")
print(f"  {t_count} thick components (T)")
print(f"  {p_count} powerlaw components (P)")
print(f"  Total: {len(model.components)} components")
if TIE_RM and rm_component_count > 1:
    print(f"  TIE_RM: Enabled (shared phi_rm across {rm_component_count} components)")

print(model.__repr__())

# =============================================================================
# ISM RM CONSTRAINT — override phi_rm prior for first constrained component
# =============================================================================
if CONSTRAIN_ISM_RM:
    is_single_t = (MODEL_TYPE == 'T')

    if is_single_t:
        print("CONSTRAIN_ISM_RM: single-T model — skipping ISM RM constraint.")
    else:
        if IONO_LOC is not None:
            iono_rm, iono_err = load_iono_rm(IONO_LOC, IONO_SOURCE, data_file=DATA_FILE_ARG)
            print(f"CONSTRAIN_ISM_RM: iono_rm = {iono_rm:.3f} ± {iono_err:.3f} rad/m²"
                  f"  (from {IONO_LOC})")
        else:
            iono_rm = 0.0
            print("CONSTRAIN_ISM_RM: IONO_LOC not set — using iono_rm = 0.0")

        phi_center = PHI_ISM + PHI_SYS + iono_rm
        if PHI_ISM_RANGE > 0:
            ism_phi_prior = ('truncated_gaussian', phi_center, PHI_ISM_SIGMA,
                             phi_center - PHI_ISM_RANGE, phi_center + PHI_ISM_RANGE)
        else:
            ism_phi_prior = ('fixed', phi_center)

        if s_count > 0:
            priors['S1']['phi_rm'] = ism_phi_prior
            target = 'S1'
        elif p_count > 0:
            if TIE_RM:
                priors['tied_phi_rm'] = ism_phi_prior
                target = 'tied_phi_rm'
            else:
                priors['P1']['phi_rm'] = ism_phi_prior
                target = 'P1'
        else:
            raise ValueError(
                f"CONSTRAIN_ISM_RM: model '{MODEL_TYPE}' has no S or P component to constrain. "
                "Only a single pure-T model is exempt from the ISM RM constraint."
            )

        print(f"CONSTRAIN_ISM_RM: {target} phi_rm prior → {ism_phi_prior}")

print(f"\nModel type: {MODEL_TYPE}")
setup = FitSetup(model, priors, data, enforce_ordering=True, tie_phi_rm=TIE_RM, lambda_sq_ref=0.05)
print(f"Setup created: {setup.ndim} parameters")

# =============================================================================
# FITTING & RESULTS  (looped N_TRIALS times)
# =============================================================================
phi_range = 250
if '4U1630' in DATA_FILE_ARG:
    phi_range = 2500

for _trial in range(N_TRIALS):
    FIT_KWARGS['sampler']['rstate'] = np.random.default_rng(_seed_seq.spawn(1)[0])

    if N_TRIALS > 1:
        _sfx = f'_{_trial:02d}'
        _results_dir = RESULTS_BASE / f"{data_id}_{MODEL_TYPE}_trials" / f"{data_id}_{MODEL_TYPE}{_sfx}"
        _plots_dir   = PLOTS_BASE   / f"{data_id}_{MODEL_TYPE}_trials" / f"{data_id}_{MODEL_TYPE}{_sfx}"
        _results_dir.mkdir(parents=True, exist_ok=True)
        _plots_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n{'='*60}\nTrial {_trial + 1}/{N_TRIALS}  ({_sfx})\n{'='*60}")
        print(f"Results: {_results_dir}")
        print(f"Plots:   {_plots_dir}")
    else:
        _sfx = ''
        _results_dir = RESULTS_DIR
        _plots_dir   = PLOTS_DIR

    pkl_path = _results_dir / f'{data_id}_{MODEL_TYPE}{_sfx}_sampler.pkl'
    if SKIP_IF_EXISTS and pkl_path.exists():
        print(f"\nFound existing results at {pkl_path}, skipping fit and loading...")
        results = load_results(pkl_path)
    else:
        print(f"\nStarting fit for {data_id} with model {MODEL_TYPE}...")
        fitter = FaradayFitter(setup, dynesty_kwargs=FIT_KWARGS, n_cores=N_CORES, sampler=SAMPLER)
        results = fitter.fit(verbose=True, progress_bar=True)

    if not (SKIP_IF_EXISTS and pkl_path.exists()):
        results.save(pkl_path)

    plot_fit_results(results, true_model=None, I_fitter=I_fitter, title=None,
                    filename=f"{data_id}_{MODEL_TYPE}{_sfx}_fit.pdf", output_dir=_plots_dir,
                    plot_fdf_posterior_samples=True)

    plot_corner(
        samples=results.samples,
        param_names=results.setup.param_names,
        param_summary=results.param_summary,
        weights=results.weights,
        setup=results.setup,
        truths=None,
        method='median',
        title=None,
        filename=f"{data_id}_{MODEL_TYPE}{_sfx}_corner_cartesian.png",
        output_dir=_plots_dir
    )

    pconfig = ProcessingConfig(
        filename_prefix=f'{MODEL_TYPE}{_sfx}',
        detect_multiple_modes=True,
        kde_bandwidth='scott',
        plot_diagnostics=True,
        plot_dir=_plots_dir / f"{data_id}_{MODEL_TYPE}{_sfx}_modes/",
        verbose=False,
    )

    results = process_posterior_modes(results, pconfig)

    results.save_json(_results_dir / f'{data_id}_{MODEL_TYPE}{_sfx}_results.json', I_fitter=I_fitter)
    results.print_summary()

    plot_corner(
        samples=results.samples_processed,
        param_names=results.param_names_processed,
        param_summary=results.param_summary,
        processed_parameters=results.processed_parameters,
        weights=results.weights,
        setup=results.setup,
        method='mode',
        title=None,
        filename=f"{data_id}_{MODEL_TYPE}{_sfx}_corner_derived.png",
        output_dir=_plots_dir
    )

    plot_convergence(
        results,
        title=None,
        filename=f"{data_id}_{MODEL_TYPE}{_sfx}_convergence.png",
        output_dir=_plots_dir
    )

    plot_thick_component_band_response(
        results,
        filename=f"{data_id}_{MODEL_TYPE}{_sfx}_thick_band.png",
        output_dir=_plots_dir
    )

    if 'P' not in MODEL_TYPE:
        plot_fdf_diagnostic(results, phi_range=phi_range * 2, n_phi=4096, title=None, model_method='representative',
                           filename=f"{data_id}_{MODEL_TYPE}{_sfx}_fdf.png", output_dir=_plots_dir,
                           plot_posterior_samples=True)

print("\nAll done!")