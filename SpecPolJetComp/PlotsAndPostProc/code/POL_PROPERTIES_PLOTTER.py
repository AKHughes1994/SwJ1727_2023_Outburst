#!/usr/bin/env python
"""
Polarisation properties evolution figure for Swift J1727.8-1613.

Six stacked panels (single MNRAS column, ~3.5" wide):
  (1) Flux density at 1.28 GHz  [log scale]
  (2) Intraband spectral index alpha_nu
  (3) Linear polarisation fraction   [log scale, %]
  (4) Polarisation angle psi_0  (deg)
  (5) Faraday depth phi_rm / phi_peak  (rad/m^2)  [spinifex-corrected]
  (6) Faraday depth dispersion sigma_phi  (rad/m^2)  [thick components only]

Marker convention
  Circles  (o) = thin Faraday components (S-type)   — black filled (main) / open (ejecta)
  Diamonds (D) = thick Faraday components (T-type)  — black open (core and ejecta)
  Filled       = main source epochs
  Open         = ejecta epochs (_ejecta in filename)

Central values: posterior mode.  Error bars: 99% HDI.
Stokes I error: sqrt(I0_err^2 + (0.03 * I)^2); I0_err taken from JSON if
present, else 0.0 (systematic only).

Spinifex ionospheric corrections are applied to the Bayesian phi_rm values
(phi_rm is the only Faraday depth quantity available in the results JSON).
"""

# =============================================================================
# IMPORTS
# =============================================================================
import re
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from astropy.time import Time

# =============================================================================
# CONFIGURATION
# =============================================================================

# --- Paths -------------------------------------------------------------------
_HERE        = Path(__file__).resolve().parent   # SpecPolJetComp/PlotsAndPostProc/code
_SPEC_POL    = _HERE.parents[1]                  # SpecPolJetComp
DATA_DIR     = _SPEC_POL / 'QUfitRMsynth/results/J1727'
SPINIFEX_DIR = _SPEC_POL / 'QUfitRMsynth/files/spinifex_iono'

# --- JSON results files (chronological) --------------------------------------
ALL_JSON_PATHS = [
    DATA_DIR / 'SwiftJ1727_WAPITI_20230904/SwiftJ1727_WAPITI_20230904_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20230906/SwiftJ1727_QU_20230906_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20230908/SwiftJ1727_QU_20230908_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20230916/SwiftJ1727_QU_20230916_SS_results.json',
    DATA_DIR / 'SwiftJ1727_WAPITI_20230923/SwiftJ1727_WAPITI_20230923_ST_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231001/SwiftJ1727_QU_20231001_SS_results.json',
    DATA_DIR / 'SwiftJ1727_WAPITI_20231006/SwiftJ1727_WAPITI_20231006_STT_results.json',
    DATA_DIR / 'SwiftJ1727_WAPITI_20231014/SwiftJ1727_WAPITI_20231014_SSTT_results.json',
    DATA_DIR / 'SwiftJ1727_WAPITI_20231016/SwiftJ1727_WAPITI_20231016_ST_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231022/SwiftJ1727_QU_20231022_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231028/SwiftJ1727_QU_20231028_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231106/SwiftJ1727_QU_20231106_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231112/SwiftJ1727_QU_20231112_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231118/SwiftJ1727_QU_20231118_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20231125/SwiftJ1727_QU_20231125_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240210/SwiftJ1727_QU_20240210_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240210_ejecta/SwiftJ1727_QU_20240210_ejecta_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240219/SwiftJ1727_QU_20240219_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240219_ejecta/SwiftJ1727_QU_20240219_ejecta_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240225/SwiftJ1727_QU_20240225_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240225_ejecta/SwiftJ1727_QU_20240225_ejecta_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240331/SwiftJ1727_QU_20240331_S_results.json',
    DATA_DIR / 'SwiftJ1727_QU_20240331_ejecta/SwiftJ1727_QU_20240331_ejecta_S_results.json',
]

# --- Spinifex ionospheric corrections ----------------------------------------
# Applied to the Bayesian phi_rm from the results JSON (the only Faraday depth
# available here; this differs from the RMclean peak phi used in other scripts).
SPINIFEX_J1727_SOURCE = 'SwiftJ1727'
SPINIFEX_OFFSET       = -0.3   # rad/m^2 additional offset after spinifex subtraction
SPINIFEX_SYS_ERR      = 0.5    # rad/m^2 systematic added in quadrature to corrected errors

# --- Time range --------------------------------------------------------------
T_START = '2023-09-01T00:00:00'
T_END   = '2024-04-10T00:00:00'

# --- Reference frequency -----------------------------------------------------
MEERKAT_FREQ_GHZ = 1.28
STOKES_I_SYS_ERR = 0.03          # 3% flux calibration systematic

# --- Weighted mean phi_rm ----------------------------------------------------
# Primary component used for the weighted mean (fall back to T1 if S1 absent).
# Ejecta epochs excluded.  Result plotted as a horizontal line on phi_rm panel.
PHI_RM_PRIMARY_COMP = 'S1'

# --- Polarisation angle reference bands --------------------------------------
PSI_REF_CENTRES = (-1.0, 89.0)   # grey ±2° bands on psi_0 panel
PSI_REF_WIDTH   = 2.0             # half-width in degrees

# --- State transitions -------------------------------------------------------
TRANSITION_MJD = [60222]
HIGHLIGHT_MJD  = Time('2023-09-20T01:00:00', format='isot').mjd
HIGHLIGHT_SPAN = ('2023-09-15T00:00:00', '2023-10-17T00:00:00')

# --- Marker styles -----------------------------------------------------------
# Thin (S-type) → black filled/open;  Thick (T-type) → black open diamonds
# Main → filled (thin) or open (thick);  Ejecta → open
MK_SIZE_THIN  = 2.75   # ms for thin components  (original 5 × 0.75)
MK_SIZE_THICK = 3.00   # ms for thick components (original 4 × 0.75)
CAPSIZE       = 2
CAPTHICK      = 0.7
ELINEWIDTH    = 0.7

def _psi_min_offsets(psi_raw, psi_lo, psi_hi):
    """Return (abs_off, abs_elo, abs_ehi, sigma_off) to the nearest PSI_REF band edge.

    Uses the same pi-ambiguity wrap as panel 4.  Distance is zero if psi lies
    inside the band, otherwise it is the gap to the nearer edge.

    Error propagation for abs_off (d is a simple difference from a fixed edge):
      psi above edge: d = psi - hi  →  elo = psi_lo, ehi = psi_hi
      psi below edge: d = lo - psi  →  elo = psi_hi, ehi = psi_lo  (signs flip)
      inside band:    d = 0         →  elo = ehi = 0

    sigma_off uses the error bar on the side facing the band edge.
    """
    psi = psi_raw + (180.0 if psi_raw < -50.0 else 0.0)
    best_d    = np.inf
    best_elo  = np.nan
    best_ehi  = np.nan
    best_nsig = np.nan
    for centre in PSI_REF_CENTRES:
        lo = centre - PSI_REF_WIDTH
        hi = centre + PSI_REF_WIDTH
        if psi > hi:
            d, sigma        = psi - hi, psi_lo
            d_elo, d_ehi    = psi_lo,   psi_hi
        elif psi < lo:
            d, sigma        = lo - psi, psi_hi
            d_elo, d_ehi    = psi_hi,   psi_lo
        else:
            d, sigma        = 0.0,  np.nan
            d_elo, d_ehi    = 0.0,  0.0
        if d < best_d:
            best_d    = d
            best_elo  = d_elo
            best_ehi  = d_ehi
            best_nsig = 0.0 if d == 0.0 else (d / sigma if sigma > 0 else np.nan)
    return best_d, best_elo, best_ehi, best_nsig


def _upper_limit(comp):
    """Return True if p0 / p0_lo_68 < 2 (mode < 2-sigma from zero)."""
    dp0 = comp['p0_lo_68']
    return (comp['p0'] / dp0 < 2.0) if dp0 > 0 else True


def _plot_p0_upper_limit(ax, dt, comp, mk):
    """Marker + downward arrow at the 99% HDI upper bound, arrow spanning 1/3 decade."""
    p0_ul = comp['p0'] + comp['p0_hi']   # 99% HDI upper bound (%)
    ax.plot(dt, p0_ul, marker=mk['marker'], ms=mk['ms'],
            mfc=mk['mfc'], mec=mk['mec'], mew=mk['mew'],
            linestyle='none', zorder=mk['zorder'])
    ax.annotate('',
                xy=(dt, p0_ul * 10 ** (-1/3)),
                xytext=(dt, p0_ul),
                arrowprops=dict(arrowstyle='->', color='black',
                                lw=mk['elinewidth'], mutation_scale=5),
                zorder=mk['zorder'])


def _mk(comp_type, is_ejecta):
    """Return dict of errorbar marker kwargs."""
    if comp_type == 'thin':
        ms  = MK_SIZE_THIN
        mfc = 'white' if is_ejecta else 'black'
        mew = 0.6 if is_ejecta else 1.0
    else:
        ms  = MK_SIZE_THICK
        mfc = 'white'   # thick components always white-faced
        mew = 0.6
    return dict(marker='o' if comp_type == 'thin' else 'D',
                ms=ms, mfc=mfc, mec='black', mew=mew,
                ecolor='black', elinewidth=ELINEWIDTH,
                capsize=CAPSIZE, capthick=CAPTHICK,
                linestyle='none', zorder=10)

# --- Figure layout -----------------------------------------------------------
FIG_SIZE      = (4.0, 7.0)
LEFT_MARGIN   = 0.20
RIGHT_MARGIN  = 0.97
TOP_MARGIN    = 0.96
BOTTOM_MARGIN = 0.07
PANEL_HSPACE  = 0.06
HEIGHT_RATIOS = [1, 1, 1, 1, 1, 1]
TICK_INTERVAL = 30    # days between major x-ticks
TICK_PAD      = 15    # days of padding beyond first/last epoch on x-axis

# --- Font sizes --------------------------------------------------------------
FONT_BASE   = 8
FONT_LABEL  = 8
FONT_TICK   = 7
FONT_LEGEND = 7

# --- Output ------------------------------------------------------------------
OUTPUT_FORMAT    = 'pdf'
SAVE_PATH        = '../plots/results/results_pol_properties'
SAVE_PATH_PSI_SIG = '../plots/results/results_psi_angle_significance'
FIG_SIZE_PSI_SIG  = (3.5, 3.5)
SHOW_PLOT        = False

# =============================================================================
# RCPARAMS
# =============================================================================
plt.rcParams.update({
    'font.family':         'serif',
    'font.serif':          ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
    'mathtext.fontset':    'dejavuserif',
    'font.size':           FONT_BASE,
    'axes.labelsize':      FONT_LABEL,
    'axes.titlesize':      FONT_BASE,
    'axes.linewidth':      0.8,
    'xtick.major.width':   0.8,
    'ytick.major.width':   0.8,
    'xtick.minor.width':   0.5,
    'ytick.minor.width':   0.5,
    'xtick.major.size':    3,
    'ytick.major.size':    3,
    'xtick.minor.size':    1.5,
    'ytick.minor.size':    1.5,
    'xtick.direction':     'in',
    'ytick.direction':     'in',
    'xtick.top':           True,
    'xtick.bottom':        True,
    'ytick.left':          True,
    'ytick.right':         True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.labelsize':     FONT_TICK,
    'ytick.labelsize':     FONT_TICK,
    'legend.fontsize':     FONT_LEGEND,
    'legend.framealpha':   1.0,
    'legend.facecolor':    'white',
    'legend.edgecolor':    'black',
    'legend.fancybox':     False,
    'lines.linewidth':     0.8,
})

# =============================================================================
# AXIS FORMATTING HELPERS  (identical to BIG_FIGURE_PLOTTER)
# =============================================================================

def plot2mjd(t):
    '''Convert from matplotlib plot date to mjd'''
    return Time(t, format="plot_date", scale='utc').mjd

def mjd2plot(mjd):
    '''Convert from mjd to matplotlib plot'''
    return Time(mjd, format="mjd", scale='utc').plot_date

def FormatAxis(ax, mjd, dt=None, interval=60, tick_fontsize=None, xlabel_fontsize=None):
    if dt is None:
        dt = TICK_PAD
    if tick_fontsize is None:
        tick_fontsize = FONT_TICK - 1
    ax[0].set_xlabel('Observing Date (UTC)', fontfamily='serif',
                     fontsize=xlabel_fontsize)
    ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    ax[0].xaxis.set_minor_locator(mdates.DayLocator(interval=max(1, interval // 5)))
    ax[0].set_xlim(Time(mjd[0] - dt, format='mjd').datetime,
                   Time(mjd[-1] + dt, format='mjd').datetime)
    ax[0].xaxis.set_label_position('top')
    plt.gcf().axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax[0].tick_params(axis='x', which='major', rotation=15,
                      labeltop=True, labelbottom=False, labelsize=tick_fontsize)
    plt.setp(ax[0].get_xticklabels(), rotation=15, ha='left', fontsize=tick_fontsize)

    mjd_ax = ax[-1].secondary_xaxis('bottom', functions=(plot2mjd, mjd2plot))
    mjd_ax.set_xlabel('Observing Date (MJD)', fontfamily='serif', fontsize=xlabel_fontsize)
    mjd_ax.tick_params(which='major', direction='in', length=0.0, width=0.0,
                       labelsize=tick_fontsize)
    plt.draw()

    mjd_ticks = []
    for lab in ax[0].get_xticklabels(which='major'):
        mjd_ticks.append(lab.get_text() + 'T00:00:00')
    mjd_ticks = Time(mjd_ticks, format='isot').mjd.astype(int)
    mjd_ax.set_xticks(mjd_ticks, labels=mjd_ticks)
    if tick_fontsize is not None:
        mjd_ax.set_xticklabels(mjd_ticks, fontsize=tick_fontsize)
    mjd_ax.tick_params(which='both', length=0.0, width=0.0)
    return mjd_ax

# =============================================================================
# SPINIFEX IONOSPHERIC CORRECTION LOADING
# =============================================================================

print('\n' + '='*60)
print('SPINIFEX IONOSPHERIC CORRECTIONS')
print('='*60)

spinifex_by_date = {}
if SPINIFEX_DIR.exists():
    sf_files = sorted(SPINIFEX_DIR.glob('*.txt'))
    print(f'  Found {len(sf_files)} spinifex file(s) in {SPINIFEX_DIR}')
    for sf in sf_files:
        try:
            rows = []
            with open(sf) as fh:
                for line in fh:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if parts[0] == SPINIFEX_J1727_SOURCE:
                        try:
                            rows.append((float(parts[-2]), float(parts[-1])))
                        except (ValueError, IndexError):
                            pass
            if not rows:
                continue
            rows    = np.array(rows)
            unix_ts = int(sf.name.split('_')[0])
            sf_date = datetime.fromtimestamp(unix_ts, tz=timezone.utc).strftime('%Y-%m-%d')
            iono_rm, iono_err = float(np.median(rows[:, 0])), float(np.median(rows[:, 1]))
            spinifex_by_date[sf_date] = (iono_rm, iono_err)
            print(f'    {sf_date}  RM={iono_rm:+.4f}  err={iono_err:.4f} rad/m^2')
        except Exception as e:
            print(f'    [WARN] {sf.name}: {e}')
    print(f'  Loaded {len(spinifex_by_date)} dates with corrections for "{SPINIFEX_J1727_SOURCE}"')
else:
    print(f'  [WARN] SPINIFEX_DIR not found: {SPINIFEX_DIR}')
    print('  Proceeding without ionospheric corrections.')

# =============================================================================
# EPOCH PARSING
# =============================================================================

def _parse_date(filepath):
    m = re.search(r'(\d{8})', Path(filepath).stem)
    if m is None:
        raise ValueError(f'No 8-digit date in: {filepath}')
    d   = m.group(1)
    iso = f'{d[:4]}-{d[4:6]}-{d[6:8]}'
    t   = Time(iso, format='iso')
    return t.mjd, t.datetime


def _mode_hdi(pp, key):
    """Return (mode, hdi_lo, hdi_hi) from a processed_parameters entry."""
    entry = pp[key]['mode_1']
    return entry['mode'], entry['hdi_99'][0], entry['hdi_99'][1]


def _comp_type(name):
    return 'thin' if name.startswith('S') else 'thick'


def parse_epoch(json_path):
    json_path = Path(json_path)
    mjd, dt   = _parse_date(json_path)
    is_ejecta = 'ejecta' in json_path.stem

    with open(json_path) as f:
        data = json.load(f)

    pp         = data['processed_parameters']
    comp_names = data['metadata']['component_names']
    model_type = data['metadata']['model_type']
    sif        = data.get('stokes_i_fit', {})

    # --- Stokes I at 1.28 GHz ------------------------------------------------
    I0        = sif.get('I0', np.nan)
    I0_err    = sif.get('I0_err', 0.0)           # statistical; 0 if absent
    freq0     = sif.get('freq0_GHz', MEERKAT_FREQ_GHZ)
    alpha_nu  = sif.get('alpha_nu',  0.0)
    I_1p28    = I0 * (MEERKAT_FREQ_GHZ / freq0) ** alpha_nu * 1e3   # → mJy
    I_err     = np.sqrt(I0_err**2 + (STOKES_I_SYS_ERR * abs(I_1p28))**2)

    # --- Per-component parameters --------------------------------------------
    components = []
    for name in comp_names:
        ctype = _comp_type(name)

        p0, p0_lo, p0_hi   = _mode_hdi(pp, f'{name}_p0')
        p0_hdi68_lo        = pp[f'{name}_p0']['mode_1']['hdi_68'][0]
        psi, psi_lo, psi_hi = _mode_hdi(pp, f'{name}_psi_0')

        phi_key = (f'{name}_phi_rm'   if f'{name}_phi_rm'   in pp else
                   f'{name}_phi_peak' if f'{name}_phi_peak' in pp else None)
        if phi_key:
            phi, phi_lo_abs, phi_hi_abs = _mode_hdi(pp, phi_key)
            phi_lo = phi - phi_lo_abs
            phi_hi = phi_hi_abs - phi
        else:
            phi = phi_lo = phi_hi = np.nan

        sig_key = f'{name}_sigma_phi' if ctype == 'thick' else None
        if sig_key and sig_key in pp:
            sphi, sphi_lo_abs, sphi_hi_abs = _mode_hdi(pp, sig_key)
            sphi_lo = sphi - sphi_lo_abs
            sphi_hi = sphi_hi_abs - sphi
        else:
            sphi = sphi_lo = sphi_hi = None

        # --- Spinifex correction for phi_rm (main source epochs only) --------
        iso_date_str = Time(mjd, format='mjd').iso[:10]
        phi_corr       = phi
        phi_lo_corr    = phi_lo
        phi_hi_corr    = phi_hi
        has_iono_corr  = False
        if not is_ejecta and not np.isnan(phi) and iso_date_str in spinifex_by_date:
            iono_rm, iono_err = spinifex_by_date[iso_date_str]
            phi_corr    = phi - iono_rm + SPINIFEX_OFFSET
            phi_lo_corr = np.sqrt(phi_lo**2 + iono_err**2 + SPINIFEX_SYS_ERR**2)
            phi_hi_corr = np.sqrt(phi_hi**2 + iono_err**2 + SPINIFEX_SYS_ERR**2)
            has_iono_corr = True

        components.append({
            'name':          name,
            'comp_type':     ctype,
            'p0':            p0  * 100.0,
            'p0_lo':         (p0  - p0_lo)  * 100.0,
            'p0_hi':         (p0_hi  - p0)  * 100.0,
            'p0_lo_68':      (p0  - p0_hdi68_lo) * 100.0,
            'psi_0':         psi,
            'psi_0_lo':      psi - psi_lo,
            'psi_0_hi':      psi_hi - psi,
            'phi':           phi_corr,
            'phi_lo':        phi_lo_corr,
            'phi_hi':        phi_hi_corr,
            'has_iono_corr': has_iono_corr,
            'sigma_phi':     sphi,
            'sigma_phi_lo':  sphi_lo,
            'sigma_phi_hi':  sphi_hi,
        })

    return {
        'mjd':          mjd,
        'dt':           dt,
        'iso_date':     Time(mjd, format='mjd').iso[:10],
        'is_ejecta':    is_ejecta,
        'model_type':   model_type,
        'stokes_i_mJy': I_1p28,
        'stokes_i_err': I_err,
        'alpha_nu':     sif.get('alpha_nu',     np.nan),
        'alpha_nu_err': sif.get('alpha_nu_err', np.nan),
        'components':   components,
    }

# =============================================================================
# LOAD ALL EPOCHS
# =============================================================================

print('\n' + '='*60)
print('EPOCH LOADING')
print('='*60)

epochs = []
n_skip = 0
for path in ALL_JSON_PATHS:
    path = Path(path)
    if not path.exists():
        print(f'  [SKIP] Not found: {path.name}')
        n_skip += 1
        continue
    try:
        ep = parse_epoch(path)
        epochs.append(ep)
        comps_str = ', '.join(
            f'{c["name"]}({c["comp_type"][0].upper()}'
            f'{"*" if c["has_iono_corr"] else ""})'
            for c in ep['components'])
        ejecta_tag = ' [EJECTA]' if ep['is_ejecta'] else ''
        print(f'  OK  {ep["iso_date"]}{ejecta_tag:9s}  '
              f'I={ep["stokes_i_mJy"]:8.2f} mJy  '
              f'alpha={ep["alpha_nu"]:+.3f}  '
              f'comps=[{comps_str}]  (* = iono corrected)')
    except Exception as e:
        print(f'  [ERR] {path.name}: {e}')
        n_skip += 1

epochs.sort(key=lambda e: e['mjd'])
print(f'\n  Loaded {len(epochs)} epochs  ({n_skip} skipped)')

# =============================================================================
# UPPER LIMIT DIAGNOSTIC
# =============================================================================

print('\n' + '='*60)
print('UPPER LIMIT CHECK  (SNR = p0 / hdi68_lo_err < 3)')
print(f'  {"Date":<12} {"Comp":<5} {"p0 mode (%)":>12} {"hdi68_lo_err":>13} {"SNR":>6}  UL?')
print('  ' + '-'*58)
for ep in epochs:
    for comp in ep['components']:
        dp0 = comp['p0_lo_68']
        snr = (comp['p0'] / dp0) if dp0 > 0 else np.nan
        is_ul = _upper_limit(comp)
        ul_tag = '<-- UL' if is_ul else ''
        print(f'  {ep["iso_date"]:<12} {comp["name"]:<5} '
              f'{comp["p0"]:>12.4f} {dp0:>13.4f} '
              f'{snr:>6.2f}  {ul_tag}')
print('='*60)

# =============================================================================
# WEIGHTED MEAN PHI_RM  (main source, S1 preferred, fall back to T1)
# =============================================================================

print('\n' + '='*60)
print('WEIGHTED MEAN phi_rm')
print('='*60)

phi_vals, phi_weights = [], []
for ep in epochs:
    if ep['is_ejecta']:
        continue
    primary = next((c for c in ep['components'] if c['name'] == PHI_RM_PRIMARY_COMP),
                   next((c for c in ep['components'] if c['name'] == 'T1'), None))
    if primary is None or np.isnan(primary['phi']):
        print(f'  {ep["iso_date"]}  no usable primary component — skipped')
        continue
    sigma = 0.5 * (primary['phi_lo'] + primary['phi_hi'])
    if sigma <= 0:
        print(f'  {ep["iso_date"]}  sigma_phi=0, skipped')
        continue
    phi_vals.append(primary['phi'])
    phi_weights.append(1.0 / sigma**2)
    print(f'  {ep["iso_date"]}  phi={primary["phi"]:+.3f}  sigma={sigma:.3f}  '
          f'comp={primary["name"]}  iono={primary["has_iono_corr"]}')

if phi_vals:
    phi_wt       = np.array(phi_weights)
    phi_mean     = float(np.average(phi_vals, weights=phi_wt))
    phi_mean_err = 1.0 / np.sqrt(np.sum(phi_wt))
    n_phi        = len(phi_vals)
    chi2         = float(np.sum(phi_wt * (np.array(phi_vals) - phi_mean)**2))
    chi2_red     = chi2 / (n_phi - 1) if n_phi > 1 else np.nan
    scale        = np.sqrt(max(chi2_red, 1.0)) if np.isfinite(chi2_red) else 1.0
    phi_mean_err_infl = phi_mean_err * scale
    print(f'\n  Weighted mean phi_rm = {phi_mean:.3f} rad/m^2  (n={n_phi})')
    print(f'  chi2_red             = {chi2_red:.3f}  (inflation factor x{scale:.3f})')
    print(f'  Raw error            : +/- {phi_mean_err:.3f} rad/m^2')
    print(f'  Inflated error       : +/- {phi_mean_err_infl:.3f} rad/m^2')
else:
    phi_mean = None
    print('  No valid phi_rm values — weighted mean not computed.')

# =============================================================================
# STATE TRANSITIONS
# =============================================================================

start_dt    = Time(T_START, format='isot').datetime
end_dt      = Time(T_END,   format='isot').datetime
trans_dt    = [Time(m, format='mjd').datetime for m in TRANSITION_MJD]
hl_dt       = Time(HIGHLIGHT_MJD, format='mjd').datetime
hl_span_dt  = (Time(HIGHLIGHT_SPAN[0], format='isot').datetime,
               Time(HIGHLIGHT_SPAN[1], format='isot').datetime)
mjd_all     = np.array([e['mjd'] for e in epochs])


def _add_transitions(ax):
    ax.axvspan(*hl_span_dt, color='C0', alpha=0.15, zorder=0)
    for tdt in trans_dt:
        if start_dt <= tdt <= end_dt:
            ax.axvline(tdt, ls='--', c='k', lw=0.7, zorder=0)
    if start_dt <= hl_dt <= end_dt:
        ax.axvline(hl_dt, ls=':', c='k', lw=0.7, zorder=0)

# =============================================================================
# FIGURE CONSTRUCTION
# =============================================================================

fig = plt.figure(figsize=FIG_SIZE)
gs  = gridspec.GridSpec(
    6, 1, figure=fig,
    hspace=PANEL_HSPACE,
    height_ratios=HEIGHT_RATIOS,
    left=LEFT_MARGIN, right=RIGHT_MARGIN,
    top=TOP_MARGIN, bottom=BOTTOM_MARGIN,
)
ax_I     = fig.add_subplot(gs[0])
ax_alpha = fig.add_subplot(gs[1], sharex=ax_I)
ax_p0    = fig.add_subplot(gs[2], sharex=ax_I)
ax_psi   = fig.add_subplot(gs[3], sharex=ax_I)
ax_phi   = fig.add_subplot(gs[4], sharex=ax_I)
ax_sig   = fig.add_subplot(gs[5], sharex=ax_I)
all_axes = [ax_I, ax_alpha, ax_p0, ax_psi, ax_phi, ax_sig]

# =============================================================================
# PANEL 1 — Stokes I  [log]
# =============================================================================

for ep in epochs:
    mfc = 'white' if ep['is_ejecta'] else 'black'
    mew = 0.6 if ep['is_ejecta'] else 1.0
    mk = dict(marker='o', ms=MK_SIZE_THIN, mfc=mfc, mec='black', mew=mew,
              ecolor='black', elinewidth=ELINEWIDTH,
              capsize=CAPSIZE, capthick=CAPTHICK, linestyle='none', zorder=10)
    ax_I.errorbar(ep['dt'], ep['stokes_i_mJy'], yerr=ep['stokes_i_err'], **mk)

ax_I.set_yscale('log')
ax_I.set_ylabel('Flux Density\n(1.28 GHz; mJy)', fontsize=FONT_LABEL)
_add_transitions(ax_I)

# =============================================================================
# PANEL 2 — Spectral index
# =============================================================================

for ep in epochs:
    if np.isnan(ep['alpha_nu']):
        continue
    mfc = 'white' if ep['is_ejecta'] else 'black'
    mew = 0.6 if ep['is_ejecta'] else 1.0
    mk = dict(marker='o', ms=MK_SIZE_THIN, mfc=mfc, mec='black', mew=mew,
              ecolor='black', elinewidth=ELINEWIDTH,
              capsize=CAPSIZE, capthick=CAPTHICK, linestyle='none', zorder=10)
    ax_alpha.errorbar(ep['dt'], ep['alpha_nu'], yerr=ep['alpha_nu_err'], **mk)

ax_alpha.axhline(0.0, ls=':', c='k', lw=0.7)
ax_alpha.set_ylabel(r'$\alpha_\nu$', fontsize=FONT_LABEL)
_add_transitions(ax_alpha)

# =============================================================================
# PANELS 3–6 — per-component polarisation properties
# =============================================================================

_legend_handles = {}   # (comp_type, is_ejecta) → handle for legend

for ep in epochs:
    for comp in ep['components']:
        ct  = comp['comp_type']
        key = (ct, ep['is_ejecta'])
        mk  = _mk(ct, ep['is_ejecta'])
        is_ul = _upper_limit(comp)

        # Panel 3: polarisation fraction (log, %)
        if is_ul:
            _plot_p0_upper_limit(ax_p0, ep['dt'], comp, mk)
        else:
            h = ax_p0.errorbar(ep['dt'], comp['p0'],
                               yerr=[[comp['p0_lo']], [comp['p0_hi']]], **mk)
            _legend_handles.setdefault(key, h)

        if is_ul:
            continue   # skip psi, phi, sigma panels for non-detections

        # Panel 4: polarisation angle
        # Wrap angles below -50 deg by adding 180 (pi-ambiguity)
        _psi = comp['psi_0'] + (180.0 if comp['psi_0'] < -50.0 else 0.0)
        ax_psi.errorbar(ep['dt'], _psi,
                        yerr=[[comp['psi_0_lo']], [comp['psi_0_hi']]], **mk)

        # Panel 5: Faraday depth (iono-corrected where available)
        if not np.isnan(comp['phi']):
            ax_phi.errorbar(ep['dt'], comp['phi'],
                            yerr=[[comp['phi_lo']], [comp['phi_hi']]], **mk)

        # Panel 6: sigma_phi (thick only)
        if ct == 'thick' and comp['sigma_phi'] is not None:
            ax_sig.errorbar(ep['dt'], comp['sigma_phi'],
                            yerr=[[comp['sigma_phi_lo']], [comp['sigma_phi_hi']]], **mk)

# --- Panel 3 -----------------------------------------------------------------
ax_p0.set_yscale('log')
ax_p0.set_ylabel(r'$p_0$ (%)', fontsize=FONT_LABEL)
_add_transitions(ax_p0)

# --- Panel 4 -----------------------------------------------------------------
for centre in PSI_REF_CENTRES:
    ax_psi.axhspan(centre - PSI_REF_WIDTH, centre + PSI_REF_WIDTH,
                   color='grey', alpha=0.25, zorder=0)
ax_psi.set_ylabel(r'$\psi_0$ (deg)', fontsize=FONT_LABEL)
ax_psi.set_ylim(-50, 120)
_add_transitions(ax_psi)

# --- Panel 5 -----------------------------------------------------------------
ax_phi.axhline(0.0, ls=':', c='k', lw=0.7, zorder=0)
if phi_mean is not None:
    ax_phi.axhline(phi_mean, ls='-', c='grey', lw=1.0, zorder=1)
    ax_phi.axhspan(phi_mean - phi_mean_err_infl, phi_mean + phi_mean_err_infl,
                   color='grey', alpha=0.25, zorder=0)
ax_phi.set_ylabel(r'$\phi_{\rm rm}$ (rad m$^{-2}$)', fontsize=FONT_LABEL)
_add_transitions(ax_phi)

# Zoomed inset: +/-5 rad/m^2, top right
ax_phi_inset = ax_phi.inset_axes([0.408, 0.555, 0.5625, 0.405])
for ep in epochs:
    for comp in ep['components']:
        if not np.isnan(comp['phi']):
            mk = _mk(comp['comp_type'], ep['is_ejecta'])
            mk['ms']        = mk['ms'] * 0.5 * 1.1 * 1.1
            mk['capsize']   = CAPSIZE * 0.5
            mk['capthick']  = CAPTHICK * 0.5
            mk['elinewidth']= ELINEWIDTH * 0.5
            ax_phi_inset.errorbar(ep['dt'], comp['phi'],
                                  yerr=[[comp['phi_lo']], [comp['phi_hi']]], **mk)
ax_phi_inset.axhline(0.0, ls=':', c='k', lw=0.7, zorder=0)
if phi_mean is not None:
    ax_phi_inset.axhline(phi_mean, ls='-', c='grey', lw=1.0, zorder=1)
    ax_phi_inset.axhspan(phi_mean - phi_mean_err_infl, phi_mean + phi_mean_err_infl,
                         color='grey', alpha=0.25, zorder=0)
ax_phi_inset.set_xlim(ax_phi.get_xlim())
ax_phi_inset.set_ylim(-5, 3)
ax_phi_inset.tick_params(axis='x', which='both',
                          bottom=False, top=False, labelbottom=False)
ax_phi_inset.tick_params(axis='y', which='both', left=True, right=True,
                          direction='in', labelsize=FONT_TICK - 1)
ax_phi_inset.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax_phi_inset.patch.set_facecolor('white')
ax_phi_inset.set_facecolor('white')
for spine in ax_phi_inset.spines.values():
    spine.set_linewidth(0.6)

# --- Panel 6 -----------------------------------------------------------------
ax_sig.set_ylabel(r'$\sigma_\phi$ (rad m$^{-2}$)', fontsize=FONT_LABEL)
_add_transitions(ax_sig)

# =============================================================================
# LEGEND  (on Stokes I panel, top right)
# =============================================================================

label_map = {
    ('thin',  False): 'Thin (core)',
    ('thin',  True):  'Thin (ejecta)',
    ('thick', False): 'Thick (core)',
    ('thick', True):  'Thick (ejecta)',
}
legend_items = [(label_map[k], v) for k, v in _legend_handles.items() if k in label_map]
if legend_items:
    ax_I.legend([h for _, h in legend_items], [l for l, _ in legend_items],
                fontsize=FONT_LEGEND, loc='upper right',
                handlelength=1.0, handletextpad=0.4)

# =============================================================================
# X-AXIS FORMATTING  (FormatAxis handles UTC top + MJD bottom)
# =============================================================================

ax_I.set_xlim(start_dt, end_dt)
for ax in all_axes[:-1]:
    ax.tick_params(labelbottom=False)

mjd_secondary_ax = FormatAxis(
    np.array(all_axes), mjd=mjd_all,
    interval=TICK_INTERVAL,
)
ax_sig.tick_params(labelbottom=False)   # FormatAxis secondary_xaxis handles bottom

fig.align_ylabels(all_axes)

# =============================================================================
# SAVE / SHOW
# =============================================================================

output_format = OUTPUT_FORMAT.strip().lower()
save_path     = Path(SAVE_PATH).with_suffix(f'.{output_format}')
save_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(save_path, dpi=300, bbox_inches='tight')
print(f'\nSaved -> {save_path}')

# =============================================================================
# FIGURE 2 — psi_0 offsets: absolute (top) and sigma (bottom)
# =============================================================================

fig2 = plt.figure(figsize=FIG_SIZE_PSI_SIG)
gs2  = gridspec.GridSpec(
    2, 1, figure=fig2, hspace=PANEL_HSPACE,
    left=LEFT_MARGIN, right=RIGHT_MARGIN,
    top=TOP_MARGIN, bottom=0.14,
)
ax2_abs = fig2.add_subplot(gs2[0])
ax2_sig = fig2.add_subplot(gs2[1], sharex=ax2_abs)
all_axes2 = [ax2_abs, ax2_sig]

_sig_handles = {}

for ep in epochs:
    for comp in ep['components']:
        if _upper_limit(comp):
            continue
        ct  = comp['comp_type']
        key = (ct, ep['is_ejecta'])
        abs_off, abs_elo, abs_ehi, sig_off = _psi_min_offsets(
            comp['psi_0'], comp['psi_0_lo'], comp['psi_0_hi'])
        mk      = _mk(ct, ep['is_ejecta'])
        is_sig  = np.isfinite(sig_off) and sig_off >= 3.0
        alpha   = 1.0 if is_sig else 0.33
        ms      = mk['ms'] + 1 + (1 if is_sig else 0)   # +1 base, +1 extra if significant
        pkw = dict(marker=mk['marker'], ms=ms, alpha=alpha,
                   mfc=mk['mfc'], mec=mk['mec'], mew=mk['mew'],
                   linestyle='none', zorder=10)
        if np.isfinite(abs_off):
            h = ax2_abs.errorbar(ep['dt'], abs_off,
                                 yerr=[[abs_elo], [abs_ehi]],
                                 ecolor='black', elinewidth=ELINEWIDTH,
                                 capsize=CAPSIZE, capthick=CAPTHICK, **pkw)
            _sig_handles.setdefault(key, h)
        if np.isfinite(sig_off):
            ax2_sig.plot(ep['dt'], sig_off, **pkw)

for ax in all_axes2:
    ax.axhline(0.0, ls=':', c='k', lw=0.7, zorder=1)
    ax.axvspan(*hl_span_dt, color='C0', alpha=0.15, zorder=0)

ax2_sig.axhline(3.0, ls='--', c='k', lw=0.7, zorder=1)

ax2_abs.set_ylim(-1.5, 35)
ax2_sig.set_ylim(bottom=-0.2)

ax2_abs.set_ylabel(r'Min $\psi_0$ abs-offset (deg)', fontsize=FONT_LABEL)
ax2_sig.set_ylabel(r'Min $\psi_0$ $\sigma$-offset', fontsize=FONT_LABEL)

sig_legend_items = [(label_map[k], v) for k, v in _sig_handles.items() if k in label_map]
if sig_legend_items:
    leg2 = ax2_abs.legend([h for _, h in sig_legend_items], [l for l, _ in sig_legend_items],
                          fontsize=FONT_LEGEND, loc='upper right',
                          handlelength=1.0, handletextpad=0.4)
    leg2.set_zorder(20)
    for lh in leg2.legend_handles:
        lh.set_alpha(1.0)

ax2_abs.tick_params(labelbottom=False)
FormatAxis(np.array(all_axes2), mjd=mjd_all, interval=TICK_INTERVAL,
           tick_fontsize=FONT_TICK - 1, xlabel_fontsize=FONT_TICK - 1)
ax2_sig.tick_params(labelbottom=False)   # FormatAxis secondary_xaxis handles bottom

fig2.align_ylabels(all_axes2)

save_path2 = Path(SAVE_PATH_PSI_SIG).with_suffix(f'.{output_format}')
save_path2.parent.mkdir(parents=True, exist_ok=True)
fig2.savefig(save_path2, dpi=300, bbox_inches='tight')
print(f'Saved -> {save_path2}')

# =============================================================================
# FIGURE 3 — Zoomed two-panel: p0 and psi_0  (2023-10-19 to 2023-11-27)
# =============================================================================

ZOOM_START     = '2023-10-19T00:00:00'
ZOOM_END       = '2023-11-27T00:00:00'
SAVE_PATH_ZOOM = '../plots/results/results_pol_zoom'
FIG_SIZE_ZOOM  = (4.0, 5.0)

zoom_start_dt = Time(ZOOM_START, format='isot').datetime
zoom_end_dt   = Time(ZOOM_END,   format='isot').datetime

fig3 = plt.figure(figsize=FIG_SIZE_ZOOM)
gs3  = gridspec.GridSpec(
    3, 1, figure=fig3,
    hspace=PANEL_HSPACE,
    left=LEFT_MARGIN, right=RIGHT_MARGIN,
    top=TOP_MARGIN, bottom=0.14,
)
ax3_I   = fig3.add_subplot(gs3[0])
ax3_p0  = fig3.add_subplot(gs3[1], sharex=ax3_I)
ax3_psi = fig3.add_subplot(gs3[2], sharex=ax3_I)
all_axes3 = [ax3_I, ax3_p0, ax3_psi]

# --- Panel 1: Stokes I -------------------------------------------------------
for ep in epochs:
    mfc = 'white' if ep['is_ejecta'] else 'black'
    mew = 0.6    if ep['is_ejecta'] else 1.0
    mk  = dict(marker='o', ms=MK_SIZE_THIN, mfc=mfc, mec='black', mew=mew,
               ecolor='black', elinewidth=ELINEWIDTH,
               capsize=CAPSIZE, capthick=CAPTHICK, linestyle='none', zorder=10)
    ax3_I.errorbar(ep['dt'], ep['stokes_i_mJy'], yerr=ep['stokes_i_err'], **mk)

ax3_I.set_yscale('log')
ax3_I.set_ylim(1, 100)
ax3_I.set_ylabel('Flux Density\n(1.28 GHz; mJy)', fontsize=FONT_LABEL)
_add_transitions(ax3_I)

# --- Panel 2: polarisation fraction ------------------------------------------
for ep in epochs:
    for comp in ep['components']:
        mk    = _mk(comp['comp_type'], ep['is_ejecta'])
        is_ul = _upper_limit(comp)
        if is_ul:
            _plot_p0_upper_limit(ax3_p0, ep['dt'], comp, mk)
        else:
            ax3_p0.errorbar(ep['dt'], comp['p0'],
                            yerr=[[comp['p0_lo']], [comp['p0_hi']]], **mk)

ax3_p0.set_yscale('log')
ax3_p0.set_ylim(0.1, 15)
ax3_p0.set_ylabel(r'$p_0$ (%)', fontsize=FONT_LABEL)
_add_transitions(ax3_p0)

# --- Panel 3: polarisation angle — only the ~0° reference band ---------------
ax3_psi.axhspan(PSI_REF_CENTRES[0] - PSI_REF_WIDTH,
                PSI_REF_CENTRES[0] + PSI_REF_WIDTH,
                color='grey', alpha=0.25, zorder=0)

for ep in epochs:
    for comp in ep['components']:
        if _upper_limit(comp):
            continue
        mk   = _mk(comp['comp_type'], ep['is_ejecta'])
        _psi = comp['psi_0'] + (180.0 if comp['psi_0'] < -50.0 else 0.0)
        ax3_psi.errorbar(ep['dt'], _psi,
                         yerr=[[comp['psi_0_lo']], [comp['psi_0_hi']]], **mk)

ax3_psi.set_ylim(-30, 30)
ax3_psi.set_ylabel(r'$\psi_0$ (deg)', fontsize=FONT_LABEL)
_add_transitions(ax3_psi)

# --- X-axis: 10-day ticks, then override xlim to zoom window -----------------
ax3_I.tick_params(labelbottom=False)
ax3_p0.tick_params(labelbottom=False)

FormatAxis(
    np.array(all_axes3), mjd=mjd_all,
    interval=10,
    tick_fontsize=FONT_TICK - 1,
    xlabel_fontsize=FONT_TICK - 1,
)
ax3_psi.tick_params(labelbottom=False)

# FormatAxis sets xlim from all epochs; override here to enforce the zoom
ax3_I.set_xlim(zoom_start_dt, zoom_end_dt)

fig3.align_ylabels(all_axes3)

save_path3 = Path(SAVE_PATH_ZOOM).with_suffix(f'.{output_format}')
save_path3.parent.mkdir(parents=True, exist_ok=True)
fig3.savefig(save_path3, dpi=300, bbox_inches='tight')
print(f'Saved -> {save_path3}')

if SHOW_PLOT:
    plt.show()
plt.close(fig)
