#!/usr/bin/env python
"""
Publication figure for Swift J1727.8-1613 spectropolarimetric analysis.

Layout
------
Portrait figure (MNRAS double-column, 6.85" × 9.21") split into two equal halves:

  Top 50%    Three-panel radio light curve:
               (1) Flux density at 1.28 GHz
               (2) Spectral index α_ν
               (3) Intrinsic linear polarisation fraction p_0
             X-axis spans LC_START to LC_END.  Letters A–G are annotated on
             the flux sub-panel at the MeerKAT data points within the labelling
             window (LABEL_WINDOW_START to LABEL_WINDOW_END).

  Bottom 50% 3×3 grid of FDF panels for the 9 selected epochs (FDF_PKL_PATHS).
             Each panel shows posterior samples, best-fit model, and EVPA dials.
             Labels: PRE_FLARE_LABEL (epoch before window), sequential letters
             A–G (epochs within window), POST_FLARE_LABEL (epoch after window).

The polarisation fraction sub-panel draws on ALL_PKL_PATHS (all available epochs)
to show the complete time series out to LC_END.

Dependencies
  numpy, matplotlib, pandas, astropy, plus local faraday_model / faraday_utils.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
import copy
import json
import re
import string
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.patches import Wedge, Circle
from pathlib import Path
from astropy.time import Time

_HERE     = Path(__file__).resolve().parent   # SpecPolJetComp/PlotsAndPostProc/code
_SPEC_POL = _HERE.parents[1]                  # SpecPolJetComp

sys.path.insert(0, str(_SPEC_POL / 'QUfitRMsynth/code'))
from faraday_model import load_results, update_model_from_params
from faraday_utils import ThinComponent

# =============================================================================
# CONFIGURATION — all plot controls are here
# =============================================================================

# --- Data paths --------------------------------------------------------------
DATA_DIR  = _SPEC_POL / 'QUfitRMsynth/results/J1727'
LC_DIR    = _HERE.parent / 'files/LC_files'
RADIO_CSV = LC_DIR / 'SW1727_Radio.csv'

# --- All epochs (sorted chronologically).  FDF_PANEL_DATES marks the 9
#     epochs that appear in the 3x3 FDF grid; all entries feed the pol LC. ----
ALL_PKL_PATHS = [
    DATA_DIR / 'SwiftJ1727_WAPITI_20230904/SwiftJ1727_WAPITI_20230904_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20230906/SwiftJ1727_QU_20230906_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20230908/SwiftJ1727_QU_20230908_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20230916/SwiftJ1727_QU_20230916_SS_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_WAPITI_20230923/SwiftJ1727_WAPITI_20230923_ST_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231001/SwiftJ1727_QU_20231001_SS_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_WAPITI_20231006/SwiftJ1727_WAPITI_20231006_STT_sampler.pkl',
    # DATA_DIR / 'SwiftJ1727_WAPITI_20231014/SwiftJ1727_WAPITI_20231014_STTT_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_WAPITI_20231014/SwiftJ1727_WAPITI_20231014_SSTT_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_WAPITI_20231016/SwiftJ1727_WAPITI_20231016_ST_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231022/SwiftJ1727_QU_20231022_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231028/SwiftJ1727_QU_20231028_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231106/SwiftJ1727_QU_20231106_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231112/SwiftJ1727_QU_20231112_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231118/SwiftJ1727_QU_20231118_S_sampler.pkl',
    DATA_DIR / 'SwiftJ1727_QU_20231125/SwiftJ1727_QU_20231125_S_sampler.pkl',
    #DATA_DIR / 'SwiftJ1727_QU_20240210/SwiftJ1727_QU_20240210_S_sampler.pkl',
    #DATA_DIR / 'SwiftJ1727_QU_20240219/SwiftJ1727_QU_20240219_S_sampler.pkl',
    #DATA_DIR / 'SwiftJ1727_QU_20240225/SwiftJ1727_QU_20240225_S_sampler.pkl',
]

# Dates (YYYYMMDD) of the 9 epochs shown in the 3x3 FDF grid, in row-major order.
FDF_PANEL_DATES = [
    '20230904',
    '20230916',
    '20230923',
    '20231001',
    '20231006',
    '20231014',
    '20231016',
    '20231022',
    '20231028',
]

# Derived — do not edit directly.
_fdf_date_set = set(FDF_PANEL_DATES)
FDF_PKL_PATHS = sorted(
    [p for p in ALL_PKL_PATHS if re.search(r'(\d{8})', p.stem).group(1) in _fdf_date_set],
    key=lambda p: re.search(r'(\d{8})', p.stem).group(1),
)

# --- Epoch labelling ---------------------------------------------------------
LABEL_WINDOW_START = '2023-09-12'   # first date (inclusive) for sequential letters
LABEL_WINDOW_END   = '2023-10-22'   # last date (inclusive) for sequential letters
PRE_FLARE_LABEL    = 'Pre-flaring'
POST_FLARE_LABEL   = 'Post-flaring'
HIGHLIGHT_DATES    = ['20231006', '20231014']   # YYYYMMDD dates whose label and FDF panel are shown in red

# --- FDF panel controls ------------------------------------------------------
PHI_RANGE        = 250            # x-axis: +/-PHI_RANGE (rad/m^2)
Y_RANGE          = (1e-3, 30)     # log y-axis limits in per cent
N_PHI            = 4096           # Faraday-depth grid resolution
N_SAMPLES        = 100            # posterior samples drawn per panel
LAMBDA_SQ_REF    = 0.05           # reference wavelength squared (m^2) for ThickComponent amplitude normalisation
FDF_THIN_MS_SAMPLE = 4            # triangle marker size for posterior sample thin components
FDF_THIN_MS_BEST   = 3            # triangle marker size for best-fit thin components

# --- Polarisation-angle dial controls ----------------------------------------
DIAL_SCALE    = 0.375       # dial diameter as fraction of axes height (portrait)
DIAL_SCALE_LS = 0.300       # dial diameter in landscape (0.8x portrait)
DIAL_PAD      = 0.00        # padding from axes edge (axes fraction)
DIAL_SEP      = 0.00        # gap between adjacent dials (axes fraction)

# --- Light-curve time range --------------------------------------------------
LC_START = '2023-09-01T00:00:00'
LC_END   = '2023-11-26T00:00:00'

# --- Reference frequency and spectral index fitting --------------------------
MEERKAT_FREQ_GHZ       = 1.28
ALPHA_TIME_WINDOW_DAYS = 1.0 / 1440.0   # 1 minute: max separation for simultaneous group
ALPHA_CLIP_DRY_RUN     = True            # True: report clips but keep all data

# --- Smoothed flux-density band ----------------------------------------------
SMOOTH_GROUP_WINDOW = 0.1                    # days: group nearby observations
SMOOTH_LOG_WIDTH    = 0.10                   # dex half-width around daily mean
SMOOTH_CUTOFF       = '2023-11-26T00:00:00'  # extend to LC_END to include all plotted data
SMOOTH_COLOR        = 'grey'
SMOOTH_ALPHA        = 0.4

# --- Light-curve marker styles -----------------------------------------------
MK_COLOR   = '#C9A7E8'   # MeerKAT fill (pastel purple)
MK_EDGE    = 'k'
MK_MARKER  = '*'
MK_SIZE    = 150         # scatter marker size in points^2; scale for figure size

OTHER_COLOR = 'k'
OTHER_ALPHA = 0.35
OTHER_SIZE  = 20

ALPHA_MFC      = 'k'
ALPHA_MEC      = 'k'
ALPHA_ALPHA    = 0.35
ALPHA_MS       = 4
ALPHA_CAPSIZE  = 2
ALPHA_CAPTHICK = 0.8

# --- Main figure layout ------------------------------------------------------
LANDSCAPE        = False          # True: LC left / FDF right (16:9); False: LC top / FDF bottom (portrait)
FIG_SIZE         = (6.85, 9.21)   # portrait (MNRAS double-column)
FIG_SIZE_LS      = (13.0, 5.5)    # landscape (16:9 approx)
OUTER_HSPACE     = 0.09           # portrait: vertical gap between LC and FDF halves
OUTER_WSPACE     = 0.12           # landscape: horizontal gap between LC and FDF halves
LEFT_MARGIN      = 0.16
RIGHT_MARGIN     = 0.97
TOP_MARGIN       = 0.96
BOTTOM_MARGIN    = 0.07

LC_HSPACE        = 0.06           # vertical gap between LC sub-panels
LC_HEIGHT_RATIOS = [2, 1, 1]      # flux density : spectral index : pol fraction
LC_TICK_INTERVAL = 10             # major x-tick spacing in days

FDF_HSPACE       = 0.06           # vertical gap between FDF rows
FDF_WSPACE       = 0.03           # horizontal gap between FDF columns

# --- Font sizes (main figure) ------------------------------------------------
FONT_BASE        = 8
FONT_LABEL       = 9
FONT_TITLE       = 8
FONT_TICK        = 8
FONT_LEGEND      = 7
FONT_FIG_TITLE   = 9

FONT_LC_YLABEL    = 10  # LC y-axis labels
FONT_LC_LEGEND    = 7   # LC legends
FONT_PANEL_LABEL  = 9   # FDF panel text label (×1.5 for legibility)
FONT_ANNOT        = 11  # LC flux-panel letter annotations
FONT_SHARED_LABEL = 10  # shared Faraday-depth / |F| fig.text labels

# --- Font sizes (standalone figures) -----------------------------------------
FONT_FDF_AXIS          = 15
FONT_FDF_TITLE         = 16
FONT_STANDALONE_YLABEL = 14
FONT_STANDALONE_LEGEND = 11
FONT_STANDALONE_TICK   = 12   # tick labels on the standalone LC figure
FONT_STANDALONE_XLABEL = 14   # x-axis labels (UTC / MJD) on the standalone LC figure

# --- Standalone figure sizes -------------------------------------------------
FIG_SIZE_LC  = (8.33, 8)
FIG_SIZE_FDF = (7, 7)

# --- Output ------------------------------------------------------------------
OUTPUT_FORMAT         = 'pdf'
SAVE_PATH             = '../plots/results/results_fullevo'
SAVE_LIGHTCURVE_PANEL = True
LIGHTCURVE_SAVE_PATH  = '../plots/results/results_stokesI'
SAVE_FDF_PANELS       = False
FDF_PANELS_DIR        = '../plots/results/fdf_panels'
SHOW_PLOT             = False

# =============================================================================
# MATPLOTLIB RCPARAMS
# =============================================================================
plt.rcParams.update({
    'font.family':           'serif',
    'font.serif':            ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
    'mathtext.fontset':      'dejavuserif',
    'font.size':             FONT_BASE,
    'axes.labelsize':        14,
    'axes.titlesize':        FONT_TITLE,
    'axes.linewidth':        1.0,
    'xtick.major.width':     1.0,
    'ytick.major.width':     1.0,
    'xtick.minor.width':     0.6,
    'ytick.minor.width':     0.6,
    'xtick.major.size':      4,
    'ytick.major.size':      4,
    'xtick.minor.size':      2,
    'ytick.minor.size':      2,
    'xtick.direction':       'in',
    'ytick.direction':       'in',
    'xtick.top':             True,
    'xtick.bottom':          True,
    'ytick.left':            True,
    'ytick.right':           True,
    'xtick.minor.visible':   True,
    'ytick.minor.visible':   True,
    'xtick.labelsize':       FONT_TICK,
    'ytick.labelsize':       FONT_TICK,
    'legend.fontsize':       FONT_LEGEND,
    'legend.framealpha':     1.0,
    'legend.facecolor':      'white',
    'legend.edgecolor':      'black',
    'legend.fancybox':       False,
    'figure.titlesize':      FONT_FIG_TITLE,
    'lines.linewidth':       1.0,
    'lines.markersize':      3,
})

# =============================================================================
# UTILITY FUNCTIONS -- axis formatting
# =============================================================================

def plot2mjd(t):
    '''Convert from matplotlib plot date to mjd'''
    return Time(t, format="plot_date", scale='utc').mjd

def mjd2plot(mjd):
    '''Convert from mjd to matplotlib plot'''
    return Time(mjd, format="mjd", scale='utc').plot_date

def FormatAxis(ax, mjd, dt=10, interval=60, tick_fontsize=None, xlabel_fontsize=None):
    ax[0].set_xlabel('Observing Date (UTC)', fontfamily='serif',
                     fontsize=xlabel_fontsize)
    ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    ax[0].set_xlim(Time(mjd[0] - dt, format='mjd').datetime, Time(mjd[-1] + dt, format='mjd').datetime)
    ax[0].xaxis.set_label_position('top')
    xformatter = mdates.DateFormatter('%Y-%m-%d')
    plt.gcf().axes[0].xaxis.set_major_formatter(xformatter)
    ax[0].tick_params(axis='x', which='major', rotation=15, labeltop=True, labelbottom=False,
                      labelsize=tick_fontsize)
    plt.setp(ax[0].get_xticklabels(), rotation=15, ha='left',
             fontsize=tick_fontsize)

    mjd_ax = ax[-1].secondary_xaxis('bottom', functions=(plot2mjd, mjd2plot))
    mjd_ax.set_xlabel('Observing Date (MJD)', fontfamily='serif',
                      fontsize=xlabel_fontsize)
    mjd_ax.tick_params(which='major', direction='in', length=0.0, width=0.0,
                       labelsize=tick_fontsize)
    plt.draw()

    mjd_ticks = []
    labels = ax[0].get_xticklabels(which='major')
    for lab in labels:
        mjd_ticks.append(lab.get_text() + 'T00:00:00')
    mjd_ticks = (Time(mjd_ticks, format='isot').mjd).astype(int)
    mjd_ax.set_xticks(mjd_ticks, labels=mjd_ticks)
    if tick_fontsize is not None:
        mjd_ax.set_xticklabels(mjd_ticks, fontsize=tick_fontsize)
    return mjd_ax


def align_axis_x(ax, ax_target):
    """Align the x-extent of *ax* with *ax_target* (fixes secondary-axis shifts)."""
    posn_old, posn_target = ax.get_position(), ax_target.get_position()
    ax.set_position([posn_target.x0, posn_old.y0, posn_target.width, posn_old.height])


# =============================================================================
# UTILITY FUNCTIONS -- JSON / epoch parsing
# =============================================================================

def pkl_to_json_path(pkl_path):
    """Infer the JSON results path from a PKL sampler path."""
    return pkl_path.parent / pkl_path.name.replace('_sampler.pkl', '_results.json')


def parse_date_from_filename(filepath):
    """Extract YYYYMMDD from *filepath* stem and return (MJD, ISO date string)."""
    match = re.search(r'(\d{8})', filepath.stem)
    if match is None:
        raise ValueError(f"No 8-digit date in filename: {filepath.name}")
    d = match.group(1)
    iso_str = f"{d[:4]}-{d[4:6]}-{d[6:8]}"
    return Time(iso_str, format='iso').mjd, iso_str


def parse_epoch(json_path):
    """
    Parse a single-epoch JSON results file.

    Returns
    -------
    dict with keys:
        mjd, iso_date, model_type, filename,
        components (list of dicts with p0/psi_0 and hdi_99 errors),
        net_p0, net_p0_err_lower, net_p0_err_upper  (quadrature sums).
    """
    with open(json_path, 'r') as f:
        data = json.load(f)

    mjd, iso_date = parse_date_from_filename(json_path)
    comp_names = data['metadata']['component_names']
    model_type = data['metadata']['model_type']
    pp = data['processed_parameters']

    components = []
    for name in comp_names:
        p0_mode = pp[f'{name}_p0']['mode_1']['mode']
        p0_hdi  = pp[f'{name}_p0']['mode_1']['hdi_99']
        psi_mode = pp[f'{name}_psi_0']['mode_1']['mode']
        psi_hdi  = pp[f'{name}_psi_0']['mode_1']['hdi_99']
        components.append({
            'name':          name,
            'p0':            p0_mode,
            'p0_err_lower':  p0_mode - p0_hdi[0],
            'p0_err_upper':  p0_hdi[1] - p0_mode,
            'psi_0':         psi_mode,
            'psi_0_hdi_99':  psi_hdi,
        })

    p0_vals = np.array([c['p0']           for c in components])
    err_lo  = np.array([c['p0_err_lower'] for c in components])
    err_hi  = np.array([c['p0_err_upper'] for c in components])

    sif = data.get('stokes_i_fit', None)

    return {
        'mjd':                    mjd,
        'iso_date':               iso_date,
        'model_type':             model_type,
        'filename':               json_path.name,
        'components':             components,
        'net_p0':                 np.sqrt(np.sum(p0_vals**2)),
        'net_p0_err_lower':       np.sqrt(np.sum(err_lo**2)),
        'net_p0_err_upper':       np.sqrt(np.sum(err_hi**2)),
        'alpha_nu_intraband':     sif['alpha_nu']     if sif else None,
        'alpha_nu_intraband_err': sif['alpha_nu_err'] if sif else None,
    }


# =============================================================================
# UTILITY FUNCTIONS -- FDF plotting helpers
# =============================================================================

def _get_psi_data(comp_name, epoch):
    """Return {mode, hdi_99} for psi_0 of a named component."""
    for c in epoch['components']:
        if c['name'] == comp_name:
            return {'mode': c['psi_0'], 'hdi_99': c['psi_0_hdi_99']}
    raise KeyError(f"Component '{comp_name}' not found in {epoch['filename']}")


def _get_phi_median(comp, epoch):
    """Return median Faraday depth (phi_rm or phi_peak) for a component."""
    json_path = pkl_to_json_path(epoch['pkl_path'])
    with open(json_path, 'r') as f:
        data = json.load(f)
    pp = data['processed_parameters']
    key = f'{comp.name}_phi_rm' if hasattr(comp, 'phi_rm') else f'{comp.name}_phi_peak'
    return pp[key]['median']


def _draw_dial(ax, epoch, comps_sorted, dial_scale, dial_pad, dial_sep):
    """
    Draw one polarisation-angle dial per Faraday component in the top-left
    corner of *ax*.

    Angle convention: 0 deg = North (up), positive = East-of-North (CCW in
    astronomical PA).  The dial shows the posterior mode EVPA as a solid line
    and the 99% HDI as a shaded wedge, both duplicated at +/-180 deg for the
    pi-ambiguity.
    """
    fig_w, fig_h = ax.get_figure().get_size_inches()
    ax_pos       = ax.get_position()
    ax_w_in      = ax_pos.width  * fig_w
    ax_h_in      = ax_pos.height * fig_h
    dial_scale_x = dial_scale * (ax_h_in / ax_w_in)

    for i, comp in enumerate(comps_sorted):
        x0 = dial_pad + i * (dial_scale_x * 0.85 + dial_sep)
        y0 = 1.0 - dial_pad - dial_scale

        ax_d = ax.inset_axes([x0, y0, dial_scale_x, dial_scale])
        ax_d.set_xlim(-1.3, 1.3)
        ax_d.set_ylim(-1.3, 1.3)
        ax_d.axis('off')
        ax_d.patch.set_visible(True)

        ax_d.add_patch(Circle((0, 0), 1.0, fill=True, facecolor='white',
                               edgecolor='black', linewidth=1.0, zorder=50))
        ax_d.plot([0, 0], [-1, 1], 'k:', linewidth=1.5, zorder=10000)

        psi_data       = _get_psi_data(comp.name, epoch)
        psi_med        = psi_data['mode']
        psi_lo, psi_hi = psi_data['hdi_99']

        mpl_t1, mpl_t2 = 90.0 + psi_lo, 90.0 + psi_hi

        ax_d.add_patch(Wedge((0, 0), 1.0, mpl_t1,         mpl_t2,
                             color='C0', alpha=0.30, zorder=1000))
        ax_d.add_patch(Wedge((0, 0), 1.0, mpl_t1 + 180.0, mpl_t2 + 180.0,
                             color='C0', alpha=0.30, zorder=1000))

        mpl_med = np.radians(90.0 + psi_med)
        cx, sx  = np.cos(mpl_med), np.sin(mpl_med)
        ax_d.plot([0,  cx], [0,  sx], color='C0', linewidth=1.2, zorder=2000)
        ax_d.plot([0, -cx], [0, -sx], color='C0', linewidth=1.2, zorder=2000)


# =============================================================================
# DATA LOADING
# =============================================================================

start_datetime = Time(LC_START, format='isot').datetime
end_datetime   = Time(LC_END,   format='isot').datetime

# --- Radio light-curve CSV ---------------------------------------------------
big_data = pd.read_csv(RADIO_CSV)
big_data['Midpoint (DT)'] = Time(big_data['Midpoint (MJD)'], format='mjd').datetime
big_data = (big_data[(big_data['Frequency (GHz)'] >= 1.0) &
                     (big_data['Frequency (GHz)'] <= 15.0)]
            .reset_index(drop=True))

# --- Inter-band spectral indices ---------------------------------------------
# Group non-MeerKAT data by (telescope, simultaneous group).  Groups with >=2
# distinct frequencies get a weighted log-space power-law fit; uncertainty via
# Monte Carlo.  Groups with outlier alpha_err (>3 sigma robust) are flagged.

N_MC = 500

row_alpha     = np.full(len(big_data), np.nan)
row_alpha_err = np.full(len(big_data), np.nan)

non_mk_tmp = big_data[big_data['Telescope'].isin(['VLA', 'ATCA', 'ATA'])].copy()

group_fits = []

for tele_g, tele_grp in non_mk_tmp.groupby('Telescope'):
    tele_grp = tele_grp.sort_values('Midpoint (MJD)')
    mjds = tele_grp['Midpoint (MJD)'].values
    sim_group = np.concatenate([[0], np.cumsum(np.diff(mjds) > ALPHA_TIME_WINDOW_DAYS)])
    tele_grp = tele_grp.copy()
    tele_grp['sim_group'] = sim_group

    for gid, grp in tele_grp.groupby('sim_group'):
        if grp['Frequency (GHz)'].nunique() < 2:
            continue
        nu    = grp['Frequency (GHz)'].values.astype(float)
        S     = grp['Integrated flux (mJy)'].values.astype(float)
        dS    = np.sqrt(grp['Error (mJy)'].values.astype(float)**2 + (0.03 * S)**2)
        ln_nu = np.log(nu)
        ln_S  = np.log(S)
        w     = S / dS

        try:
            a_fit = float(np.polyfit(ln_nu, ln_S, 1, w=w)[0])
        except (np.linalg.LinAlgError, ValueError):
            continue

        rng = np.random.default_rng(seed=int(grp['Midpoint (MJD)'].iloc[0] * 1440) % (2**31))
        mc_alphas = []
        for _ in range(N_MC):
            S_mc = S + dS * rng.standard_normal(len(S))
            if np.any(S_mc <= 0):
                continue
            try:
                mc_alphas.append(float(np.polyfit(ln_nu, np.log(S_mc), 1, w=w)[0]))
            except (np.linalg.LinAlgError, ValueError):
                continue
        if len(mc_alphas) < N_MC // 2:
            continue
        a_err = float(np.std(mc_alphas))

        group_fits.append((grp.index, grp['Midpoint (MJD)'].mean(),
                           tele_g, a_fit, a_err, len(grp)))

all_aerrs = np.array([r[4] for r in group_fits]) if group_fits else np.array([])
if len(all_aerrs) > 3:
    med_ae   = np.median(all_aerrs)
    mad_ae   = np.median(np.abs(all_aerrs - med_ae))
    aerr_thr = med_ae + 3.0 * (1.4826 * mad_ae)
else:
    aerr_thr = np.inf

alpha_records = []
for idx, mean_mjd, tele_g, a_fit, a_err, n_obs in group_fits:
    iso_g = Time(mean_mjd, format='mjd').isot
    if a_err > aerr_thr:
        tag = '[DRY-CLIP]' if ALPHA_CLIP_DRY_RUN else '[CLIP]'
        print(f'  {tag}  {iso_g}  {tele_g:<8}  alpha={a_fit:.4f}  '
              f'alpha_err={a_err:.4f}  threshold={aerr_thr:.4f}  n={n_obs}')
        if not ALPHA_CLIP_DRY_RUN:
            continue
    row_alpha[idx]     = a_fit
    row_alpha_err[idx] = a_err
    alpha_records.append((mean_mjd, a_fit, a_err))

alpha_records.sort(key=lambda x: x[0])
if alpha_records:
    alpha_mjds_plot = np.array([r[0] for r in alpha_records])
    alpha_vals_plot = np.array([r[1] for r in alpha_records])
    alpha_errs_plot = np.array([r[2] for r in alpha_records])
    alpha_dt_plot   = Time(alpha_mjds_plot, format='mjd').datetime
else:
    alpha_mjds_plot = alpha_vals_plot = alpha_errs_plot = np.array([])
    alpha_dt_plot   = np.array([])

# --- Extrapolate all fluxes to MEERKAT_FREQ_GHZ (1.28 GHz) ------------------
is_mk_1p28 = ((big_data['Telescope'] == 'MeerKAT') &
              (big_data['Frequency (GHz)'] == MEERKAT_FREQ_GHZ)).values
has_alpha  = ~np.isnan(row_alpha)
valid      = is_mk_1p28 | has_alpha

alpha_for_scale    = np.where(is_mk_1p28, 0.0, row_alpha)
scale              = (MEERKAT_FREQ_GHZ / big_data['Frequency (GHz)'].values) ** alpha_for_scale
ln_ratio           = np.log(MEERKAT_FREQ_GHZ / big_data['Frequency (GHz)'].values)
alpha_err_for_prop = np.where(is_mk_1p28, 0.0, row_alpha_err)
S_obs              = big_data['Integrated flux (mJy)'].values
dS_obs             = big_data['Error (mJy)'].values

flux_1p28 = S_obs * scale
err_1p28  = np.abs(scale) * np.sqrt(dS_obs**2 + (S_obs * ln_ratio * alpha_err_for_prop)**2)

big_data['flux_1p28'] = np.where(valid, flux_1p28, np.nan)
big_data['err_1p28']  = np.where(valid, err_1p28,  np.nan)

radio_teles = ['VLA', 'MeerKAT', 'ATCA', 'ATA']

# --- Write extrapolation diagnostic table ------------------------------------
_extrap_rows = []
for _i, _row in big_data.iterrows():
    if not (start_datetime <= _row['Midpoint (DT)'] <= end_datetime):
        continue
    if _row['Telescope'] not in radio_teles:
        continue
    if pd.isna(_row['flux_1p28']):
        continue
    _pos   = big_data.index.get_loc(_i)
    _iso   = Time(_row['Midpoint (MJD)'], format='mjd').isot
    _tele  = _row['Telescope']
    _f_obs = _row['Integrated flux (mJy)']
    _f_scl = _row['flux_1p28']
    _freq  = _row['Frequency (GHz)']
    _alpha = row_alpha[_pos]
    _aerr  = row_alpha_err[_pos]
    _extrap_rows.append(
        f'{_iso}  {_tele:<10}  {_f_obs:12.4f}  {_f_scl:12.4f}  '
        f'{_freq:8.4f}  {_alpha:9.5f}  {_aerr:10.4f}')

_extrap_header = (
    f'# Flux extrapolation to {MEERKAT_FREQ_GHZ} GHz   S(1.28) = S(f) * (1.28/f)^alpha\n'
    f'# NaN scaled flux means <2 frequencies on same day (non-MeerKAT) or MeerKAT non-ref freq\n'
    f'# alpha and alpha_err from same-day scipy.curve_fit in log-space\n'
    f'#\n'
    f'# {"Date (ISOT)":<24}  {"Telescope":<10}  {"S_obs (mJy)":>12}  '
    f'{"S_1p28 (mJy)":>12}  {"Freq (GHz)":>10}  {"alpha":>9}  {"alpha_err":>10}'
)
_extrap_path = LC_DIR / 'extrapolation_diagnostic.txt'
with open(_extrap_path, 'w') as _ef:
    _ef.write(_extrap_header + '\n')
    _ef.write('\n'.join(_extrap_rows) + '\n')
print(f'Wrote extrapolation table -> {_extrap_path}')

# --- State-transition dates --------------------------------------------------
transition_datetime      = Time([60389, 60222], format='mjd').datetime
highlight_sep20_datetime = Time('2023-09-20T01:00:00', format='isot').datetime

# --- Smoothed radio band (inverse-variance weighted daily averages) ----------
time_mask  = ((big_data['Midpoint (DT)'] >= start_datetime) &
              (big_data['Midpoint (DT)'] <= end_datetime))
mk_non_ref = ((big_data['Telescope'] == 'MeerKAT') &
              (big_data['Frequency (GHz)'] != MEERKAT_FREQ_GHZ))
smooth_pool = (big_data[big_data['Telescope'].isin(radio_teles) & time_mask & ~mk_non_ref]
               .sort_values('Midpoint (MJD)').copy())

cutoff_dt = Time(SMOOTH_CUTOFF, format='isot').datetime

smoothed_times, smoothed_flux = [], []
used_mask = np.zeros(len(smooth_pool), dtype=bool)
for i in range(len(smooth_pool)):
    if used_mask[i] or smooth_pool['Midpoint (DT)'].iloc[i] > cutoff_dt:
        continue
    current_mjd = smooth_pool.iloc[i]['Midpoint (MJD)']
    group_mask  = np.abs(smooth_pool['Midpoint (MJD)'] - current_mjd) <= SMOOTH_GROUP_WINDOW
    group       = smooth_pool[group_mask].dropna(subset=['flux_1p28', 'err_1p28'])
    if len(group) == 0:
        used_mask[group_mask.values] = True
        continue
    weights = 1.0 / (group['err_1p28']**2 + 1e-10)
    smoothed_flux.append(np.average(group['flux_1p28'], weights=weights))
    smoothed_times.append(
        Time(np.average(group['Midpoint (MJD)'], weights=weights), format='mjd').datetime)
    used_mask[group_mask.values] = True
smoothed_times = np.array(smoothed_times)
smoothed_flux  = np.array(smoothed_flux)

mjd_range = np.array([Time(start_datetime).mjd, Time(end_datetime).mjd])

output_format = OUTPUT_FORMAT.strip().lower()
if output_format not in {'pdf', 'png'}:
    raise ValueError("OUTPUT_FORMAT must be 'pdf' or 'png'.")


def _path_with_output_format(path_like):
    """Return a path with the correct output-format suffix."""
    return Path(path_like).with_suffix(f'.{output_format}')


# --- Parse epochs ------------------------------------------------------------
# all_epochs: all available JSONs (ALL_PKL_PATHS) for the pol fraction LC panel
# fdf_epochs: the 9 selected FDF panel epochs, as a sorted subset of all_epochs

all_epochs = []
for pkl_path in ALL_PKL_PATHS:
    json_path = pkl_to_json_path(pkl_path)
    epoch = parse_epoch(json_path)
    epoch['pkl_path'] = pkl_path
    all_epochs.append(epoch)
all_epochs.sort(key=lambda e: e['mjd'])

fdf_epoch_stems = {p.stem for p in FDF_PKL_PATHS}
fdf_epochs = [e for e in all_epochs if e['pkl_path'].stem in fdf_epoch_stems]
fdf_epochs.sort(key=lambda e: e['mjd'])

# --- Assign panel labels to FDF epochs ---------------------------------------
label_window_start_dt = Time(LABEL_WINDOW_START, format='iso').datetime
label_window_end_dt   = Time(LABEL_WINDOW_END,   format='iso').datetime

letter_idx = 0
for epoch in fdf_epochs:
    epoch_dt = Time(epoch['mjd'], format='mjd').datetime
    if label_window_start_dt <= epoch_dt <= label_window_end_dt:
        epoch['panel_label'] = string.ascii_uppercase[letter_idx]
        letter_idx += 1
    elif epoch_dt < label_window_start_dt:
        epoch['panel_label'] = PRE_FLARE_LABEL
    else:
        epoch['panel_label'] = POST_FLARE_LABEL

n_fdf = len(fdf_epochs)
print(f'Loaded {len(all_epochs)} total epochs, {n_fdf} FDF panel epochs.')

# =============================================================================
# STANDALONE LIGHT-CURVE FIGURE
# =============================================================================

if SAVE_LIGHTCURVE_PANEL and LIGHTCURVE_SAVE_PATH is not None:
    lc_save_path = _path_with_output_format(LIGHTCURVE_SAVE_PATH)
    lc_save_path.parent.mkdir(parents=True, exist_ok=True)

    lc_fig = plt.figure(figsize=FIG_SIZE_LC)
    lc_gs  = gridspec.GridSpec(2, 1, figure=lc_fig, hspace=0.02, height_ratios=[2, 1],
                               left=0.12, right=0.95, top=0.92, bottom=0.12)
    lc_ax_new = [lc_fig.add_subplot(lc_gs[0])]
    lc_ax_new.append(lc_fig.add_subplot(lc_gs[1], sharex=lc_ax_new[0]))

    if len(smoothed_flux) > 1:
        lc_ax_new[0].fill_between(
            smoothed_times,
            smoothed_flux / 10**SMOOTH_LOG_WIDTH,
            smoothed_flux * 10**SMOOTH_LOG_WIDTH,
            color=SMOOTH_COLOR, alpha=SMOOTH_ALPHA, zorder=1,
            label=f'Daily avg (\u00b1{SMOOTH_LOG_WIDTH} dex)')

    _other = big_data[big_data['Telescope'].isin(['VLA', 'ATCA', 'ATA']) & time_mask]
    lc_ax_new[0].scatter(_other['Midpoint (DT)'], _other['flux_1p28'],
                         marker='o', s=OTHER_SIZE * 2, color=OTHER_COLOR,
                         alpha=OTHER_ALPHA, ec='none', zorder=500)
    lc_ax_new[0].errorbar(_other['Midpoint (DT)'], _other['flux_1p28'],
                          _other['err_1p28'],
                          fmt='none', ecolor='k', alpha=OTHER_ALPHA, zorder=400)
    _mk = big_data[
        (big_data['Telescope'] == 'MeerKAT') &
        (big_data['Frequency (GHz)'] == MEERKAT_FREQ_GHZ) &
        time_mask]
    lc_ax_new[0].scatter(_mk['Midpoint (DT)'], _mk['flux_1p28'],
                         marker=MK_MARKER, s=MK_SIZE * 2, color=MK_COLOR,
                         ec=MK_EDGE, zorder=10000, label='MeerKAT 1.28 GHz')
    lc_ax_new[0].errorbar(_mk['Midpoint (DT)'], _mk['flux_1p28'], _mk['err_1p28'],
                          fmt='--', ecolor='k', zorder=5000, color='k')
    lc_ax_new[0].set_yscale('log')
    lc_ax_new[0].set_ylabel('Flux Density (1.28 GHz; mJy)',
                             fontsize=FONT_STANDALONE_YLABEL)
    lc_ax_new[0].set_ylim(3, 1.4e3)
    lc_ax_new[0].legend(fontsize=FONT_STANDALONE_LEGEND, loc='upper right')

    if len(alpha_dt_plot):
        lc_ax_new[1].errorbar(alpha_dt_plot, alpha_vals_plot, alpha_errs_plot,
                              fmt='o', mfc=ALPHA_MFC, mec=ALPHA_MEC, ecolor='k',
                              ms=ALPHA_MS * 1.5, capsize=ALPHA_CAPSIZE * 1.5,
                              capthick=ALPHA_CAPTHICK,
                              linestyle='none', alpha=ALPHA_ALPHA, zorder=10,
                              label='Inter-band')
    _intra_lc = [(Time(e['mjd'], format='mjd').datetime,
                  e['alpha_nu_intraband'],
                  e['alpha_nu_intraband_err'])
                 for e in all_epochs if e['alpha_nu_intraband'] is not None]
    if _intra_lc:
        _idt_lc, _ialpha_lc, _ierr_lc = zip(*_intra_lc)
        lc_ax_new[1].errorbar(list(_idt_lc), list(_ialpha_lc), list(_ierr_lc),
                              fmt=MK_MARKER, mfc=MK_COLOR, mec=MK_EDGE, ecolor=MK_EDGE,
                              ms=np.sqrt(MK_SIZE * 2), capsize=ALPHA_CAPSIZE * 1.5,
                              capthick=ALPHA_CAPTHICK,
                              linestyle='none', zorder=20, label='MeerKAT intra-band')
        lc_ax_new[1].legend(fontsize=FONT_STANDALONE_LEGEND, loc='upper right')
    lc_ax_new[1].axhline(0.0, ls=':', c='k')
    lc_ax_new[1].set_ylabel(r'$\alpha_\nu$', fontsize=FONT_STANDALONE_YLABEL)
    lc_ax_new[1].set_ylim(-1.7, 1.2)

    lc_ax_new[0].set_xlim(start_datetime, end_datetime)
    for a in lc_ax_new:
        if start_datetime <= transition_datetime[0] <= end_datetime:
            a.axvline(transition_datetime[0], ls='--', c='k', lw=1, zorder=0)
        if start_datetime <= transition_datetime[1] <= end_datetime:
            a.axvline(transition_datetime[1], ls=':', c='k', lw=1)
        if start_datetime <= highlight_sep20_datetime <= end_datetime:
            a.axvline(highlight_sep20_datetime, ls='--', c='k', lw=1, zorder=0)
    for a in lc_ax_new:
        a.tick_params(labelbottom=False)
    FormatAxis(np.array(lc_ax_new), mjd=mjd_range, interval=10, dt=0,
               tick_fontsize=FONT_STANDALONE_TICK, xlabel_fontsize=FONT_STANDALONE_XLABEL)
    for _a in lc_ax_new:
        _a.tick_params(axis='y', labelsize=FONT_STANDALONE_TICK)
    lc_ax_new[-1].tick_params(labelbottom=False)
    align_axis_x(lc_ax_new[1], lc_ax_new[0])

    lc_fig.savefig(lc_save_path, dpi=300, bbox_inches='tight')
    print(f'Saved LC panel -> {lc_save_path}')
    plt.close(lc_fig)


# =============================================================================
# MAIN FIGURE CONSTRUCTION
# =============================================================================

fig = plt.figure(figsize=FIG_SIZE_LS if LANDSCAPE else FIG_SIZE)

if LANDSCAPE:
    # Left half = LC, right half = FDF 3x3 grid
    outer_gs = gridspec.GridSpec(
        1, 2, figure=fig,
        wspace=OUTER_WSPACE,
        left=LEFT_MARGIN, right=RIGHT_MARGIN,
        top=TOP_MARGIN, bottom=BOTTOM_MARGIN,
        width_ratios=[2, 3],
    )
    lc_slot  = outer_gs[0]
    fdf_slot = outer_gs[1]
else:
    # Top half = LC, bottom half = FDF 3x3 grid
    outer_gs = gridspec.GridSpec(
        2, 1, figure=fig,
        hspace=OUTER_HSPACE,
        left=LEFT_MARGIN, right=RIGHT_MARGIN,
        top=TOP_MARGIN, bottom=BOTTOM_MARGIN,
        height_ratios=[1, 1],
    )
    lc_slot  = outer_gs[0]
    fdf_slot = outer_gs[1]

# LC axes created first so plt.gcf().axes[0] == ax_lc[0] when FormatAxis runs
lc_gs = lc_slot.subgridspec(3, 1, hspace=LC_HSPACE, height_ratios=LC_HEIGHT_RATIOS)
ax_lc = [fig.add_subplot(lc_gs[0])]
ax_lc.append(fig.add_subplot(lc_gs[1], sharex=ax_lc[0]))
ax_lc.append(fig.add_subplot(lc_gs[2], sharex=ax_lc[0]))

# FDF 3x3 grid
fdf_gs      = fdf_slot.subgridspec(3, 3, hspace=FDF_HSPACE, wspace=FDF_WSPACE)
ax_fdf_list = [fig.add_subplot(fdf_gs[r, c]) for r in range(3) for c in range(3)]


# =============================================================================
# FDF PANELS -- posterior samples + best-fit model + EVPA dials
# =============================================================================

phi_grid        = np.linspace(-PHI_RANGE, PHI_RANGE, N_PHI)
_active_dial_scale = DIAL_SCALE_LS if LANDSCAPE else DIAL_SCALE

for idx, epoch in enumerate(fdf_epochs):
    ax       = ax_fdf_list[idx]
    pkl_path = epoch['pkl_path']

    results    = load_results(pkl_path)
    rep_params = results.get_representative_sample(percentile=1.0, method='median')
    best_model = copy.deepcopy(results.setup.model)
    update_model_from_params(best_model, rep_params, results.setup.param_names,
                             lambda_sq_ref=LAMBDA_SQ_REF)

    comps_sorted = sorted(best_model.components,
                          key=lambda c: _get_phi_median(c, epoch))

    weights_norm = results.weights / results.weights.sum()
    sample_idx   = np.random.choice(len(results.samples),
                                    size=min(N_SAMPLES, len(results.samples)),
                                    replace=False, p=weights_norm)

    for si in sample_idx:
        sample_model = copy.deepcopy(results.setup.model)
        update_model_from_params(sample_model, results.samples[si],
                                 results.setup.param_names,
                                 lambda_sq_ref=LAMBDA_SQ_REF)
        for comp in sample_model.components:
            amp = np.abs(comp.compute_fdf(phi_grid)) * 100.0
            amp = np.clip(amp, None, Y_RANGE[1])
            if isinstance(comp, ThinComponent):
                ax.plot(phi_grid, amp, color='C0', lw=0.5, alpha=0.05, zorder=1)
                ax.plot(comp.phi_rm, np.clip(np.max(amp), Y_RANGE[0], Y_RANGE[1]),
                        '^', color='C0', ms=FDF_THIN_MS_SAMPLE, alpha=0.05, zorder=2)
            else:
                ax.plot(phi_grid, amp, color='C0', lw=0.5, alpha=0.05, zorder=1)

    for comp in best_model.components:
        comp_amp = np.abs(comp.compute_fdf(phi_grid)) * 100.0
        ax.plot(phi_grid, comp_amp, 'k-', lw=1.5, zorder=5)
        if isinstance(comp, ThinComponent):
            ax.plot(comp.phi_rm, np.max(comp_amp), '^',
                    color='black', ms=FDF_THIN_MS_BEST, zorder=6)

    ax.set_yscale('log')
    ax.set_xlim(-PHI_RANGE, PHI_RANGE)
    ax.set_ylim(*Y_RANGE)
    ax.axvline(0, color='black', ls=':', lw=1.0)

    _hl_color = 'red' if epoch['iso_date'].replace('-', '') in HIGHLIGHT_DATES else 'black'

    _draw_dial(ax, epoch, comps_sorted, _active_dial_scale, DIAL_PAD, DIAL_SEP)

    ax.text(0.96, 0.96, epoch['panel_label'], transform=ax.transAxes,
            fontsize=FONT_PANEL_LABEL, fontweight='bold', ha='right', va='top',
            color=_hl_color,
            bbox=dict(facecolor='white', edgecolor=_hl_color, boxstyle='round,pad=0.2'))

    print(f'  [{idx+1}/{n_fdf}] Plotted FDF: {pkl_path.stem}')

    # --- Save individual FDF panel as standalone figure ----------------------
    if SAVE_FDF_PANELS:
        fdf_fig, fdf_ax = plt.subplots(figsize=FIG_SIZE_FDF)
        fdf_fig.patch.set_facecolor('white')

        for si in sample_idx:
            sample_model = copy.deepcopy(results.setup.model)
            update_model_from_params(sample_model, results.samples[si],
                                     results.setup.param_names,
                                     lambda_sq_ref=LAMBDA_SQ_REF)
            for comp in sample_model.components:
                amp = np.abs(comp.compute_fdf(phi_grid)) * 100.0
                amp = np.clip(amp, None, Y_RANGE[1])
                if isinstance(comp, ThinComponent):
                    fdf_ax.plot(phi_grid, amp, color='C0', lw=0.5, alpha=0.05, zorder=1)
                    fdf_ax.plot(comp.phi_rm,
                                np.clip(np.max(amp), Y_RANGE[0], Y_RANGE[1]),
                                '^', color='C0', ms=FDF_THIN_MS_SAMPLE, alpha=0.05, zorder=2)
                else:
                    fdf_ax.plot(phi_grid, amp, color='C0', lw=0.5, alpha=0.05, zorder=1)

        for comp in best_model.components:
            comp_amp = np.abs(comp.compute_fdf(phi_grid)) * 100.0
            fdf_ax.plot(phi_grid, comp_amp, 'k-', lw=2.0, zorder=5)
            if isinstance(comp, ThinComponent):
                fdf_ax.plot(comp.phi_rm, np.max(comp_amp), '^',
                            color='black', ms=FDF_THIN_MS_BEST, zorder=6)

        fdf_ax.set_yscale('log')
        fdf_ax.set_xlim(-PHI_RANGE, PHI_RANGE)
        fdf_ax.set_ylim(*Y_RANGE)
        fdf_ax.axvline(0, color='black', ls=':', lw=1.5)
        fdf_ax.set_xlabel(r'Faraday Depth $\phi_f\;(\mathrm{rad\,m^{-2}})$',
                          fontsize=FONT_FDF_AXIS)
        fdf_ax.set_ylabel(r'$|F(\phi_f)|$ (%) $(\mathrm{rad\,m^{-2}})^{-1}$', fontsize=FONT_FDF_AXIS)
        fdf_ax.set_title(f'FDF: {epoch["iso_date"]} (MJD {epoch["mjd"]:.1f})',
                         fontsize=FONT_FDF_TITLE, fontweight='bold')

        fdf_fig.canvas.draw()
        _draw_dial(fdf_ax, epoch, comps_sorted, DIAL_SCALE, DIAL_PAD, DIAL_SEP)

        fdf_dir = Path(FDF_PANELS_DIR)
        fdf_dir.mkdir(parents=True, exist_ok=True)
        label_clean  = epoch['panel_label'].replace('/', '_').replace(' ', '_')
        fdf_filename = f'fdf_{label_clean}_{epoch["iso_date"]}.{output_format}'
        fdf_path     = fdf_dir / fdf_filename
        fdf_fig.savefig(fdf_path, dpi=300, bbox_inches='tight')
        print(f'      Saved FDF panel -> {fdf_path}')
        plt.close(fdf_fig)


# =============================================================================
# FDF AXIS CLEANUP -- shared labels, tick-label visibility
# =============================================================================

for idx in range(n_fdf):
    ax  = ax_fdf_list[idx]
    row = idx // 3
    col = idx  % 3
    if row != 2:
        ax.tick_params(labelbottom=False)
    if col != 0:
        ax.tick_params(labelleft=False)

fig.canvas.draw()
pos_bl = ax_fdf_list[6].get_position()   # bottom-left panel
pos_tl = ax_fdf_list[0].get_position()   # top-left panel
pos_br = ax_fdf_list[8].get_position()   # bottom-right panel

_ylabel_offset = 0.03 if LANDSCAPE else 0.08   # tighter in landscape

fig.text((pos_bl.x0 + pos_br.x1) / 2.0, pos_bl.y0 - 0.025,
         r'Faraday Depth $\phi_f\;(\mathrm{rad\,m^{-2}})$',
         ha='center', va='top', fontsize=FONT_SHARED_LABEL)
fig.text(pos_bl.x0 - _ylabel_offset, (pos_bl.y0 + pos_tl.y1) / 2.0,
         r'$|F(\phi_f)|$ (%) $(\mathrm{rad\,m^{-2}})^{-1}$',
         ha='center', va='center', rotation='vertical', fontsize=FONT_SHARED_LABEL)


# =============================================================================
# LC PANEL 1 -- Flux density
# =============================================================================

if len(smoothed_flux) > 1:
    ax_lc[0].fill_between(
        smoothed_times,
        smoothed_flux / 10**SMOOTH_LOG_WIDTH,
        smoothed_flux * 10**SMOOTH_LOG_WIDTH,
        color=SMOOTH_COLOR, alpha=SMOOTH_ALPHA, zorder=1,
        label=f'Daily avg (\u00b1{SMOOTH_LOG_WIDTH} dex)')

other_data = big_data[big_data['Telescope'].isin(['VLA', 'ATCA', 'ATA']) & time_mask]
ax_lc[0].scatter(other_data['Midpoint (DT)'], other_data['flux_1p28'],
                 marker='o', s=OTHER_SIZE, color=OTHER_COLOR,
                 alpha=OTHER_ALPHA, ec='none', zorder=500)
ax_lc[0].errorbar(other_data['Midpoint (DT)'], other_data['flux_1p28'],
                  other_data['err_1p28'],
                  fmt='none', ecolor='k', alpha=OTHER_ALPHA, zorder=400)

mk_data = big_data[
    (big_data['Telescope'] == 'MeerKAT') &
    (big_data['Frequency (GHz)'] == MEERKAT_FREQ_GHZ) &
    time_mask]
ax_lc[0].scatter(mk_data['Midpoint (DT)'], mk_data['flux_1p28'],
                 marker=MK_MARKER, s=MK_SIZE, color=MK_COLOR,
                 ec=MK_EDGE, zorder=10000, label='MeerKAT 1.28 GHz')
ax_lc[0].errorbar(mk_data['Midpoint (DT)'], mk_data['flux_1p28'],
                  mk_data['err_1p28'],
                  fmt='--', ecolor='k', zorder=5000, color='k')

ax_lc[0].set_yscale('log')
ax_lc[0].set_ylabel('Flux Density\n(1.28 GHz; mJy)', fontsize=FONT_LC_YLABEL)
ax_lc[0].set_ylim(3, 1.4e3)
ax_lc[0].legend(fontsize=FONT_LC_LEGEND, loc='upper right')

# Letter annotations (A-G only) on flux panel
mk_mjd_arr  = np.array(mk_data['Midpoint (MJD)'])
mk_flux_arr = np.array(mk_data['flux_1p28'])

for epoch in fdf_epochs:
    label = epoch['panel_label']
    if not (len(label) == 1 and label.isalpha()):
        continue
    epoch_dt = Time(epoch['mjd'], format='mjd').datetime
    if not (start_datetime <= epoch_dt <= end_datetime):
        continue
    dists   = np.abs(mk_mjd_arr - epoch['mjd'])
    nearest = np.argmin(dists)
    if dists[nearest] > 1.0:
        continue
    x_ann = mk_data['Midpoint (DT)'].iloc[nearest]
    y_ann = mk_flux_arr[nearest]
    _ann_color = 'red' if epoch['iso_date'].replace('-', '') in HIGHLIGHT_DATES else 'black'
    ax_lc[0].annotate(
        label, (x_ann, y_ann),
        xytext=(2, 8), textcoords='offset points',
        fontsize=FONT_ANNOT, fontweight='bold', ha='left', va='bottom',
        color=_ann_color,
        zorder=20000,
        bbox=dict(facecolor='white', edgecolor='none', pad=1.0))


# =============================================================================
# LC PANEL 2 -- Spectral index
# =============================================================================

if len(alpha_dt_plot):
    ax_lc[1].errorbar(alpha_dt_plot, alpha_vals_plot, alpha_errs_plot,
                      fmt='o', mfc=ALPHA_MFC, mec=ALPHA_MEC, ecolor='k',
                      ms=ALPHA_MS, capsize=ALPHA_CAPSIZE, capthick=ALPHA_CAPTHICK,
                      linestyle='none', alpha=ALPHA_ALPHA, zorder=10,
                      label='Inter-band')

_intra = [(Time(e['mjd'], format='mjd').datetime,
           e['alpha_nu_intraband'],
           e['alpha_nu_intraband_err'])
          for e in all_epochs if e['alpha_nu_intraband'] is not None]
if _intra:
    _idt, _ialpha, _ierr = zip(*_intra)
    ax_lc[1].errorbar(list(_idt), list(_ialpha), list(_ierr),
                      fmt=MK_MARKER, mfc=MK_COLOR, mec=MK_EDGE, ecolor=MK_EDGE,
                      ms=np.sqrt(MK_SIZE), capsize=ALPHA_CAPSIZE, capthick=ALPHA_CAPTHICK,
                      linestyle='none', zorder=20, label='MeerKAT intra-band')
    ax_lc[1].legend(fontsize=FONT_LC_LEGEND, loc='upper right')

ax_lc[1].axhline(0.0, ls=':', c='k')
ax_lc[1].set_ylabel(r'$\alpha_\nu$', fontsize=FONT_LC_YLABEL)
ax_lc[1].set_ylim(-1.7, 1.2)


# =============================================================================
# LC PANEL 3 -- Intrinsic polarisation fraction (all epochs)
# =============================================================================

pol_mjd    = np.array([e['mjd']              for e in all_epochs])
pol_dt     = Time(pol_mjd, format='mjd').datetime
pol_p0     = np.array([e['net_p0']           for e in all_epochs])
pol_err_lo = np.array([e['net_p0_err_lower'] for e in all_epochs])
pol_err_hi = np.array([e['net_p0_err_upper'] for e in all_epochs])

in_range = (pol_dt >= start_datetime) & (pol_dt <= end_datetime)

ax_lc[2].errorbar(
    pol_dt[in_range], pol_p0[in_range] * 100.0,
    yerr=[pol_err_lo[in_range] * 100.0, pol_err_hi[in_range] * 100.0],
    fmt='o', mfc=MK_COLOR, mec=MK_EDGE, ecolor='k',
    ms=6, capsize=ALPHA_CAPSIZE, capthick=ALPHA_CAPTHICK,
    linestyle='none', zorder=10, label=r'Net $p_0$')

for ei, epoch in enumerate(all_epochs):
    epoch_dt = Time(epoch['mjd'], format='mjd').datetime
    if not (start_datetime <= epoch_dt <= end_datetime):
        continue
    for ci, comp in enumerate(epoch['components']):
        ax_lc[2].errorbar(
            epoch_dt, comp['p0'] * 100.0,
            yerr=[[comp['p0_err_lower'] * 100.0], [comp['p0_err_upper'] * 100.0]],
            fmt='s', mfc=MK_COLOR, mec=MK_EDGE, ecolor='k',
            ms=4.5, capsize=1.5, capthick=0.6,
            linestyle='none', alpha=0.5, zorder=5,
            label=r'Component $p_0$' if (ei == 0 and ci == 0) else None)

ax_lc[2].set_yscale('log')
ax_lc[2].set_ylabel(r'$p_0$ (%)', fontsize=FONT_LC_YLABEL)


# =============================================================================
# LC -- shared x-axis, state transitions, formatting
# =============================================================================

ax_lc[0].set_xlim(start_datetime, end_datetime)
for a in ax_lc:
    if start_datetime <= transition_datetime[0] <= end_datetime:
        a.axvline(transition_datetime[0], ls='--', c='k', lw=0.8, zorder=0)
    if start_datetime <= transition_datetime[1] <= end_datetime:
        a.axvline(transition_datetime[1], ls=':', c='k', lw=0.8)
    if start_datetime <= highlight_sep20_datetime <= end_datetime:
        a.axvline(highlight_sep20_datetime, ls='--', c='k', lw=0.8, zorder=0)

for a in ax_lc:
    a.tick_params(labelbottom=False)

mjd_secondary_ax = FormatAxis(np.array(ax_lc), mjd=mjd_range,
                               interval=LC_TICK_INTERVAL, dt=0)
ax_lc[-1].tick_params(labelbottom=False)

align_axis_x(ax_lc[1], ax_lc[0])
align_axis_x(ax_lc[2], ax_lc[0])


# =============================================================================
# SAVE / SHOW
# =============================================================================

if SAVE_PATH is not None:
    save_path = _path_with_output_format(SAVE_PATH)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f'\nSaved -> {save_path}')

if SHOW_PLOT:
    plt.show()

plt.close(fig)

# --- Summary table -----------------------------------------------------------
print('\nEpoch summary:')
for e in all_epochs:
    comp_str = ', '.join(f"{c['name']}: p0={c['p0']*100:.3f}%" for c in e['components'])
    print(f"  {e['iso_date']}  MJD={e['mjd']:.1f}  {e['model_type']:>3s}  "
          f"net_p0={e['net_p0']*100:.3f}% "
          f"(+{e['net_p0_err_upper']*100:.3f}/-{e['net_p0_err_lower']*100:.3f})  "
          f"[{comp_str}]")
