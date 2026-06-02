#!/usr/bin/env python
"""
Publication summary figure for Swift J1727.8-1613.

Layout
------
3 × 2 grid of epoch panels (row-major, left→right, top→bottom):
  2023-09-04  2023-09-23  2023-10-06
  2023-10-14  2023-10-16  2023-10-22

Each panel contains three stacked sub-plots:
  Top    : fractional polarisation p(λ²) in %  — posterior sample lines (C0) +
                                                  68% CI fill (C0) + median model (black)
           λ² axis label and tick labels on the TOP of this panel;
           shown for ALL outer rows. All six p panels share the same λ² x-range.
  Middle : polarisation angle ψ(λ²) (deg) — same style.
           x-axis tick labels shown at BOTTOM of ψ panel (shared x-range).
           ψ panels in the same outer row share the same y-axis range.
  Bottom : FDF in % — posterior samples (C0, log y-axis) + CLEANed FDF (grey,
                 linear-mapped to log axis, foreground zorder > model) + model (black)

Sized for MNRAS double-column page  (174 mm × 234 mm ≈ 6.85" × 9.21").
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
import copy
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

_HERE     = Path(__file__).resolve().parent   # SpecPolJetComp/PlotsAndPostProc/code
_SPEC_POL = _HERE.parents[1]                  # SpecPolJetComp

sys.path.insert(0, str(_SPEC_POL / 'QUfitRMsynth/code'))
from faraday_model import load_results, update_model_from_params
from faraday_utils import (ThinComponent, ThinPowerLawComponent,
                           compute_polarization_quantities)

# =============================================================================
# CONFIGURATION — edit paths & controls here
# =============================================================================

PHI_RANGE           = 250           # FDF x-axis ± (rad m⁻²)
FDF_Y_RANGE         = (1e-3, 5e1)  # model FDF log y-axis limits — percentage units
N_PHI               = 4096          # Faraday-depth grid points
N_POSTERIOR_SAMPLES = 200           # posterior draws for CI bands on all panels
N_MC_OBS            = 500           # MC draws for observed p / ψ errors
LSQ_PAD_FRAC        = 0.10          # fractional padding on each side of global λ² range

DATA_DIR    = _SPEC_POL / 'QUfitRMsynth/results/J1727'
RMSYNTH_DIR = _SPEC_POL / 'QUfitRMsynth/rmsynth/J1727'

# Each entry: pkl path (QU-fit posterior) + rmsynth directory (for CLEANed FDF)
EPOCHS = [
    dict(pkl=DATA_DIR / 'SwiftJ1727_WAPITI_20230904/SwiftJ1727_WAPITI_20230904_S_sampler.pkl',
         rmsynth=RMSYNTH_DIR / 'SwiftJ1727_WAPITI_20230904'),
    dict(pkl=DATA_DIR / 'SwiftJ1727_WAPITI_20230923/SwiftJ1727_WAPITI_20230923_ST_sampler.pkl',
         rmsynth=RMSYNTH_DIR / 'SwiftJ1727_WAPITI_20230923'),
    dict(pkl=DATA_DIR / 'SwiftJ1727_WAPITI_20231006/SwiftJ1727_WAPITI_20231006_STT_sampler.pkl',
         rmsynth=RMSYNTH_DIR / 'SwiftJ1727_WAPITI_20231006'),
    dict(pkl=DATA_DIR / 'SwiftJ1727_WAPITI_20231014/SwiftJ1727_WAPITI_20231014_SSTT_sampler.pkl',
         rmsynth=RMSYNTH_DIR / 'SwiftJ1727_WAPITI_20231014'),
    dict(pkl=DATA_DIR / 'SwiftJ1727_WAPITI_20231016/SwiftJ1727_WAPITI_20231016_ST_sampler.pkl',
         rmsynth=RMSYNTH_DIR / 'SwiftJ1727_WAPITI_20231016'),
    dict(pkl=DATA_DIR / 'SwiftJ1727_QU_20231022/SwiftJ1727_QU_20231022_S_sampler.pkl',
         rmsynth=RMSYNTH_DIR / 'SwiftJ1727_QU_20231022'),
]

# Grid
NROWS, NCOLS  = 2, 3
HSPACE_OUTER  = 0.08
WSPACE_OUTER  = 0.10
HSPACE_INNER  = 0.12   # slightly larger to accommodate x tick labels at bottom of ψ panel
HEIGHT_RATIOS = [1, 1, 1.2]  # p : ψ : FDF

# Output
OUTPUT_FORMAT = 'pdf'   # 'pdf' or 'png'
SAVE_PATH     = '../plots/results/results_summary_panel'
SHOW_PLOT     = False

# =============================================================================
# MATPLOTLIB RCPARAMS — MNRAS publication style
# =============================================================================
FONT_BASE  = 9
FONT_LABEL = 10
FONT_TICK  = 9
FONT_DATE  = 9

plt.rcParams.update({
    'font.family':         'serif',
    'font.serif':          ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
    'mathtext.fontset':    'dejavuserif',
    'font.size':           FONT_BASE,
    'axes.labelsize':      FONT_LABEL,
    'axes.titlesize':      FONT_LABEL,
    'axes.linewidth':      0.8,
    'xtick.major.width':   0.8,
    'ytick.major.width':   0.8,
    'xtick.minor.width':   0.5,
    'ytick.minor.width':   0.5,
    'xtick.major.size':    3.0,
    'ytick.major.size':    3.0,
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
    'legend.fontsize':     FONT_TICK,
    'legend.framealpha':   1.0,
    'legend.facecolor':    'white',
    'legend.edgecolor':    'black',
    'legend.fancybox':     False,
    'lines.linewidth':     1.0,
    'lines.markersize':    3,
})

# =============================================================================
# PRE-PASS — global λ² range shared by all p / ψ panels
# =============================================================================

output_fmt = OUTPUT_FORMAT.strip().lower()
if output_fmt not in {'pdf', 'png'}:
    raise ValueError("OUTPUT_FORMAT must be 'pdf' or 'png'.")

print('Pre-computing global λ² range ...')
lsq_global_min, lsq_global_max = np.inf, -np.inf
for ec in EPOCHS:
    r = load_results(ec['pkl'])
    lsq_global_min = min(lsq_global_min, r.setup.data.lambda_sq.min())
    lsq_global_max = max(lsq_global_max, r.setup.data.lambda_sq.max())

lsq_pad      = LSQ_PAD_FRAC * (lsq_global_max - lsq_global_min)
lsq_xlim     = (lsq_global_min - lsq_pad, lsq_global_max + lsq_pad)
lsq_fine     = np.linspace(*lsq_xlim, 400)   # shared fine grid for all model curves
print(f'  λ² range: [{lsq_xlim[0]:.4f}, {lsq_xlim[1]:.4f}] m²')

phi_grid = np.linspace(-PHI_RANGE, PHI_RANGE, N_PHI)

# Precompute log-scale mapping constants for linear CLEANed FDF overlay
_log_ymin = np.log10(FDF_Y_RANGE[0])
_log_ymax = np.log10(FDF_Y_RANGE[1])

# =============================================================================
# FIGURE CONSTRUCTION
# =============================================================================

fig = plt.figure(figsize=(6.85, 9.21))
outer_gs = gridspec.GridSpec(
    NROWS, NCOLS, figure=fig,
    hspace=HSPACE_OUTER, wspace=WSPACE_OUTER,
    left=0.11, right=0.97, top=0.95, bottom=0.05,
)

# Collect ψ axes/ylim extents per outer row to share y-axis after the loop
row_psi_axes  = {r: [] for r in range(NROWS)}
row_psi_ylims = {r: [] for r in range(NROWS)}

for idx, epoch_cfg in enumerate(EPOCHS):
    row, col = divmod(idx, NCOLS)
    is_bottom_row = (row == NROWS - 1)

    inner_gs = outer_gs[row, col].subgridspec(
        3, 1, hspace=HSPACE_INNER, height_ratios=HEIGHT_RATIOS)
    ax_p   = fig.add_subplot(inner_gs[0])
    ax_psi = fig.add_subplot(inner_gs[1], sharex=ax_p)
    ax_fdf = fig.add_subplot(inner_gs[2])

    pkl_path  = epoch_cfg['pkl']
    rmsynth_d = epoch_cfg['rmsynth']
    stem      = rmsynth_d.name

    dm = re.search(r'(\d{8})', pkl_path.stem)
    d  = dm.group(1)
    iso_date = f'{d[:4]}-{d[4:6]}-{d[6:]}'

    # ----------------------------------------------------------------
    # Load posterior & draw samples
    # ----------------------------------------------------------------
    results = load_results(pkl_path)
    data    = results.setup.data
    lsq     = data.lambda_sq
    lsq_ref = getattr(results.setup, 'lambda_sq_ref', None)

    n_tot    = len(results.samples)
    w_norm   = results.weights / results.weights.sum()
    samp_idx = np.random.choice(n_tot, size=min(N_POSTERIOR_SAMPLES, n_tot),
                                replace=False, p=w_norm)

    # Representative model (needed for peak caps before the sample loop)
    rep_params = results.get_representative_sample(percentile=1.0, method='median')
    best_model = copy.deepcopy(results.setup.model)
    update_model_from_params(best_model, rep_params,
                             results.setup.param_names, lambda_sq_ref=lsq_ref)
    # Peak caps in fractional units (clipping before × 100 scale)
    best_peaks = [np.max(np.abs(comp.compute_fdf(phi_grid)))
                  for comp in best_model.components]

    n_comp = len(best_model.components)
    print(f'  [{idx+1}/6] {iso_date}  n_components={n_comp}')

    p_fine_s   = []
    psi_fine_s = []
    for si in samp_idx:
        ms = copy.deepcopy(results.setup.model)
        update_model_from_params(ms, results.samples[si],
                                 results.setup.param_names, lambda_sq_ref=lsq_ref)

        # p(λ²) and ψ(λ²): thin sample lines + accumulate for percentile band
        P     = ms.compute_polarization(lsq_fine)
        p_s   = np.abs(P)
        psi_s = np.degrees(np.angle(P) / 2.0)
        p_fine_s.append(p_s)
        psi_fine_s.append(psi_s)
        ax_p.plot(lsq_fine,  p_s * 100,   color='C0', lw=0.3, alpha=0.08, zorder=2)
        ax_psi.plot(lsq_fine, psi_s, color='C0', lw=0.3, alpha=0.08, zorder=2)

        # FDF posterior samples in % (clip before scaling)
        for comp, peak_cap in zip(ms.components, best_peaks):
            amp = np.abs(comp.compute_fdf(phi_grid))
            if isinstance(comp, ThinComponent):
                amp = np.clip(amp, None, peak_cap)
            ax_fdf.plot(phi_grid, amp * 100, color='C0', lw=0.5, alpha=0.12, zorder=2)

    p_arr   = np.array(p_fine_s)
    psi_arr = np.array(psi_fine_s)
    p_16,   p_50,   p_84   = np.percentile(p_arr,   [16, 50, 84], axis=0)
    psi_16, psi_50, psi_84 = np.percentile(psi_arr, [16, 50, 84], axis=0)

    # Observed p and ψ
    p_obs, p_err_mc, chi_obs, chi_err_mc = compute_polarization_quantities(
        data.Q, data.U, data.Q_err, data.U_err,
        method='mc', n_samples=N_MC_OBS, debias=True)
    psi_obs_deg = np.degrees(chi_obs)
    p_err_lo,   p_err_hi   = p_err_mc
    psi_err_lo, psi_err_hi = np.degrees(chi_err_mc[0]), np.degrees(chi_err_mc[1])

    # ----------------------------------------------------------------
    # p(λ²) panel — values in %, log y-scale
    # ----------------------------------------------------------------
    ax_p.set_yscale('log')
    ax_p.fill_between(lsq_fine, p_16 * 100, p_84 * 100,
                      color='C0', alpha=0.35, lw=0, zorder=3)
    ax_p.errorbar(lsq, p_obs * 100,
                  yerr=[p_err_lo * 100, p_err_hi * 100],
                  fmt='o', mfc='lightgrey', mec='k', ecolor='k',
                  ms=3.0, mew=0.02, lw=0.4, capsize=1.0, capthick=0.4, zorder=5)
    ax_p.plot(lsq_fine, p_50 * 100, color='k', lw=1.2, zorder=10)

    ax_p.set_ylim(1e-1, 10)
    ax_p.set_xlim(*lsq_xlim)

    ax_p.set_ylabel(r'$p\,(\lambda^2)$ (%) $(\mathrm{m^{2}})^{-1}$')

    # λ² axis on top — top outer row only
    if row == 0:
        ax_p.tick_params(axis='x', which='both',
                         top=True, labeltop=True,
                         bottom=True, labelbottom=False)
        ax_p.xaxis.set_label_position('top')
        ax_p.set_xlabel(r'$\lambda^2\;(\mathrm{m}^2)$')
    else:
        ax_p.tick_params(axis='x', which='both', labeltop=False, labelbottom=False)
        ax_p.set_xlabel('')

    # ----------------------------------------------------------------
    # ψ(λ²) panel — values in degrees, x tick labels at bottom
    # ----------------------------------------------------------------
    ax_psi.fill_between(lsq_fine, psi_16, psi_84, color='C0', alpha=0.35, lw=0, zorder=3)
    ax_psi.errorbar(lsq, psi_obs_deg, yerr=[psi_err_lo, psi_err_hi],
                    fmt='o', mfc='lightgrey', mec='k', ecolor='k',
                    ms=3.0, mew=0.02, lw=0.4, capsize=1.0, capthick=0.4, zorder=5)
    ax_psi.plot(lsq_fine, psi_50, color='k', lw=1.2, zorder=10)

    psi_err_med = np.median(np.concatenate([psi_err_lo, psi_err_hi]))
    psi_ymin = np.percentile(psi_obs_deg, 2)  - 3 * psi_err_med
    psi_ymax = np.percentile(psi_obs_deg, 98) + 3 * psi_err_med
    row_psi_axes[row].append(ax_psi)
    row_psi_ylims[row].append((psi_ymin, psi_ymax))

    ax_psi.set_ylabel(r'$\psi\,(\lambda^2)$ (deg) $(\mathrm{m^{2}})^{-1}$')
    ax_psi.tick_params(axis='x', which='both', labelbottom=False, labeltop=False)
    ax_psi.set_xlabel('')

    # ----------------------------------------------------------------
    # FDF panel — amplitudes in %
    # ----------------------------------------------------------------
    ax_fdf.set_yscale('log')
    ax_fdf.set_xlim(-PHI_RANGE, PHI_RANGE)
    ax_fdf.set_ylim(*FDF_Y_RANGE)
    ax_fdf.axvline(0, color='k', ls=':', lw=0.6, zorder=0)
    ax_fdf.set_ylabel(r'$|F(\phi_f)|$ (%) $(\mathrm{rad\,m^{-2}})^{-1}$')

    if is_bottom_row:
        ax_fdf.set_xlabel(r'$\phi_f\;(\mathrm{rad\,m^{-2}})$')
    else:
        ax_fdf.set_xlabel('')
        ax_fdf.tick_params(axis='x', labelbottom=False)

    # CLEANed FDF — mapped linearly into the log y-axis space (% units via FDF_Y_RANGE).
    # zorder=10000: foreground above both posterior samples (2) and model (15/16).
    fdf_file = rmsynth_d / f'{stem}_rmsynth_FDFclean.dat'
    if fdf_file.exists():
        fdf_dat  = np.loadtxt(fdf_file)
        phi_c    = fdf_dat[:, 0]
        amp_c    = np.sqrt(fdf_dat[:, 1]**2 + fdf_dat[:, 2]**2)
        mask_c   = np.abs(phi_c) <= PHI_RANGE
        amp_peak = amp_c[mask_c].max()
        if amp_peak > 0:
            amp_norm   = amp_c[mask_c] / amp_peak          # 0 → 1 (linear)
            amp_mapped = 10 ** (_log_ymin + amp_norm * (_log_ymax - _log_ymin) * 0.90)
            ax_fdf.plot(phi_c[mask_c], amp_mapped,
                        color='grey', alpha=1.0, lw=1.0, zorder=10000)
    else:
        print(f'  WARNING: missing FDFclean for {stem}')

    # Representative model — zorder 15/16, below CLEANed FDF (20)
    for comp in best_model.components:
        comp_amp = np.abs(comp.compute_fdf(phi_grid)) * 100   # fractional → %
        ax_fdf.plot(phi_grid, comp_amp, 'k-', lw=1.2, zorder=15)
        if isinstance(comp, (ThinComponent, ThinPowerLawComponent)):
            ax_fdf.plot(comp.phi_rm, np.max(comp_amp), '^',
                        color='k', ms=4, zorder=16)

    # Date label top-right of FDF panel
    ax_fdf.text(0.97, 0.96, iso_date, transform=ax_fdf.transAxes,
                fontsize=FONT_DATE, ha='right', va='top',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=1.0))

    # ----------------------------------------------------------------
    # Suppress y-labels/tick-labels on non-left-column panels
    # ----------------------------------------------------------------
    if col != 0:
        ax_p.set_ylabel('', labelpad=0)
        ax_psi.set_ylabel('', labelpad=0)
        ax_fdf.set_ylabel('', labelpad=0)
        ax_p.tick_params(labelleft=False)
        ax_psi.tick_params(labelleft=False)
        ax_fdf.tick_params(labelleft=False)

    print(f'  [{idx+1}/6] {iso_date}  done')

# =============================================================================
# SHARE ψ Y-AXIS WITHIN EACH OUTER ROW
# =============================================================================
for r in range(NROWS):
    if row_psi_axes[r]:
        shared_ymin = min(y[0] for y in row_psi_ylims[r])
        shared_ymax = max(y[1] for y in row_psi_ylims[r])
        for ax in row_psi_axes[r]:
            ax.set_ylim(shared_ymin, shared_ymax)

# =============================================================================
# SAVE / SHOW
# =============================================================================
save_path = Path(f'{SAVE_PATH}.{output_fmt}')
save_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(save_path, dpi=300, bbox_inches='tight')
print(f'\nSaved → {save_path}')

if SHOW_PLOT:
    plt.show()
plt.close(fig)
