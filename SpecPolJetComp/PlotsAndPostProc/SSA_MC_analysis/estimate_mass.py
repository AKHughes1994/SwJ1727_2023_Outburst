#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimate magnetic field geometry, RM-inferred particle density, and rotating
mass from the joint SSA + RM posterior samples.

Quantities derived
------------------
B_perp      : SSA equipartition B-field perpendicular to the LOS (mG)
B_parallel  : LOS component of B, inferred via depolarisation + inclination (mG)
ne_rel      : relativistic electron density from SSA (cm^-3)
ne_rot      : Faraday-rotating particle density inferred from RM (cm^-3)
ne_frac     : relativistic electron fraction  ne_rel / (ne_rel + ne_rot)
M_rot       : mass of Faraday-rotating plasma (g, M_sun)
f_rot       : M_rot / M_accrete  (dimensionless mass fraction)

Plotting convention
-------------------
Skewed / wide-dynamic-range posteriors (B_parallel, ne_frac, f_rot) are
plotted in log-space on the x-axis using _posterior_plot_logx, which computes
the histogram, KDE, mode, and HDI entirely in log10-space before converting
back to linear for display.  All other posteriors use linear or log10 samples
directly with _posterior_plot.

@author: andrew KENNETH hughes
"""

import numpy as np
import matplotlib.pyplot as plt
import arviz as az
from pathlib import Path
from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar

# ── Matplotlib style (MNRAS-compatible) ──────────────────────────────────────
plt.rcParams.update({
    'font.family':          'serif',
    'font.serif':           ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
    'mathtext.fontset':     'dejavuserif',
    'font.size':            8,
    'axes.labelsize':       8,
    'axes.titlesize':       8,
    'axes.linewidth':       0.8,
    'xtick.major.width':    0.8,
    'ytick.major.width':    0.8,
    'xtick.minor.width':    0.5,
    'ytick.minor.width':    0.5,
    'xtick.major.size':     3,
    'ytick.major.size':     3,
    'xtick.minor.size':     1.5,
    'ytick.minor.size':     1.5,
    'xtick.direction':      'in',
    'ytick.direction':      'in',
    'xtick.top':            True,
    'xtick.bottom':         True,
    'ytick.left':           True,
    'ytick.right':          True,
    'xtick.minor.visible':  True,
    'ytick.minor.visible':  True,
    'xtick.labelsize':      7,
    'ytick.labelsize':      7,
    'legend.fontsize':      7,
    'legend.framealpha':    1.0,
    'legend.facecolor':     'white',
    'legend.edgecolor':     'black',
    'legend.fancybox':      False,
    'lines.linewidth':      0.8,
})

# ── Paths and constants ───────────────────────────────────────────────────────
SAMPLES_DIR = Path("samples")
PLOTS_DIR   = Path("plots_inferred")
PLOTS_DIR.mkdir(exist_ok=True)

M_SUN  = 1.989e33  # g
F_V    = 1.0       # volume filling / geometric factor (fixed to unity)
W_RM99 = 100.0     # RM width at 99% level (rad m^-2); sets ne_rot scale
BINS   = 180       # histogram bin count used by all posterior plots

# ── Load SSA posterior samples ────────────────────────────────────────────────
# All SSA output arrays are stored as log10 of the physical quantity.
print(f"\n{'=' * 62}")
print(f"  estimate_mass  —  loading SSA posteriors")
print(f"  samples_dir : {SAMPLES_DIR.resolve()}")
print(f"  plots_dir   : {PLOTS_DIR.resolve()}")
print(f"  W_rm99      : {W_RM99} rad m^-2  |  f_V = {F_V}")
print(f"{'=' * 62}")

_B_log10  = np.load(SAMPLES_DIR / "mc_out_B_energy_form_samples.npy")
_R_log10  = np.load(SAMPLES_DIR / "mc_out_R_energy_form_samples.npy")
_ne_log10 = np.load(SAMPLES_DIR / "mc_out_ne_energy_form_samples.npy")
D_kpc     = np.load(SAMPLES_DIR / "mc_out_dist_kpc_samples.npy")

# Convert SSA samples to linear CGS
B_gauss = 10.0 ** _B_log10   # G  — perpendicular component from SSA equipartition
R_cm    = 10.0 ** _R_log10   # cm — source radius
ne_rel  = 10.0 ** _ne_log10  # cm^-3 — relativistic electron number density

B_perp_mG = B_gauss * 1e3    # G → mG

n = len(B_gauss)
print(f"\n  Loaded {n:,} samples for B, R, n_e, D")

# ── Draw nuisance parameters ──────────────────────────────────────────────────
# Inclination: uniform in cos(i) over [10°, 80°] — avoids face-on / edge-on extremes
cos_i = np.random.uniform(np.cos(np.deg2rad(80)), np.cos(np.deg2rad(10)), size=n)
i_rad = np.arccos(cos_i)

# Depolarisation fraction: uniform over plausible compact-jet range
depol = np.random.uniform(0.01, 0.1, size=n)

# ── Derived quantities ────────────────────────────────────────────────────────

# LOS magnetic field component via depolarisation and inclination geometry
B_parallel_mG    = np.sqrt(depol) * B_gauss / np.tan(i_rad) * 1e3   # G → mG
log10_B_parallel = np.log10(B_parallel_mG)

# Path length through the source: SSA diameter
l_cm = R_cm * 2

# Jet path length: l ~ prefactor × 10^14 × (D/5.5 kpc) × (φ_j/5°) cm
# Prefactor and opening angle φ_j drawn from physically motivated priors
_prefactor = np.random.uniform(0.7,  1.4,  size=n)
_phi_j_deg = np.random.uniform(1.0, 10.0,  size=n)
l_jet_cm   = _prefactor * 1e14 * (D_kpc / 5.5) * (_phi_j_deg / 5.0)

# RM-inferred Faraday-rotating particle density (eq. 1)
ne_rot = 3.8e3 * (W_RM99 / 100.0) * B_parallel_mG**(-1) * (l_cm / 1e14)**(-1)

# Rotating plasma mass (eq. 3), f_V = 1
M_rot_g    = 6.4e21 * F_V * (W_RM99 / 100.0) * B_parallel_mG**(-1) * (l_cm / 1e14)**2
M_rot_Msun = M_rot_g / M_SUN

# Relativistic electron fraction (bounded 0–1); computed in log-space
ne_ratio      = ne_rel / (ne_rel + ne_rot)
log10_ne_frac = np.log10(ne_ratio)

# Accreted mass budget
C_CGS    = 2.99792458e10  # cm s^-1
ETA      = 0.05           # accretion radiative efficiency
T_FLARE  = 12.0 * 3600.0  # s  (12 hr flare duration)
F_XRAY   = 2e-7           # erg cm^-2 s^-1  (2–20 keV flux, WebPIMMS, Γ=2.0, N_H=2.68e21)
BOL_CORR = 2.0            # bolometric correction factor

D_cm        = np.median(D_kpc) * 3.0857e21
L_bol       = 4.0 * np.pi * D_cm**2 * F_XRAY * BOL_CORR
M_accrete_g = L_bol * T_FLARE / (ETA * C_CGS**2)

# Mass fraction; both derived in log-space to preserve dynamic range
log10_M_rot_g = np.log10(M_rot_g)
log10_f_rot   = log10_M_rot_g - np.log10(M_accrete_g)
f_rot         = M_rot_g / M_accrete_g

# ── Statistics helpers ────────────────────────────────────────────────────────
def _kde_mode(x):
    """KDE-based mode of samples x, searched between the 1st and 99th percentiles."""
    kde = gaussian_kde(x)
    lo, hi = np.percentile(x, [1, 99])
    result = minimize_scalar(lambda v: -kde(v).item(), bounds=(lo, hi), method="bounded")
    return float(result.x)

def _linear_hdi(mode_log, hdi):
    """Convert log10-space mode and HDI bounds to linear asymmetric uncertainties."""
    val = 10.0 ** mode_log
    neg = val - 10.0 ** float(hdi[0])
    pos = 10.0 ** float(hdi[1]) - val
    return val, neg, pos

def _print_summary(label, samples, unit, is_log10=False):
    """Print mode, mean, median, std, and HDI for a posterior sample array.

    If is_log10=True the samples are assumed to be log10 of the physical
    quantity and linear-space equivalents are printed as well.
    """
    mode   = _kde_mode(samples)
    mean   = float(np.mean(samples))
    median = float(np.median(samples))
    std    = float(np.std(samples, ddof=1))
    hdi68  = az.hdi(samples, hdi_prob=0.68)
    hdi95  = az.hdi(samples, hdi_prob=0.95)
    hdi997 = az.hdi(samples, hdi_prob=0.997)
    print(f"\n  Results — {label}")
    print(f"  {'─' * 56}")
    print(f"  Mode (KDE)    : {mode:.4g} {unit}")
    print(f"  Mean          : {mean:.4g} {unit}")
    print(f"  Median        : {median:.4g} {unit}")
    print(f"  Std dev       : {std:.4g} {unit}")
    print(f"  68% HDI       : [{hdi68[0]:.4g},  {hdi68[1]:.4g}] {unit}")
    print(f"  95% HDI       : [{hdi95[0]:.4g},  {hdi95[1]:.4g}] {unit}")
    print(f"  99.7% HDI     : [{hdi997[0]:.4g},  {hdi997[1]:.4g}] {unit}")
    if is_log10:
        val, n68, p68   = _linear_hdi(mode, hdi68)
        _,   n95, p95   = _linear_hdi(mode, hdi95)
        _,   n997, p997 = _linear_hdi(mode, hdi997)
        print(f"  ·· linear ··")
        print(f"  Mode (linear) : {val:.4g}")
        print(f"  68% HDI (lin) : {val:.4g}  (-{n68:.4g} / +{p68:.4g})")
        print(f"  95% HDI (lin) : {val:.4g}  (-{n95:.4g} / +{p95:.4g})")
        print(f"  99.7%   (lin) : {val:.4g}  (-{n997:.4g} / +{p997:.4g})")
    print(f"{'─' * 62}")
    return dict(mode=mode, mean=mean, median=median, std=std,
                hdi68=hdi68, hdi95=hdi95, hdi997=hdi997)

# ── Plotting helpers ──────────────────────────────────────────────────────────
def _posterior_plot(samples, xlabel, fname, stats, xscale="linear", legend_loc="upper right"):
    """Standard posterior plot: histogram + KDE with 68% HDI and mode.

    Mode and HDI are always computed from the samples in their natural
    (linear) space.  xscale='log' switches to log-spaced bins and a log
    x-axis for display only — statistics are unchanged.
    """
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    if xscale == "log":
        bin_edges = np.logspace(np.log10(samples[samples > 0].min()),
                                np.log10(samples.max()), BINS + 1)
        xkde = np.logspace(np.log10(samples[samples > 0].min()),
                           np.log10(samples.max()), 600)
    else:
        bin_edges = BINS
        xkde = np.linspace(samples.min(), samples.max(), 600)

    ax.hist(samples, bins=bin_edges, density=True, color="C0", alpha=0.5, zorder=2)
    ax.plot(xkde, gaussian_kde(samples)(xkde), color="C0", lw=1.0, zorder=3)

    ax.axvspan(*stats["hdi68"], alpha=0.5, color="grey", label="68% HDI", zorder=1)
    ax.axvline(stats["mode"], color="black", lw=1.0, ls="--", zorder=5, label="Mode")

    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend(loc=legend_loc)

    fig.tight_layout()
    path = PLOTS_DIR / fname
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved posterior plot:        {path}")
    plt.close(fig)

def _posterior_plot_logx(log10_samples, xlabel, fname, legend_loc="upper right"):
    """Log-space posterior plot for skewed / wide-dynamic-range distributions.

    Everything — histogram bins, KDE, mode, HDI — is computed in log10-space
    then converted to linear for display on a log x-axis.  Avoids the dynamic-
    range compression that arises when a right-skewed distribution is treated
    as if it were symmetric in linear space.  Density is normalised to peak = 1.
    """
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    log_s = log10_samples

    mode_lin  = 10.0 ** _kde_mode(log_s)
    hdi68_lin = 10.0 ** az.hdi(log_s, hdi_prob=0.68)

    counts, edges = np.histogram(log_s, bins=BINS)
    log_grid = np.linspace(log_s.min(), log_s.max(), 600)
    y_kde    = gaussian_kde(log_s)(log_grid)

    ax.stairs(counts / counts.max(), 10.0 ** edges,
              fill=True, color="C0", alpha=0.5, zorder=2)
    ax.plot(10.0 ** log_grid, y_kde / y_kde.max(), color="C0", lw=1.0, zorder=3)
    ax.axvspan(*hdi68_lin, alpha=0.5, color="grey", label="68% HDI", zorder=1)
    ax.axvline(mode_lin, color="black", lw=1.0, ls="--", zorder=5, label="Mode")

    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalised density")
    ax.legend(loc=legend_loc)

    fig.tight_layout()
    path = PLOTS_DIR / fname
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved posterior plot:        {path}")
    plt.close(fig)

# ── Print summaries ───────────────────────────────────────────────────────────
# B_parallel and B_perp: B_parallel is log-space because depol × tan(i) sampling
# produces a highly skewed distribution; B_perp is roughly symmetric.
stats_l_jet   = _print_summary("l_jet (jet diam)",          np.log10(l_jet_cm),    "log10(cm)",    is_log10=True)
stats_Bperp   = _print_summary("B_perp (SSA)",              B_perp_mG,             "mG")
stats_Bpar    = _print_summary("B_parallel",                log10_B_parallel,      "log10(mG)",    is_log10=True)
stats_ne_rel  = _print_summary("ne_rel (SSA)",              _ne_log10,             "log10(cm^-3)", is_log10=True)
stats_ne_rot  = _print_summary("ne_rot (RM)",               np.log10(ne_rot),      "log10(cm^-3)", is_log10=True)
stats_ne_frac = _print_summary("ne_rel / (ne_rel+ne_rot)",  ne_ratio,              "(linear fraction)")
ll99_ne_frac  = float(np.percentile(ne_ratio, 1))
print(f"  99% lower limit (ne_ratio) : {ll99_ne_frac:.4f}  ({ll99_ne_frac*100:.2f}%)")
print(f"{'─' * 62}")
stats_M_g     = _print_summary("M_rot",                     log10_M_rot_g,         "log10(g)",     is_log10=True)
stats_M_sun   = _print_summary("M_rot",                     np.log10(M_rot_Msun),  "log10(M_sun)", is_log10=True)

# f_rot stats computed in log10-space to recover the geometric mode
hdi68_f  = 10.0 ** az.hdi(log10_f_rot, hdi_prob=0.68)
hdi95_f  = 10.0 ** az.hdi(log10_f_rot, hdi_prob=0.95)
ul99_f   = float(10.0 ** np.percentile(log10_f_rot, 99))
mode_f   = float(10.0 ** _kde_mode(log10_f_rot))
stats_f_rot = dict(
    mode=mode_f, mean=float(np.mean(f_rot)), median=float(np.median(f_rot)),
    std=float(np.std(f_rot, ddof=1)),
    hdi68=hdi68_f, hdi95=hdi95_f, hdi997=10.0 ** az.hdi(log10_f_rot, hdi_prob=0.997),
)

_M_acc_55 = 4.0 * np.pi * (5.5 * 3.0857e21)**2 * F_XRAY * BOL_CORR * T_FLARE / (ETA * C_CGS**2)
print(f"\n  Results — M_rot / M_accrete  (η={ETA}, t_flare={T_FLARE/3600:.0f} hr)")
print(f"  {'─' * 56}")
print(f"  M_accrete   : {M_accrete_g:.3g} g  ({M_accrete_g/M_SUN:.3g} M_sun)")
print(f"  M_acc scale : {_M_acc_55:.3g} × (D/5.5 kpc)² × ({ETA}/η) × (t_flare/{T_FLARE/3600:.0f} hr)  g")
print(f"  Mode        : {mode_f:.3g}")
print(f"  68% HDI     : [{hdi68_f[0]:.3g},  {hdi68_f[1]:.3g}]")
print(f"  95% HDI     : [{hdi95_f[0]:.3g},  {hdi95_f[1]:.3g}]")
print(f"  99% UL      : < {ul99_f:.3g}")
print(f"  f scaling   : f ~ {mode_f:.3g}  ×  (η/{ETA})  ×  ({T_FLARE/3600:.0f} hr/t_flare)")
print(f"{'─' * 62}\n")

# ── Save individual posterior plots ──────────────────────────────────────────
_posterior_plot(
    B_perp_mG,
    r"$B_\perp\ (\mathrm{mG})$",
    "B_perp_posterior.pdf",
    stats_Bperp,
)

# B_parallel is skewed (depol × inclination geometry) — use log-space plot
_posterior_plot_logx(
    log10_B_parallel,
    r"$B_\parallel\ (\mathrm{mG})$",
    "B_parallel_posterior.pdf",
)

_posterior_plot(
    _ne_log10,
    r"$\log_{10}(n_{e,\rm rel}\ /\ \mathrm{cm}^{-3})$",
    "ne_rel_posterior.pdf",
    stats_ne_rel,
)

_posterior_plot(
    np.log10(ne_rot),
    r"$\log_{10}(n_{e,\rm rot}\ /\ \mathrm{cm}^{-3})$",
    "ne_rot_posterior.pdf",
    stats_ne_rot,
)

_posterior_plot(
    ne_ratio,
    r"$n_{e,\rm rel}\ /\ (n_{e,\rm rel} + n_{e,\rm rot})$",
    "ne_ratio_posterior.pdf",
    stats_ne_frac,
    xscale="log",
    legend_loc="upper left",
)

_posterior_plot(
    log10_M_rot_g,
    r"$\log_{10}(M_{\rm rot}\ /\ \mathrm{g})$",
    "M_rot_g_posterior.pdf",
    stats_M_g,
)

_posterior_plot(
    np.log10(M_rot_Msun),
    r"$\log_{10}(M_{\rm rot}\ /\ M_\odot)$",
    "M_rot_Msun_posterior.pdf",
    stats_M_sun,
)

# f_rot is right-skewed over several decades — use log-space plot
_posterior_plot_logx(
    log10_f_rot,
    r"$f_{\rm rot}$",
    "f_rot_posterior.pdf",
)

_posterior_plot(
    np.log10(l_jet_cm),
    r"$\log_{10}(l_{\rm jet}\ /\ \mathrm{cm})$",
    "l_jet_posterior.pdf",
    stats_l_jet,
)

# ── Two-panel: M_rot (g) and f_rot on log x-axes ─────────────────────────────
def _two_panel_log_x_plot(log10_top, xlabel_top, log10_bot, xlabel_bot, fname):
    """Two-panel figure with log x-axes.  Both panels computed entirely in
    log10-space (uniform bins, KDE, mode, HDI) then displayed in linear space
    on a log x-axis.  Legend shown on top panel only."""
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(3.5, 5.0))
    for i, (ax, log_s, xlabel) in enumerate((
        (ax_top, log10_top, xlabel_top),
        (ax_bot, log10_bot, xlabel_bot),
    )):
        mode_lin  = 10.0 ** _kde_mode(log_s)
        hdi68_lin = 10.0 ** az.hdi(log_s, hdi_prob=0.68)

        counts, edges = np.histogram(log_s, bins=BINS)
        log_grid = np.linspace(log_s.min(), log_s.max(), 600)
        y_kde    = gaussian_kde(log_s)(log_grid)

        ax.stairs(counts / counts.max(), 10.0 ** edges,
                  fill=True, color="C0", alpha=0.5, zorder=2)
        ax.plot(10.0 ** log_grid, y_kde / y_kde.max(), color="C0", lw=1.0, zorder=3)
        ax.axvspan(*hdi68_lin, alpha=0.5, color="grey", label="68% HDI", zorder=1)
        ax.axvline(mode_lin, color="black", lw=1.0, ls="--", zorder=5, label="Mode")
        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalised density")
        if i == 0:
            ax.legend(loc="upper right")

    fig.tight_layout()
    path = PLOTS_DIR / fname
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved two-panel plot:        {path}")
    plt.close(fig)

_two_panel_log_x_plot(
    log10_M_rot_g, r"$M_{\rm rot}\ (\mathrm{g})$",
    log10_f_rot,   r"$f_{\rm rot}$",
    "M_rot_f_rot_two_panel.pdf",
)

# ── Diagnostic: ne_rel and ne_rot vs distance ─────────────────────────────────
_panels = [
    (np.log10(ne_rel), r"$\log_{10}(n_{e,\rm rel}\ /\ \mathrm{cm}^{-3})$",
     r"$n_{e,\rm rel}$ vs distance (SSA)",    "Blues"),
    (np.log10(ne_rot), r"$\log_{10}(n_{e,\rm rot}\ /\ \mathrm{cm}^{-3})$",
     r"$n_{e,\rm rot}$ vs distance (RM)",     "Oranges"),
    (ne_ratio,         r"$n_{e,\rm rel}\ /\ (n_{e,\rm rel} + n_{e,\rm rot})$",
     "Relativistic fraction vs distance",     "Greens"),
]

fig, axes = plt.subplots(2, 3, figsize=(11, 7),
                          gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05})

for col, (ydata, ylabel, title, cmap) in enumerate(_panels):
    ax_top = axes[0, col]
    ax_bot = axes[1, col]

    hb = ax_top.hexbin(D_kpc, ydata, gridsize=60, cmap=cmap, mincnt=1, linewidths=0.1)
    ax_top.set_ylabel(ylabel)
    ax_top.set_title(title)
    ax_top.set_xticklabels([])
    fig.colorbar(hb, ax=ax_top, label="count")

    # Bottom panel: posterior histogram matching _posterior_plot style
    colour = hb.cmap(0.6)
    mode_val  = _kde_mode(ydata)
    hdi68_val = az.hdi(ydata, hdi_prob=0.68)

    ax_bot.hist(ydata, bins=BINS, density=True, color=colour, alpha=0.5, zorder=2)
    xkde = np.linspace(ydata.min(), ydata.max(), 600)
    ax_bot.plot(xkde, gaussian_kde(ydata)(xkde), color=colour, lw=1.0, zorder=3)
    ax_bot.axvspan(*hdi68_val, alpha=0.4, color="grey", label="68% HDI", zorder=1)
    ax_bot.axvline(mode_val, color="black", lw=1.0, ls="--", zorder=5, label="Mode")
    ax_bot.set_xlim(ax_top.get_ylim())
    ax_bot.set_xlabel(ylabel)
    ax_bot.set_ylabel("Density")
    ax_bot.legend(fontsize=6, loc="upper right")

fig.tight_layout()
path = PLOTS_DIR / "ne_vs_distance_diagnostic.pdf"
fig.savefig(path, dpi=300, bbox_inches="tight")
print(f"Saved diagnostic plot:       {path}")
plt.close(fig)
