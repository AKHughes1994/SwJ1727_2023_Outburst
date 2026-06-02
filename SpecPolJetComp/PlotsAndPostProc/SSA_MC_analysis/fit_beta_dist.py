"""
fit_beta.py

Fits a shifted/scaled Beta distribution to a user-specified mode, median,
and equal-tailed interval (ETI), using deterministic analytical objectives
and scipy differential evolution.

Outputs the best-fit parameters and saves a diagnostic plot.
"""

import os
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import differential_evolution

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

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================

# Output directory for plots. Override this to save elsewhere.
PLOT_DIR = os.path.join(os.getcwd(), "plots")

# Plot filename (no extension)
PLOT_NAME = "fit_beta_result"

# X-axis range for the plot. Set to None to use automatic ranging.
PLOT_XRANGE = (2, 12)

# Whether to call plt.show() and display the plot interactively.
SHOW_PLOT = False

# =============================================================================
# TARGET CONSTRAINTS  —  edit these
# =============================================================================

MODE_TARGET    = 5.1   # Target mode
MEDIAN_TARGET  = 5.5   # Target median
Q_LOWER_TARGET = MEDIAN_TARGET - 1.1   # Value at lower quantile
Q_UPPER_TARGET = MEDIAN_TARGET + 1.4   # Value at upper quantile
P_LOWER        = 0.16  # Lower quantile probability (e.g. 0.16 = 16th percentile)
P_UPPER        = 0.84  # Upper quantile probability (e.g. 0.84 = 84th percentile)

# Optional far-tail constraints. Set either target to None to disable that side.
Q_FAR_LEFT_TARGET  = 3.2   # e.g. 3.0 — value at the far-left quantile
P_FAR_LEFT         = 0.01   # Far-left quantile probability  (default 1st percentile)
Q_FAR_RIGHT_TARGET = 9.8   # e.g. 9.0 — value at the far-right quantile
P_FAR_RIGHT        = 0.99   # Far-right quantile probability (default 99th percentile)

# Optional initial guess for [a, b, loc, scale]. Set to None to let DE start fresh.
INIT_PARAMS = None

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def analytical_mode(a, b, loc, scale):
    """Closed-form mode of a shifted/scaled Beta(a, b, loc, scale). Requires a, b > 1."""
    return loc + scale * (a - 1) / (a + b - 2)


def objective(params, mode_target, median_target, q_lower_target, q_upper_target,
              p_lower, p_upper,
              q_far_left_target, p_far_left,
              q_far_right_target, p_far_right,
              weight_far_left = 1.0, weight_far_right = 0.10):
    """
    Sum of squared differences between target and analytical Beta statistics.

    Inputs:
    - params (list)              : [a, b, loc, scale] Beta distribution parameters.
    - mode_target (float)        : Target mode.
    - median_target (float)      : Target median.
    - q_lower_target (float)     : Target value at lower ETI quantile.
    - q_upper_target (float)     : Target value at upper ETI quantile.
    - p_lower (float)            : Lower quantile probability (e.g. 0.16).
    - p_upper (float)            : Upper quantile probability (e.g. 0.84).
    - q_far_left_target (float)  : Target value at far-left quantile (or None to skip).
    - p_far_left (float)         : Far-left quantile probability (e.g. 0.01).
    - q_far_right_target (float) : Target value at far-right quantile (or None to skip).
    - p_far_right (float)        : Far-right quantile probability (e.g. 0.99).

    Output:
    - float: Objective value (lower = better fit).
    """
    a, b, loc, scale = params

    # a, b > 1 required for a well-defined interior mode
    if a <= 1.0 or b <= 1.0 or scale <= 0.0:
        return np.inf

    try:
        dist = ss.beta(a, b, loc=loc, scale=scale)

        mode_val   = analytical_mode(a, b, loc, scale)
        median_val = dist.median()
        q_lower    = dist.ppf(p_lower)
        q_upper    = dist.ppf(p_upper)

        loss = (
            (mode_target    - mode_val)   ** 2 +
            (median_target  - median_val) ** 2 +
            (q_lower_target - q_lower)    ** 2 +
            (q_upper_target - q_upper)    ** 2
        )

        if q_far_left_target is not None:
            loss += weight_far_left * (q_far_left_target - dist.ppf(p_far_left)) ** 2

        if q_far_right_target is not None:
            loss += weight_far_right * (q_far_right_target - dist.ppf(p_far_right)) ** 2

        return loss

    except Exception as e:
        print(f"  [objective] Exception: {e} | params: {params}")
        return np.inf


def fit_beta(mode_target, median_target, q_lower_target, q_upper_target,
             p_lower=0.16, p_upper=0.84,
             q_far_left_target=None, p_far_left=0.01,
             q_far_right_target=None, p_far_right=0.99,
             init_params=None, bounds=None, tol=1e-8, maxiter=5000, seed=42):
    """
    Optimise Beta distribution parameters to match a target mode, median, and ETI.

    Inputs:
    - mode_target (float)        : Target mode.
    - median_target (float)      : Target median.
    - q_lower_target (float)     : Target value at lower quantile.
    - q_upper_target (float)     : Target value at upper quantile.
    - p_lower (float)            : Lower quantile probability (default 0.16).
    - p_upper (float)            : Upper quantile probability (default 0.84).
    - q_far_left_target (float)  : Target value at far-left quantile (or None to skip).
    - p_far_left (float)         : Far-left quantile probability (default 0.01).
    - q_far_right_target (float) : Target value at far-right quantile (or None to skip).
    - p_far_right (float)        : Far-right quantile probability (default 0.99).
    - init_params (list)         : Optional initial guess for [a, b, loc, scale].
    - bounds (list)              : Parameter bounds (defaults to broad uninformative ranges).
    - tol (float)                : Convergence tolerance (default 1e-8).
    - maxiter (int)              : Maximum DE iterations (default 5000).
    - seed (int)                 : Random seed for reproducibility (default 42).

    Output:
    - np.ndarray: Optimised [a, b, loc, scale], or None if optimisation fails.
    """
    if bounds is None:
        bounds = [(1.01, 30.0), (1.01, 30.0), (-100.0, 100.0), (0.01, 100.0)]

    # Seed population with initial guess if provided
    if init_params is not None:
        rng = np.random.default_rng(seed)
        extra = [
            np.array([rng.uniform(lo, hi) for (lo, hi) in bounds])
            for _ in range(14)  # DE default popsize=15, so 14 extras + init
        ]
        init_population = [np.asarray(init_params)] + extra
    else:
        init_population = "latinhypercube"

    result = differential_evolution(
        objective,
        bounds,
        strategy="best1bin",
        tol=tol,
        maxiter=maxiter,
        workers=-1,
        args=(mode_target, median_target, q_lower_target, q_upper_target,
              p_lower, p_upper,
              q_far_left_target, p_far_left,
              q_far_right_target, p_far_right),
        init=init_population,
        seed=seed,
        polish=True,   # Final L-BFGS-B polish pass
    )

    if not result.success:
        print(f"  [fit_beta] Optimisation warning: {result.message}")

    return result.x if result.fun < 1e6 else None


# =============================================================================
# PLOTTING
# =============================================================================

def make_plot(a, b, loc, scale,
              mode_target, median_target, q_lower_target, q_upper_target,
              p_lower, p_upper, plot_dir, plot_name,
              plot_xrange=None, show_plot=False):
    """
    Plot the fitted Beta distribution alongside the target constraints.
    Saves to plot_dir/plot_name.pdf and plot_dir/plot_name.png.
    """
    os.makedirs(plot_dir, exist_ok=True)

    dist = ss.beta(a, b, loc=loc, scale=scale)

    median_fit = dist.median()
    q_lo_fit   = dist.ppf(p_lower)
    q_hi_fit   = dist.ppf(p_upper)

    if plot_xrange is not None:
        x_lo, x_hi = plot_xrange
    else:
        x_lo = min(q_lower_target, q_lo_fit) - 2.0
        x_hi = max(q_upper_target, q_hi_fit) + 2.0
    x   = np.linspace(x_lo, x_hi, 2000)
    pdf = dist.pdf(x)

    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    # ETI shading — translucent
    ax.axvspan(q_lower_target, q_upper_target, color="C0",  alpha=0.15, zorder=1)
    ax.axvspan(q_lo_fit,       q_hi_fit,       color="red", alpha=0.15, zorder=1)

    # Median lines — dashed
    ax.axvline(median_target, color="C0",  lw=0.8, ls="--", zorder=4)
    ax.axvline(median_fit,    color="red", lw=0.8, ls="--", zorder=4)

    # Fitted distribution — black
    ax.plot(x, pdf, color="black", lw=0.8, zorder=3)

    ax.axhline(0, color="black", lw=0.5, zorder=0)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel("Distance (kpc)")
    ax.set_ylabel("Probability density")

    handles = [
        mpatches.Patch(facecolor="C0",  edgecolor="C0",  label="Target"),
        mpatches.Patch(facecolor="red", edgecolor="red", label="Approximated"),
    ]
    ax.legend(handles=handles, loc="upper right")

    fig.tight_layout()

    for ext in ("pdf", "png"):
        out = os.path.join(plot_dir, f"{plot_name}.{ext}")
        fig.savefig(out, dpi=150)
        print(f"  Saved: {out}")

    if show_plot:
        plt.show()
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("  Beta distribution fitter")
    print("=" * 60)
    print(f"  Targets:")
    print(f"    Mode   : {MODE_TARGET}")
    print(f"    Median : {MEDIAN_TARGET}")
    if Q_FAR_LEFT_TARGET is not None:
        print(f"    q{P_FAR_LEFT:.0%}   : {Q_FAR_LEFT_TARGET}")
    print(f"    q{P_LOWER:.0%}   : {Q_LOWER_TARGET}")
    print(f"    q{P_UPPER:.0%}   : {Q_UPPER_TARGET}")
    if Q_FAR_RIGHT_TARGET is not None:
        print(f"    q{P_FAR_RIGHT:.0%}  : {Q_FAR_RIGHT_TARGET}")
    print()

    params = fit_beta(
        MODE_TARGET, MEDIAN_TARGET, Q_LOWER_TARGET, Q_UPPER_TARGET,
        p_lower=P_LOWER, p_upper=P_UPPER,
        q_far_left_target=Q_FAR_LEFT_TARGET, p_far_left=P_FAR_LEFT,
        q_far_right_target=Q_FAR_RIGHT_TARGET, p_far_right=P_FAR_RIGHT,
        init_params=INIT_PARAMS,
    )

    if params is None:
        print("  Optimisation failed — no valid solution found.")
        return

    a, b, loc, scale = params
    dist = ss.beta(a, b, loc=loc, scale=scale)

    print("  Best-fit parameters:")
    print(f"  alpha: {a}")
    print(f"  beta: {b}")
    print(f"  loc: {loc}")
    print(f"  scale: {scale}")
    print()
    print("  Recovered statistics vs. targets:")
    print(f"    {'Statistic':<12} {'Target':>10} {'Fitted':>10} {'Residual':>10}")
    print(f"    {'-'*44}")

    stats = [
        ("Mode",              MODE_TARGET,    analytical_mode(a, b, loc, scale)),
        ("Median",            MEDIAN_TARGET,  dist.median()),
        (f"q{P_LOWER:.0%}",  Q_LOWER_TARGET, dist.ppf(P_LOWER)),
        (f"q{P_UPPER:.0%}",  Q_UPPER_TARGET, dist.ppf(P_UPPER)),
    ]
    if Q_FAR_LEFT_TARGET is not None:
        stats.insert(2, (f"q{P_FAR_LEFT:.0%}", Q_FAR_LEFT_TARGET, dist.ppf(P_FAR_LEFT)))
    if Q_FAR_RIGHT_TARGET is not None:
        stats.append((f"q{P_FAR_RIGHT:.0%}", Q_FAR_RIGHT_TARGET, dist.ppf(P_FAR_RIGHT)))
    for name, target, fitted in stats:
        print(f"    {name:<12} {target:>10.4f} {fitted:>10.4f} {fitted - target:>+10.4f}")

    print()
    make_plot(
        a, b, loc, scale,
        MODE_TARGET, MEDIAN_TARGET, Q_LOWER_TARGET, Q_UPPER_TARGET,
        P_LOWER, P_UPPER,
        PLOT_DIR, PLOT_NAME,
        plot_xrange=PLOT_XRANGE,
        show_plot=SHOW_PLOT,
    )


if __name__ == "__main__":
    main()
