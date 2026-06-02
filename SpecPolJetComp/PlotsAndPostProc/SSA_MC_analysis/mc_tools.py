#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:54:19 2026

@author: cowie
"""

import os
import time
import numpy as np
from typing import Callable, Dict, Tuple, Optional, Sequence, List, Union, Any
from dataclasses import dataclass
import math
import matplotlib.pyplot as plt
import warnings
from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar
import arviz as az

# --- Types & results ---------------------------------------------------------
ParamSpec = Union[Tuple[float, float], Dict[str, Any]]

@dataclass
class MCResult:
    mean: float
    std: float
    median: float
    mode: float
    ci: Tuple[float, float]            # central equal-tailed interval
    hdi: Tuple[float, float]           # highest density interval (arviz)
    samples: np.ndarray                # posterior samples of the output quantity
    param_samples: Dict[str, np.ndarray]  # raw draws per input parameter
    func_name: str = ""                # name of the sampled function (set automatically)

# --- Spec handling & sampling -----------------------------------------------
def _to_lognormal_log_params(mean: float, sigma: float) -> Tuple[float, float]:
    if mean <= 0:
        raise ValueError("Lognormal 'mean' must be > 0 when space='linear'.")
    if sigma < 0:
        raise ValueError("Lognormal 'sigma' must be >= 0.")
    if sigma == 0.0:
        return float(np.log(mean)), 0.0
    sigma2 = sigma**2
    mu_log = np.log(mean**2 / np.sqrt(sigma2 + mean**2))
    sigma_log = np.sqrt(np.log(1.0 + sigma2 / (mean**2)))
    return float(mu_log), float(sigma_log)

def _kde_mode(samples: np.ndarray) -> float:
    """Estimate the mode of a continuous sample via KDE maximum."""
    kde = gaussian_kde(samples)
    lo, hi = np.percentile(samples, [1, 99])
    res = minimize_scalar(lambda x: -float(kde(x)), bounds=(lo, hi), method="bounded")
    return float(res.x)

# --- Spec handling & sampling -----------------------------------------------
def _resolve_spec(spec: ParamSpec) -> Dict[str, Any]:
    if isinstance(spec, tuple):
        m, s = spec
        return {
            "dist": "normal",
            "sampler_params": (float(m), float(s)),
            "space": "normal",
            "bounds": None,
            "nominal": float(m),
            "strict_positive": False,
            "eps": 1e-300,
        }

    if not isinstance(spec, dict) or "dist" not in spec:
        raise ValueError("Spec must be tuple or dict with a 'dist' key.")

    dist = spec["dist"].lower()
    bounds = spec.get("bounds", None)
    strict_positive = bool(spec.get("strict_positive", False))
    eps = float(spec.get("eps", 1e-300))

    if dist == "uniform":
        low = float(spec["low"])
        high = float(spec["high"])
        if not (high > low):
            raise ValueError("For uniform, require high > low.")
        space = spec.get("space", "linear").lower()
        if space == "linear":
            nominal = 0.5 * (low + high)
        elif space == "log":
            if low <= 0:
                raise ValueError("Log-uniform requires low > 0.")
            nominal = np.sqrt(low * high)  # geometric mean
        else:
            raise ValueError("space must be 'linear' or 'log' for uniform.")
        return {
            "dist": "uniform",
            "sampler_params": (low, high),
            "space": space,
            "bounds": bounds,
            "nominal": nominal,
            "strict_positive": strict_positive,
            "eps": eps,
        }

    if dist == "normal":
        m = float(spec["mean"])
        s = float(spec["sigma"]) if "sigma" in spec else float(abs(m) * spec["sigma_frac"])
        return {"dist": "normal", "sampler_params": (m, s), "space": "normal", "bounds": bounds,
                "nominal": m, "strict_positive": strict_positive, "eps": eps}

    if dist == "lognormal":
        space = spec.get("space", "linear").lower()
        if space == "linear":
            m = float(spec["mean"])
            s = float(spec["sigma"]) if "sigma" in spec else float(abs(m) * spec["sigma_frac"])
            mu_log, sigma_log = _to_lognormal_log_params(m, s)
            nominal = m
        elif space == "log":
            mu_log = float(spec["mean"])
            sigma_log = float(spec["sigma"])
            nominal = float(np.exp(mu_log))
        else:
            raise ValueError("space must be 'linear' or 'log' for lognormal.")
        return {"dist": "lognormal", "sampler_params": (mu_log, sigma_log), "space": "log", "bounds": bounds,
                "nominal": nominal, "strict_positive": True, "eps": eps}

    # --- New: asymmetric (split) normal ---
    if dist in ("asymm_normal", "asymmetric_normal", "split_normal", "two_piece_normal"):
        m = float(spec["mean"])
        # Prefer explicit sigma_up/down if provided, otherwise use sigma_frac if present.
        if "sigma_up" in spec or "sigma_down" in spec:
            if "sigma_up" not in spec or "sigma_down" not in spec:
                raise ValueError("asymm_normal requires both 'sigma_up' and 'sigma_down' if one is provided.")
            sigma_up = float(spec["sigma_up"])
            sigma_down = float(spec["sigma_down"])
        else:
            if "sigma_frac" not in spec:
                raise ValueError("asymm_normal requires 'sigma_up' and 'sigma_down', or 'sigma_frac' to set both.")
            sigma_up = float(abs(m) * spec["sigma_frac"])
            sigma_down = float(abs(m) * spec["sigma_frac"])
        if sigma_up < 0 or sigma_down < 0:
            raise ValueError("Asymmetric normal sigmas must be >= 0.")
        # Store as (mean, sigma_down, sigma_up) to be consistent when sampling
        return {"dist": "asymm_normal", "sampler_params": (m, sigma_down, sigma_up), "space": "normal",
                "bounds": bounds, "nominal": m, "strict_positive": strict_positive, "eps": eps}

    if dist == "beta":
        a     = float(spec["alpha"])
        b     = float(spec["beta"])
        if a <= 0 or b <= 0:
            raise ValueError("Beta distribution requires alpha > 0 and beta > 0.")
        loc   = float(spec.get("loc",   0.0))
        scale = float(spec.get("scale", 1.0))
        if scale <= 0:
            raise ValueError("Beta distribution requires scale > 0.")
        nominal = loc + (a / (a + b)) * scale   # distribution mean
        return {
            "dist": "beta",
            "sampler_params": (a, b, loc, scale),
            "space": "linear",
            "bounds": bounds,
            "nominal": nominal,
            "strict_positive": strict_positive,
            "eps": eps,
        }

    raise ValueError("Unsupported 'dist'.")

def _sample_param(rng: np.random.Generator, resolved: Dict[str, Any], n: int) -> np.ndarray:
    dist = resolved["dist"]
    lo_hi = resolved.get("bounds", None)
    eps = float(resolved.get("eps", 1e-300))

    if dist == "normal":
        mu, sigma = resolved["sampler_params"]
        if sigma < 0:
            raise ValueError("Sigma must be >= 0.")
        x = rng.normal(mu, sigma, size=n)

    elif dist == "lognormal":
        mu_log, sigma_log = resolved["sampler_params"]
        if sigma_log < 0:
            raise ValueError("Sigma must be >= 0.")
        x = rng.lognormal(mean=mu_log, sigma=sigma_log, size=n)

    elif dist == "uniform":
        low, high = resolved["sampler_params"]
        if resolved["space"] == "linear":
            x = rng.uniform(low, high, size=n)
        elif resolved["space"] == "log":
            loglow, loghigh = np.log10(low), np.log10(high)
            logx = rng.uniform(loglow, loghigh, size=n)
            x = np.power(10.0, logx)
        else:
            raise RuntimeError("Unexpected space for uniform.")

    elif dist == "asymm_normal":
        # Split/Two-piece normal: draw standard normal and scale positive/negative sides differently
        mu, sigma_down, sigma_up = resolved["sampler_params"]
        if (sigma_down < 0) or (sigma_up < 0):
            raise ValueError("Asymmetric normal sigmas must be >= 0.")
        z = rng.normal(0.0, 1.0, size=n)
        # For z >= 0 use sigma_up, for z < 0 use sigma_down
        scales = np.where(z >= 0.0, sigma_up, sigma_down)
        x = mu + z * scales

    elif dist == "beta":
        a, b, loc, scale = resolved["sampler_params"]
        x = rng.beta(a, b, size=n) * scale + loc

    else:
        raise RuntimeError("Unexpected dist.")

    if lo_hi is not None:
        lo, hi = lo_hi
        if lo is not None:
            x = np.maximum(x, lo)
        if hi is not None:
            x = np.minimum(x, hi)

    if resolved.get("strict_positive", False):
        x = np.where(x <= eps, eps, x)

    return x

def draw_params(
    params: Dict[str, ParamSpec],
    n: int,
    seed: Optional[int] = None,
    subset: Optional[Sequence[str]] = None,
) -> Dict[str, np.ndarray]:
    """Draw n samples for all (or a named subset of) parameters."""
    rng = np.random.default_rng(seed)
    keys = list(subset) if subset is not None else list(params.keys())
    resolved = {k: _resolve_spec(params[k]) for k in keys}
    return {k: _sample_param(rng, resolved[k], n) for k in keys}

# --- Diagnostic helpers ------------------------------------------------------
def _param_summary_line(name: str, spec: Dict[str, Any]) -> str:
    dist = spec["dist"]
    nominal = spec["nominal"]

    if dist == "normal":
        mu, sigma = spec["sampler_params"]
        desc = f"normal(μ={mu:.4g}, σ={sigma:.4g})"
    elif dist == "lognormal":
        mu_log, sigma_log = spec["sampler_params"]
        desc = f"lognormal(μ_log={mu_log:.4g}, σ_log={sigma_log:.4g})  nominal={nominal:.4g}"
    elif dist == "uniform":
        low, high = spec["sampler_params"]
        space = spec.get("space", "linear")
        desc = f"uniform({low:.4g} → {high:.4g})  space={space}"
    elif dist == "asymm_normal":
        mu, sigma_down, sigma_up = spec["sampler_params"]
        desc = f"asymm_normal(μ={mu:.4g}, σ↓={sigma_down:.4g}, σ↑={sigma_up:.4g})"
    elif dist == "beta":
        a, b, loc, scale = spec["sampler_params"]
        desc = f"beta(α={a:.4g}, β={b:.4g}, loc={loc:.4g}, scale={scale:.4g})  support=[{loc:.4g}, {loc + scale:.4g}]"
    else:
        desc = dist

    extras = []
    bounds = spec.get("bounds")
    if bounds is not None:
        lo_b, hi_b = bounds
        lo_s = str(lo_b) if lo_b is not None else "−∞"
        hi_s = str(hi_b) if hi_b is not None else "+∞"
        extras.append(f"bounds=[{lo_s}, {hi_s}]")
    if spec.get("strict_positive"):
        extras.append("strict_positive")

    suffix = "  " + ", ".join(extras) if extras else ""
    return f"    {name:<25} {desc}{suffix}"


# --- Monte Carlo core --------------------------------------------------------
def mc_uncertainty(
    f: Callable[..., float],
    params: Dict[str, ParamSpec],
    n: int = 100_000,
    seed: Optional[int] = None,
    ci_level: float = 0.68,
    fixed_kwargs: Optional[Dict[str, Any]] = None,
    fixed_draws: Optional[Dict[str, np.ndarray]] = None,
) -> MCResult:
    """
    Monte Carlo error propagation for y = f(**params).

    Pass fixed_draws (a dict of pre-drawn arrays) to reuse a shared parameter
    draw across multiple function evaluations (joint sampling mode).  When
    fixed_draws is provided, n and seed are ignored for sampling.
    """
    if fixed_kwargs is None:
        fixed_kwargs = {}

    resolved = {k: _resolve_spec(v) for k, v in params.items()}

    if fixed_draws is not None:
        draws = {k: v.copy() for k, v in fixed_draws.items()}
        n = len(next(iter(draws.values())))
        seed_label = "joint (fixed draws)"
    else:
        rng = np.random.default_rng(seed)
        draws = {k: _sample_param(rng, resolved[k], n) for k in resolved.keys()}
        seed_label = str(seed)

    # --- Run header ---
    print(f"\n{'─' * 62}")
    print(f"  mc_uncertainty  →  {f.__name__}")
    print(f"  n = {n:,}  |  seed = {seed_label}  |  CI = {ci_level * 100:.0f}%")
    if fixed_kwargs:
        print(f"  fixed_kwargs : {fixed_kwargs}")

    draw_label = "Fixed draws" if fixed_draws is not None else "Sampling"
    print(f"  {'─' * 58}")
    print(f"  {draw_label} — {len(resolved)} parameter(s):")
    for name, spec in resolved.items():
        print(_param_summary_line(name, spec))
    print(f"{'─' * 62}")

    # --- Vectorised evaluation (fast path) -----------------------------------
    # Try calling f with the full draw arrays.  Functions whose internals are
    # all numpy-compatible (tau_m_interp, K_E_func, etc.) will complete in one
    # call.  Functions with a per-sample scalar solver (e.g. gamma_min_constraint)
    # will raise an exception and fall through to the sample-by-sample loop.
    t0 = time.perf_counter()
    try:
        y_raw = f(**draws, **fixed_kwargs)
        y = np.asarray(y_raw, dtype=np.complex128)
        if y.shape != (n,):
            raise ValueError(f"Expected output length {n}, got shape {y.shape}.")
        elapsed_total = time.perf_counter() - t0
        print(f"  Vectorised evaluation done in {elapsed_total:.2f}s")

    except Exception as exc:
        # --- Sample-by-sample fallback (for functions with scalar solvers) ---
        print(f"  Vectorised call failed ({type(exc).__name__}); "
              f"falling back to sample-by-sample loop ...")
        y = np.empty(n, dtype=np.complex128)
        report_every = max(1, n // 20)
        for i in range(n):
            if i % report_every == 0:
                pct = i / n * 100
                elapsed = time.perf_counter() - t0
                print(f"  Evaluating ...  {pct:5.1f}%   ({elapsed:.1f}s elapsed)",
                      end="\r", flush=True)
            kwargs_i = {k: v[i] for k, v in draws.items()}
            y[i] = f(**kwargs_i, **fixed_kwargs)
        elapsed_total = time.perf_counter() - t0
        print(f"  Evaluating ...  100.0%  ({elapsed_total:.1f}s for {n:,} samples)    ")

    # Warn if any complex values were produced
    imag = np.abs(np.imag(y))
    n_complex = int(np.sum(imag > 0))
    if n_complex > 0:
        warnings.warn(
            f"mc_uncertainty: {n_complex} samples returned complex values; "
            f"their imaginary parts will be discarded."
        )
        print(f"  [warn] {n_complex:,} complex samples — imaginary parts discarded.")

    # Keep only the real part for downstream statistics
    y = np.real(y)

    # Remove invalid values
    mask = np.isfinite(y)
    n_valid = int(np.sum(mask))
    n_dropped = n - n_valid
    if n_dropped > 0:
        warnings.warn(
            f"mc_uncertainty: {n_dropped} invalid samples were discarded; "
            f"effective sample size = {n_valid}."
        )
        print(f"  [warn] {n_dropped:,} non-finite samples discarded  "
              f"(effective n = {n_valid:,}, {n_valid / n * 100:.2f}%).")

    y = y[mask]
    draws = {k: v[mask] for k, v in draws.items()}

    if len(y) == 0:
        raise RuntimeError("All Monte Carlo samples were invalid.")

    mean   = float(np.mean(y))
    median = float(np.median(y))
    std    = float(np.std(y, ddof=1))
    mode   = _kde_mode(y)
    lo_ci, hi_ci = np.quantile(y, [(1 - ci_level) / 2, 1 - (1 - ci_level) / 2])
    lo_hdi, hi_hdi = az.hdi(y, hdi_prob=ci_level)

    # --- Results summary ---
    print(f"\n  Results — {f.__name__}")
    print(f"  {'─' * 56}")
    valid_pct = n_valid / n * 100
    print(f"  Valid samples : {n_valid:,} / {n:,}  ({valid_pct:.1f}%)")
    print(f"  Mean          : {mean:.6g}")
    print(f"  Median        : {median:.6g}")
    print(f"  Mode (KDE)    : {mode:.6g}")
    print(f"  Std dev       : {std:.6g}")
    print(f"  {ci_level * 100:.0f}% CI (equal-tail) : [{float(lo_ci):.6g},  {float(hi_ci):.6g}]")
    print(f"  {ci_level * 100:.0f}% HDI             : [{float(lo_hdi):.6g},  {float(hi_hdi):.6g}]")
    print(f"{'─' * 62}\n")

    return MCResult(
        mean=mean,
        std=std,
        median=median,
        mode=mode,
        ci=(float(lo_ci), float(hi_ci)),
        hdi=(float(lo_hdi), float(hi_hdi)),
        samples=y,
        param_samples=draws,
        func_name=f.__name__,
    )

# --- Display label map -------------------------------------------------------
# Maps function names to clean LaTeX axis labels.  Falls back to the raw
# function name (underscores stripped) for anything not listed here.
_FUNC_LABELS: Dict[str, str] = {
    "E_energy_form":        r"$\log_{10}(E\ /\ \mathrm{erg})$",
    "E_frequency_form":     r"$\log_{10}(E\ /\ \mathrm{erg})$",
    "R_energy_form":        r"$\log_{10}(R\ /\ \mathrm{cm})$",
    "R_frequency_form":     r"$\log_{10}(R\ /\ \mathrm{cm})$",
    "B_energy_form":        r"$\log_{10}(B\ /\ \mathrm{G})$",
    "B_frequency_form":     r"$\log_{10}(B\ /\ \mathrm{G})$",
    "ne_energy_form":       r"$\log_{10}(n_e\ /\ \mathrm{cm}^{-3})$",
    "ne_frequency_form":    r"$\log_{10}(n_e\ /\ \mathrm{cm}^{-3})$",
    "N0_energy_form":       r"$\log_{10}(N_0\ /\ \mathrm{cm}^{-3})$",
    "N0_frequency_form":    r"$\log_{10}(N_0\ /\ \mathrm{cm}^{-3})$",
    "TB_energy_form":       r"$\log_{10}(T_B\ /\ \mathrm{K})$",
    "TB_frequency_form":    r"$\log_{10}(T_B\ /\ \mathrm{K})$",
    "gamma_min_constraint": r"$\log_{10}(\gamma_{\min})$",
}

# --- Plotting utilities ------------------------------------------------------
def plot_mc_diagnostics(
    result,
    params,
    ci_level: float = 0.68,
    bins: int = 60,
    max_cols: int = 3,
    plots_dir: str = "",
    output_stem: str = "mc_out",
    show_plot: bool = False,
    show_params: bool = True,
):
    """
    Save two diagnostic PDFs and optionally display them:
      1. <stem>_<func>_params.pdf     — grid of input-parameter histograms
      2. <stem>_<func>_posterior.pdf  — posterior histogram of the output quantity

    The function name is taken from result.func_name (set automatically by
    mc_uncertainty) and is used as both the output-file suffix and the x-axis
    label on the posterior plot.

    Parameters
    ----------
    plots_dir   : directory to write files into; created if it does not exist.
    output_stem : base filename stem before the function name is appended.
    show_plot   : call plt.show() after saving (set False for batch/headless runs).
    """
    func_name = result.func_name or "output"
    x_label = _FUNC_LABELS.get(func_name, func_name.replace("_", " "))

    # Build output paths; create directory if needed
    dir_ = plots_dir.strip() if (plots_dir and plots_dir.strip()) else os.getcwd()
    os.makedirs(dir_, exist_ok=True)
    base = f"{output_stem}_{func_name}"
    path_params = os.path.join(dir_, f"{base}_params.pdf")
    path_post   = os.path.join(dir_, f"{base}_posterior.pdf")

    # Resolve distributions so we know which parameters are log-space
    resolved = {k: _resolve_spec(v) for k, v in params.items()}
    names = list(result.param_samples.keys())
    K = len(names)

    # --- 1. Parameter samples grid (optional) ---
    fig_params = None
    if show_params:
        rows = math.ceil(K / max_cols)
        fig_params, axes = plt.subplots(rows, min(K, max_cols),
                                        figsize=(5 * min(K, max_cols), 3.6 * rows))
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        axes = axes.reshape(rows, -1)

        for idx, name in enumerate(names):
            r, c = divmod(idx, max_cols)
            ax = axes[r, c]
            x = result.param_samples[name]
            spec = resolved[name]

            if spec["dist"] == "lognormal" or (spec["dist"] == "uniform" and spec["space"] == "log"):
                logx = np.log10(x)
                ax.hist(logx, bins=bins, color="C0", alpha=0.7)
                ax.set_xlabel("log10(" + name + ")")
                ax.axvline(np.log10(spec["nominal"]), color="k")
            else:
                ax.hist(x, bins=bins, color="C0", alpha=0.7)
                ax.set_xlabel(name)
                ax.axvline(spec["nominal"], color="k")

            ax.set_title(name)
            ax.set_ylabel("count")

        # Hide any unused axes
        for idx in range(K, rows * min(K, max_cols)):
            r, c = divmod(idx, max_cols)
            axes[r, c].set_visible(False)

        fig_params.suptitle(f"Input parameters — {x_label}", y=1.01)
        fig_params.tight_layout()
        fig_params.savefig(path_params, dpi=300, bbox_inches="tight")
        print(f"Saved parameter diagnostics: {path_params}")

    # --- 2. Posterior histogram ---
    y = result.samples
    lo_hdi, hi_hdi = result.hdi
    lo_ci,  hi_ci  = result.ci
    lo_95,  hi_95  = az.hdi(y, hdi_prob=0.95)
    lo_997, hi_997 = az.hdi(y, hdi_prob=0.997)

    fig_post, axp = plt.subplots(figsize=(10, 6))

    # Normalised histogram so the KDE overlay is on the same scale
    axp.hist(y, bins=bins, density=True, color="C1", alpha=0.45, zorder=2)

    # Smooth KDE overlay
    kde  = gaussian_kde(y)
    xkde = np.linspace(y.min(), y.max(), 600)
    axp.plot(xkde, kde(xkde), color="C1", lw=1.8, zorder=3)

    # Shaded credible intervals (back to front so narrower sits on top)
    axp.axvspan(lo_997, hi_997, alpha=0.10, color="steelblue", label="99.7% HDI")
    axp.axvspan(lo_95,  hi_95,  alpha=0.15, color="steelblue", label="95% HDI")
    axp.axvspan(lo_hdi, hi_hdi, alpha=0.25, color="steelblue",
                label=f"{int(ci_level * 100)}% HDI")

    # Point estimates
    axp.axvline(result.mode,   color="C3", lw=1.8, ls="-",  zorder=5, label="mode")
    axp.axvline(result.mean,   color="k",  lw=1.5, ls="--", zorder=5, label="mean")
    axp.axvline(result.median, color="C0", lw=1.5, ls=":",  zorder=5, label="median")

    axp.set_xlabel(x_label, fontsize=13)
    axp.set_ylabel("Density", fontsize=12)
    axp.set_title(x_label, fontsize=14)
    axp.legend(fontsize=10, loc="upper left")

    # Stats box — right side to avoid the histogram peak
    stats_lines = [
        f"Mode    {result.mode:.4g}",
        f"Mean    {result.mean:.4g}",
        f"Median  {result.median:.4g}",
        f"σ       {result.std:.4g}",
        "",
        f"{int(ci_level*100)}% HDI  [{lo_hdi:.4g}, {hi_hdi:.4g}]",
        f"95% HDI  [{lo_95:.4g}, {hi_95:.4g}]",
        f"99.7% HDI  [{lo_997:.4g}, {hi_997:.4g}]",
    ]
    axp.text(0.975, 0.97, "\n".join(stats_lines),
             transform=axp.transAxes, va="top", ha="right",
             fontsize=9, family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                       edgecolor="0.75", alpha=0.9))

    fig_post.tight_layout()
    fig_post.savefig(path_post, dpi=300, bbox_inches="tight")
    print(f"Saved posterior histogram:   {path_post}")

    if show_plot:
        plt.show()
    else:
        plt.close("all")

    out = {"posterior": fig_post}
    if fig_params is not None:
        out["params"] = fig_params
    return out


def plot_mc_input_params(
    draws: Dict[str, np.ndarray],
    params: Dict[str, ParamSpec],
    plots_dir: str = "",
    output_stem: str = "mc_out",
    show_plot: bool = False,
    bins: int = 60,
    max_cols: int = 3,
    title: str = "Input parameters — joint draw",
):
    """Save a parameter-grid PDF from a pre-drawn samples dict."""
    dir_ = plots_dir.strip() if (plots_dir and plots_dir.strip()) else os.getcwd()
    os.makedirs(dir_, exist_ok=True)
    path = os.path.join(dir_, f"{output_stem}_joint_params.pdf")

    resolved = {k: _resolve_spec(v) for k, v in params.items() if k in draws}
    names = [k for k in draws if k in resolved]
    K = len(names)
    rows = math.ceil(K / max_cols)

    fig, axes = plt.subplots(rows, min(K, max_cols),
                              figsize=(5 * min(K, max_cols), 3.6 * rows))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.reshape(rows, -1)

    for idx, name in enumerate(names):
        r, c = divmod(idx, max_cols)
        ax = axes[r, c]
        x = draws[name]
        spec = resolved[name]

        if spec["dist"] == "lognormal" or (spec["dist"] == "uniform" and spec["space"] == "log"):
            logx = np.log10(x)
            ax.hist(logx, bins=bins, color="C0", alpha=0.7)
            ax.set_xlabel("log10(" + name + ")")
            ax.axvline(np.log10(spec["nominal"]), color="k")
        else:
            ax.hist(x, bins=bins, color="C0", alpha=0.7)
            ax.set_xlabel(name)
            ax.axvline(spec["nominal"], color="k")

        ax.set_title(name)
        ax.set_ylabel("count")

    for idx in range(K, rows * min(K, max_cols)):
        r, c = divmod(idx, max_cols)
        axes[r, c].set_visible(False)

    fig.suptitle(title, y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved joint parameter diagnostics: {path}")

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    return fig


def plot_mc_corner(
    chain: np.ndarray,
    labels: List[str],
    plots_dir: str = "",
    output_stem: str = "mc_out",
    show_plot: bool = False,
):
    """
    Corner plot of jointly-sampled output quantities using the corner library.

    chain  : (n, d) float array — one column per output quantity.
    labels : d LaTeX strings for axis labels.
    Truths (KDE mode) and ranges (99.7% HDI) are computed automatically.
    """
    import corner as corner_lib

    _, d = chain.shape

    truths = [_kde_mode(chain[:, i]) for i in range(d)]
    ranges = [tuple(az.hdi(chain[:, i], hdi_prob=0.997)) for i in range(d)]

    corner_kwargs = {
        "labels":          labels,
        "label_kwargs":    {"fontsize": 14},
        "title_kwargs":    {"fontsize": 10},
        "quantiles":       [0.16, 0.5, 0.84],
        "show_titles":     False,
        "smooth":          2.5,
        "bins":            50,
        "hist_kwargs":     {"alpha": 0.8, "color": "steelblue"},
        "color":           "steelblue",
        "plot_datapoints": True,
        "data_kwargs":     {"alpha": 0.1},
        "plot_density":    True,
        "contour_kwargs":  {"colors": "darkblue", "linewidths": 1.2},
        "contourf_kwargs": {"colors": ["lightblue", "steelblue"], "alpha": 0.6},
        "range":           ranges,
        "truths":          truths,
        "truth_color":     "red",
        "truth_kwargs":    {"linewidth": 2, "alpha": 0.8},
    }

    fig = corner_lib.corner(chain, **corner_kwargs)

    dir_ = plots_dir.strip() if (plots_dir and plots_dir.strip()) else os.getcwd()
    os.makedirs(dir_, exist_ok=True)
    path = os.path.join(dir_, f"{output_stem}_corner.pdf")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved corner plot:           {path}")

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    return fig


def save_mc_samples(result, samples_dir: str = "", output_stem: str = "mc_out"):
    """Save the posterior output samples as a plain .npy file.

    File name: <stem>_<func_name>_samples.npy
    Load with: np.load(path)
    """
    dir_ = samples_dir.strip() if (samples_dir and samples_dir.strip()) else os.getcwd()
    os.makedirs(dir_, exist_ok=True)
    func_name = result.func_name or "output"
    path = os.path.join(dir_, f"{output_stem}_{func_name}_samples.npy")
    np.save(path, result.samples)
    print(f"Saved MC samples:            {path}")