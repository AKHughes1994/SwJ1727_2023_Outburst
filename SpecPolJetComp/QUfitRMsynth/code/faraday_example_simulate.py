#!/usr/bin/env python3
"""
Example: Simulate and Fit Broadband Polarization Data

Generate synthetic polarized data from configurable multi-component models,
fit with Bayesian nested sampling (dynesty), and produce a full diagnostic
plotting suite.

Workflow:
1. Define true model components (thin, power-law, thick)
2. Generate synthetic data with power-law Stokes I spectrum
3. Build fit model and specify priors
4. Run nested sampling with dynesty
5. Post-process posteriors (mode detection, derived parameters)
6. Diagnostic plots: fit results, corner plots, convergence, FDF

All outputs saved to: ../simulated/
"""

import matplotlib
matplotlib.use('Agg')

import copy
import subprocess
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from typing import Optional

from faraday_data import DataLoader, PolarizationData
from faraday_utils import (ThinComponent, ThickComponent, ThinPowerLawComponent, CompositeModel,
                           compute_polarization_quantities, compute_rmtf, rm_synthesis, C_LIGHT)
from faraday_model import FitSetup, FaradayFitter, update_model_from_params
from faraday_plot import (plot_data_diagnostic, plot_fit_results, plot_corner,
                         plot_convergence, plot_fdf_diagnostic,
                         setup_publication_quality, add_ticks_all_sides)
from faraday_processing import process_posterior_modes, ProcessingConfig


# =============================================================================
# PATHS
# =============================================================================

output_dir = Path(__file__).parent.parent / "simulated"
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Output directory: {output_dir.absolute()}")

# =============================================================================
# COMPONENT CONFIGURATION
# =============================================================================
# Define all components here using arrays - easily add/remove components
# Components are automatically labeled S1, S2, ... (thin screens),
# P1, P2, ... (powerlaw), T1, T2, ... (thick screens)

# THIN SCREEN COMPONENTS (S1, S2, ...)
THIN_COMPONENTS = [
    {'p0': 0.02, 'psi_0_deg': -2.6, 'phi_rm': -5.76},  # S1 (e.g., ISM)
    #{'p0': 0.01, 'psi_0_deg': 37, 'phi_rm': 0.0},     # S2 (e.g., Jet)
]

# POWERLAW COMPONENTS (P1, P2, ...)
POWERLAW_COMPONENTS = [
#    {'p0': 0.035, 'psi_0_deg': 0.0, 'phi_rm': 30.0, 'beta': -1.5, 'lambda_sq_0': 0.07},  # P1
#    {'p0': 0.025, 'psi_0_deg': -33.0, 'phi_rm': 30.0, 'beta': 4.0, 'lambda_sq_0': 0.07},  # P2
]

# THICK SCREEN COMPONENTS (T1, T2, ...)
THICK_COMPONENTS = [
#   {'p0': 0.006, 'psi_0_deg': 62.4, 'phi_peak': -88.0, 'sigma_phi': 20.8, 'N': 15.6},  # T1 (e.g., Source)
    {'p0': 0.05, 'psi_0_deg': -35.0, 'phi_peak': 53.0, 'sigma_phi':35.0, 'N': 40},  # T2 (e.g., Jet)
]

# MODEL TO FIT
# Specify which components to fit: 'S' (thin), 'P' (powerlaw), 'T' (thick)
# Examples: 'S' (one thin), 'ST' (one thin + one thick), 'SPP' (one thin + two powerlaw)
# Default: fit the same model as the truth (automatically determined)
# Set to None to auto-fit all components in truth model
MODEL_TO_FIT = None  # e.g., 'P', 'ST', 'SPP', 'SST', etc. or None for auto

# If True, run only: simulate → fit → presentation plot (skip corners, convergence, FDF)
QUICK_MODE = False

# Uncomment for reproducible results
# np.random.seed(42)


# =============================================================================
# HELPER FUNCTION: BUILD MODEL FROM CONFIGURATION
# =============================================================================

def build_model_from_config(thin_comps, powerlaw_comps, thick_comps):
    """
    Build a CompositeModel from configuration arrays.

    Args:
        thin_comps: List of dicts with thin component parameters
        powerlaw_comps: List of dicts with powerlaw component parameters
        thick_comps: List of dicts with thick component parameters

    Returns:
        model: CompositeModel
        component_info: List of dicts with component metadata
    """
    model = CompositeModel()
    component_info = []

    # Add thin components (S1, S2, ...)
    for i, comp_params in enumerate(thin_comps, start=1):
        name = f'S{i}'
        p0 = comp_params['p0']
        psi_0 = np.radians(comp_params['psi_0_deg'])
        phi_rm = comp_params['phi_rm']
        a = p0 * np.cos(2 * psi_0)
        b = p0 * np.sin(2 * psi_0)
        model.add_component(ThinComponent(a=a, b=b, phi_rm=phi_rm, name=name))
        component_info.append({
            'label': name,
            'type': 'ThinComponent',
            'params': {'p0': p0, 'psi_0': psi_0, 'psi_0_deg': comp_params['psi_0_deg'],
                      'phi_rm': phi_rm, 'a': a, 'b': b}
        })

    # Add powerlaw components (P1, P2, ...)
    for i, comp_params in enumerate(powerlaw_comps, start=1):
        name = f'P{i}'
        p0 = comp_params['p0']
        psi_0 = np.radians(comp_params['psi_0_deg'])
        phi_rm = comp_params['phi_rm']
        beta = comp_params['beta']
        lambda_sq_0 = comp_params['lambda_sq_0']
        a = p0 * np.cos(2 * psi_0)
        b = p0 * np.sin(2 * psi_0)
        model.add_component(ThinPowerLawComponent(
            a=a, b=b, phi_rm=phi_rm, beta=beta, lambda_sq_0=lambda_sq_0, name=name))
        component_info.append({
            'label': name,
            'type': 'ThinPowerLawComponent',
            'params': {'p0': p0, 'psi_0': psi_0, 'psi_0_deg': comp_params['psi_0_deg'],
                      'phi_rm': phi_rm, 'beta': beta, 'lambda_sq_0': lambda_sq_0,
                      'a': a, 'b': b}
        })

    # Add thick components (T1, T2, ...)
    for i, comp_params in enumerate(thick_comps, start=1):
        name = f'T{i}'
        p0 = comp_params['p0']
        psi_0 = np.radians(comp_params['psi_0_deg'])
        phi_peak = comp_params['phi_peak']
        sigma_phi = comp_params['sigma_phi']
        N = comp_params['N']
        a = p0 * np.cos(2 * psi_0)
        b = p0 * np.sin(2 * psi_0)
        model.add_component(ThickComponent(
            a=a, b=b, phi_peak=phi_peak, sigma_phi=sigma_phi, N=N, name=name))
        component_info.append({
            'label': name,
            'type': 'ThickComponent',
            'params': {'p0': p0, 'psi_0': psi_0, 'psi_0_deg': comp_params['psi_0_deg'],
                      'phi_peak': phi_peak, 'sigma_phi': sigma_phi, 'N': N,
                      'a': a, 'b': b}
        })

    return model, component_info


def parse_model_spec(model_spec, thin_comps, powerlaw_comps, thick_comps):
    """
    Parse model specification string (e.g., 'SPT' or 'SST') to determine which
    components to include.

    Args:
        model_spec: String like 'S', 'ST', 'SPP', 'SST', etc.
        thin_comps, powerlaw_comps, thick_comps: Configuration arrays

    Returns:
        thin_to_use: List of thin component configs
        powerlaw_to_use: List of powerlaw component configs
        thick_to_use: List of thick component configs
    """
    if model_spec is None:
        # Use all components
        return thin_comps, powerlaw_comps, thick_comps

    # Count how many of each type
    n_thin = model_spec.count('S')
    n_powerlaw = model_spec.count('P')
    n_thick = model_spec.count('T')

    # Check we have enough components defined
    if n_thin > len(thin_comps):
        raise ValueError(f"Model spec requests {n_thin} thin components but only {len(thin_comps)} defined")
    if n_powerlaw > len(powerlaw_comps):
        raise ValueError(f"Model spec requests {n_powerlaw} powerlaw components but only {len(powerlaw_comps)} defined")
    if n_thick > len(thick_comps):
        raise ValueError(f"Model spec requests {n_thick} thick components but only {len(thick_comps)} defined")

    return thin_comps[:n_thin], powerlaw_comps[:n_powerlaw], thick_comps[:n_thick]


def match_components(fit_component_info, true_component_info, param_summary):
    """
    Match fitted components to true components based on parameter similarity.

    For each fitted component, find the best matching true component by comparing
    the median fitted parameters to the true parameter values.

    Args:
        fit_component_info: List of dicts with fitted component metadata
        true_component_info: List of dicts with true component metadata
        param_summary: Dict mapping param names to summary statistics (must have 'median')

    Returns:
        matches: Dict mapping fit_idx -> (true_label, true_idx)
                 e.g., {0: ('P2', 1), 1: ('P1', 0)} means fit comp 0 matches true P2
    """
    from scipy.optimize import linear_sum_assignment

    # Build cost matrix: cost[i, j] = distance between fit component i and true component j
    n_fit = len(fit_component_info)
    n_true = len(true_component_info)
    cost_matrix = np.full((n_fit, n_true), np.inf)

    for i, fit_info in enumerate(fit_component_info):
        fit_type = fit_info['type']
        fit_label = fit_info['label']

        # Extract median fitted parameters - handle tied/fixed RM case
        if fit_type == 'ThinComponent':
            fit_params = {
                'p0': param_summary[f'{fit_label}_p0']['median'],
                'psi_0': param_summary[f'{fit_label}_psi_0']['median'],
            }
            param_weights = {'p0': 1.0, 'psi_0': 0.5}  # Relative importance

            # Only include phi_rm if this component has one
            if f'{fit_label}_phi_rm' in param_summary:
                fit_params['phi_rm'] = param_summary[f'{fit_label}_phi_rm']['median']
                param_weights['phi_rm'] = 0.1

        elif fit_type == 'ThinPowerLawComponent':
            fit_params = {
                'p0': param_summary[f'{fit_label}_p0']['median'],
                'psi_0': param_summary[f'{fit_label}_psi_0']['median'],
                'beta': param_summary[f'{fit_label}_beta']['median']
            }
            param_weights = {'p0': 1.0, 'psi_0': 0.5, 'beta': 0.3}

            # Only include phi_rm if this component has one
            if f'{fit_label}_phi_rm' in param_summary:
                fit_params['phi_rm'] = param_summary[f'{fit_label}_phi_rm']['median']
                param_weights['phi_rm'] = 0.1

        elif fit_type == 'ThickComponent':
            fit_params = {
                'p0': param_summary[f'{fit_label}_p0']['median'],
                'psi_0': param_summary[f'{fit_label}_psi_0']['median'],
                'phi_peak': param_summary[f'{fit_label}_phi_peak']['median'],
                'sigma_phi': param_summary[f'{fit_label}_sigma_phi']['median'],
                'N': param_summary[f'{fit_label}_N']['median']
            }
            param_weights = {'p0': 1.0, 'psi_0': 0.5, 'phi_peak': 0.1, 'sigma_phi': 0.2, 'N': 0.2}

        # Compare to each true component
        for j, true_info in enumerate(true_component_info):
            if true_info['type'] != fit_type:
                continue  # Can only match same type

            true_params = true_info['params']

            # Calculate weighted distance
            distance = 0.0
            for param_name in fit_params.keys():
                if param_name == 'psi_0':
                    # Handle angle wrapping: circular distance on [-π/2, π/2]
                    # Map to [0, π] space and compute circular distance
                    angle1 = fit_params[param_name] % np.pi
                    angle2 = true_params[param_name] % np.pi
                    diff = np.abs(angle1 - angle2)
                    diff = np.minimum(diff, np.pi - diff)  # Circular distance on [0, π]
                else:
                    # Normalize by true value (or 1 if true is ~0)
                    scale = np.abs(true_params[param_name]) if np.abs(true_params[param_name]) > 1e-3 else 1.0
                    diff = np.abs((fit_params[param_name] - true_params[param_name]) / scale)

                distance += param_weights[param_name] * diff**2

            cost_matrix[i, j] = np.sqrt(distance)

    # Use Hungarian algorithm to find optimal matching
    if n_fit == n_true:
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        matches = {i: (true_component_info[j]['label'], j) for i, j in zip(row_ind, col_ind)}
    else:
        # Greedy matching if different number of components
        matches = {}
        used_true = set()
        for i in range(n_fit):
            # Find best match among unused true components
            best_j = None
            best_cost = np.inf
            for j in range(n_true):
                if j not in used_true and cost_matrix[i, j] < best_cost:
                    best_j = j
                    best_cost = cost_matrix[i, j]
            if best_j is not None:
                matches[i] = (true_component_info[best_j]['label'], best_j)
                used_true.add(best_j)

    return matches


# =============================================================================
# HELPER FUNCTION: GENERATE SYNTHETIC DATA WITH POWER-LAW SPECTRUM
# =============================================================================

def generate_synthetic_data(model, output_file, freq_ghz_range=(0.9, 1.7), n_channels=128,
                           spectral_index=-0.35, I_ref=1.0, snr=25,
                           oscillation=False, only_I=True):
    """
    Generate synthetic polarized data and save directly to file.

    Args:
        model: CompositeModel defining the polarized emission
        output_file: Path to save the data file
        freq_ghz_range: Frequency range in GHz (min, max)
        n_channels: Number of frequency channels
        spectral_index: Spectral index alpha (I ∝ ν^alpha)
        I_ref: Reference total intensity at 1.4 GHz
        snr: Signal-to-noise ratio for noise calculation
        oscillation: If True, add sinusoidal oscillation (3 periods over λ band)
        only_I: If True, fractional pol constant; if False, absolute pol constant
    """
    # Generate frequency grid
    freq_ghz = np.linspace(freq_ghz_range[0], freq_ghz_range[1], n_channels)
    lambda_m = 0.2998 / freq_ghz  # Wavelength in m
    lambda_sq = lambda_m**2  # Wavelength squared in m^2
    freq_hz = freq_ghz * 1e9  # Frequency in Hz

    # Generate power-law total intensity spectrum: I(ν) = I_ref * (ν / ν_ref)^alpha
    nu_ref = 1.3  # Reference frequency in GHz
    I = I_ref * (freq_ghz / nu_ref)**spectral_index

    # Compute complex polarization P = Q + iU from model (fractional)
    P = model.compute_polarization(lambda_sq)
    Q_frac = P.real
    U_frac = P.imag

    # Convert to absolute units (mJy)
    Q_abs = Q_frac * I
    U_abs = U_frac * I

    # Add oscillation if requested (3 periods over λ bandwidth)
    if oscillation:
        lambda_norm = (lambda_m - lambda_m.min()) / (lambda_m.max() - lambda_m.min())
        osc_factor = 1.0 + 0.02 * np.sin(2 * np.pi * 3 * lambda_norm)
        if only_I:
            # Fractional pol constant → absolute Q, U oscillate with I
            I = I * osc_factor
            Q_abs = Q_abs * osc_factor
            U_abs = U_abs * osc_factor
        else:
            # Absolute pol constant → fractional pol oscillates inversely with I
            I = I * osc_factor
            # Q_abs, U_abs unchanged

    # Add noise based on SNR (defined relative to peak absolute polarization)
    p_peak = np.max(np.sqrt(Q_abs**2 + U_abs**2))
    noise_level = p_peak / snr

    Q_noise = np.random.normal(0, noise_level, size=len(lambda_sq))
    U_noise = np.random.normal(0, noise_level, size=len(lambda_sq))
    I_noise = np.random.normal(0, noise_level, size=len(lambda_sq))

    I_obs = I + I_noise
    Q_obs = Q_abs + Q_noise
    U_obs = U_abs + U_noise

    # Uncertainties (homoscedastic)
    I_err = np.full_like(I, noise_level)
    Q_err = np.full_like(Q_abs, noise_level)
    U_err = np.full_like(U_abs, noise_level)

    # Save to file
    with open(output_file, 'w') as f:
        f.write("# Synthetic polarization data\n")
        f.write("# Freq[Hz]  I[mJy]  Q[mJy]  U[mJy]  err_I[mJy]  err_Q[mJy]  err_U[mJy]\n")
        for i in range(len(lambda_sq)):
            f.write(f"{freq_hz[i]:.6e}  {I_obs[i]:.6e}  {Q_obs[i]:.6e}  {U_obs[i]:.6e}  ")
            f.write(f"{I_err[i]:.6e}  {Q_err[i]:.6e}  {U_err[i]:.6e}\n")

    print(f"Saved synthetic data to: {output_file}")


# =============================================================================
# PRESENTATION SUMMARY PLOT
# =============================================================================

def plot_presentation_summary(results: 'FitResults',
                              true_model: Optional['CompositeModel'] = None,
                              n_posterior_samples: int = 100,
                              phi_range: Optional[float] = None,
                              n_phi: int = 2048,
                              weight_type: str = 'variance',
                              model_method: str = 'representative',
                              title: Optional[str] = None,
                              filename: Optional[str] = None,
                              output_dir: Optional[Path] = None,
                              show: bool = False) -> plt.Figure:
    """
    Create a 4-panel summary plot sized for landscape presentations (16:9).

    Panels:
    - Top-left: Fractional polarization p vs λ² with posterior credible bands
    - Bottom-left: EVPA ψ vs λ² with posterior credible bands
    - Top-right: Observed FDF from RM synthesis of data (linear scale)
    - Bottom-right: Intrinsic model component FDFs with posterior sample bands

    Args:
        results: FitResults object from fitting
        true_model: Optional true model for overplotting (synthetic data)
        n_posterior_samples: Number of posterior draws for bands
        phi_range: Faraday depth range (±phi_range). Auto-computed if None.
        n_phi: Number of Faraday depth grid points
        weight_type: Weighting for RM synthesis ('variance' or 'uniform')
        model_method: Method for FDF model parameters ('representative', 'median',
                     'mean', 'map', 'mode')
        title: Optional figure title
        filename: If provided, save figure to disk
        output_dir: Directory for saving
        show: If True, display plot interactively

    Returns:
        Figure object
    """
    setup_publication_quality()

    fig, axes = plt.subplots(2, 2, figsize=(13.3, 7.5))
    ax_p, ax_fdf_obs = axes[0]
    ax_psi, ax_fdf_model = axes[1]

    data = results.setup.data
    lambda_sq = data.lambda_sq
    lambda_sq_fine = np.linspace(lambda_sq.min(), lambda_sq.max(), 2000)

    # --- Sample posterior for credible bands ---
    posterior_fine = []
    if n_posterior_samples > 0:
        n_total = len(results.samples)
        indices = np.random.choice(n_total,
                                   size=min(n_posterior_samples, n_total),
                                   replace=False, p=results.weights)
        for idx in indices:
            model_sample = copy.deepcopy(results.setup.model)
            update_model_from_params(model_sample, results.samples[idx],
                                     results.setup.param_names,
                                     lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
            posterior_fine.append(model_sample.compute_polarization(lambda_sq_fine))

    # Compute percentile bands on fine grid
    if posterior_fine:
        P_abs_fine = np.array([np.abs(P) for P in posterior_fine])
        psi_fine = np.array([np.degrees(np.angle(P) / 2) for P in posterior_fine])
        p_pct = np.percentile(P_abs_fine, [2.5, 16, 50, 84, 97.5], axis=0)
        psi_pct = np.percentile(psi_fine, [2.5, 16, 50, 84, 97.5], axis=0)
    else:
        p_pct = psi_pct = None

    # Compute observed p and ψ with errors (Ricean debiased)
    p_data, p_err, chi_data, chi_err = compute_polarization_quantities(
        data.Q, data.U, data.Q_err, data.U_err,
        method='mc', n_samples=1000, debias=True
    )
    psi_data_deg = np.degrees(chi_data)

    # ========================================================================
    # Top-left: Fractional polarization p vs λ²
    # ========================================================================
    if p_pct is not None:
        ax_p.fill_between(lambda_sq_fine, p_pct[1], p_pct[3],
                          color='C0', alpha=0.25, zorder=3, label='68% CI')
        ax_p.plot(lambda_sq_fine, p_pct[2], '-', color='C0',
                  linewidth=1.5, alpha=0.8, zorder=4, label='Median')

    p_lower, p_upper = p_err
    ax_p.errorbar(lambda_sq, p_data, yerr=[p_lower, p_upper],
                  fmt='o', markerfacecolor='lightgrey', markeredgecolor='black',
                  ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                  zorder=5, label='Data')

    if true_model is not None:
        P_true_fine = true_model.compute_polarization(lambda_sq_fine)
        ax_p.plot(lambda_sq_fine, np.abs(P_true_fine), '-', color='black',
                  linewidth=2.0, label='True model', zorder=6)

    ax_p.set_xlabel(r'$\lambda^2$ (m$^2$)')
    ax_p.set_ylabel('Fractional Polarisation')
    ax_p.legend(loc='best', framealpha=0.9, fontsize=8)
    ax_p.grid(True, alpha=0.3)
    add_ticks_all_sides(ax_p)

    # ========================================================================
    # Bottom-left: EVPA ψ vs λ²
    # ========================================================================
    if psi_pct is not None:
        ax_psi.fill_between(lambda_sq_fine, psi_pct[1], psi_pct[3],
                            color='C0', alpha=0.25, zorder=3)
        ax_psi.plot(lambda_sq_fine, psi_pct[2], '-', color='C0',
                    linewidth=1.5, alpha=0.8, zorder=4)

    chi_lower, chi_upper = chi_err
    ax_psi.errorbar(lambda_sq, psi_data_deg,
                    yerr=[np.degrees(chi_lower), np.degrees(chi_upper)],
                    fmt='o', markerfacecolor='lightgrey', markeredgecolor='black',
                    ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                    zorder=5)

    if true_model is not None:
        P_true_fine = true_model.compute_polarization(lambda_sq_fine)
        psi_true = np.degrees(np.angle(P_true_fine) / 2)
        ax_psi.plot(lambda_sq_fine, psi_true, '-', color='black',
                    linewidth=2.0, zorder=6)

    ax_psi.set_xlabel(r'$\lambda^2$ (m$^2$)')
    ax_psi.set_ylabel('EVPA (deg)')
    ax_psi.grid(True, alpha=0.3)
    ax_psi.set_ylim(-90, 90)
    add_ticks_all_sides(ax_psi)

    # ========================================================================
    # FDF setup (shared by both right panels)
    # ========================================================================
    if model_method == 'representative':
        rep_params = results.get_representative_sample(percentile=1.0, method='median')
        fdf_model = copy.deepcopy(results.setup.model)
        update_model_from_params(fdf_model, rep_params, results.setup.param_names,
                                 lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
    else:
        fdf_model = results.get_best_fit_model(method=model_method)

    # Auto-compute phi_range from RMTF and fitted RMs
    if phi_range is None:
        _, rmtf_props = compute_rmtf(lambda_sq, np.array([0.0]),
                                      weight_type=weight_type,
                                      Q_err=data.Q_err, U_err=data.U_err)
        fwhm = rmtf_props['expected_fwhm_nominal']
        rm_vals = []
        for comp in fdf_model.components:
            if hasattr(comp, 'phi_rm'):
                rm_vals.append(comp.phi_rm)
            elif hasattr(comp, 'phi_peak'):
                rm_vals.append(comp.phi_peak)
        if rm_vals:
            rm_extent = max(abs(max(rm_vals) + 2*fwhm), abs(min(rm_vals) - 2*fwhm))
        else:
            rm_extent = 0.0
        phi_range = max(10.0 * fwhm, rm_extent, 100.0)

    phi_grid = np.linspace(-phi_range, phi_range, n_phi)
    phi_grid_wide = np.linspace(-5*phi_range, 5*phi_range, 5*n_phi)
    n_start = 2 * n_phi
    n_end = 3 * n_phi

    # Observed FDF from data
    observed_fdf = rm_synthesis(data, phi_grid_wide, weight_type=weight_type)[n_start:n_end]

    # ========================================================================
    # Top-right: Observed FDF (linear scale)
    # ========================================================================
    ax_fdf_obs.plot(phi_grid, np.abs(observed_fdf), color='black', linewidth=2.0,
                    label='|F(φ)|', zorder=10)
    ax_fdf_obs.plot(phi_grid, np.real(observed_fdf), 'b--', linewidth=1.0, alpha=0.7,
                    label='Re')
    ax_fdf_obs.plot(phi_grid, np.imag(observed_fdf), 'r--', linewidth=1.0, alpha=0.7,
                    label='Im')
    ax_fdf_obs.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax_fdf_obs.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax_fdf_obs.set_xlabel(r'Faraday Depth $\phi$ (rad/m$^2$)')
    ax_fdf_obs.set_ylabel('F(φ) (fractional)')
    ax_fdf_obs.legend(loc='best', framealpha=0.9, fontsize=8)
    ax_fdf_obs.grid(True, alpha=0.3)
    add_ticks_all_sides(ax_fdf_obs)

    # ========================================================================
    # Bottom-right: Intrinsic model component FDFs
    # ========================================================================
    # Posterior sample bands (intrinsic FDFs per component)
    if n_posterior_samples > 0:
        weights_norm = results.weights / np.sum(results.weights)
        n_fdf_samples = min(100, len(results.samples))
        fdf_indices = np.random.choice(len(results.samples), size=n_fdf_samples,
                                        replace=False, p=weights_norm)
        for idx in fdf_indices:
            model_s = copy.deepcopy(results.setup.model)
            update_model_from_params(model_s, results.samples[idx],
                                      results.setup.param_names,
                                      lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
            for comp in model_s.components:
                comp_fdf_s = comp.compute_fdf(phi_grid)
                ax_fdf_model.plot(phi_grid, np.abs(comp_fdf_s), 'C0-',
                                  linewidth=0.5, alpha=0.15, zorder=1)

    # Best-fit intrinsic component FDFs
    max_amp = 0.0
    for comp in fdf_model.components:
        comp_fdf = comp.compute_fdf(phi_grid)
        comp_amp = np.abs(comp_fdf)
        max_amp = max(max_amp, np.max(comp_amp))
        ax_fdf_model.plot(phi_grid, comp_amp, 'k-', linewidth=2.0,
                          alpha=1.0, zorder=12)
        if isinstance(comp, (ThinComponent, ThinPowerLawComponent)):
            ax_fdf_model.plot(comp.phi_rm, np.max(comp_amp),
                              marker='^', markersize=10, color='black', zorder=15)

    ax_fdf_model.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax_fdf_model.set_yscale('log')
    if max_amp > 0:
        ax_fdf_model.set_ylim(1e-5, min(10.0 * max_amp, 1.0))
    else:
        ax_fdf_model.set_ylim(1e-5, 1e-3)
    ax_fdf_model.set_xlabel(r'Faraday Depth $\phi$ (rad/m$^2$)')
    ax_fdf_model.set_ylabel('|F(φ)| (fractional)')
    ax_fdf_model.grid(True, alpha=0.3, which='both')
    add_ticks_all_sides(ax_fdf_model)

    # ========================================================================
    # Final layout
    # ========================================================================
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 0.96] if title else [0, 0, 1, 1])

    if filename:
        out = Path(output_dir) if output_dir else Path.cwd()
        out.mkdir(parents=True, exist_ok=True)
        filepath = out / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved presentation plot: {filepath}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


# =============================================================================
# PAPER EXAMPLE PLOT
# =============================================================================

def plot_paper_example(results: 'FitResults',
                       true_model: Optional['CompositeModel'] = None,
                       n_posterior_samples: int = 100,
                       phi_range: Optional[float] = None,
                       n_phi: int = 2048,
                       rmsynth_peak_threshold: str = '-9',
                       rmsynth_finer_threshold: str = '-4',
                       filename: Optional[str] = None,
                       output_dir: Optional[Path] = None,
                       show: bool = False) -> plt.Figure:
    """
    Three-panel publication figure: p(λ²), ψ(λ²), and the Faraday depth spectrum.

    The FDF panel shows posterior sample component FDFs (C0), the representative
    model (black), the true FDF per component (hard black, if true_model provided),
    and the RM-CLEANed FDF linearly mapped onto the log axis (grey) as a guide.

    RM synthesis and CLEAN are run via the rmsynth1d/rmclean1d CLI tools on a
    temporary file written from the data attached to results.

    Args:
        results: FitResults object from fitting
        true_model: Optional true model for overplotting
        n_posterior_samples: Number of posterior draws for bands
        phi_range: Faraday depth range (±phi_range). Auto-computed if None.
        n_phi: Number of Faraday depth grid points
        rmsynth_peak_threshold: rmclean1d -c argument (log10 peak fraction)
        rmsynth_finer_threshold: rmclean1d -w argument (log10 finer mask fraction)
        filename: If provided, save figure (PDF)
        output_dir: Directory for saving
        show: If True, display interactively
    """
    setup_publication_quality()

    data       = results.setup.data
    lambda_sq  = data.lambda_sq
    lsq_ref    = getattr(results.setup, 'lambda_sq_ref', None)
    lambda_sq_fine = np.linspace(lambda_sq.min(), lambda_sq.max(), 2000)

    # =========================================================================
    # RM synthesis and CLEAN on the simulated data
    # =========================================================================
    freq_hz = C_LIGHT / np.sqrt(lambda_sq)

    I_vals   = data.I   if data.I   is not None else np.ones_like(freq_hz)
    I_err    = data.I_err if data.I_err is not None else np.full_like(freq_hz, 1e-3)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        rmsynth_txt = tmpdir / 'sim_rmsynth.txt'
        np.savetxt(rmsynth_txt,
                   np.column_stack([freq_hz, I_vals, data.Q, data.U,
                                    I_err, data.Q_err, data.U_err]))

        subprocess.run([
            'rmsynth1d', str(rmsynth_txt),
            '-S', '-o', '4', '-l', '2500', '--super-resolution'
        ], check=False, capture_output=True)

        subprocess.run([
            'rmclean1d', str(rmsynth_txt),
            '-S', '-c', rmsynth_peak_threshold, '-w', rmsynth_finer_threshold
        ], check=False, capture_output=True)

        clean_file = tmpdir / 'sim_rmsynth_FDFclean.dat'
        cleaned_phi = cleaned_amp = None
        if clean_file.exists():
            fdf_dat     = np.loadtxt(clean_file)
            cleaned_phi = fdf_dat[:, 0]
            cleaned_amp = np.sqrt(fdf_dat[:, 1]**2 + fdf_dat[:, 2]**2)
        else:
            print('  WARNING: rmclean output not found — CLEANed FDF will not be shown')

    # =========================================================================
    # Posterior samples for p, psi, and FDF panels
    # =========================================================================
    posterior_fine = []
    posterior_fdf  = []
    if n_posterior_samples > 0:
        n_total      = len(results.samples)
        weights_norm = results.weights / np.sum(results.weights)
        indices      = np.random.choice(n_total,
                                        size=min(n_posterior_samples, n_total),
                                        replace=False, p=weights_norm)
        for idx in indices:
            m = copy.deepcopy(results.setup.model)
            update_model_from_params(m, results.samples[idx],
                                     results.setup.param_names,
                                     lambda_sq_ref=lsq_ref)
            posterior_fine.append(m.compute_polarization(lambda_sq_fine))
            posterior_fdf.append(copy.deepcopy(m))

    # Auto phi_range using true model if available, else first posterior sample
    if phi_range is None:
        ref_model = true_model if true_model is not None else (
            posterior_fdf[0] if posterior_fdf else None)
        _, rmtf_props = compute_rmtf(lambda_sq, np.array([0.0]),
                                     Q_err=data.Q_err, U_err=data.U_err)
        fwhm = rmtf_props['expected_fwhm_nominal']
        rm_vals = []
        if ref_model is not None:
            for comp in ref_model.components:
                if hasattr(comp, 'phi_rm'):
                    rm_vals.append(comp.phi_rm)
                elif hasattr(comp, 'phi_peak'):
                    rm_vals.append(comp.phi_peak)
        rm_extent = max(abs(max(rm_vals) + 2*fwhm), abs(min(rm_vals) - 2*fwhm)) if rm_vals else 0.0
        phi_range = max(10.0 * fwhm, rm_extent, 100.0)

    phi_grid = np.linspace(-phi_range, phi_range, n_phi)

    # =========================================================================
    # Observed p and psi with errors
    # =========================================================================
    p_data, p_err, chi_data, chi_err = compute_polarization_quantities(
        data.Q, data.U, data.Q_err, data.U_err,
        method='mc', n_samples=1000, debias=True
    )
    psi_data_deg = np.degrees(chi_data)

    # =========================================================================
    # Figure layout: 3 vertical panels
    # =========================================================================
    fig, (ax_p, ax_psi, ax_fdf) = plt.subplots(3, 1, figsize=(6.0, 9.0),
                                                 gridspec_kw={'hspace': 0.35})

    # -------------------------------------------------------------------------
    # Panel 1: p(lambda^2) in percent
    # -------------------------------------------------------------------------
    for P_samp in posterior_fine:
        ax_p.plot(lambda_sq_fine, np.abs(P_samp) * 100.0,
                  '-', color='C0', linewidth=1.0, alpha=0.15, zorder=5)

    p_lower, p_upper = p_err
    ax_p.errorbar(lambda_sq, p_data * 100.0,
                  yerr=[p_lower * 100.0, p_upper * 100.0],
                  fmt='o', markerfacecolor='lightgrey', markeredgecolor='darkgrey',
                  ecolor='black', alpha=1.0, markersize=4, markeredgewidth=0.5,
                  zorder=1)

    if true_model is not None:
        ax_p.plot(lambda_sq_fine,
                  np.abs(true_model.compute_polarization(lambda_sq_fine)) * 100.0,
                  '-', color='black', linewidth=0.75, zorder=10)

    ax_p.set_xlabel(r'$\lambda^2\;(\mathrm{m}^2)$')
    ax_p.set_ylabel(r'$p(\lambda^2)$ (%)')
    add_ticks_all_sides(ax_p)

    # -------------------------------------------------------------------------
    # Panel 2: psi(lambda^2)
    # -------------------------------------------------------------------------
    for P_samp in posterior_fine:
        ax_psi.plot(lambda_sq_fine, np.degrees(np.angle(P_samp) / 2),
                    '-', color='C0', linewidth=1.0, alpha=0.15, zorder=5)

    chi_lower, chi_upper = chi_err
    ax_psi.errorbar(lambda_sq, psi_data_deg,
                    yerr=[np.degrees(chi_lower), np.degrees(chi_upper)],
                    fmt='o', markerfacecolor='lightgrey', markeredgecolor='darkgrey',
                    ecolor='black', alpha=1.0, markersize=4, markeredgewidth=0.5,
                    zorder=1)

    if true_model is not None:
        P_true_fine = true_model.compute_polarization(lambda_sq_fine)
        ax_psi.plot(lambda_sq_fine, np.degrees(np.angle(P_true_fine) / 2),
                    '-', color='black', linewidth=0.75, zorder=10)

    ax_psi.set_xlabel(r'$\lambda^2\;(\mathrm{m}^2)$')
    ax_psi.set_ylabel(r'$\psi(\lambda^2)$ (deg)')
    add_ticks_all_sides(ax_psi)

    # -------------------------------------------------------------------------
    # Panel 3: FDF -- log scale, samples in C0, true in black, CLEANed in grey
    # -------------------------------------------------------------------------
    ax_fdf.set_yscale('log')

    # y-limits from true model if available, else envelope of posterior samples
    if true_model is not None:
        true_amps = [np.abs(c.compute_fdf(phi_grid)) for c in true_model.components]
        max_amp = max(a.max() for a in true_amps) if true_amps else 1e-3
    elif posterior_fdf:
        samp_amps = [np.abs(c.compute_fdf(phi_grid))
                     for m in posterior_fdf[:10] for c in m.components]
        max_amp = max(a.max() for a in samp_amps) if samp_amps else 1e-3
    else:
        max_amp = 1e-3

    fdf_ymin = max(1e-7, max_amp * 1e-4) * 100.0
    fdf_ymax = min(100.0, max_amp * 10.0) * 100.0
    ax_fdf.set_ylim(fdf_ymin, fdf_ymax)

    _log_ymin = np.log10(fdf_ymin)
    _log_ymax = np.log10(fdf_ymax)

    def _draw_arrow(ax, phi_rm, amp, color, alpha, lw, zorder):
        """Draw a narrow upward arrow at phi_rm to mark a point-source component."""
        ax.annotate('',
                    xy=(phi_rm, amp),
                    xytext=(phi_rm, amp * 0.3),
                    arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                   mutation_scale=3),
                    alpha=alpha, zorder=zorder)

    # Posterior FDF samples
    for m in posterior_fdf:
        for comp in m.components:
            comp_amp = np.abs(comp.compute_fdf(phi_grid)) * 100.0
            ax_fdf.plot(phi_grid, comp_amp,
                        'C0-', linewidth=0.5, alpha=0.15, zorder=1)
            if isinstance(comp, (ThinComponent, ThinPowerLawComponent)):
                _draw_arrow(ax_fdf, comp.phi_rm, comp_amp.max(),
                            color='C0', alpha=0.15, lw=1.0, zorder=2)

    # True FDF per component
    if true_model is not None:
        for comp in true_model.components:
            comp_amp = np.abs(comp.compute_fdf(phi_grid)) * 100.0
            ax_fdf.plot(phi_grid, comp_amp, 'k-', linewidth=2.0, zorder=12)
            if isinstance(comp, (ThinComponent, ThinPowerLawComponent)):
                _draw_arrow(ax_fdf, comp.phi_rm, comp_amp.max(),
                            color='black', alpha=1.0, lw=0.75, zorder=13)
                ax_fdf.plot(comp.phi_rm, comp_amp.max(), marker='^',
                            markersize=12, color='black', zorder=14)

    # CLEANed FDF linearly mapped onto log axis
    if cleaned_phi is not None:
        mask_c   = np.abs(cleaned_phi) <= phi_range
        amp_peak = cleaned_amp[mask_c].max()
        if amp_peak > 0:
            amp_norm   = cleaned_amp[mask_c] / amp_peak
            amp_mapped = 10 ** (_log_ymin + amp_norm * (_log_ymax - _log_ymin) * 0.90)
            ax_fdf.plot(cleaned_phi[mask_c], amp_mapped,
                        color='grey', alpha=1.0, linewidth=1.5, zorder=20)

    ax_fdf.axvline(0, color='gray', linestyle=':', linewidth=0.6, zorder=0)
    ax_fdf.set_xlim(-150, 150)
    ax_fdf.set_xlabel(r'$\phi_f$ (rad m$^{-2}$)')
    ax_fdf.set_ylabel(r'$|F(\phi_f)|$ (%) (rad m$^{-2}$)$^{-1}$')
    add_ticks_all_sides(ax_fdf)

    # =========================================================================
    # Save
    # =========================================================================
    if filename:
        out = Path(output_dir) if output_dir else Path.cwd()
        out.mkdir(parents=True, exist_ok=True)
        filepath = out / filename
        fig.savefig(filepath, format='pdf', bbox_inches='tight')
        print(f"Saved paper example plot: {filepath}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


# =============================================================================
# CONFIGURABLE TEST SCENARIO
# =============================================================================

def run_test_scenario():
    """
    Flexible test scenario - easily add/remove components and modify settings.

    This function is partitioned into clear sections for easy modification:
    - TRUE MODEL CONFIGURATION
    - DATA GENERATION
    - FIT MODEL SETUP
    - PRIOR SPECIFICATION
    - FITTING CONFIGURATION
    - PLOTTING
    """

    print("\n" + "="*80)
    print("FLEXIBLE TEST SCENARIO")
    print("="*80)

    # =========================================================================
    # SECTION 1: TRUE MODEL CONFIGURATION
    # =========================================================================
    # Build the "true" model from configured components at top of file

    print("\n[1] Configuring TRUE MODEL from component arrays...")

    true_model, true_component_info = build_model_from_config(
        THIN_COMPONENTS, POWERLAW_COMPONENTS, THICK_COMPONENTS
    )

    print("\nTrue model components:")
    for info in true_component_info:
        print(f"  {info['label']} ({info['type']}): {info['params']}")
    print(f"Total: {len(true_model.components)} components")

    # =========================================================================
    # SECTION 2: DATA GENERATION AND LOADING
    # =========================================================================

    print("\n[2] Generating SYNTHETIC DATA and saving to file...")

    # Data generation parameters
    FREQ_RANGE = (0.9, 1.65)      # Frequency range in GHz
    N_CHANNELS = 250               # Number of frequency channels
    SPECTRAL_INDEX = -0.35         # Spectral index (I ∝ ν^α)
    SNR = 10.0                      # Signal-to-noise ratio
    OSCILLATION = False           # Add sinusoidal oscillation (3 periods over λ²)
    ONLY_I = False                 # True: frac pol constant, False: abs pol constant

    data_file = output_dir / "synthetic_data.txt"

    generate_synthetic_data(
        true_model,
        output_file=data_file,
        freq_ghz_range=FREQ_RANGE,
        n_channels=N_CHANNELS,
        spectral_index=SPECTRAL_INDEX,
        snr=SNR,
        oscillation=OSCILLATION,
        only_I=ONLY_I
    )

    print(f"Generated {N_CHANNELS} channels from {FREQ_RANGE[0]}-{FREQ_RANGE[1]} GHz, SNR={SNR}")

    # Load data with total intensity fit
    print("\n[2B] Loading data and fitting TOTAL INTENSITY...")
    central_freq_ghz = (FREQ_RANGE[0] + FREQ_RANGE[1]) / 2
    lambda_sq_ref = (0.2998 / central_freq_ghz) ** 2
    data, I_fitter = DataLoader.load_physical_with_fit(
        str(data_file),
        I_model='power-law',
        poly_degree=4,
        central_freq=central_freq_ghz,
        flagging=None,
        load_stokes_v=False
    )

    print(f"Power-law fit: I₀ = {I_fitter.fit_params[0]:.4f} mJy, α = {I_fitter.fit_params[1]:.4f}")
    print(f"True values: I₀ = 1.0 mJy, α = {SPECTRAL_INDEX}")

    # Plot data diagnostic with I model fit
    if not QUICK_MODE:
        plot_data_diagnostic(
            data,
            I_fitter=I_fitter,
            title=None,
            filename="test_data.png",
            output_dir=output_dir
        )

    # =========================================================================
    # SECTION 3: FIT MODEL SETUP
    # =========================================================================
    # Build fit model from MODEL_TO_FIT specification

    print("\n[3] Setting up FIT MODEL...")

    # Parse model specification
    thin_to_fit, powerlaw_to_fit, thick_to_fit = parse_model_spec(
        MODEL_TO_FIT, THIN_COMPONENTS, POWERLAW_COMPONENTS, THICK_COMPONENTS
    )

    # Build fit model (with zero initial parameters)
    fit_model, fit_component_info = build_model_from_config(
        thin_to_fit, powerlaw_to_fit, thick_to_fit
    )

    # Set all parameters to zero for fitting
    for comp in fit_model.components:
        if isinstance(comp, ThinComponent):
            comp.a, comp.b, comp.phi_rm = 0.0, 0.0, 0.0
        elif isinstance(comp, ThinPowerLawComponent):
            comp.a, comp.b, comp.phi_rm = 0.0, 0.0, 0.0
            # beta and lambda_sq_0 keep their values
        elif isinstance(comp, ThickComponent):
            comp.a, comp.b = 0.0, 0.0
            comp.phi_peak = 0.0
            # sigma_phi and N keep their values

    if MODEL_TO_FIT is None:
        print(f"Fit model: AUTO (same as truth) - {len(fit_model.components)} components")
    else:
        print(f"Fit model: '{MODEL_TO_FIT}' - {len(fit_model.components)} components")

    print("Fit model components:")
    for info in fit_component_info:
        print(f"  {info['label']} ({info['type']})")

    # =========================================================================
    # SECTION 4: PRIOR SPECIFICATION
    # =========================================================================
    # Define priors for each component
    # Choose between Cartesian (a, b) or polar_conversion (p0, psi_0) styles

    print("\n[4] Specifying PRIORS...")

    # --- PRIOR STYLE OPTIONS ---
    USE_POLAR_CONVERSION = True   # True: use (p0, psi_0), False: use (a, b)
    USE_LOG_UNIFORM_P0 = False     # Only applies if USE_POLAR_CONVERSION=True

    # --- DEFINE BASE PRIORS (shared across components of same type) ---

    # Base prior for ALL Thin components
    if USE_POLAR_CONVERSION:
        BASE_THIN_PRIOR = {
            'type': 'ThinComponent',
            'phi_rm': ('uniform', -200.0, 200.0),
            'polar_conversion': {
                'p0_bounds': (0.00, 0.5, 'radial_uniform'),
                'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked')
            }
        }
    else:
        BASE_THIN_PRIOR = {
            'type': 'ThinComponent',
            'a': ('uniform', -0.1, 0.1),
            'b': ('uniform', -0.1, 0.1)
        }

    # Base prior for ALL Thick components
    if USE_POLAR_CONVERSION:
        BASE_THICK_PRIOR = {
            'type': 'ThickComponent',
            'polar_conversion': {
                'p0_bounds': (0.0, 0.5, 'radial_uniform'),
                'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked')
            },
            'sigma_phi': ('log_uniform', 1e-1, 50.0),
            'phi_peak': ('uniform', -200.0, 200.0),
            'N': ('log_uniform', 2.0, 50.0),
            'effective_width_mode': False  # Set True to use effective width interpretation
        }
    else:
        BASE_THICK_PRIOR = {
            'type': 'ThickComponent',
            'a': ('uniform', -0.1, 0.1),
            'b': ('uniform', -0.1, 0.1),
            'sigma_phi': ('log_uniform', 1e-1, 50.0),
            'phi_peak': ('uniform', -200.0, 200.0),
            'N': ('log_uniform', 2.0, 50.0),
            'effective_width_mode': False  # Set True to use effective width interpretation
        }

    # Base prior for ALL Thin Power-Law components
    if USE_POLAR_CONVERSION:
        BASE_PL_PRIOR = {
            'type': 'ThinPowerLawComponent',
            'phi_rm': ('uniform', -200.0, 200.0),
            'polar_conversion': {
                'p0_bounds': (0.00, 0.5, 'radial_uniform'),
                'psi_0_bounds': (-np.pi/2, np.pi/2, 'peaked')
            },
            'beta': ('uniform', -5.0, 5.0),
        }
    else:
        BASE_PL_PRIOR = {
            'type': 'ThinPowerLawComponent',
            'a': ('uniform', -0.1, 0.1),
            'b': ('uniform', -0.1, 0.1),
            'beta': ('uniform', -3.0, 3.0),
        }

    # --- APPLY PRIORS TO EACH COMPONENT ---

    priors = {'tied_phi_rm': ('uniform', -200, 200)}

    # Loop through fit model components and assign priors
    for i, comp in enumerate(fit_model.components):
        comp_name = comp.name if comp.name else f'component_{i}'

        if isinstance(comp, ThinComponent):
            # Copy base thin prior
            priors[comp_name] = BASE_THIN_PRIOR.copy()

        elif isinstance(comp, ThickComponent):
            # Copy base thick prior
            priors[comp_name] = BASE_THICK_PRIOR.copy()

        elif isinstance(comp, ThinPowerLawComponent):
            # Copy base thin power-law prior
            priors[comp_name] = BASE_PL_PRIOR.copy()
        else:
            sys.exit(f"ERROR: Unknown component type for {comp_name}")

    # --- OPTIONAL: SYSTEMATIC NOISE ---
    # priors['systematic_noise'] = {
    #     'sigma_sys_Q': ('uniform', 1e-4, 3e-3),
    #     'sigma_sys_U': ('uniform', 1e-4, 3e-3)
    # }

    print(f"Priors specified for {len([k for k in priors.keys() if k != 'systematic_noise'])} components")
    print(f"Using {'polar conversion (p0, ψ₀)' if USE_POLAR_CONVERSION else 'Cartesian (a, b)'} priors")

    # =========================================================================
    # SECTION 5: FITTING CONFIGURATION
    # =========================================================================

    print("\n[5] Configuring BAYESIAN FIT...")

    # --- FIT SETTINGS ---
    ENFORCE_ORDERING = True       # Enforce RM ordering between components
    P0_PHYSICAL_MAX = 0.7         # Maximum allowed polarization fraction
    N_LIVE = 900                 # Number of live points for nested sampling
    SAMPLER = 'dynamic'             # Sampling method: 'rwalk' or 'rslice'
    N_CORES = 15
    WALKS = 150
    SLICES = 3 + 9
    FIX_RM = False

    setup = FitSetup(
        model=fit_model,
        priors=priors,
        data=data,
        enforce_ordering=ENFORCE_ORDERING,
        p0_physical_max=P0_PHYSICAL_MAX,
        tie_phi_rm=FIX_RM,  # Tie RM across components (if True, use 'tied_phi_rm' prior)
        lambda_sq_ref=lambda_sq_ref
    )

    print(f"\nFit configuration:")
    print(f"  Parameters: {setup.ndim}")
    print(f"  Ordering: {'Enabled' if ENFORCE_ORDERING else 'Disabled'}")
    print(f"  p0_max: {P0_PHYSICAL_MAX}")

    print("\n[6] Running NESTED SAMPLING...")
    if SAMPLER == 'dynamic':
        FIT_KWARGS = {
            'sampler': {'bound': 'multi', 'sample': 'rwalk', 'walks': WALKS},#, 'slices': SLICES},
            'run': {'dlogz_init': 0.01, 'nlive_init': int(N_LIVE / 2.0), 'nlive_batch': int(N_LIVE/6.0), 'n_effective': 10_000, 'maxiter': None, 'maxcall': None}
        }
    else:
        FIT_KWARGS = {
            'sampler': {'nlive': N_LIVE, 'bound': 'multi', 'sample': 'rwalk', 'walks': WALKS},
            'run': {'dlogz': 0.1, 'maxiter': None, 'maxcall': None}
        }
    fitter = FaradayFitter(setup, dynesty_kwargs=FIT_KWARGS, n_cores=N_CORES, sampler=SAMPLER)
    results = fitter.fit()

    # Save results (pickle immediately; JSON after post-processing)
    results.save(output_dir / "test_results.pkl")

    # =========================================================================
    # SECTION 6: PLOTTING
    # =========================================================================

    print("\n[7] Creating PLOTS...")

    if QUICK_MODE:
        # Quick mode: only the presentation summary
        plot_presentation_summary(
            results,
            true_model=true_model,
            title=None,
            filename="test_presentation.png",
            output_dir=output_dir
        )
        plot_paper_example(
            results,
            true_model=true_model,
            filename="test_paper_example.pdf",
            output_dir=output_dir
        )
        results.save_json(output_dir / "test_results.json")

    else:
        # --- Determine if we should plot truths or medians ---
        # If true model and fit model have the same structure, plot matched truths
        models_match = (len(true_component_info) == len(fit_component_info) and
                       all(t['type'] == f['type'] for t, f in zip(true_component_info, fit_component_info)))

        if models_match:
            print("\nTrue and fit models have matching structure")

            # Match components based on fitted parameter values
            matches = match_components(fit_component_info, true_component_info, results.param_summary)

            print("Component matching (fit → truth):")
            for fit_idx, (true_label, true_idx) in matches.items():
                fit_label = fit_component_info[fit_idx]['label']
                true_params = true_component_info[true_idx]['params']
                print(f"  {fit_label} → {true_label}: p0={true_params['p0']:.4f}, psi_0={np.degrees(true_params['psi_0']):.1f}°", end='')
                if 'beta' in true_params:
                    print(f", beta={true_params['beta']:.2f}")
                elif 'phi_peak' in true_params:
                    print(f", phi_peak={true_params['phi_peak']:.1f}")
                else:
                    print(f", phi_rm={true_params['phi_rm']:.1f}")

            plot_method = 'mode'  # Show truths

            # Build truth arrays in the order of fitted components (using matches)
            truths_cartesian = []
            truths_derived = []

            # Get the shared phi_rm value if RM is fixed (from the first matched component)
            if FIX_RM and 0 in matches:
                shared_phi_rm = true_component_info[matches[0][1]]['params']['phi_rm']
            else:
                shared_phi_rm = None

            for fit_idx in range(len(fit_component_info)):
                # Get corresponding true component via matching
                if fit_idx in matches:
                    true_label, true_idx = matches[fit_idx]
                    true_info = true_component_info[true_idx]
                    params = true_info['params']
                    comp_type = true_info['type']  # Use TRUE component type, not fit
                else:
                    # No match found, use None (will not plot truth line)
                    params = None
                    comp_type = fit_component_info[fit_idx]['type']

                if params is not None:
                    if comp_type == 'ThinComponent':
                        truths_cartesian.extend([params['a'], params['b']])
                        truths_derived.extend([params['p0'], np.degrees(params['psi_0'])])
                        # Only first component adds phi_rm if FIX_RM
                        if not FIX_RM or fit_idx == 0:
                            truths_cartesian.append(params['phi_rm'])
                            truths_derived.append(params['phi_rm'])

                    elif comp_type == 'ThinPowerLawComponent':
                        truths_cartesian.extend([params['a'], params['b']])
                        truths_derived.extend([params['p0'], np.degrees(params['psi_0'])])
                        # Only first component adds phi_rm if FIX_RM
                        if not FIX_RM or fit_idx == 0:
                            truths_cartesian.append(params['phi_rm'])
                            truths_derived.append(params['phi_rm'])
                        # Always add beta
                        truths_cartesian.append(params['beta'])
                        truths_derived.append(params['beta'])

                    elif comp_type == 'ThickComponent':
                        truths_cartesian.extend([params['a'], params['b'], params['phi_peak'],
                                                params['sigma_phi'], params['N']])
                        truths_derived.extend([params['p0'], np.degrees(params['psi_0']), params['phi_peak'],
                                              params['sigma_phi'], params['N']])
                else:
                    # Fill with None for unmatched components
                    if comp_type == 'ThinComponent':
                        truths_cartesian.extend([None, None])
                        truths_derived.extend([None, None])
                        if not FIX_RM or fit_idx == 0:
                            truths_cartesian.append(None)
                            truths_derived.append(None)
                    elif comp_type == 'ThinPowerLawComponent':
                        truths_cartesian.extend([None, None])
                        truths_derived.extend([None, None])
                        if not FIX_RM or fit_idx == 0:
                            truths_cartesian.append(None)
                            truths_derived.append(None)
                        truths_cartesian.append(None)  # beta
                        truths_derived.append(None)    # beta
                    elif comp_type == 'ThickComponent':
                        truths_cartesian.extend([None, None, None, None, None])
                        truths_derived.extend([None, None, None, None, None])
        else:
            print("\nTrue and fit models have different structure - will plot medians")
            plot_method = 'median'
            truths_cartesian = None
            truths_derived = None

        # --- Fit results plot ---
        plot_fit_results(
            results,
            I_fitter=I_fitter,
            title=None,
            filename="test_fit.png",
            output_dir=output_dir,
            plot_fdf_posterior_samples=True
        )

        # --- Cartesian corner plot ---
        plot_corner(
            samples=results.samples,
            param_names=results.setup.param_names,
            param_summary=results.param_summary,
            weights=results.weights,
            setup=results.setup,
            method=plot_method,
            truths=truths_cartesian,
            title=None,
            filename="test_corner_cartesian.png",
            output_dir=output_dir
        )

        # --- Post-processing: mode detection and derived parameters ---
        pconfig = ProcessingConfig(
            filename_prefix='test',
            detect_multiple_modes=True,
            kde_bandwidth='scott',
            plot_diagnostics=True,
            plot_dir=output_dir / "test_modes/",
            verbose=False,
        )
        results = process_posterior_modes(results, pconfig)
        results.save_json(output_dir / "test_results.json")
        results.print_summary()

        # --- Derived (polar) corner plot ---
        plot_corner(
            samples=results.samples_processed,
            param_names=results.param_names_processed,
            param_summary=results.param_summary,
            processed_parameters=results.processed_parameters,
            weights=results.weights,
            setup=results.setup,
            method=plot_method,
            title=None,
            filename="test_corner_polar.png",
            truths=truths_derived,
            output_dir=output_dir
        )

        # --- Convergence diagnostic ---
        plot_convergence(
            results,
            title=None,
            filename="test_convergence.png",
            output_dir=output_dir
        )

        # --- FDF diagnostic ---
        if not any(isinstance(c, ThinPowerLawComponent) for c in fit_model.components):
            plot_fdf_diagnostic(
                results,
                title=None,
                model_method='representative',
                filename="test_fdf.png",
                output_dir=output_dir,
                plot_posterior_samples=True
            )

        # --- Presentation summary (landscape 4-panel) ---
        plot_presentation_summary(
            results,
            true_model=true_model,
            title=None,
            filename="test_presentation.png",
            output_dir=output_dir,
            phi_range=40
        )

        # --- Paper example (3-panel publication figure) ---
        plot_paper_example(
            results,
            true_model=true_model,
            filename="test_paper_example.pdf",
            output_dir=output_dir
        )

    print(f"\n{'='*80}")
    print("TEST COMPLETE!")
    print(f"{'='*80}")
    print(f"All outputs saved to: {output_dir}/")
    print(f"{'='*80}\n")

    return results


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":

    # Run the flexible test scenario
    results = run_test_scenario()
