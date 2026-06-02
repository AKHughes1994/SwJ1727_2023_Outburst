# faraday_plot.py
"""
Publication-quality plotting functions for Faraday rotation analysis.

This module provides visualization tools with consistent formatting:
- Thick axes, major/minor ticks on all sides
- Proper font sizes and spacing
- Clean, professional appearance
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path
from typing import Optional, Tuple, List, Dict
import corner
from dynesty import utils as dyfunc
from faraday_data import PolarizationData, IFitter
from faraday_utils import (compute_rmtf, rm_synthesis, compute_polarization_quantities,
                           ThinComponent, ThinPowerLawComponent, ThickComponent, C_LIGHT)
from faraday_model import update_model_from_params


# Publication quality matplotlib settings
def setup_publication_quality():
    """
    Configure matplotlib for publication-quality plots.
    
    Sets consistent styling across all plots:
    - Thicker axes and ticks
    - Proper tick sizes (major and minor)
    - Ticks on all four sides (top, bottom, left, right)
    - Reasonable font sizes
    - Serif font family for professional appearance
    """
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
        'mathtext.fontset': 'dejavuserif',
        'axes.linewidth': 1.5,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'xtick.minor.width': 1.0,
        'ytick.minor.width': 1.0,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.minor.size': 3,
        'ytick.minor.size': 3,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,           # Show ticks on top
        'xtick.bottom': True,        # Show ticks on bottom
        'ytick.left': True,          # Show ticks on left
        'ytick.right': True,         # Show ticks on right
        'xtick.minor.visible': True, # Show minor ticks
        'ytick.minor.visible': True, # Show minor ticks
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
    })


def add_ticks_all_sides(ax):
    """
    Add minor tick locators to axis.
    
    Tick visibility on all 4 sides is now handled globally by 
    setup_publication_quality(). This function just adds the 
    automatic minor tick positioning.
    
    Skips minor tick locators for logarithmic scales (AutoMinorLocator
    doesn't work on log scales).
    
    Args:
        ax: Matplotlib axis object
    """
    # Only add AutoMinorLocator for linear scales
    if ax.get_xscale() == 'linear':
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if ax.get_yscale() == 'linear':
        ax.yaxis.set_minor_locator(AutoMinorLocator())


def lambda_sq_to_freq_ghz(lambda_sq: np.ndarray) -> np.ndarray:
    """Convert wavelength squared (m²) to frequency (GHz)."""
    return (C_LIGHT / np.sqrt(lambda_sq + 1e-8) / 1e9)


def plot_data_diagnostic(data: PolarizationData, 
                         I_fitter: Optional[IFitter] = None,
                         title: str = "Polarization Data Diagnostic",
                         filename: Optional[str] = None,
                         output_dir: Optional[Path] = None,
                         show: bool = False,
                         error_method: str = 'mc',
                         n_mc_samples: int = 100) -> plt.Figure:
    """
    Create comprehensive 2×2 diagnostic plot for polarization data.
    
    Panels:
    - Top-left: Stokes I vs λ² (with optional fit), frequency on top axis
    - Top-right: Q, U, P vs λ²
    - Bottom-left: Polarization angle ψ vs λ²
    - Bottom-right: Q-U plane
    
    Args:
        data: PolarizationData object
        I_fitter: Optional IFitter object (if I was fitted)
        title: Overall plot title
        filename: If provided, save figure to output_dir/filename
        output_dir: Directory for saving (default: current directory)
        show: If True, display plot with plt.show()
        error_method: Error computation method - 'mc' (Monte Carlo) or 'analytic'
        n_mc_samples: Number of Monte Carlo samples (only if error_method='mc')
        
    Returns:
        Figure object
    """
    setup_publication_quality()
    
    # Compute P and chi with errors using helper function
    P, P_err, chi, chi_err = compute_polarization_quantities(
        data.Q, data.U, data.Q_err, data.U_err, 
        method=error_method, n_samples=n_mc_samples
    )
    chi_deg = np.degrees(chi)  # Convert to degrees for plotting
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Convert lambda_sq to frequency for top axis
    freq_ghz = lambda_sq_to_freq_ghz(data.lambda_sq)
    
    # ========================================================================
    # Top-left: Stokes I vs λ²
    # ========================================================================
    ax = axes[0, 0]
    
    if data.I is not None:
        # Plot observed I
        ax.errorbar(data.lambda_sq, data.I, yerr=data.I_err,
                   fmt='o', color='black', alpha=0.6, markersize=4,
                   label='Observed I', zorder=1)
        
        # If I was fitted, show the fit
        if I_fitter is not None and I_fitter.fit_params is not None:
            # Get fitted I and error from stored values
            if I_fitter.I_fit is not None:
                I_fit = I_fitter.I_fit
                I_fit_err = I_fitter.I_fit_err
                lsq_plot = I_fitter.lambda_sq
            else:
                # Fallback: recompute from params (no error band)
                if I_fitter.model == 'power-law':
                    # 2-parameter power-law with fixed reference
                    I0, alpha = I_fitter.fit_params
                    lsq0 = I_fitter.lambda_sq_0 if hasattr(I_fitter, 'lambda_sq_0') else np.median(data.lambda_sq)
                    I_fit = I0 * (data.lambda_sq / lsq0)**alpha
                else:
                    I_fit = np.polyval(I_fitter.fit_params, data.lambda_sq)
                I_fit_err = None
                lsq_plot = data.lambda_sq
            
            # Plot fit line
            if I_fitter.model == 'power-law':
                alpha = I_fitter.fit_params[1]  # Alpha is now at index 1
                label = f'Fit: power-law (α={alpha:.2f})'
            else:
                label = f'Fit: polynomial (deg={I_fitter.poly_degree})'
            
            ax.plot(lsq_plot, I_fit, '-', color='red', 
                   linewidth=2, label=label, zorder=2)
            
            # Plot uncertainty band if available
            if I_fit_err is not None and np.any(I_fit_err > 0):
                ax.fill_between(lsq_plot, I_fit - I_fit_err, I_fit + I_fit_err,
                               color='red', alpha=0.2, zorder=1,
                               label='Fit uncertainty')
        
        ax.set_ylabel('Stokes I (Jy)', fontsize=12)
        ax.legend(loc='best', framealpha=0.9)
    else:
        ax.text(0.5, 0.5, 'No Stokes I data\n(fractional mode)', 
               ha='center', va='center', transform=ax.transAxes,
               fontsize=12, color='gray')
        ax.set_ylabel('Stokes I (Jy)', fontsize=12)
    
    ax.set_xlabel('λ² (m²)', fontsize=12)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)
    
    # Add frequency labels to top axis
    if data.I is not None:
        # Get tick positions
        xticks = ax.get_xticks()
        freq_labels = [f'{lambda_sq_to_freq_ghz(np.array([t]))[0]:.2f}' for t in xticks]
        
        # Create secondary top axis with frequency labels
        ax_top = ax.secondary_xaxis('top')
        ax_top.set_xticks(xticks)
        ax_top.set_xticklabels(freq_labels)
        ax_top.set_xlabel('Frequency (GHz)', fontsize=12)
        ax_top.tick_params(which='both', direction='in')
        ax_top.xaxis.set_minor_locator(AutoMinorLocator())
    
    # ========================================================================
    # Top-right: Q, U, P vs λ²
    # ========================================================================
    ax = axes[0, 1]
    
    ax.errorbar(data.lambda_sq, data.Q, yerr=data.Q_err,
               fmt='o', color='blue', alpha=0.6, markersize=4,
               label='Q', zorder=2)
    ax.errorbar(data.lambda_sq, data.U, yerr=data.U_err,
               fmt='s', color='red', alpha=0.6, markersize=4,
               label='U', zorder=2)
    
    # Plot P with asymmetric errors (if MC) or symmetric (if analytic)
    if error_method == 'mc':
        P_lower, P_upper = P_err
        ax.errorbar(data.lambda_sq, P, yerr=[P_lower, P_upper],
                   fmt='o', linestyle='', color='purple', alpha=0.7, markersize=4,
                   label='P = √(Q²+U²)', zorder=3)
    else:
        ax.errorbar(data.lambda_sq, P, yerr=P_err,
                   fmt='o', linestyle='', color='purple', alpha=0.7, markersize=4,
                   label='P = √(Q²+U²)', zorder=3)
    
    ax.axhline(0, color='gray', linestyle=':', linewidth=1, zorder=0)
    ax.set_xlabel('λ² (m²)', fontsize=12)
    ax.set_ylabel('Fractional Stokes Q, U, P', fontsize=12)
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Bottom-left: Polarization angle ψ vs λ²
    # ========================================================================
    ax = axes[1, 0]
    
    # Plot angle with asymmetric errors (if MC) or symmetric (if analytic)
    if error_method == 'mc':
        chi_lower, chi_upper = chi_err
        chi_lower_deg = np.degrees(chi_lower)
        chi_upper_deg = np.degrees(chi_upper)
        ax.errorbar(data.lambda_sq, chi_deg, yerr=[chi_lower_deg, chi_upper_deg],
                   fmt='o', linestyle='', color='purple', alpha=0.7, markersize=4)
    else:
        chi_err_deg = np.degrees(chi_err)
        ax.errorbar(data.lambda_sq, chi_deg, yerr=chi_err_deg,
                   fmt='o', linestyle='', color='purple', alpha=0.7, markersize=4)
    
    ax.set_xlabel('λ² (m²)', fontsize=12)
    ax.set_ylabel('Polarization Angle ψ (deg)', fontsize=12)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Bottom-right: Q-U plane
    # ========================================================================
    ax = axes[1, 1]
    
    ax.errorbar(data.Q, data.U, xerr=data.Q_err, yerr=data.U_err,
               fmt='o', color='darkgreen', alpha=0.6, markersize=4,
               elinewidth=1)
    
    # Add origin lines
    ax.axhline(0, color='gray', linestyle=':', linewidth=1)
    ax.axvline(0, color='gray', linestyle=':', linewidth=1)
    
    # Equal aspect ratio for Q-U plane
    ax.set_aspect('equal', adjustable='box')
    
    ax.set_xlabel('Fractional Q', fontsize=12)
    ax.set_ylabel('Fractional U', fontsize=12)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Overall formatting
    # ========================================================================
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save if requested
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig


def plot_stokes_comparison(data_list: list, 
                           labels: list,
                           title: str = "Stokes Q, U Comparison",
                           filename: Optional[str] = None,
                           output_dir: Optional[Path] = None,
                           show: bool = False) -> plt.Figure:
    """
    Compare Q and U from multiple datasets (useful for comparing loading modes).
    
    Args:
        data_list: List of PolarizationData objects
        labels: List of labels for each dataset
        title: Plot title
        filename: If provided, save figure
        output_dir: Directory for saving
        show: If True, display plot
        
    Returns:
        Figure object
    """
    setup_publication_quality()
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(data_list)))
    
    for data, label, color in zip(data_list, labels, colors):
        # Q panel
        axes[0].errorbar(data.lambda_sq, data.Q, yerr=data.Q_err,
                        fmt='o-', alpha=0.7, label=label, color=color,
                        markersize=4)
        
        # U panel
        axes[1].errorbar(data.lambda_sq, data.U, yerr=data.U_err,
                        fmt='o-', alpha=0.7, label=label, color=color,
                        markersize=4)
    
    axes[0].axhline(0, color='gray', linestyle=':', linewidth=1)
    axes[0].set_ylabel('Fractional Q', fontsize=12)
    axes[0].legend(loc='best', framealpha=0.9)
    axes[0].grid(True, alpha=0.3)
    add_ticks_all_sides(axes[0])
    
    axes[1].axhline(0, color='gray', linestyle=':', linewidth=1)
    axes[1].set_xlabel('λ² (m²)', fontsize=12)
    axes[1].set_ylabel('Fractional U', fontsize=12)
    axes[1].legend(loc='best', framealpha=0.9)
    axes[1].grid(True, alpha=0.3)
    add_ticks_all_sides(axes[1])
    
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig


def plot_fit_results(results: 'FitResults',
                    true_model: Optional['CompositeModel'] = None,
                    I_fitter: Optional[IFitter] = None,
                    title: str = "Fit Results",
                    filename: Optional[str] = None,
                    output_dir: Optional[Path] = None,
                    show: bool = False,
                    error_method: str = 'mc',
                    n_mc_samples: int = 100,
                    n_posterior_samples: int = 100,
                    padding: float = 1.25,
                    phi_range: Optional[float] = None,
                    n_phi: int = 2048,
                    weight_type: str = 'variance',
                    model_method: str = 'representative',
                    plot_fdf_posterior_samples: bool = False) -> plt.Figure:
    """
    Create diagnostic plot for fit results.
    
    Shows posterior percentile bands (68% and 95% credible intervals) overlaid on data.
    Residual panels show χ_Q and χ_U computed from posterior median.
    
    Nine-panel plot:
    - Top: Stokes I vs λ² with fit (spanning both columns, frequency axis on top)
    - Left: Q vs λ² with posterior bands
    - Q residual: χ_Q vs λ² (below Q panel, shared x-axis)
    - Right: U vs λ² with posterior bands
    - U residual: χ_U vs λ² (below U panel, shared x-axis)
    - Middle-left: Q-U plane
    - Middle-right: Polarization angle ψ vs λ²
    - Lower: Fractional polarization p vs λ² (spanning both columns)
    - Bottom: FDF diagnostic (spanning both columns)
    
    Args:
        results: FitResults object from fitting
        true_model: Optional true model for comparison (synthetic data)
        I_fitter: Optional IFitter object for Stokes I panel. If None, panel
                  shows placeholder text.
        title: Plot title
        filename: If provided, save figure
        output_dir: Directory for saving
        show: If True, display plot
        error_method: Error computation method - 'mc' (Monte Carlo) or 'analytic'
        n_mc_samples: Number of Monte Carlo samples (only if error_method='mc')
        n_posterior_samples: Number of posterior realizations to show (0 to disable)
        padding: λ² range multiplier for model evaluation and plotting (default=1.0)
                 Values > 1.0 extend range beyond data for debugging.
        phi_range: Faraday depth range (±phi_range) for FDF panel. If None, 
                   auto-computed from RMTF FWHM and fitted RM values.
        n_phi: Number of points in Faraday depth grid for FDF panel
        weight_type: Weighting for RM synthesis ('variance' or 'uniform')
        model_method: Method for FDF model parameters ('representative', 'median',
                     'mean', 'map', 'mode')
        plot_fdf_posterior_samples: If True, plot posterior samples on FDF panel
        
    Returns:
        Figure object
    """
    setup_publication_quality()
    
    # Widened for better visibility
    fig = plt.figure(figsize=(12, 17))
    gs_outer = fig.add_gridspec(2, 1,
                                height_ratios=[0.6 + 1 + 0.25, 1 + 1 + 1],
                                hspace=0.09)
    gs_top = gs_outer[0].subgridspec(3, 2,
                                     height_ratios=[0.6, 1, 0.25],
                                     hspace=0.15, wspace=0.3)
    gs_bot = gs_outer[1].subgridspec(3, 2,
                                     height_ratios=[1, 1, 1],
                                     hspace=0.35, wspace=0.3)
    class _GS:
        def __getitem__(self, key):
            row, col = key
            return gs_top[row, col] if row <= 2 else gs_bot[row - 3, col]
    gs = _GS()
    
    # Get data
    data = results.setup.data
    lambda_sq = data.lambda_sq
    
    # Apply padding to plotting range
    lambda_sq_min = lambda_sq.min()
    lambda_sq_max = lambda_sq.max()
    delta = lambda_sq_max - lambda_sq_min
    extra = (padding - 1.0) * delta / 2.0
    lambda_sq_plot_min = lambda_sq_min - extra
    lambda_sq_plot_max = lambda_sq_max + extra
    
    # Create fine grid for smooth model curves
    lambda_sq_fine = np.linspace(lambda_sq_plot_min, lambda_sq_plot_max, 500)
    
    
    # Sample posterior realizations for uncertainty visualization
    posterior_models_data = []  # On data points (for residuals)
    posterior_models_fine = []  # On fine grid (for smooth plotting)
    if n_posterior_samples > 0:
        # Sample parameter sets from weighted posterior
        n_samples_total = len(results.samples)
        sample_indices = np.random.choice(n_samples_total, 
                                         size=min(n_posterior_samples, n_samples_total),
                                         replace=False,
                                         p=results.weights)
        
        # Compute model for each sample on both grids
        for idx in sample_indices:
            params = results.samples[idx]
            model_sample = copy.deepcopy(results.setup.model)
            update_model_from_params(model_sample, params, results.setup.param_names, lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
            P_sample_data = model_sample.compute_polarization(lambda_sq)
            P_sample_fine = model_sample.compute_polarization(lambda_sq_fine)
            posterior_models_data.append(P_sample_data)
            posterior_models_fine.append(P_sample_fine)
    
    # Compute percentile bands from posterior samples
    if len(posterior_models_data) > 0:
        # Stack samples on data grid: shape (n_samples, n_lambda)
        Q_samples_data = np.array([P.real for P in posterior_models_data])
        U_samples_data = np.array([P.imag for P in posterior_models_data])
        
        # Stack samples on fine grid
        Q_samples_fine = np.array([P.real for P in posterior_models_fine])
        U_samples_fine = np.array([P.imag for P in posterior_models_fine])
        P_abs_samples_fine = np.array([np.abs(P) for P in posterior_models_fine])
        psi_samples_fine = np.array([np.angle(P)/2 for P in posterior_models_fine])
        
        # Compute percentiles at data points (for residuals)
        Q_percentiles_data = np.percentile(Q_samples_data, [2.5, 16, 50, 84, 97.5], axis=0)
        U_percentiles_data = np.percentile(U_samples_data, [2.5, 16, 50, 84, 97.5], axis=0)
        
        # Compute percentiles on fine grid (for smooth plotting)
        Q_percentiles_fine = np.percentile(Q_samples_fine, [2.5, 16, 50, 84, 97.5], axis=0)
        U_percentiles_fine = np.percentile(U_samples_fine, [2.5, 16, 50, 84, 97.5], axis=0)
        P_percentiles_fine = np.percentile(P_abs_samples_fine, [2.5, 16, 50, 84, 97.5], axis=0)
        psi_percentiles_fine = np.percentile(psi_samples_fine, [2.5, 16, 50, 84, 97.5], axis=0)
        
        # Compute chi residuals from median at data points
        Q_median = Q_percentiles_data[2]
        U_median = U_percentiles_data[2]
        
        # Chi with statistical errors only
        chi_Q_stat = (data.Q - Q_median) / data.Q_err
        chi_U_stat = (data.U - U_median) / data.U_err
        
    else:
        Q_percentiles_data = U_percentiles_data = None
        Q_percentiles_fine = U_percentiles_fine = P_percentiles_fine = psi_percentiles_fine = None
        chi_Q_stat = chi_U_stat = None
    
    # Compute P and chi with errors for observed data
    # Ricean debiasing is on by default (Hales 2012): affects P panel only, not EVPA or Q/U.
    # Compute with statistical errors only
    p_data_stat, p_data_err_stat, chi_data_stat, chi_data_err_stat = compute_polarization_quantities(
        data.Q, data.U, data.Q_err, data.U_err,
        method=error_method, n_samples=n_mc_samples, debias=True
    )
    psi_data_deg_stat = np.degrees(chi_data_stat)
    psi_data_err_deg_stat = np.degrees(chi_data_err_stat)
    
    
    # Get true model if provided
    if true_model is not None:
        P_true = true_model.compute_polarization(lambda_sq)
    
    # Complex data
    P_data = data.Q + 1j * data.U
    
    # ========================================================================
    # Top: Stokes I vs λ² (full width, frequency axis on top)
    # ========================================================================
    ax_I = fig.add_subplot(gs[0, :])
    
    if I_fitter is not None and data.I is not None:
        # Plot observed I
        ax_I.errorbar(data.lambda_sq, data.I, yerr=data.I_err,
                     fmt='o', markerfacecolor='lightgrey', markeredgecolor='black',
                     ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                     zorder=1)
        
        # Overlay Stokes I fit
        if I_fitter.fit_params is not None:
            if I_fitter.I_fit is not None:
                I_fit = I_fitter.I_fit
                I_fit_err = I_fitter.I_fit_err
                lsq_plot = I_fitter.lambda_sq
            else:
                if I_fitter.model == 'power-law':
                    I0, alpha_I = I_fitter.fit_params
                    lsq0 = I_fitter.lambda_sq_0 if hasattr(I_fitter, 'lambda_sq_0') else np.median(data.lambda_sq)
                    I_fit = I0 * (data.lambda_sq / lsq0)**alpha_I
                else:
                    I_fit = np.polyval(I_fitter.fit_params, data.lambda_sq)
                I_fit_err = None
                lsq_plot = data.lambda_sq
            
            if I_fitter.model == 'power-law':
                alpha_freq = -2.0 * I_fitter.fit_params[1]
                fit_label = f'Fit: power-law (α={alpha_freq:.2f})'
            else:
                fit_label = f'Fit: polynomial (deg={I_fitter.poly_degree})'
            
            ax_I.plot(lsq_plot, I_fit, '-', color='C0', 
                     linewidth=2, label=fit_label, zorder=2)
            
            if I_fit_err is not None and np.any(I_fit_err > 0):
                ax_I.fill_between(lsq_plot, I_fit - I_fit_err, I_fit + I_fit_err,
                                 color='C0', alpha=0.2, zorder=1)
        
        ax_I.set_ylabel('Stokes I (Jy)', fontsize=13)
        ax_I.legend(loc='best', framealpha=0.9)
    else:
        ax_I.text(0.5, 0.5, 'No Stokes I data\n(fractional mode)', 
                 ha='center', va='center', transform=ax_I.transAxes,
                 fontsize=12, color='gray')
        ax_I.set_ylabel('Stokes I (Jy)', fontsize=13)
    
    ax_I.set_xlim(lambda_sq_plot_min, lambda_sq_plot_max)
    ax_I.tick_params(axis='both', labelsize=11, labelbottom=False)
    ax_I.grid(True, alpha=0.3)
    add_ticks_all_sides(ax_I)
    
    # Frequency axis on top
    ax_I_freq = ax_I.secondary_xaxis('top')
    xticks_I = ax_I.get_xticks()
    # Filter ticks to those within the plot range
    xticks_I = xticks_I[(xticks_I >= lambda_sq_plot_min) & (xticks_I <= lambda_sq_plot_max)]
    freq_labels_I = [f'{lambda_sq_to_freq_ghz(np.array([t]))[0]:.2f}' for t in xticks_I]
    ax_I_freq.set_xticks(xticks_I)
    ax_I_freq.set_xticklabels(freq_labels_I)
    ax_I_freq.set_xlabel('Frequency (GHz)', fontsize=13)
    ax_I_freq.tick_params(which='both', direction='in', labelsize=11)
    ax_I_freq.xaxis.set_minor_locator(AutoMinorLocator())
    
    # ========================================================================
    # Q vs λ²
    # ========================================================================
    ax_q = fig.add_subplot(gs[1, 0])
    
    # Plot posterior percentile bands on fine grid
    if Q_percentiles_fine is not None:
        ax_q.fill_between(lambda_sq_fine, Q_percentiles_fine[1], Q_percentiles_fine[3],
                          color='C0', alpha=0.3, zorder=3, label='68% credible')
        ax_q.plot(lambda_sq_fine, Q_percentiles_fine[2], '-', color='C0',
                 linewidth=1.5, alpha=0.8, zorder=4, label='Median')
    
    # Plot data with error bars
    ax_q.errorbar(lambda_sq, data.Q, yerr=data.Q_err, fmt='o',
                markerfacecolor='lightgrey', markeredgecolor='black',
                ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                label='Data', zorder=1)
    
    if true_model is not None:
        ax_q.plot(lambda_sq, P_true.real, '--', color='green', 
               linewidth=2.5, alpha=0.7, label='True model', zorder=5)
    
    ax_q.set_ylabel(r'$q\,(\lambda^2)$', fontsize=13)
    ax_q.grid(True, alpha=0.3)
    ax_q.tick_params(axis='both', labelsize=11, labelbottom=False)
    ax_q.set_xlim(lambda_sq_plot_min, lambda_sq_plot_max)
    # Range set by 2nd-98th percentile of data ± 2 × median(err); outliers poke outside
    q_err_median = np.median(data.Q_err)
    q_ylim_min = np.percentile(data.Q, 2) - 3 * q_err_median
    q_ylim_max = np.percentile(data.Q, 98) + 3 * q_err_median
    ax_q.set_ylim(q_ylim_min, q_ylim_max)
    add_ticks_all_sides(ax_q)
    
    # Q chi residuals (below Q panel, shared x-axis)
    ax_q_res = fig.add_subplot(gs[2, 0], sharex=ax_q)
    if chi_Q_stat is not None:
        ax_q_res.scatter(lambda_sq, chi_Q_stat, c='black', s=20, alpha=1.0, edgecolors='none')
    ax_q_res.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8, zorder=10)
    ax_q_res.set_xlabel('λ² (m²)', fontsize=13)
    ax_q_res.set_ylabel(r'$\chi_q$', fontsize=13)
    ax_q_res.tick_params(axis='both', labelsize=11)
    ax_q_res.grid(True, alpha=0.3)
    add_ticks_all_sides(ax_q_res)
    
    # Top-right: U vs λ²
    ax_u = fig.add_subplot(gs[1, 1])
    
    # Plot posterior percentile bands on fine grid
    if U_percentiles_fine is not None:
        ax_u.fill_between(lambda_sq_fine, U_percentiles_fine[1], U_percentiles_fine[3],
                          color='C0', alpha=0.3, zorder=3, label='68% credible')
        ax_u.plot(lambda_sq_fine, U_percentiles_fine[2], '-', color='C0',
                 linewidth=1.5, alpha=0.8, zorder=4, label='Median')
    
    # Plot data with error bars
    ax_u.errorbar(lambda_sq, data.U, yerr=data.U_err, fmt='o',
                markerfacecolor='lightgrey', markeredgecolor='black',
                ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                label='Data', zorder=1)
    
    if true_model is not None:
        ax_u.plot(lambda_sq, P_true.imag, '--', color='green', 
               linewidth=2.5, alpha=0.7, label='True model', zorder=5)
    
    ax_u.set_ylabel(r'$u\,(\lambda^2)$', fontsize=13)
    ax_u.grid(True, alpha=0.3)
    ax_u.tick_params(axis='both', labelsize=11, labelbottom=False)
    ax_u.set_xlim(lambda_sq_plot_min, lambda_sq_plot_max)
    # Range set by 2nd-98th percentile of data ± 2 × median(err); outliers poke outside
    u_err_median = np.median(data.U_err)
    u_ylim_min = np.percentile(data.U, 2) - 3 * u_err_median
    u_ylim_max = np.percentile(data.U, 98) + 3 * u_err_median
    ax_u.set_ylim(u_ylim_min, u_ylim_max)
    add_ticks_all_sides(ax_u)
    
    # U chi residuals (below U panel, shared x-axis)
    ax_u_res = fig.add_subplot(gs[2, 1], sharex=ax_u)
    if chi_U_stat is not None:
        ax_u_res.scatter(lambda_sq, chi_U_stat, c='black', s=20, alpha=1.0, edgecolors='none')
    ax_u_res.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8, zorder=10)
    ax_u_res.set_xlabel('λ² (m²)', fontsize=13)
    ax_u_res.set_ylabel(r'$\chi_u$', fontsize=13)
    ax_u_res.tick_params(axis='both', labelsize=11)
    ax_u_res.grid(True, alpha=0.3)
    add_ticks_all_sides(ax_u_res)
    
    # Middle-left: Q-U plane
    ax = fig.add_subplot(gs[3, 0])
    
    # Plot posterior percentile bands on fine grid
    if Q_percentiles_fine is not None and U_percentiles_fine is not None:
        vertices_68 = np.column_stack([
            np.concatenate([Q_percentiles_fine[1], Q_percentiles_fine[3][::-1]]),
            np.concatenate([U_percentiles_fine[1], U_percentiles_fine[3][::-1]])
        ])
        poly_68 = Polygon(vertices_68, facecolor='C0', alpha=0.3,
                         edgecolor='none', zorder=3, label='68% credible')
        ax.add_patch(poly_68)
        ax.plot(Q_percentiles_fine[2], U_percentiles_fine[2], '-', color='C0',
               linewidth=1.5, alpha=0.8, zorder=4, label='Median')
    
    ax.errorbar(data.Q, data.U, xerr=data.Q_err, yerr=data.U_err,
               fmt='o', markerfacecolor='lightgrey', markeredgecolor='black',
               ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
               elinewidth=1, label='Data', zorder=1)
    
    if true_model is not None:
        ax.plot(P_true.real, P_true.imag, '--', color='green',
               linewidth=2.5, alpha=0.7, label='True model', zorder=5)
    
    ax.axhline(0, color='gray', linestyle=':', linewidth=1)
    ax.axvline(0, color='gray', linestyle=':', linewidth=1)
    ax.set_xlabel(r'$q\,(\lambda^2)$', fontsize=13)
    ax.set_ylabel(r'$u\,(\lambda^2)$', fontsize=13)
    ax.tick_params(axis='both', labelsize=11)
    ax.tick_params(axis='x', rotation=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(q_ylim_min, q_ylim_max)
    ax.set_ylim(u_ylim_min, u_ylim_max)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Middle-right: Polarization angle ψ vs λ²
    # ========================================================================
    ax = fig.add_subplot(gs[3, 1])
    
    # Plot posterior percentile bands on fine grid
    if psi_percentiles_fine is not None:
        psi_deg = np.degrees(psi_percentiles_fine)
        ax.fill_between(lambda_sq_fine, psi_deg[1], psi_deg[3],
                        color='C0', alpha=0.3, zorder=3, label='68% credible')
        # Median
        ax.plot(lambda_sq_fine, psi_deg[2], '-', color='C0',
               linewidth=1.5, alpha=0.8, zorder=4, label='Median')
    
    # Plot data with asymmetric errors (if MC) or symmetric (if analytic)
    if error_method == 'mc':
        chi_lower, chi_upper = chi_data_err_stat
        psi_lower_deg = np.degrees(chi_lower)
        psi_upper_deg = np.degrees(chi_upper)
        ax.errorbar(lambda_sq, psi_data_deg_stat, yerr=[psi_lower_deg, psi_upper_deg],
                   fmt='o', linestyle='', markerfacecolor='lightgrey', markeredgecolor='black',
                   ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                   label='Data', zorder=1)
    else:
        psi_err_deg = np.degrees(chi_data_err_stat)
        ax.errorbar(lambda_sq, psi_data_deg_stat, yerr=psi_err_deg,
                   fmt='o', linestyle='', markerfacecolor='lightgrey', markeredgecolor='black',
                   ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                   label='Data', zorder=1)
    
    if true_model is not None:
        psi_true = np.angle(P_true) / 2
        ax.plot(lambda_sq, np.degrees(psi_true), '--', color='green',
               linewidth=2.5, alpha=0.7, label='True model', zorder=5)
    
    ax.set_xlabel('λ² (m²)', fontsize=13)
    ax.set_ylabel(r'$\psi\,(\lambda^2)$ (deg)', fontsize=13)
    ax.tick_params(axis='both', labelsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(lambda_sq_plot_min, lambda_sq_plot_max)
    # Range set by 2nd-98th percentile of data ± 2 × median(err); outliers poke outside
    psi_data_use = psi_data_deg_stat
    psi_err_raw = psi_data_err_deg_stat
    psi_err_median = np.median(np.concatenate(psi_err_raw) if isinstance(psi_err_raw, (tuple, list)) else psi_err_raw)
    psi_ylim_min = np.percentile(psi_data_use, 2) - 3 * psi_err_median
    psi_ylim_max = np.percentile(psi_data_use, 98) + 3 * psi_err_median
    ax.set_ylim(psi_ylim_min, psi_ylim_max)
    add_ticks_all_sides(ax)
    
    # Bottom: Fractional polarization p vs λ² (spanning both columns)
    ax = fig.add_subplot(gs[4, :])
    
    # Plot posterior percentile bands on fine grid
    if P_percentiles_fine is not None:
        ax.fill_between(lambda_sq_fine, P_percentiles_fine[1], P_percentiles_fine[3],
                        color='C0', alpha=0.3, zorder=3, label='68% credible')
        # Median
        ax.plot(lambda_sq_fine, P_percentiles_fine[2], '-', color='C0',
               linewidth=1.5, alpha=0.8, zorder=4, label='Median')
    
    # Plot data with asymmetric errors (if MC) or symmetric (if analytic)
    if error_method == 'mc':
        p_lower, p_upper = p_data_err_stat
        ax.errorbar(lambda_sq, p_data_stat, yerr=[p_lower, p_upper],
                   fmt='o', linestyle='', markerfacecolor='lightgrey', markeredgecolor='black',
                   ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                   zorder=1)
    else:
        ax.errorbar(lambda_sq, p_data_stat, yerr=p_data_err_stat,
                   fmt='o', linestyle='', markerfacecolor='lightgrey', markeredgecolor='black',
                   ecolor='black', alpha=1.0, markersize=4, markeredgewidth=1,
                   zorder=1)
    
    if true_model is not None:
        p_true = np.abs(P_true)
        ax.plot(lambda_sq, p_true, '--', color='green',
               linewidth=2.5, alpha=0.7, label='True model', zorder=5)
    
    ax.set_xlabel('λ² (m²)', fontsize=13)
    ax.set_ylabel(r'$p\,(\lambda^2)$', fontsize=13)
    ax.tick_params(axis='both', labelsize=11)
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(lambda_sq_plot_min, lambda_sq_plot_max)
    # Range set by 2nd-98th percentile of data ± 2 × median(err); outliers poke outside
    p_data_use = p_data_stat
    p_err_raw = p_data_err_stat
    p_err_median = np.median(np.concatenate(p_err_raw) if isinstance(p_err_raw, (tuple, list)) else p_err_raw)
    p_ylim_min = np.percentile(p_data_use, 2) - 3 * p_err_median
    p_ylim_max = np.percentile(p_data_use, 98) + 3 * p_err_median
    ax.set_ylim(p_ylim_min, p_ylim_max)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # FDF diagnostic (spanning both columns)
    # ========================================================================
    ax_fdf = fig.add_subplot(gs[5, :])
    
    # Get model for FDF computation
    if model_method == 'representative':
        representative_params = results.get_representative_sample(percentile=1.0, method='median')
        fdf_model = copy.deepcopy(results.setup.model)
        update_model_from_params(fdf_model, representative_params, results.setup.param_names, lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
    else:
        fdf_model = results.get_best_fit_model(method=model_method)
    
    # Auto-compute phi_range if not provided
    if phi_range is None:
        _, rmtf_props_temp = compute_rmtf(lambda_sq, np.array([0.0]),
                                          weight_type=weight_type,
                                          Q_err=data.Q_err, U_err=data.U_err)
        fwhm = rmtf_props_temp['expected_fwhm_nominal']
        
        rm_values = []
        for comp in fdf_model.components:
            if hasattr(comp, 'phi_rm'):
                rm_values.append(comp.phi_rm)
            elif hasattr(comp, 'phi_peak'):
                rm_values.append(comp.phi_peak)
        
        if rm_values:
            rm_max = np.max(rm_values)
            rm_min = np.min(rm_values)
            phi_range_comps = max(abs(rm_max + 1.0 * fwhm), abs(rm_min - 1.0 * fwhm))
        else:
            phi_range_comps = 0.0
        
        phi_range = max(5.0 * fwhm, phi_range_comps, 100.0)
    
    # Faraday depth grids (wide for computation, central for display)
    phi_grid = np.linspace(-phi_range, phi_range, n_phi)
    phi_grid_wide = np.linspace(-5*phi_range, 5*phi_range, 5*n_phi)
    n_start = 2*n_phi
    n_end = 3*n_phi
    
    
    # Per-component intrinsic FDFs with optional posterior samples
    if plot_fdf_posterior_samples:
        n_samples_total = len(results.samples)
        weights_normalized = results.weights / np.sum(results.weights)
        n_samples_plot = min(100, n_samples_total)
        sample_indices = np.random.choice(n_samples_total, size=n_samples_plot,
                                         replace=False, p=weights_normalized)
        
        for idx in sample_indices:
            params = results.samples[idx]
            model_sample = copy.deepcopy(results.setup.model)
            update_model_from_params(model_sample, params, results.setup.param_names, lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
            for ci, comp in enumerate(model_sample.components):
                comp_fdf_sample = comp.compute_fdf(phi_grid)
                ax_fdf.plot(phi_grid, np.abs(comp_fdf_sample), '-', color='C0', linewidth=0.5,
                           alpha=0.15, zorder=1)
    
    # Plot individual component FDFs (representative model)
    for i, comp in enumerate(fdf_model.components):
        comp_fdf = comp.compute_fdf(phi_grid)
        comp_amp = np.abs(comp_fdf)
        lbl = 'Model' if i == 0 else '_nolegend_'
        ax_fdf.plot(phi_grid, comp_amp, 'k-', linewidth=2.0, alpha=1.0, zorder=12,
                   label=lbl)
        
        if isinstance(comp, (ThinComponent, ThinPowerLawComponent)):
            peak_amp = np.max(comp_amp)
            ax_fdf.plot(comp.phi_rm, peak_amp, marker='^',
                       markersize=10, color='black', zorder=15)
    
    ax_fdf.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax_fdf.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax_fdf.set_yscale('log')
    comp_fdfs = [comp.compute_fdf(phi_grid) for comp in fdf_model.components]
    max_fdf_amp = max(np.max(np.abs(f)) for f in comp_fdfs) if comp_fdfs else 0.0
    if max_fdf_amp > 0:
        y_upper = min(max_fdf_amp * 10**0.5, 1.0)
        ax_fdf.set_ylim(1e-5, y_upper)
    else:
        ax_fdf.set_ylim(1e-5, 1e-3)
    ax_fdf.set_xlabel(r'$\phi_{\rm f}$ (rad m$^{-2}$)', fontsize=13)
    ax_fdf.set_ylabel(r'$|F(\phi_{\rm f})|$ (fractional)', fontsize=13)
    ax_fdf.tick_params(axis='both', labelsize=11)
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color='k', linewidth=2.0, label='Model'),
        Line2D([0], [0], color='C0', linewidth=1.5, label='Posterior samples'),
    ]
    if not plot_fdf_posterior_samples:
        legend_handles = legend_handles[:1]
    ax_fdf.legend(handles=legend_handles, loc='best', framealpha=0.9, fontsize=13)
    ax_fdf.grid(True, alpha=0.3, which='both')
    add_ticks_all_sides(ax_fdf)
    
    # ========================================================================
    # Overall formatting
    # ========================================================================
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig


def format_parameter_label(param_name: str, setup: Optional['FitSetup'] = None) -> str:
    r"""
    Convert parameter name to publication-quality LaTeX label.
    
    Maps internal parameter names (e.g., 'Source_phi_rm') to 
    formatted LaTeX labels (e.g., r'$\phi_{RM}$ (Source)').
    
    Supports both Cartesian (a, b) and derived (p0, psi_0) parameters:
    - Cartesian: a, b (intrinsic polarization in Stokes space)
    - Derived: p0 (polarization fraction), psi_0 (polarization angle)
    
    Special handling for p0:
    - ThinComponent: $p_{0,\mathrm{thin}}$
    - ThickComponent: $p_{0,\mathrm{thick}}$
    
    If mapping fails for any reason, returns the original parameter name.
    
    Args:
        param_name: Internal parameter name (e.g., 'Component_p0' or 'Component_a')
        setup: Optional FitSetup object to determine component types
        
    Returns:
        Formatted LaTeX string for plotting (or original name if mapping fails)
    """
    try:
        # Define known parameter names (checked as suffixes)
        # Order matters for matching: check longer names first to avoid partial matches
        known_params = ['phi_rm', 'sigma_phi', 'phi_peak', 'delta_phi',
                       'phi_fg', 'psi_0', 'p0', 'N', 'beta', 'a', 'b']

        # Map parameter names to LaTeX symbols
        param_map = {
            'phi_rm': r'$\phi_{\rm rm}$',
            'psi_0': r'$\psi_0$',
            'phi_peak': r'$\phi_{peak}$',
            'sigma_phi': r'$\sigma_\phi$',
            'N': r'$N$',
            'beta': r'$\beta$',
            'delta_phi': r'$\Delta\phi$',
            'phi_fg': r'$\phi_{FG}$',
            'a': r'$a$',
            'b': r'$b$',
        }
        
        # Check if param_name ends with any known parameter
        comp_name = None
        param = None
        
        for known_param in known_params:
            if param_name.endswith('_' + known_param):
                # Extract component name (everything before _param)
                comp_name = param_name[:-(len(known_param) + 1)]
                param = known_param
                break
            elif param_name == known_param:
                # Just the parameter, no component
                comp_name = None
                param = known_param
                break
        
        # If no match found, return original name
        if param is None:
            return param_name
        
        # Special handling for p0 - depends on component type
        if param == 'p0':
            if setup is not None and comp_name is not None:
                # Determine component type from setup
                
                # Find the component in the model
                comp_idx = None
                for i, comp in enumerate(setup.model.components):
                    comp_check_name = comp.name if comp.name else f'component_{i}'
                    if comp_check_name == comp_name or f'component_{i}' == comp_name:
                        comp_idx = i
                        break
                
                if comp_idx is not None:
                    comp = setup.model.components[comp_idx]
                    if isinstance(comp, ThinPowerLawComponent):
                        param_label = r'$p_{0,\mathrm{PL}}$'
                    elif isinstance(comp, ThinComponent):
                        param_label = r'$p_{0,\mathrm{thin}}$'
                    elif isinstance(comp, ThickComponent):
                        param_label = r'$p_{0,\mathrm{thick}}$'
                    else:
                        param_label = r'$p_0$'  # Fallback
                else:
                    param_label = r'$p_0$'  # Fallback if component not found
            else:
                # No setup provided, use generic label
                param_label = r'$p_0$'
        else:
            # Get formatted parameter name from map
            param_label = param_map.get(param, param)
        
        # Add component name if present
        if comp_name:
            return f"{param_label} ({comp_name})"
        else:
            return param_label
    
    except Exception:
        # If anything goes wrong, just return the original name
        return param_name


def adaptive_decimal_precision(median: float, lower_err: float, upper_err: float) -> int:
    """
    Determine decimal places based on error size.
    Shows 1 more decimal than the first non-zero digit in the smaller error.
    
    Examples:
        error=0.234  → first non-zero at 10^-1 → show 2 decimals
        error=0.0234 → first non-zero at 10^-2 → show 3 decimals  
        error=2.34   → first non-zero at 10^0  → show 1 decimal
        error=23.4   → first non-zero at 10^1  → show 0 decimals
    
    Args:
        median: Central value
        lower_err: Lower error (positive value)
        upper_err: Upper error (positive value)
        
    Returns:
        Number of decimal places to display
    """
    # Use smaller error to determine precision
    min_err = min(abs(lower_err), abs(upper_err))
    
    if min_err == 0 or not np.isfinite(min_err):
        return 3  # Safe default
    
    # Find position of first non-zero digit
    magnitude = np.floor(np.log10(min_err))
    
    # Show 1 more decimal place than that position
    decimals = int(-magnitude + 1)
    
    # Clamp to reasonable range [0, 10]
    return max(0, min(decimals, 10))


def get_parameter_units(param_name: str) -> str:
    """
    Get the appropriate units for a parameter.
    
    Args:
        param_name: Parameter name (e.g., 'ISM_psi_0', 'component_0_phi_rm')
        
    Returns:
        Unit string (e.g., ' deg', ' rad m^{-2}') or empty string if no units
        Note: Returns LaTeX-ready string for use inside math mode (no $ delimiters)
    """
    param_lower = param_name.lower()
    
    # Angle parameters - degrees
    if 'psi_0' in param_lower or 'psi0' in param_lower:
        return r' \, \mathrm{deg}'
    
    # Faraday depth parameters - rad m^{-2}
    if any(x in param_lower for x in ['phi_rm', 'phi_peak', 'sigma_phi']):
        return r' \, \mathrm{rad} \, \mathrm{m}^{-2}'
    
    # Spectral index - dimensionless (could add if desired)
    if 'alpha' in param_lower:
        return ''
    
    # No units for a, b, p0, N
    return ''


def plot_corner(samples: np.ndarray,
                param_names: List[str],
                param_summary: Dict,
                weights: np.ndarray,
                setup: 'FitSetup',
                processed_parameters: Optional[Dict] = None,
                truths: Optional[List[float]] = None,
                method: str = 'median',
                title: Optional[str] = None,
                filename: Optional[str] = None,
                output_dir: Optional[Path] = None,
                show: bool = False) -> plt.Figure:
    """
    Create publication-quality corner plot for posterior distribution.
    
    Generic function that plots whatever samples/parameters you provide.
    
    Recommended workflows:
    
    1. Basic fit (no processing):
       ```python
       plot_corner(
           samples=results.samples_derived,
           param_names=results.param_names_derived,
           param_summary=results.param_summary,
           weights=results.weights,
           setup=results.setup,
           method='median'
       )
       ```
    
    2. After mode processing (recommended):
       ```python
       results = process_posterior_modes(results, config)
       plot_corner(
           samples=results.samples_processed,  # Unwrapped samples
           param_names=results.param_names_processed,
           param_summary=results.param_summary,
           processed_parameters=results.processed_parameters,  # Contains mode_1
           weights=results.weights,
           setup=results.setup,
           method='mode'  # Uses highest likelihood mode
       )
       ```
    
    Shows only free parameters (fixed parameters excluded from plot).
    Uses the corner package with improved styling:
    - Formatted LaTeX labels for parameters
    - Density contours with filled regions
    - Best-fit values from mode_1 (highest likelihood mode) when method='mode'
    - Clean, professional appearance
    
    Args:
        samples: Posterior samples to plot (n_samples, n_params)
                 Recommend: samples_processed (unwrapped) after processing
        param_names: List of parameter names matching samples columns
        param_summary: Dictionary of parameter summaries (mean, median, quantiles)
        weights: Importance weights for resampling (from dynesty)
        setup: FitSetup object for context (fixed params, component names, etc.)
        processed_parameters: Optional dict from process_posterior_modes with mode_1 stats
                             Contains mode_1, mean, median, eti_68, eti_99 from unwrapped samples
        truths: Optional list of true parameter values (for synthetic data)
                If None, uses best-fit values based on method
        method: Method for determining best-fit when truths=None
                'median': Use median from param_summary or processed_parameters
                'mean': Use mean from param_summary or processed_parameters
                'mode': Use mode from processed_parameters['param']['mode_1']['mode']
                Default: 'median'
        title: Optional plot title (auto-generated if None)
        filename: If provided, save figure
        output_dir: Directory for saving
        show: If True, display plot
        
    Returns:
        Figure object
    """
    # Extract free parameters only (exclude fixed parameters from plot)
    # Note: If using samples_processed from process_posterior_modes, fixed parameters
    # are already filtered out. This logic mainly applies to Cartesian samples.
    if len(setup.fixed_params) > 0 and param_names == setup.param_names:
        free_idx = setup.free_param_indices
        param_names_plot = setup.free_param_names
        n_free = len(free_idx)
        n_fixed = len(setup.fixed_params)
        title_suffix = f" ({n_free} free, {n_fixed} fixed)"
    else:
        # Using derived/processed parameters (fixed already filtered) or no fixed params
        free_idx = list(range(len(param_names)))
        param_names_plot = param_names
        
        # Check if any fixed params exist but aren't in param_names (already filtered)
        if len(setup.fixed_params) > 0:
            title_suffix = f" ({len(param_names)} free, {len(setup.fixed_params)} fixed)"
        else:
            title_suffix = ""
    
    # Create formatted parameter labels
    labels = []
    for name in param_names_plot:
        label = format_parameter_label(name, setup=setup)
        labels.append(label)
    
    # If no truths provided, use best-fit values
    if truths is None:
        truths = []
        mode_fallback_warned = False
        for i in free_idx:
            pname = param_names[i]
            
            # Priority 1: Use processed_parameters if available
            if processed_parameters is not None and pname in processed_parameters:
                proc_param = processed_parameters[pname]
                if method == 'mode' and 'mode_1' in proc_param:
                    truths.append(proc_param['mode_1']['mode'])
                elif method == 'mean' and 'mean' in proc_param:
                    truths.append(proc_param['mean'])
                elif method == 'median' and 'median' in proc_param:
                    truths.append(proc_param['median'])
                else:
                    # Fallback to median in processed_parameters
                    truths.append(proc_param.get('median', proc_param.get('mean', 0.0)))
            
            # Priority 2: Use param_summary
            elif pname in param_summary:
                if method == 'mode' and 'mode' in param_summary[pname]:
                    # Mode available in param_summary (backward compatibility)
                    truths.append(param_summary[pname]['mode'])
                elif method in param_summary[pname]:
                    truths.append(param_summary[pname][method])
                elif method == 'mode':
                    # Mode requested but not available
                    if not mode_fallback_warned:
                        print("\nWarning: 'mode' not available. Run process_posterior_modes() first.")
                        print("Falling back to 'median' for corner plot.")
                        mode_fallback_warned = True
                    truths.append(param_summary[pname]['median'])
                else:
                    # Fallback to median
                    truths.append(param_summary[pname]['median'])
            else:
                # Parameter not found anywhere - use 0
                truths.append(0.0)
        
        truth_label = f"Best Fit ({method.capitalize()})"
    else:
        truths = [truths[i] for i in free_idx]
        truth_label = "True Values"
    
    # Resample to get equal-weighted samples
    print(f"   Resampling posterior for corner plot...")
    samples_full = dyfunc.resample_equal(samples, weights)
    samples_eq = samples_full[:, free_idx]
    print(f"   Resampled to {len(samples_eq)} equal-weighted samples")
    
    # Compute quantile-based ranges with padding, clipped to actual sample range
    ranges = []
    for i in range(samples_eq.shape[1]):
        q_low, q_high = np.percentile(samples_eq[:, i], [1, 99])
        span = q_high - q_low
        padding = 0.1 * span
        
        range_low = q_low - padding
        range_high = q_high + padding
        
        # Clip to actual sample range
        sample_min = np.min(samples_eq[:, i])
        sample_max = np.max(samples_eq[:, i])
        range_low = max(range_low, sample_min)
        range_high = min(range_high, sample_max)
        
        ranges.append((range_low, range_high))
    
    # Corner plot configuration
    corner_kwargs = {
        'labels': labels,
        'label_kwargs': {'fontsize': 14},  # Increased for publication
        'title_kwargs': {'fontsize': 10},   # Increased for publication
        'quantiles': [0.16, 0.5, 0.84],
        'show_titles': False,
        'smooth': 2.5,
        'bins': 50,
        'hist_kwargs': {'alpha': 0.8, 'color': 'steelblue'},
        'color': 'steelblue',
        'plot_datapoints': True,
        'data_kwargs': {'alpha': 0.1},
        'plot_density': True,
        'contour_kwargs': {'colors': 'darkblue', 'linewidths': 1.2},
        'contourf_kwargs': {'colors': ['lightblue', 'steelblue'], 'alpha': 0.6},
        'range': ranges,
        'truths': truths,
        'truth_color': 'red',
        'truth_kwargs': {'linewidth': 2, 'alpha': 0.8}
    }
    
    # Create corner plot
    fig = corner.corner(samples_eq, **corner_kwargs)
    
    # Add confidence interval shading and custom titles to 1D histograms
    ndim = len(labels)
    axes = np.array(fig.axes).reshape((ndim, ndim))
    
    for i in range(ndim):
        ax = axes[i, i]
        
        # Get parameter name
        param_name = param_names_plot[i]
        
        # Get stats from processed_parameters if available, otherwise param_summary
        if processed_parameters is not None and param_name in processed_parameters:
            proc_param = processed_parameters[param_name]
            median = proc_param['median']
            eti_68 = proc_param['eti_68']
        elif param_name in param_summary:
            summary = param_summary[param_name]
            median = summary['median']
            eti_68 = summary['eti_68']
        else:
            continue  # Skip if no stats available
        
        lower = eti_68[0]
        upper = eti_68[1]
        
        # Compute errors
        lower_err = median - lower
        upper_err = upper - median
        
        # Determine adaptive precision
        precision = adaptive_decimal_precision(median, lower_err, upper_err)
        
        # Get units for this parameter
        units = get_parameter_units(param_name)
        
        # Format title with parameter label, adaptive precision, and units
        param_label = labels[i].replace('$', '')  # Remove outer $ for embedding
        # title_str = f"${param_label} = {median:.{precision}f}_{{-{lower_err:.{precision}f}}}^{{+{upper_err:.{precision}f}}}{units}$"
        title_str = f"${param_label}$\n${median:.{precision}f}_{{-{lower_err:.{precision}f}}}^{{+{upper_err:.{precision}f}}}{units}$"
        ax.set_title(title_str, fontsize=11)  # Decreased for better spacing
        
        # Shade the confidence interval
        ylim = ax.get_ylim()
        ax.axvspan(lower, upper, alpha=0.2, color='lightblue', zorder=-1)
    
    # Create informative title
    if title is not None:
        # Use provided title with suffix
        title = f"{title}{title_suffix}"
        fig.suptitle(title, fontsize=14, y=0.98)
    
    # Adjust subplot spacing for better readability
    fig.subplots_adjust(top=0.93, bottom=0.08, left=0.08, right=0.95, 
                       wspace=0.05, hspace=0.05)
    
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig


def plot_convergence(results: 'FitResults',
                    title: str = "Convergence Diagnostics",
                    filename: Optional[str] = None,
                    output_dir: Optional[Path] = None,
                    show: bool = False) -> plt.Figure:
    """
    Create convergence diagnostic plots for nested sampling.
    
    Four-panel plot:
    - Top-left: Normalized log-likelihood (Δln(L) with inset showing final 25%)
    - Top-right: Log-evidence (ln(Z) with inset showing final 25%)
    - Bottom-left: Chi-squared evolution (log scale)
    - Bottom-right: Effective sample size evolution
    
    Main plots show full convergence trajectory on linear scale.
    Insets zoom to final 25% of iterations to show convergence detail.
    This handles large dynamic range without log scale complications.
    
    Args:
        results: FitResults object from fitting
        title: Plot title
        filename: If provided, save figure
        output_dir: Directory for saving
        show: If True, display plot
        
    Returns:
        Figure object
    """
    setup_publication_quality()
    
    # Verify results has required chi-squared statistics
    required_attrs = ['chi_squared_mean', 'chi_squared_std', 
                     'chi_squared_mean_history', 'chi_squared_std_history']
    missing = [attr for attr in required_attrs if not hasattr(results, attr) or getattr(results, attr) is None]
    if missing:
        raise AttributeError(
            f"FaradayFitter object missing required attributes: {missing}. "
            f"Run compute_summaries() first."
        )
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Get iteration numbers and filter to last 95%
    iterations_all = np.arange(len(results.sampler.results.logl))
    cutoff_idx = int(0.05 * len(iterations_all))  # Skip first 5%
    iterations = iterations_all[cutoff_idx:]
    
    # ========================================================================
    # Top-left: Normalized log-likelihood with inset (last 95%)
    # ========================================================================
    ax = axes[0, 0]
    logl = results.sampler.results.logl[cutoff_idx:]
    
    # Δln(L) = ln(L) - max(ln(L)) (normalized log-likelihood)
    logl_max = np.max(logl)
    delta_logl = logl - logl_max  # Negative values, 0 at maximum
    
    # Main plot: Full range
    ax.plot(iterations, delta_logl, 'k-', linewidth=1, alpha=0.7)
    ax.axhline(0.0, color='darkgreen', linestyle='--', linewidth=1, 
               alpha=0.5, label='Maximum')
    
    ax.set_ylabel('Δln(L) = ln(L) - ln(L$_{max}$)', fontsize=12)
    ax.set_xlabel('Posterior sample index', fontsize=12)
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title('Normalized Log-Likelihood', fontsize=13)
    add_ticks_all_sides(ax)
    
    # Inset: Final 25% of iterations
    ax_inset = inset_axes(ax, width="40%", height="40%", loc='lower right',
                          borderpad=2.0)
    cutoff_inset = int(0.75 * len(iterations))
    iterations_inset = iterations[cutoff_inset:]
    delta_logl_inset = delta_logl[cutoff_inset:]
    
    ax_inset.plot(iterations_inset, delta_logl_inset, 'k-', linewidth=1, alpha=0.7)
    ax_inset.axhline(0.0, color='darkgreen', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_inset.set_title('Final 25%', fontsize=9)
    ax_inset.grid(True, alpha=0.3)
    ax_inset.tick_params(labelsize=8)
    
    # ========================================================================
    # Top-right: Log-evidence with inset (last 95%)
    # ========================================================================
    ax = axes[0, 1]
    logz = results.sampler.results.logz[cutoff_idx:]
    logzerr = results.sampler.results.logzerr[cutoff_idx:]
    
    # Main plot: ln(Z) evolution
    ax.plot(iterations, logz, 'b-', linewidth=2, label='ln(Z)')
    ax.fill_between(iterations, logz - logzerr, logz + logzerr,
                    alpha=0.3, color='blue', label='±1σ')
    
    ax.set_xlabel('NS iteration index', fontsize=12)
    ax.set_ylabel('Log-Evidence ln(Z)', fontsize=12)
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title('Log-Evidence Evolution', fontsize=13)
    add_ticks_all_sides(ax)
    
    # Inset: Final 25% of iterations
    ax_inset = inset_axes(ax, width="40%", height="40%", loc='lower right',
                          borderpad=2.0)
    cutoff_inset = int(0.75 * len(iterations))
    iterations_inset = iterations[cutoff_inset:]
    logz_inset = logz[cutoff_inset:]
    logzerr_inset = logzerr[cutoff_inset:]
    
    ax_inset.plot(iterations_inset, logz_inset, 'b-', linewidth=2)
    ax_inset.fill_between(iterations_inset, 
                          logz_inset - logzerr_inset, 
                          logz_inset + logzerr_inset,
                          alpha=0.3, color='blue')
    ax_inset.set_title('Final 25%', fontsize=9)
    ax_inset.grid(True, alpha=0.3)
    ax_inset.tick_params(labelsize=8)
    
    # ========================================================================
    # Bottom-left: Chi-squared evolution
    # ========================================================================
    ax = axes[1, 0]
    
    # Plot chi-squared history (already thinned during computation for speed)
    # Need to compute corresponding iteration indices for thinned samples
    n_samples_total = len(results.samples)
    thin_factor = max(1, n_samples_total // 1000)
    iter_thinned = iterations_all[::thin_factor]
    
    # Handle case where thinning might create length mismatch
    n_plot = min(len(iter_thinned), len(results.chi_squared_history))
    
    # Plot importance-weighted mean (bold line)
    ax.plot(iter_thinned[:n_plot], results.chi_squared_mean_history[:n_plot],
            'darkred', linewidth=2.5, label=f'Weighted mean: χ² = {results.chi_squared_mean:.1f}')
    
    # Shade ±1σ region
    mean = results.chi_squared_mean_history[:n_plot]
    std = results.chi_squared_std_history[:n_plot]
    ax.fill_between(iter_thinned[:n_plot], mean - std, mean + std,
                    alpha=0.3, color='red', label=f'±1σ: {results.chi_squared_std:.1f}')
    
    # Plot best sample (horizontal line)
    ax.axhline(results.chi_squared_best, color='blue', linestyle='--', 
               linewidth=2, label=f'Best sample: χ² = {results.chi_squared_best:.1f}')
    
    # Mark expected value (degrees of freedom)
    dof = results.n_data - results.n_params
    ax.axhline(dof, color='green', linestyle=':', linewidth=2,
              label=f'Expected: χ² = {dof} (dof)')
    
    # Use log scale for y-axis
    ax.set_yscale('log')
    
    ax.set_xlabel('Posterior sample index', fontsize=12)
    ax.set_ylabel('Chi-Squared χ² (log scale)', fontsize=12)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Chi-Squared Evolution', fontsize=13)
    add_ticks_all_sides(ax)
    
    # Inset: Final 25% of iterations
    ax_inset = inset_axes(ax, width="40%", height="40%", loc='lower left',
                          borderpad=2.5)
    ax_inset.set_zorder(10)
    
    cutoff_inset = int(0.75 * n_plot)
    iter_inset = iter_thinned[cutoff_inset:n_plot]
    mean_inset = results.chi_squared_mean_history[cutoff_inset:n_plot]
    std_inset = results.chi_squared_std_history[cutoff_inset:n_plot]
    
    # Plot weighted mean and ±1σ band in inset
    ax_inset.plot(iter_inset, mean_inset, 'darkred', linewidth=2)
    ax_inset.fill_between(iter_inset, mean_inset - std_inset, mean_inset + std_inset,
                          alpha=0.3, color='red')
    
    # Plot best sample and DOF lines in inset
    ax_inset.axhline(results.chi_squared_best, color='blue', linestyle='--', 
                     linewidth=1.5, alpha=0.7)
    ax_inset.axhline(dof, color='green', linestyle=':', linewidth=1, alpha=0.5)
    
    # Linear scale with ±10% padding around data range
    # Include both mean±std and best sample in range calculation
    # Filter NaN values (from early overflow samples)
    lower_band = mean_inset - std_inset
    upper_band = mean_inset + std_inset
    
    # Get finite range from bands
    finite_lower = lower_band[np.isfinite(lower_band)]
    finite_upper = upper_band[np.isfinite(upper_band)]
    
    values_for_range = []
    if len(finite_lower) > 0:
        values_for_range.append(np.min(finite_lower))
    if len(finite_upper) > 0:
        values_for_range.append(np.max(finite_upper))
    if np.isfinite(results.chi_squared_best):
        values_for_range.append(results.chi_squared_best)
    
    if len(values_for_range) >= 2:
        data_min = min(values_for_range)
        data_max = max(values_for_range)
        data_range = data_max - data_min
        ax_inset.set_ylim(data_min - 0.1*data_range, data_max + 0.1*data_range)
    else:
        # Fallback to DOF range
        ax_inset.set_ylim(dof * 0.8, dof * 1.2)
    
    # Move y-axis labels to right side
    ax_inset.yaxis.tick_right()
    ax_inset.yaxis.set_label_position('right')
    
    ax_inset.set_title('Final 25%', fontsize=9)
    ax_inset.grid(True, alpha=0.3)
    ax_inset.tick_params(labelsize=8)
    
    # ========================================================================
    # Bottom-right: Effective sample size (Kish ESS)
    # ========================================================================
    ax = axes[1, 1]
    
    # Compute running Kish ESS
    weights = results.weights
    n_samples = len(weights)
    ess_running = np.zeros(n_samples)
    
    W = 0.0
    Q = 0.0
    for i in range(n_samples):
        w = weights[i]
        if np.isfinite(w) and w > 0:
            W += w
            Q += w * w
        
        if Q > 0:
            ess_running[i] = (W * W) / Q
        else:
            ess_running[i] = 0.0
    
    # Posterior sample indices
    sample_indices = np.arange(n_samples)
    
    # Plot running ESS
    ax.plot(sample_indices, ess_running, 'g-', linewidth=2, label='Running ESS')
    
    # Add horizontal line at true dynesty n_effective
    ax.axhline(results.n_eff, color='darkgreen', linestyle='--', linewidth=2,
               label=f'Dynesty n_effective: {results.n_eff:.1f}')
    
    ax.set_xlabel('Posterior sample index', fontsize=12)
    ax.set_ylabel('Approximate effective sample size', fontsize=12)
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_title('Effective Sample Size Evolution', fontsize=13)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Overall formatting
    # ========================================================================
    fig.suptitle(title, fontsize=14, fontweight='bold')
    fig.subplots_adjust(hspace=0.35, wspace=0.30, top=0.92, bottom=0.08, left=0.08, right=0.95)
    
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig


def plot_fdf_diagnostic(results: 'FitResults',
                       phi_range: Optional[float] = None,
                       n_phi: int = 2048,
                       weight_type: str = 'variance',
                       model_method: str = 'representative',
                       plot_posterior_samples: bool = False,
                       title: str = "FDF Diagnostic",
                       filename: Optional[str] = None,
                       output_dir: Optional[Path] = None,
                       show: bool = False) -> plt.Figure:
    """
    Create FDF diagnostic plot showing RMTF and observed/model FDFs.
    
    Three-panel plot:
    - Top: RMTF (|R(φ)|, Re, Im components) with FWHM annotation
    - Middle: Observed FDF from data and model FDF from RM synthesis
    - Bottom: Intrinsic model component amplitudes (zoomed)
    
    The observed FDF is computed from actual Q,U data using RM synthesis.
    The model FDF is computed by applying RM synthesis to model Q,U predictions.
    
    Args:
        results: FitResults object from fitting
        phi_range: Range in Faraday depth (±phi_range). If None, auto-computed as
                   max(10×FWHM, max_RM+2×FWHM) to cover full parameter space
        n_phi: Number of points in Faraday depth grid
        weight_type: Weighting for RM synthesis ('variance' or 'uniform')
        model_method: Method for model parameters ('representative', 'median', 'mean', 'map', 'mode')
                     'representative' (default): High-likelihood sample near median (accounts for covariances)
                     'mode': Falls back to median if mode not available
                     Other methods: Use summary statistic (may not correspond to high-likelihood sample)
        plot_posterior_samples: If True, plot 100 random samples from high-likelihood region
                               with low alpha to show posterior uncertainty
        title: Plot title
        filename: If provided, save figure
        output_dir: Directory for saving
        show: If True, display plot
        
    Returns:
        Figure object
    """
    setup_publication_quality()
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.05)  # Reduced spacing between subplots
    
    # Get data and model based on method
    data = results.setup.data
    lambda_sq = data.lambda_sq
    
    if model_method == 'representative':
        # Get representative sample from top 1% likelihood near median
        representative_params = results.get_representative_sample(percentile=1.0, method='median')
        best_model = copy.deepcopy(results.setup.model)
        update_model_from_params(best_model, representative_params, results.setup.param_names, lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))
        method_label = "Median-like sample"
    else:
        # Use summary statistic method
        best_model = results.get_best_fit_model(method=model_method)
        method_label = f"{model_method.capitalize()} Parameters"
    
    # Auto-compute phi_range if not provided
    if phi_range is None:
        # Get FWHM from RMTF properties
        _, rmtf_props_temp = compute_rmtf(lambda_sq, np.array([0.0]), 
                                          weight_type=weight_type,
                                          Q_err=data.Q_err, U_err=data.U_err)
        fwhm = rmtf_props_temp['expected_fwhm_nominal']
        
        # Option 1: 10x FWHM
        phi_range_fwhm = 10.0 * fwhm
        
        # Option 2: Extract RM values from fitted model
        rm_values = []
        for comp in best_model.components:
            if hasattr(comp, 'phi_rm'):
                rm_values.append(comp.phi_rm)
            elif hasattr(comp, 'phi_peak'):
                rm_values.append(comp.phi_peak)
        
        if rm_values:
            rm_max = np.max(rm_values)
            rm_min = np.min(rm_values)
            phi_range_comps = max(abs(rm_max + 2.0 * fwhm), abs(rm_min - 2.0 * fwhm))
        else:
            phi_range_comps = 0.0
        
        phi_range = max(phi_range_fwhm, phi_range_comps, 100.0)
        print(f"   Auto-computed phi_range = {phi_range:.1f} rad/m² (FWHM = {fwhm:.1f})")
    
    # Create Faraday depth grid for display
    phi_grid = np.linspace(-phi_range, phi_range, n_phi)
    
    # Use broader grid for computation to avoid edge effects
    # Compute on 5x wider range, then extract central region
    phi_grid_wide = np.linspace(-5*phi_range, 5*phi_range, 5*n_phi)
    
    # Compute RMTF on wide grid
    rmtf_wide, rmtf_props = compute_rmtf(lambda_sq, phi_grid_wide, 
                                         weight_type=weight_type,
                                         Q_err=data.Q_err, U_err=data.U_err)
    
    # Compute intrinsic model FDF on WIDE grid (infinite bandwidth limit)
    model_fdf_wide = best_model.compute_fdf(phi_grid_wide)
    
    # Extract central regions for display
    n_start = 2*n_phi  # Start of central region
    n_end = 3*n_phi    # End of central region
    rmtf = rmtf_wide[n_start:n_end]
    model_fdf = model_fdf_wide[n_start:n_end]
    
    # Get FWHM for annotation (use full resolution FWHM from Rudnick & Cotton 2023)
    fwhm = rmtf_props['expected_fwhm_full']
    
    # Compute observed FDF from data using RM synthesis (Rudnick & Cotton 2023, λ₀²=0)
    observed_fdf_data_wide = rm_synthesis(data, phi_grid_wide, weight_type=weight_type)
    observed_fdf_data = observed_fdf_data_wide[n_start:n_end]
    
    # ========================================================================
    # Top panel: RMTF
    # ========================================================================
    ax = axes[0]
    
    # Actual RMTF
    ax.plot(phi_grid, np.abs(rmtf), 'k-', linewidth=2, label='|R(φ)| Actual')
    ax.plot(phi_grid, np.real(rmtf), 'b--', linewidth=1.5, label='Re[R(φ)]', alpha=0.7)
    ax.plot(phi_grid, np.imag(rmtf), 'r--', linewidth=1.5, label='Im[R(φ)]', alpha=0.7)
    
    ax.axhline(0, color='gray', linestyle='--', linewidth=1.0)
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    
    # FWHM annotation removed per user request
    
    ax.set_ylabel('RMTF', fontsize=12)
    ax.tick_params(axis='both', labelsize=12, labelbottom=False)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Middle panel: FDF (observed from data, model from RM synthesis)
    # ========================================================================
    ax = axes[1]
    
    # Compute model FDF via RM synthesis
    P_model = best_model.compute_polarization(lambda_sq)
    
    model_data_temp = PolarizationData(lambda_sq=lambda_sq, Q=np.real(P_model), 
                                       U=np.imag(P_model), Q_err=data.Q_err, U_err=data.U_err)
    observed_fdf_model_wide = rm_synthesis(model_data_temp, phi_grid_wide, weight_type=weight_type)
    observed_fdf_model = observed_fdf_model_wide[n_start:n_end]
    
    # Plot observed FDF from data
    ax.plot(phi_grid, np.abs(observed_fdf_data), 'k-', linewidth=2.5, 
           label='|F(φ)| Observed (from data)', zorder=10, alpha=0.8)
    
    # Plot model FDF from RM synthesis (representative sample)
    ax.plot(phi_grid, np.abs(observed_fdf_model), color='purple', linewidth=2.5, 
           label=f'|F(φ)| Model FDF ({method_label})', zorder=11, alpha=0.8)
    
    # Real/Imaginary components of model FDF
    ax.plot(phi_grid, np.real(observed_fdf_model), 'b--', linewidth=1.5, 
           label='Re[F(φ)] Model FDF', alpha=0.4)
    ax.plot(phi_grid, np.imag(observed_fdf_model), 'r--', linewidth=1.5, 
           label='Im[F(φ)] Model FDF', alpha=0.4)
    
    ax.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    
    ax.set_ylabel('|F(φ)| (fractional)', fontsize=12)
    ax.tick_params(axis='both', labelsize=12, labelbottom=False)
    ax.legend(loc='best', framealpha=0.9, fontsize=9)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)
    
    # ========================================================================
    # Bottom panel: Individual component FDFs on custom grids
    # ========================================================================
    ax_bottom = axes[2]
    
    # Plot posterior samples if requested (intrinsic model FDFs at back)
    if plot_posterior_samples:
        # Get top weighted samples (highest likelihoods)
        n_samples_total = len(results.samples)
        weights_normalized = results.weights / np.sum(results.weights)
        
        # Sample 100 points from weighted posterior
        n_samples_plot = min(100, n_samples_total)
        sample_indices = np.random.choice(n_samples_total, 
                                         size=n_samples_plot,
                                         replace=False,
                                         p=weights_normalized)
        
        # Plot each sample's intrinsic FDF per component with low alpha
        for idx in sample_indices:
            params = results.samples[idx]
            model_sample = copy.deepcopy(results.setup.model)
            update_model_from_params(model_sample, params, results.setup.param_names, lambda_sq_ref=getattr(results.setup, "lambda_sq_ref", None))

            for comp in model_sample.components:
                comp_fdf_sample = comp.compute_fdf(phi_grid)
                ax_bottom.plot(phi_grid, np.abs(comp_fdf_sample), 'C0-', linewidth=0.5,
                               alpha=0.15, zorder=1)
    
    # Plot each component on the same grid as the samples
    max_amp_components = 0.0
    
    for i, comp in enumerate(best_model.components):
        # Compute FDF for this component on the same grid as samples (phi_grid)
        comp_fdf = comp.compute_fdf(phi_grid)
        comp_amp = np.abs(comp_fdf)
        
        # Track maximum amplitude across all components
        max_amp_components = max(max_amp_components, np.max(comp_amp))
        
        # Plot component FDF as solid black line
        ax_bottom.plot(phi_grid, comp_amp, 'k-', linewidth=2.0, alpha=1.0, zorder=12)
        
        # Add arrow marker at peak for ThinComponents
        if isinstance(comp, ThinComponent):
            peak_amp = np.max(comp_amp)
            ax_bottom.plot(comp.phi_rm, peak_amp, marker='^', 
                          markersize=10, color='black', zorder=15)
    
    # Match x-axis to middle panel
    ax_bottom.set_xlim(axes[1].get_xlim())
    
    # Set y-axis limits
    # Upper: 10x largest component amplitude, capped at 1.0
    # Lower: 1e-5 (log scale minimum)
    ax_bottom.set_yscale('log')
    if max_amp_components > 0:
        y_upper = min(10.0 * max_amp_components, 1.0)
        ax_bottom.set_ylim(1e-5, y_upper)
    else:
        ax_bottom.set_ylim(1e-5, 1e-3)
    
    ax_bottom.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    
    ax_bottom.set_xlabel('Faraday Depth φ (rad/m²)', fontsize=12)
    ax_bottom.set_ylabel('|F(φ)| (fractional)', fontsize=12)
    ax_bottom.tick_params(axis='both', labelsize=12)
    # Legend removed for cleaner appearance
    ax_bottom.grid(True, alpha=0.3, which='both')
    add_ticks_all_sides(ax_bottom)
    
    # ========================================================================
    # Overall formatting
    # ========================================================================
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig


def plot_flagging_diagnostic(raw_data: Dict, 
                            flagging_info: Dict,
                            title: str = "Channel Flagging Diagnostic",
                            filename: Optional[str] = None,
                            output_dir: Optional[Path] = None,
                            show: bool = False) -> plt.Figure:
    """
    Diagnostic plot showing raw vs flagged data.
    
    Creates multi-panel plot with I, Q, U (+ V if present) showing:
    - Raw data (all channels, grayed out, 'x' markers)
    - Kept data (unflagged, normal colors, 'o' markers)  
    - Explicitly flagged frequency regions (RFI bands) highlighted with vertical bands
    - Other individual flags shown as 'x' markers only (no shading)
    
    Parameters
    ----------
    raw_data : dict
        Original data before flagging. Must contain:
        - 'freq': Frequencies in Hz
        - 'I', 'Q', 'U': Stokes parameters
        - Optional: 'V'
    flagging_info : dict
        From flag_channels() containing:
        - 'mask': Boolean array (True = flagged)
        - 'kept_indices': Indices of kept channels
        - 'flag_counts': Flagging statistics
        - 'n_total', 'n_flagged', 'n_kept'
        - Optional: 'flagged_regions': List of explicitly flagged frequency regions
          Each region is a dict with either:
            * 'lambda_sq_range': (min, max) in m²
            * 'freq_range': (min, max) in Hz
    title : str
        Plot title
    filename : str or None
        If provided, save figure
    output_dir : Path or None
        Directory for saving
    show : bool
        If True, display plot
        
    Returns
    -------
    fig : matplotlib Figure
    """
    setup_publication_quality()
    
    # Extract data
    freq_hz = np.asarray(raw_data['freq'])
    freq_ghz = freq_hz / 1e9
    
    # Compute lambda² (m²) for bottom x-axis
    lambda_sq = (C_LIGHT / freq_hz)**2
    
    mask = flagging_info['mask']
    kept_idx = flagging_info['kept_indices']
    
    # Determine which Stokes parameters to plot
    has_V = 'V' in raw_data
    stokes_params = ['I', 'Q', 'U', 'V'] if has_V else ['I', 'Q', 'U']
    n_panels = len(stokes_params)
    
    # Create figure
    fig, axes = plt.subplots(n_panels, 1, figsize=(14, 3.5 * n_panels))
    if n_panels == 1:
        axes = [axes]
    
    # Color scheme
    colors = {'I': 'black', 'Q': 'blue', 'U': 'red', 'V': 'green'}
    
    for idx, stokes in enumerate(stokes_params):
        ax = axes[idx]
        
        data = np.asarray(raw_data[stokes])
        
        # Get error data if available
        err_key = f'd{stokes}'
        if err_key in raw_data:
            data_err = np.asarray(raw_data[err_key])
        else:
            data_err = None
        
        # Plot all raw data (grayed out, 'x' markers, with error bars)
        if data_err is not None:
            ax.errorbar(lambda_sq, data, yerr=data_err,
                       marker='x', linestyle='', capsize=0,
                       color='gray', alpha=0.85, markersize=6, elinewidth=1,
                       label='Flagged', zorder=1)
        else:
            ax.plot(lambda_sq, data, marker='x', linestyle='', 
                   color='gray', alpha=0.85, markersize=6, 
                   label='Flagged', zorder=1)
        
        # Plot kept data (normal colors, 'o' markers, with error bars)
        if data_err is not None:
            ax.errorbar(lambda_sq[kept_idx], data[kept_idx], yerr=data_err[kept_idx],
                       marker='o', linestyle='', capsize=2,
                       color=colors[stokes], alpha=0.8, markersize=5, elinewidth=1.5,
                       label='Kept', zorder=3)
        else:
            ax.plot(lambda_sq[kept_idx], data[kept_idx], marker='o', linestyle='',
                   color=colors[stokes], alpha=0.8, markersize=5,
                   label='Kept', zorder=3)
        
        # Plot spectral baseline polynomial fit (for any Stokes parameter that was fit)
        if 'Spectral baseline clipping' in flagging_info.get('flag_counts', {}):
            spectral_clip_info = flagging_info['flag_counts']['Spectral baseline clipping']
            if 'fit' in spectral_clip_info and spectral_clip_info['fit'] is not None:
                fits_dict = spectral_clip_info['fit']
                
                # Check if this Stokes parameter has a fit
                if stokes in fits_dict:
                    fit = fits_dict[stokes]
                    # Evaluate polynomial: X = poly(freq_normalized)
                    freq_norm = (freq_hz - fit['freq_mean']) / fit['freq_std']
                    X_model = np.polyval(fit['coef'], freq_norm)
                    ax.plot(lambda_sq, X_model, '-', color='orange', linewidth=2.5, 
                           alpha=0.8, label=f'Spectral baseline ({stokes})', zorder=2)
        
        # Highlight explicitly flagged frequency regions (RFI bands specified by user)
        # Individual channel flags are shown as 'x' markers only, no shading
        if 'flagged_regions' in flagging_info and flagging_info['flagged_regions'] is not None:
            # Only shade explicitly defined frequency/wavelength regions
            for region in flagging_info['flagged_regions']:
                if 'lambda_sq_range' in region:
                    # Region specified in lambda_sq space
                    lsq_min, lsq_max = region['lambda_sq_range']
                elif 'freq_range' in region:
                    # Region specified in frequency space - convert to lambda_sq
                    freq_min, freq_max = region['freq_range']
                    lsq_max = (C_LIGHT / freq_min)**2  # Higher freq → lower lambda_sq
                    lsq_min = (C_LIGHT / freq_max)**2  # Lower freq → higher lambda_sq
                else:
                    continue  # Skip if no range specified
                
                ax.axvspan(lsq_min, lsq_max, alpha=0.15, color='red', 
                          zorder=0, linewidth=0)
        
        # Formatting
        ax.set_ylabel(f'Stokes {stokes} (Jy)', fontsize=12)
        # Only show legend in top panel (Stokes I)
        if idx == 0:
            ax.legend(loc='best', framealpha=0.9, fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Set y-axis limits to zoom on unflagged data with 25% padding
        kept_data = data[kept_idx]
        data_min = kept_data.min()
        data_max = kept_data.max()
        data_range = data_max - data_min
        ax.set_ylim(data_min - 0.25*data_range, data_max + 0.25*data_range)
        
        # Add horizontal line at 0 for Q, U, V (not I)
        if stokes in ['Q', 'U', 'V']:
            ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=2)
        
        # Set xlabel and xlim with 5% padding on each side
        ax.set_xlabel('λ² (m²)', fontsize=12)
        lsq_range = lambda_sq.max() - lambda_sq.min()
        ax.set_xlim(lambda_sq.min() - 0.05*lsq_range, lambda_sq.max() + 0.05*lsq_range)
        add_ticks_all_sides(ax)
        
        # Add frequency labels to TOP ticks (only for top panel)
        if idx == 0:
            # Get tick positions
            xticks = ax.get_xticks()
            freq_labels = [f'{lambda_sq_to_freq_ghz(np.array([t]))[0]:.2f}' for t in xticks]
            
            # Create secondary top axis with frequency labels
            ax_top = ax.secondary_xaxis('top')
            ax_top.set_xticks(xticks)
            ax_top.set_xticklabels(freq_labels)
            ax_top.set_xlabel('Frequency (GHz)', fontsize=12)
            ax_top.tick_params(which='both', direction='in')
            ax_top.xaxis.set_minor_locator(AutoMinorLocator())
        
        # Hide xlabel on non-bottom panels
        if idx != n_panels - 1:
            ax.xaxis.label.set_visible(False)
    
    # Overall title with flagging summary
    if title is not None:
        n_total = flagging_info['n_total']
        n_flagged = flagging_info['n_flagged']
        n_kept = flagging_info['n_kept']
        pct_flagged = 100 * n_flagged / n_total if n_total > 0 else 0
        
        summary = (f"{title}\n"
                  f"Total: {n_total} channels  |  "
                  f"Flagged: {n_flagged} ({pct_flagged:.1f}%)  |  "
                  f"Kept: {n_kept} ({100-pct_flagged:.1f}%)")
        
        fig.suptitle(summary, fontsize=14, fontweight='bold')
    
    # Note: tight_layout is not used as it can cause axis misalignment
    
    if filename:
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig





def plot_thick_component_band_response(results: 'FitResults',
                                       n_samples: int = 1000,
                                       filename: Optional[str] = None,
                                       output_dir: Optional[Path] = None,
                                       show: bool = False) -> Optional[plt.Figure]:
    """
    Plot posterior uncertainty on per-component polarization amplitude across
    the observed band for ThickComponents.

    For each ThickComponent, draws weighted posterior samples and computes the
    median and 68% credible interval of |P(λ²)| as a function of λ². A vertical
    line marks λ²_ref so you can verify it sits in a region of reasonable signal
    for every component.

    Returns None without plotting if the model contains no ThickComponents.

    Args:
        results: FitResults object from fitting
        n_samples: Maximum number of posterior samples to draw (default: 1000)
        filename: If provided, save figure to this filename
        output_dir: Directory for saved figure (default: current directory)
        show: If True, display figure interactively

    Returns:
        matplotlib Figure, or None if no ThickComponents present
    """
    from faraday_utils import ThickComponent

    model = results.setup.model
    thick_components = [(i, comp) for i, comp in enumerate(model.components)
                        if isinstance(comp, ThickComponent)]

    if not thick_components:
        print("No ThickComponents in model — skipping band response plot.")
        return None

    data = results.setup.data
    lambda_sq = data.lambda_sq
    lambda_sq_ref = getattr(results.setup, 'lambda_sq_ref', None)

    lambda_sq_fine = np.linspace(lambda_sq.min(), lambda_sq.max(), 500)

    n_total = len(results.samples)
    weights_norm = results.weights / np.sum(results.weights)
    n_draw = min(n_samples, n_total)
    sample_indices = np.random.choice(n_total, size=n_draw, replace=False, p=weights_norm)

    n_thick = len(thick_components)
    amp_arrays = [np.zeros((n_draw, 500)) for _ in range(n_thick)]

    model_sample = copy.deepcopy(model)
    for draw_idx, sample_idx in enumerate(sample_indices):
        params = results.samples[sample_idx]
        update_model_from_params(model_sample, params, results.setup.param_names,
                                 lambda_sq_ref=lambda_sq_ref)
        for arr_idx, (comp_idx, _) in enumerate(thick_components):
            comp = model_sample.components[comp_idx]
            amp_arrays[arr_idx][draw_idx] = np.abs(comp.compute_polarization(lambda_sq_fine))

    colours = [f'C{i}' for i in range(n_thick)]
    fig, ax = plt.subplots(figsize=(10, 5))

    for arr_idx, (comp_idx, ref_comp) in enumerate(thick_components):
        colour = colours[arr_idx]
        amps = amp_arrays[arr_idx]
        median = np.median(amps, axis=0)
        lo = np.percentile(amps, 16, axis=0)
        hi = np.percentile(amps, 84, axis=0)
        comp_name = ref_comp.name if ref_comp.name else f'component_{comp_idx}'
        ax.plot(lambda_sq_fine, median, color=colour, linewidth=2.0, label=comp_name, zorder=5)
        ax.fill_between(lambda_sq_fine, lo, hi, color=colour, alpha=0.25, zorder=3)

    if lambda_sq_ref is not None:
        ax.axvline(lambda_sq_ref, color='black', linestyle='--', linewidth=1.5,
                   label=f'λ²_ref = {lambda_sq_ref:.3f} m²', zorder=10)

    ax.set_xlabel('λ² (m²)', fontsize=12)
    ax.set_ylabel('|P (λ²)| (fractional)', fontsize=12)
    ax.tick_params(axis='both', labelsize=11)
    ax.legend(loc='best', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3)
    add_ticks_all_sides(ax)

    fig.tight_layout()

    if filename:
        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {filepath}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig
