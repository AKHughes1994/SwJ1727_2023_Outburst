# faraday_processing.py
"""
Posterior mode detection and processing for Faraday rotation analysis.

This module processes posterior samples in Cartesian (a, b) parameterization
to compute modes, credible intervals, and derived polar coordinates. Handles
circular/axial angle topology correctly via mode-conditioned unwrapping.

Physical conventions:
- Complex polarization: P = Q + iU = p*I*exp(2i*psi)
- Cartesian intrinsic: (a, b) where P₀ = a + ib
- Derived: p₀ = √(a² + b²), φ = atan2(b, a), ψ₀ = φ/2
- φ is circular (2π-periodic), ψ₀ is axial (π-periodic)

Mode detection:
- Linear parameters (p₀, φ_rm, σ_φ, N): Gaussian KDE
- Axial angles (ψ₀): Detected via circular φ-space, then converted

Multi-modal detection:
- Finds all local maxima in KDE above height threshold
- Assigns samples to nearest mode for linear parameters (Voronoi)
- Merges nearby modes within distance threshold
- Filters modes by posterior mass threshold
- Ranks modes by posterior mass or likelihood
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Union, Literal
from pathlib import Path
from scipy.stats import gaussian_kde
from scipy.special import i0, i0e
from scipy.signal import find_peaks, peak_prominences
from dynesty.utils import resample_equal, quantile as dyfunc_quantile
from dynesty import utils as dyfunc
from faraday_utils import ThinComponent, ThickComponent, ThinPowerLawComponent


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class ProcessingConfig:
    """
    Configuration for posterior mode processing.
    
    Attributes:
        detect_multiple_modes: Enable multi-modal detection (finds all significant modes)
        mode_height_threshold: Minimum peak height for initial detection (fraction of max peak, 0-1)
        mode_mass_threshold: Minimum posterior mass to qualify as mode (fraction, 0-1)
        mode_merge_fraction: Merge distance as fraction of 99% sample range (default 0.05 = 5%)
        mode_merge_strategy: How to merge nearby modes ('keep_most_mass' or 'keep_highest')
        mode_valley_ratio: Ratio for valley depth checking (valley must be 1/ratio of smaller peak height)
        mode_ranking: How to rank modes ('posterior_mass' or 'likelihood')
        max_modes: Maximum number of modes to report (None = report all above threshold)
        kde_bandwidth: KDE bandwidth selection ('scott', 'silverman', 'loo', or float)
        eti_levels: Credible interval levels to compute (default: 68%, 99%)
        plot_diagnostics: Generate diagnostic plots showing KDE, modes, HDI
        plot_dir: Directory for diagnostic plots
        filename_prefix: Optional prefix for plot filenames
        grid_n: Number of grid points for linear parameter KDE evaluation
        circular_grid_n: Number of grid points for circular/axial parameter KDE
        verbose: Print detailed diagnostic output for mode detection workflow
    """
    detect_multiple_modes: bool = True
    mode_height_threshold: float = 0.01
    mode_mass_threshold: float = 0.1
    mode_merge_fraction: float = 0.05
    mode_merge_strategy: Literal['keep_most_mass', 'keep_highest'] = 'keep_highest'
    mode_valley_ratio: float = 2.0
    mode_ranking: Literal['posterior_mass', 'likelihood'] = 'posterior_mass'
    max_modes: Optional[int] = None
    kde_bandwidth: Union[str, float] = 'scott'
    eti_levels: List[float] = field(default_factory=lambda: [0.68, 0.99])
    plot_diagnostics: bool = True
    plot_dir: Path = Path('plots/processing')
    filename_prefix: str = ''
    grid_n: int = 4000
    circular_grid_n: int = 4000
    verbose: bool = False



# =============================================================================
# WEIGHTED STATISTICS UTILITIES
# =============================================================================

def normalize_weights(weights: np.ndarray) -> np.ndarray:
    """
    Normalize weights to sum to 1.
    
    Args:
        weights: Array of importance weights
        
    Returns:
        Normalized weights
        
    Raises:
        ValueError: If all weights are zero or negative
    """
    weights = np.asarray(weights, dtype=float)
    weights = np.clip(weights, 0.0, np.inf)
    total = weights.sum()
    
    if total <= 0:
        raise ValueError("All weights are zero or negative")
    
    return weights / total


def weighted_mean_std(x: np.ndarray, weights: np.ndarray) -> Tuple[float, float]:
    """
    Compute weighted mean and standard deviation.
    
    Args:
        x: Sample values
        weights: Normalized weights
        
    Returns:
        (mean, std)
    """
    mean = np.average(x, weights=weights)
    variance = np.average((x - mean)**2, weights=weights)
    std = np.sqrt(variance)
    return float(mean), float(std)


def weighted_median(x: np.ndarray, weights: np.ndarray) -> float:
    """
    Compute weighted median.
    
    Args:
        x: Sample values
        weights: Weights (will be normalized internally)
        
    Returns:
        Weighted median value
    """
    x = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    
    # Normalize weights internally
    weight_sum = np.sum(weights)
    if weight_sum <= 0:
        raise ValueError("weights must sum to positive value")
    weights = weights / weight_sum
    
    # Sort by x
    idx = np.argsort(x)
    x_sorted = x[idx]
    w_sorted = weights[idx]
    
    # Cumulative distribution
    cdf = np.cumsum(w_sorted)
    
    # Find median (50th percentile)
    median_idx = np.searchsorted(cdf, 0.5)
    
    return float(x_sorted[median_idx])


def weighted_hdi(x: np.ndarray, weights: np.ndarray, 
                 cred_mass: float = 0.68) -> Tuple[float, float]:
    """
    Compute weighted Highest Density Interval (shortest credible interval).
    
    Finds the shortest interval [x_low, x_high] that contains the specified
    credible mass of the posterior. For unimodal distributions, this is the
    narrowest interval with the desired probability mass.
    
    Args:
        x: Sample values (1D array)
        weights: Weights (will be normalized internally)
        cred_mass: Credible mass (default: 0.68 for ~1σ)
        
    Returns:
        (lower_bound, upper_bound)
        
    Raises:
        ValueError: If arrays have incompatible shapes
    """
    x = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    
    if x.ndim != 1 or weights.ndim != 1:
        raise ValueError("x and weights must be 1D arrays")
    if x.size != weights.size:
        raise ValueError("x and weights must have same length")
    
    # Normalize weights internally
    weight_sum = np.sum(weights)
    if weight_sum <= 0:
        raise ValueError("weights must sum to positive value")
    weights = weights / weight_sum
    
    # Sort by x values
    idx = np.argsort(x)
    x_sorted = x[idx]
    w_sorted = weights[idx]
    
    # Cumulative distribution
    cdf = np.cumsum(w_sorted)
    
    # Find shortest interval containing cred_mass
    best_width = np.inf
    best_i, best_j = 0, 0
    
    for i in range(len(x_sorted)):
        mass_before = cdf[i - 1] if i > 0 else 0.0
        target_mass = mass_before + cred_mass
        
        if target_mass > 1.0:
            break
        
        # Find right endpoint
        j = np.searchsorted(cdf, target_mass, side='left')
        if j >= len(x_sorted):
            j = len(x_sorted) - 1
        
        # Check width
        width = x_sorted[j] - x_sorted[i]
        if width < best_width:
            best_width = width
            best_i = i
            best_j = j
    
    return float(x_sorted[best_i]), float(x_sorted[best_j])


# =============================================================================
# LINEAR PARAMETER MODE DETECTION
# =============================================================================

def find_mode_linear(x: np.ndarray, weights: np.ndarray,
                    bandwidth: Union[str, float] = 'scott',
                    grid_n: int = 4000) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Find mode of weighted samples using Gaussian KDE.
    
    Uses scipy.stats.gaussian_kde with weighted samples. Mode is defined
    as the maximum of the KDE on a dense grid.
    
    Args:
        x: Sample values (1D)
        weights: Normalized weights
        bandwidth: 'scott', 'silverman', or float
        grid_n: Number of grid points for KDE evaluation
        
    Returns:
        mode: Position of KDE maximum
        grid: Grid points where KDE was evaluated
        pdf: KDE values on grid
    """
    x = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    
    # Check if parameter has zero variance (should not happen since fixed params are excluded)
    x_std = np.std(x)
    x_range = np.ptp(x)
    
    if x_range < 1e-12 or x_std < 1e-12:
        # Return delta function
        mode = float(np.mean(x))
        grid = np.array([mode - 1e-6, mode, mode + 1e-6])
        pdf = np.array([0.0, 1.0, 0.0])
        return mode, grid, pdf
    
    # Create weighted KDE using scipy
    if isinstance(bandwidth, str):
        kde = gaussian_kde(x, bw_method=bandwidth, weights=weights)
    else:
        kde = gaussian_kde(x, bw_method=float(bandwidth), weights=weights)
    
    # Create evaluation grid using weighted quantiles to avoid outliers
    x_min, x_max = dyfunc_quantile(x, [0.001, 0.999], weights=weights)
    span = x_max - x_min
    padding = 0.15 * span  # 15% padding for smoother edges
    grid = np.linspace(x_min - padding, x_max + padding, grid_n)
    
    # Evaluate KDE
    pdf = kde(grid)
    
    # Find mode (maximum)
    mode_idx = np.argmax(pdf)
    mode = float(grid[mode_idx])
    
    return mode, grid, pdf


# =============================================================================
# MULTI-MODAL DETECTION UTILITIES
# =============================================================================

def find_kde_peaks(grid: np.ndarray, pdf: np.ndarray, 
                   height_threshold: float = 0.05,
                   valley_ratio: float = 3.0,
                   skip_valley_check: bool = False,
                   verbose: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Find significant local maxima in KDE probability density function.
    
    Requires peaks to have both:
    - Height above threshold (fraction of global maximum)
    - Isolation from other peaks (valley depth >= peak_height / valley_ratio)
    
    Args:
        grid: Parameter grid values
        pdf: KDE probability density values
        height_threshold: Minimum peak height as fraction of maximum (0-1)
        valley_ratio: Valley depth requirement (valley must be peak_height/valley_ratio deep)
        skip_valley_check: If True, skip valley isolation filtering (return all peaks above threshold)
        verbose: Print detailed diagnostic output
        
    Returns:
        peak_indices: Array indices where peaks occur
        peak_values: Parameter values at peaks
        peak_heights: PDF values at peaks
    """
    max_pdf = np.max(pdf)
    min_height = height_threshold * max_pdf
    
    # Find all peaks above height threshold
    peak_indices, properties = find_peaks(pdf, height=min_height)
    
    if len(peak_indices) == 0:
        return peak_indices, np.array([]), np.array([])
    
    peak_heights = pdf[peak_indices]
    
    # Print all detected peaks
    if verbose:
        print(f"  Found {len(peak_indices)} peaks above height threshold ({height_threshold:.1%} of max)")
        for i, (idx, height) in enumerate(zip(peak_indices, peak_heights)):
            print(f"    Peak {i+1}: value={grid[idx]:.6f}, height={height:.6f}")
    
    # Valley-depth isolation checking
    if not skip_valley_check:
        # For each pair of peaks, check if there's a sufficient valley between them
        # A peak is isolated if: valley_depth >= peak_height / valley_ratio
        # where valley_depth = peak_height - min(pdf between peaks)
        
        # Sort peaks by height (descending) - always keep tallest
        sorted_indices = np.argsort(peak_heights)[::-1]
        keep_mask = np.ones(len(peak_indices), dtype=bool)
        rejection_reasons = {}  # Store why peaks were rejected
        
        for i in range(len(sorted_indices)):
            if not keep_mask[sorted_indices[i]]:
                continue  # Already rejected
            
            smaller_idx = sorted_indices[i]
            smaller_peak_idx = peak_indices[smaller_idx]
            smaller_height = peak_heights[smaller_idx]
            
            # Check against all taller peaks
            for j in range(i):
                if not keep_mask[sorted_indices[j]]:
                    continue
                
                taller_idx = sorted_indices[j]
                taller_peak_idx = peak_indices[taller_idx]
                taller_height = peak_heights[taller_idx]
                
                # Find valley between these two peaks
                left = min(smaller_peak_idx, taller_peak_idx)
                right = max(smaller_peak_idx, taller_peak_idx)
                valley_height = np.min(pdf[left:right+1])
                
                # Valley depth measured from the shorter peak (which is the smaller peak)
                shorter_peak_height = smaller_height
                valley_depth = shorter_peak_height - valley_height
                required_depth = shorter_peak_height * (valley_ratio - 1) / valley_ratio
                
                if valley_depth < required_depth:
                    # Not isolated - reject smaller peak
                    keep_mask[smaller_idx] = False
                    rejection_reasons[smaller_idx] = {
                        'valley_depth': valley_depth,
                        'required_depth': required_depth,
                        'taller_peak': taller_idx
                    }
                    break
        
        # Print valley isolation details with proper context
        if verbose:
            print(f"\nValley Isolation Check")
            print(f"  Requirement: valley height ≤ 1/{valley_ratio:.0f} of shorter peak")
            
            tallest_idx = sorted_indices[0]
            print(f"  Peak {tallest_idx+1} (height={peak_heights[tallest_idx]:.6f}): TALLEST, automatically KEPT")
            
            for i in range(1, len(sorted_indices)):
                smaller_idx = sorted_indices[i]
                smaller_peak_idx = peak_indices[smaller_idx]
                smaller_height = peak_heights[smaller_idx]
                
                min_valley_height = np.inf
                limiting_taller_idx = None
                
                for j in range(i):
                    taller_idx = sorted_indices[j]
                    taller_peak_idx = peak_indices[taller_idx]
                    
                    left = min(smaller_peak_idx, taller_peak_idx)
                    right = max(smaller_peak_idx, taller_peak_idx)
                    valley_height = np.min(pdf[left:right+1])
                    
                    if valley_height < min_valley_height:
                        min_valley_height = valley_height
                        limiting_taller_idx = taller_idx
                
                shorter_peak_height = smaller_height
                max_valley_height = shorter_peak_height / valley_ratio
                status = "KEPT" if keep_mask[smaller_idx] else "REJECTED"
                
                taller_height = peak_heights[limiting_taller_idx] if limiting_taller_idx is not None else 0
                print(f"  Peak {smaller_idx+1} (h={smaller_height:.4f}) vs Peak {limiting_taller_idx+1} (h={taller_height:.4f}): " +
                      f"valley_h={min_valley_height:.4f}, max_allowed={max_valley_height:.4f} → {status}")
            
            n_kept = np.sum(keep_mask)
            n_rejected = np.sum(~keep_mask)
            print(f"→ Kept {n_kept}, rejected {n_rejected}")
        
        peak_indices = peak_indices[keep_mask]
        peak_values = grid[peak_indices]
        peak_heights = pdf[peak_indices]
    else:
        # Skip valley checking - keep all peaks above height threshold
        peak_values = grid[peak_indices]
        peak_heights = pdf[peak_indices]
    
    return peak_indices, peak_values, peak_heights


def check_psi_valley_isolation(candidate_modes: List[Dict], psi_raw: np.ndarray,
                                weights: np.ndarray, valley_ratio: float,
                                grid_n: int = 4000, verbose: bool = False) -> List[Dict]:
    """
    Check valley isolation for axial parameter modes in ψ-space using periodic KDE.
    
    Uses von Mises KDE with proper periodic boundary conditions to handle
    ψ ∈ [-π/2, π/2] where ψ=-π/2 and ψ=+π/2 are the same physical state.
    
    Args:
        candidate_modes: List of mode dictionaries with 'value' (in radians) and 'peak_height'
        psi_raw: Raw ψ samples in radians
        weights: Normalized sample weights
        valley_ratio: Valley depth requirement (valley_depth >= peak_height / valley_ratio)
        grid_n: Number of grid points for KDE
        verbose: Print diagnostic output
        
    Returns:
        Filtered list of modes that pass valley isolation check
    """
    if len(candidate_modes) <= 1:
        return candidate_modes
    
    # Build periodic KDE in ψ-space using von Mises approach
    # Map ψ ∈ [-π/2, π/2] to θ = 2ψ ∈ [-π, π] for circular statistics
    theta_raw = 2.0 * psi_raw
    
    # Estimate kappa for von Mises
    grid_spacing = 2.0 * np.pi / grid_n
    kappa = estimate_kappa_simple(theta_raw, weights, grid_spacing)
    
    # Build von Mises KDE on θ-grid
    theta_grid = np.linspace(-np.pi, np.pi, grid_n, endpoint=False)
    theta_pdf = von_mises_kde(theta_raw, weights, theta_grid, kappa)
    
    # Transform back to ψ-space
    # ψ = θ/2, and Jacobian dθ/dψ = 2, so pdf_ψ = 2 * pdf_θ
    psi_grid = theta_grid / 2.0
    psi_pdf = 2.0 * theta_pdf
    
    # Sort grid for easier indexing
    sort_idx = np.argsort(psi_grid)
    psi_grid = psi_grid[sort_idx]
    psi_pdf = psi_pdf[sort_idx]
    
    # Find grid indices closest to each mode
    mode_indices = []
    mode_heights_psi = []
    for mode in candidate_modes:
        idx = np.argmin(np.abs(psi_grid - mode['value']))
        mode_indices.append(idx)
        mode_heights_psi.append(psi_pdf[idx])
    
    # Sort by height
    sorted_order = np.argsort(mode_heights_psi)[::-1]
    keep_mask = np.ones(len(candidate_modes), dtype=bool)
    
    # Check valley isolation with periodic boundary handling
    for i in range(len(sorted_order)):
        if not keep_mask[sorted_order[i]]:
            continue
        
        smaller_idx = sorted_order[i]
        smaller_peak_idx = mode_indices[smaller_idx]
        smaller_height = mode_heights_psi[smaller_idx]
        
        # Check against all taller peaks
        for j in range(i):
            if not keep_mask[sorted_order[j]]:
                continue
            
            taller_idx = sorted_order[j]
            taller_peak_idx = mode_indices[taller_idx]
            taller_height = mode_heights_psi[taller_idx]
            
            # Find valley between these two peaks (check both directions due to periodicity)
            left = min(smaller_peak_idx, taller_peak_idx)
            right = max(smaller_peak_idx, taller_peak_idx)
            
            valley_height_direct = np.min(psi_pdf[left:right+1])
            wrapped_segment = np.concatenate([psi_pdf[right:], psi_pdf[:left+1]])
            valley_height_wrapped = np.min(wrapped_segment)
            valley_height = max(valley_height_direct, valley_height_wrapped)
            
            # Valley depth measured from the shorter peak (which is the smaller peak)
            shorter_peak_height = smaller_height
            valley_depth = shorter_peak_height - valley_height
            required_depth = shorter_peak_height * (valley_ratio - 1) / valley_ratio
            
            if valley_depth < required_depth:
                keep_mask[smaller_idx] = False
                break
    
    # Verbose output
    if verbose:
        print(f"\nValley Isolation Check in ψ-space (periodic KDE)")
        print(f"  Requirement: valley height ≤ 1/{valley_ratio:.0f} of shorter peak")
        
        tallest_order_idx = sorted_order[0]
        tallest_mode = candidate_modes[tallest_order_idx]
        print(f"  Mode {tallest_order_idx+1} (ψ={tallest_mode['value']:.3f} rad, {tallest_mode['value_degrees']:.1f}°, " +
              f"h={mode_heights_psi[tallest_order_idx]:.4f}): TALLEST, automatically KEPT")
        
        for i in range(1, len(sorted_order)):
            smaller_idx = sorted_order[i]
            smaller_mode = candidate_modes[smaller_idx]
            smaller_height = mode_heights_psi[smaller_idx]
            
            min_valley_height = np.inf
            limiting_taller_idx = None
            
            for j in range(i):
                taller_idx = sorted_order[j]
                
                left = min(mode_indices[smaller_idx], mode_indices[taller_idx])
                right = max(mode_indices[smaller_idx], mode_indices[taller_idx])
                
                valley_height_direct = np.min(psi_pdf[left:right+1])
                wrapped_segment = np.concatenate([psi_pdf[right:], psi_pdf[:left+1]])
                valley_height_wrapped = np.min(wrapped_segment)
                valley_height = max(valley_height_direct, valley_height_wrapped)
                
                if valley_height < min_valley_height:
                    min_valley_height = valley_height
                    limiting_taller_idx = taller_idx
            
            shorter_peak_height = smaller_height
            max_valley_height = shorter_peak_height / valley_ratio
            status = "KEPT" if keep_mask[smaller_idx] else "REJECTED"
            
            taller_mode = candidate_modes[limiting_taller_idx] if limiting_taller_idx is not None else None
            taller_height = mode_heights_psi[limiting_taller_idx] if limiting_taller_idx is not None else 0
            
            print(f"  Mode {smaller_idx+1} (ψ={smaller_mode['value']:.3f} rad, {smaller_mode['value_degrees']:.1f}°, h={smaller_height:.4f}) " +
                  f"vs Mode {limiting_taller_idx+1} (h={taller_height:.4f}): " +
                  f"valley_h={min_valley_height:.4f}, max_allowed={max_valley_height:.4f} → {status}")
        
        n_kept = np.sum(keep_mask)
        print(f"→ Kept {n_kept}, rejected {len(candidate_modes) - n_kept}")
    
    # Filter modes
    filtered_modes = [mode for i, mode in enumerate(candidate_modes) if keep_mask[i]]
    return filtered_modes


def compute_mode_statistics(mode_value: float, samples: np.ndarray, 
                            weights: np.ndarray, logl: Optional[np.ndarray],
                            config: ProcessingConfig,
                            all_mode_values: Optional[List[float]] = None) -> Dict:
    """
    Compute HDI, posterior mass, and likelihood statistics around a mode.
    
    For multi-modal distributions (when all_mode_values is provided), samples are
    assigned to their nearest mode before computing statistics. This ensures each
    mode reports statistics for its own region of the posterior.
    
    Args:
        mode_value: Mode location on parameter grid
        samples: Parameter samples (already unwrapped for axial parameters)
        weights: Sample weights (will be normalized internally)
        logl: Log-likelihood values per sample (optional)
        config: ProcessingConfig with HDI levels
        all_mode_values: List of all detected mode values (for mode assignment)
        
    Returns:
        Dictionary containing:
            - hdi_68, hdi_99: Highest density intervals
            - posterior_mass: Fraction of posterior weight in this mode's region
            - likelihood_max: Maximum log-likelihood in HDI_68
            - likelihood_std: Standard deviation of log-likelihood in HDI_68
            - n_samples: Number of samples assigned to this mode
            - n_in_hdi_68: Number of samples in HDI_68 region
    """
    if all_mode_values is not None and len(all_mode_values) > 1:
        distances = np.abs(samples[:, np.newaxis] - np.array(all_mode_values))
        nearest_mode_idx = np.argmin(distances, axis=1)
        target_mode_idx = np.argmin(np.abs(np.array(all_mode_values) - mode_value))
        mode_mask = nearest_mode_idx == target_mode_idx
        
        if np.sum(mode_mask) == 0:
            return {
                'hdi_68': (np.nan, np.nan),
                'hdi_99': (np.nan, np.nan),
                'posterior_mass': 0.0,
                'likelihood_max': np.nan,
                'likelihood_std': np.nan,
                'n_samples': 0,
                'n_in_hdi_68': 0
            }
        
        samples_mode = samples[mode_mask]
        weights_mode = weights[mode_mask]
        posterior_mass = float(np.sum(weights_mode))
        logl_mode = logl[mode_mask] if logl is not None else None
    else:
        samples_mode = samples
        weights_mode = weights
        posterior_mass = float(np.sum(weights))
        logl_mode = logl
        mode_mask = np.ones(len(samples), dtype=bool)
    
    # HDI computation (weights normalized internally by weighted_hdi)
    hdi_68 = weighted_hdi(samples_mode, weights_mode, cred_mass=0.68)
    hdi_99 = weighted_hdi(samples_mode, weights_mode, cred_mass=0.99)
    
    # Count samples in HDI_68
    mask_68 = (samples_mode >= hdi_68[0]) & (samples_mode <= hdi_68[1])
    n_in_hdi_68 = int(np.sum(mask_68))
    n_samples_assigned = int(np.sum(mode_mask))
    
    # Likelihood statistics within HDI_68
    if logl_mode is not None and np.sum(mask_68) > 0:
        logl_in_hdi = logl_mode[mask_68]
        
        if len(logl_in_hdi) > 0:
            # Get top 100 samples by likelihood
            logl_sorted = np.sort(logl_in_hdi)[::-1]  # Sort descending
            n_top = min(100, len(logl_sorted))
            logl_top = logl_sorted[:n_top]
            
            likelihood_max = float(logl_top[0])
            likelihood_std = float(np.nanstd(logl_top))
        else:
            likelihood_max = np.nan
            likelihood_std = np.nan
    else:
        likelihood_max = np.nan
        likelihood_std = np.nan
    
    return {
        'hdi_68': hdi_68,
        'hdi_99': hdi_99,
        'posterior_mass': posterior_mass,
        'likelihood_max': likelihood_max,
        'likelihood_std': likelihood_std,
        'n_samples': n_samples_assigned,
        'n_in_hdi_68': n_in_hdi_68
    }


def merge_nearby_modes(modes: List[Dict], merge_distance: float,
                       strategy: str = 'keep_most_mass',
                       param_type: str = 'linear',
                       verbose: bool = False) -> List[Dict]:
    """
    Merge modes that are closer than merge_distance threshold.
    
    Args:
        modes: List of mode dictionaries with 'value', 'posterior_mass', 'peak_height'
        merge_distance: Distance threshold for merging
        strategy: 'keep_most_mass' or 'keep_highest'
        param_type: 'linear' or 'axial' (for ψ₀ with 180° periodicity)
        verbose: Print detailed diagnostic output
        
    Returns:
        Filtered list of modes with nearby peaks merged
    """
    if len(modes) <= 1:
        return modes
    
    sorted_modes = sorted(modes, key=lambda m: m['value'])
    
    merged = []
    i = 0
    while i < len(sorted_modes):
        current = sorted_modes[i]
        
        group = [current]
        j = i + 1
        while j < len(sorted_modes):
            if param_type == 'axial':
                # Axial distance: d = min(|Δ|, 180 - |Δ|) in degrees
                diff = abs(sorted_modes[j]['value'] - current['value'])
                dist = min(diff, 180.0 - diff)
            else:
                # Linear distance
                dist = abs(sorted_modes[j]['value'] - current['value'])
            
            if dist <= merge_distance:
                group.append(sorted_modes[j])
                j += 1
            else:
                break
        
        if strategy == 'keep_most_mass':
            best = max(group, key=lambda m: m['posterior_mass'])
        elif strategy == 'keep_highest':
            best = max(group, key=lambda m: m['peak_height'])
        else:
            raise ValueError(f"Unknown merge strategy: {strategy}")
        
        merged.append(best)
        i = j
    
    return merged


def rank_and_filter_modes(modes: List[Dict], mass_threshold: float,
                          ranking: str = 'posterior_mass',
                          max_modes: Optional[int] = None,
                          verbose: bool = False) -> List[Dict]:
    """
    Filter modes by mass threshold and rank by importance.
    
    Args:
        modes: List of mode dictionaries
        mass_threshold: Minimum posterior mass to qualify as significant mode
        ranking: 'posterior_mass' or 'likelihood' for ranking criterion
        max_modes: Maximum number of modes to keep (None = keep all)
        verbose: Print detailed diagnostic output
        
    Returns:
        Filtered and sorted list of significant modes
    """
    passed = []
    rejected = []
    
    for m in modes:
        if m['posterior_mass'] >= mass_threshold:
            passed.append(m)
        else:
            rejected.append(m)
    
    if len(passed) == 0:
        return []
    
    # CRITICAL: When selecting exactly 1 mode, pick the HIGHEST PEAK (most visually prominent)
    # For multiple modes, use configured ranking (mass or likelihood gives importance)
    if max_modes == 1 or len(passed) == 1:
        ranked = sorted(passed, key=lambda m: m['peak_height'], reverse=True)
    else:
        # Multiple modes: use configured ranking
        if ranking == 'posterior_mass':
            ranked = sorted(passed, key=lambda m: m['posterior_mass'], reverse=True)
        elif ranking == 'likelihood':
            ranked = sorted(passed, key=lambda m: m['likelihood_max'], reverse=True)
        else:
            raise ValueError(f"Unknown ranking criterion: {ranking}")
    
    if max_modes is not None and len(ranked) > max_modes:
        ranked = ranked[:max_modes]
    
    return ranked


def process_linear_parameter(x: np.ndarray, weights: np.ndarray,
                             param_name: str,
                             config: ProcessingConfig,
                             logl: Optional[np.ndarray] = None) -> List[Dict]:
    """
    Process linear parameter with multi-modal detection.
    
    Workflow:
    1. PHASE 1: Detect candidate modes via KDE
    2. PHASE 2: Merge nearby modes and filter by significance
    3. PHASE 3: Assign ALL samples to FINAL modes only
    4. PHASE 5: Compute statistics per mode
    5. PHASE 6: Store results consistently
    
    Args:
        x: Sample values
        weights: Normalized weights
        param_name: Parameter name
        config: ProcessingConfig with multi-modal settings
        logl: Optional log-likelihood values per sample
        
    Returns:
        List of mode dictionaries with complete assignments and statistics
    """
    # =======================================================================
    # PHASE 1: MODE DETECTION
    # =======================================================================
    mode, grid, pdf = find_mode_linear(x, weights, 
                                       bandwidth=config.kde_bandwidth,
                                       grid_n=config.grid_n)
    
    # Find peaks
    if not config.detect_multiple_modes:
        # Force single mode
        peak_idx = np.argmax(pdf)
        peak_values = [grid[peak_idx]]
        peak_heights = [pdf[peak_idx]]
    else:
        peak_indices, peak_values, peak_heights = find_kde_peaks(
            grid, pdf, config.mode_height_threshold, 
            valley_ratio=config.mode_valley_ratio,
            verbose=config.verbose
        )
        
        if len(peak_indices) == 0:
            # No peaks found, use global maximum
            peak_idx = np.argmax(pdf)
            peak_values = [grid[peak_idx]]
            peak_heights = [pdf[peak_idx]]
    
    # Build candidate modes
    candidate_modes = []
    for peak_val, peak_height in zip(peak_values, peak_heights):
        candidate_modes.append({
            'value': float(peak_val),
            'peak_height': float(peak_height)
        })
    
    if config.verbose:
        # Format parameter name with proper subscripts
        display_name = param_name.replace('p0', 'p₀').replace('phi_rm', 'ϕ_RM').replace('sigma_phi', 'σ_ϕ').replace('beta', 'β')
        print(f"\n{'~'*70}")
        print(f"MULTI-MODAL PROCESSING: {display_name}")
        print(f"{'~'*70}")
        print(f"\nPHASE 1 - Peak Detection")
        print(f"  Detected {len(candidate_modes)} candidate peak(s):")
        for i, m in enumerate(candidate_modes):
            print(f"    Peak {i+1}: value={m['value']:.6f}, height={m['peak_height']:.6f}")
    
    # =======================================================================
    # PHASE 2: MERGE AND FILTER
    # =======================================================================
    if len(candidate_modes) > 1:
        # Compute merge distance as fraction of 99% sample range
        # Use weighted quantiles to avoid spurious outliers
        q_low = dyfunc.quantile(x, [0.005], weights=weights)[0]
        q_high = dyfunc.quantile(x, [0.995], weights=weights)[0]
        range_99 = q_high - q_low
        merge_dist = config.mode_merge_fraction * range_99
        
        if config.verbose:
            print(f"\nPHASE 2 - Mode Merging")
            print(f"  Sample range (99%): [{q_low:.6f}, {q_high:.6f}] (span = {range_99:.6f})")
            print(f"  Merge threshold: {config.mode_merge_fraction:.1%} × span = {merge_dist:.6f}")
        
        # Merge nearby modes (using linear distance)
        modes_for_merge = []
        for m in candidate_modes:
            modes_for_merge.append({
                'value': m['value'],
                'posterior_mass': 1.0 / len(candidate_modes),  # Placeholder
                'peak_height': m['peak_height'],
                'original': m
            })
        
        merged = merge_nearby_modes(modes_for_merge, merge_dist, 
                                   config.mode_merge_strategy, param_type='linear',
                                   verbose=config.verbose)
        candidate_modes = [m['original'] for m in merged]
        
        if config.verbose:
            print(f"  Result: {len(candidate_modes)} mode(s) after merging")
            for i, m in enumerate(candidate_modes):
                print(f"    Mode {i+1}: value={m['value']:.6f}")
    
    # If filtering would leave us with nothing, keep at least one mode
    if len(candidate_modes) == 0:
        peak_idx = np.argmax(pdf)
        candidate_modes = [{
            'value': float(grid[peak_idx]),
            'peak_height': float(pdf[peak_idx])
        }]
    
    # Extract final mode values
    final_mode_values = [m['value'] for m in candidate_modes]
    
    # =======================================================================
    # PHASE 3: ASSIGN ALL SAMPLES TO FINAL MODES
    # =======================================================================
    mode_assignment = assign_samples_to_modes_linear(x, final_mode_values)
    
    # =======================================================================
    # PHASE 5: COMPUTE STATISTICS PER MODE
    # =======================================================================
    final_modes = []
    
    for mode_idx, mode_data in enumerate(candidate_modes):
        mask = (mode_assignment == mode_idx)
        n_assigned = np.sum(mask)
        
        if n_assigned == 0:
            # This shouldn't happen with proper assignment, but handle gracefully
            continue
        
        # Get assigned samples
        x_mode = x[mask]
        weights_mode = weights[mask]
        logl_mode = logl[mask] if logl is not None else None
        
        # Posterior mass = sum of assigned weights
        posterior_mass = float(np.sum(weights_mode))
        
        # Compute statistics
        stats = compute_mode_statistics(
            mode_data['value'],
            x_mode,
            weights_mode,
            logl_mode,
            config,
            None
        )
        
        # =======================================================================
        # PHASE 6: STORE RESULTS
        # =======================================================================
        mode_dict = {
            'value': mode_data['value'],
            'peak_height': mode_data['peak_height'],
            'posterior_mass': posterior_mass,
            'hdi_68': stats['hdi_68'],
            'hdi_99': stats['hdi_99'],
            'likelihood_max': stats['likelihood_max'],
            'likelihood_std': stats['likelihood_std'],
            'n_samples': stats['n_samples'],
            'n_in_hdi_68': stats.get('n_in_hdi_68', stats['n_samples']),
            # Store full arrays for all modes
            'unwrapped_samples': x,
            'mode_assignment': mode_assignment,
            # Plotting data
            'grid': grid,
            'pdf': pdf,
            'type': 'linear',
            'param_name': param_name
        }
        
        final_modes.append(mode_dict)
    
    # =======================================================================
    # FINAL FILTER BY MASS THRESHOLD
    # =======================================================================
    if len(final_modes) > 1:
        modes_for_filter = []
        for m in final_modes:
            modes_for_filter.append({
                'value': m['value'],
                'posterior_mass': m['posterior_mass'],
                'likelihood_max': m['likelihood_max'],
                'peak_height': m.get('peak_height', 0.0),  # Include for ranking
                'original': m
            })
        
        filtered = rank_and_filter_modes(
            modes_for_filter,
            config.mode_mass_threshold,
            config.mode_ranking,
            config.max_modes,
            verbose=config.verbose
        )
        
        if len(filtered) > 0:
            if config.verbose:
                print(f"\nPHASE 3 - Final Filtering & Reassignment")
                print(f"  Mass threshold: {config.mode_mass_threshold:.1%}")
                print(f"  Kept {len(filtered)}/{len(final_modes)} mode(s):")
                for i, m in enumerate(filtered):
                    print(f"    Mode {i+1}: value={m['value']:.6f}, mass={m['posterior_mass']:.3f}")
            
            # CRITICAL: Reassign samples if modes were filtered out
            filtered_modes = [m['original'] for m in filtered]
            filtered_mode_values = [m['value'] for m in filtered_modes]
            
            # Reassign all samples to remaining modes
            mode_assignment = assign_samples_to_modes_linear(x, filtered_mode_values)
            
            # Recompute statistics with new assignments
            for mode_idx, mode_dict in enumerate(filtered_modes):
                mask = (mode_assignment == mode_idx)
                x_mode = x[mask]
                weights_mode = weights[mask]
                logl_mode = logl[mask] if logl is not None else None
                
                posterior_mass = float(np.sum(weights_mode))
                
                stats = compute_mode_statistics(
                    mode_dict['value'],
                    x_mode,
                    weights_mode,
                    logl_mode,
                    config,
                    None
                )
                
                # Update mode dictionary
                mode_dict['posterior_mass'] = posterior_mass
                mode_dict['hdi_68'] = stats['hdi_68']
                mode_dict['hdi_99'] = stats['hdi_99']
                mode_dict['likelihood_max'] = stats['likelihood_max']
                mode_dict['likelihood_std'] = stats['likelihood_std']
                mode_dict['n_samples'] = stats['n_samples']
                mode_dict['n_in_hdi_68'] = stats.get('n_in_hdi_68', stats['n_samples'])
                mode_dict['mode_assignment'] = mode_assignment
            
            return filtered_modes
        else:
            # All modes filtered out - return single mode at global peak
            peak_idx = np.argmax(pdf)
            peak_value = float(grid[peak_idx])
            peak_height = float(pdf[peak_idx])
            mode_assignment = np.zeros(len(x), dtype=int)
            
            stats = compute_mode_statistics(peak_value, x, weights, logl, config, None)
            
            mode_dict = {
                'value': peak_value,
                'peak_height': peak_height,
                'posterior_mass': 1.0,
                'hdi_68': stats['hdi_68'],
                'hdi_99': stats['hdi_99'],
                'likelihood_max': stats['likelihood_max'],
                'likelihood_std': stats['likelihood_std'],
                'n_samples': stats['n_samples'],
                'n_in_hdi_68': stats.get('n_in_hdi_68', stats['n_samples']),
                'unwrapped_samples': x,
                'mode_assignment': mode_assignment,
                'grid': grid,
                'pdf': pdf,
                'type': 'linear',
                'param_name': param_name
            }
            
            return [mode_dict]
    
    return final_modes


# =============================================================================
# PARAMETER FORMATTING
# =============================================================================

def format_param_label(param_name: str) -> str:
    """
    Format parameter name for plots with LaTeX notation.
    
    Args:
        param_name: Raw parameter name (e.g., 'S1_phi_rm', 'comp0_p0')
        
    Returns:
        Formatted LaTeX label
    """
    # Extract component name if present
    parts = param_name.split('_')
    
    if '_p0' in param_name or param_name.endswith('_p0'):
        comp = param_name.replace('_p0', '')
        if comp:
            return f'$p_0$ ({comp})'
        else:
            return r'$p_0$'
    
    elif '_psi_0' in param_name or 'psi_0' in param_name:
        comp = param_name.replace('_psi_0', '').replace('psi_0', '')
        if comp:
            return f'$\\psi_0$ ({comp})'
        else:
            return r'$\psi_0$'
    
    elif 'phi_rm' in param_name:
        # Extract component prefix (e.g., 'S1' from 'S1_phi_rm')
        comp = param_name.replace('_phi_rm', '').replace('phi_rm', '')
        if comp:
            return f'$\\phi_{{\\rm{{RM}}}}$ ({comp})'
        else:
            return r'$\phi_{\rm{RM}}$'
    
    elif 'sigma_phi' in param_name:
        comp = param_name.replace('_sigma_phi', '').replace('sigma_phi', '')
        if comp:
            return f'$\\sigma_\\phi$ ({comp})'
        else:
            return r'$\sigma_\phi$'
    
    elif param_name.endswith('_N'):
        comp = param_name.replace('_N', '')
        if comp:
            return f'$N$ ({comp})'
        else:
            return r'$N$'
    
    elif 'beta' in param_name:
        comp = param_name.replace('_beta', '').replace('beta', '')
        if comp:
            return f'$\\beta$ ({comp})'
        else:
            return r'$\beta$'
    
    # Fallback: return as-is
    return param_name


# =============================================================================
# CIRCULAR/AXIAL PARAMETER MODE DETECTION (ψ₀ via φ)
# =============================================================================

def von_mises_kde(phi_samples: np.ndarray, weights: np.ndarray,
                  phi_grid: np.ndarray, kappa: float) -> np.ndarray:
    """
    Weighted von Mises kernel density estimate for circular data.
    
    Uses numerically stable formulation with i0e (exponentially scaled Bessel)
    to prevent overflow when kappa is large.
    
    f(φ) = Σᵢ wᵢ · K(φ - φᵢ; κ)
    where K(θ; κ) = exp(κ·cos(θ)) / (2π·I₀(κ))
                  = exp(κ·(cos(θ)-1)) / (2π·i0e(κ))  [stable form]
    
    Args:
        phi_samples: Circular sample values in [-π, π]
        weights: Normalized weights
        phi_grid: Grid of evaluation points
        kappa: Concentration parameter (larger = narrower kernel)
        
    Returns:
        pdf: KDE values on grid
    """
    # Normalization using exponentially scaled Bessel function
    # i0e(kappa) = I0(kappa) * exp(-kappa)
    norm = 1.0 / (2.0 * np.pi * i0e(kappa))
    
    # Vectorized computation: (n_grid, n_samples)
    # Broadcast: phi_grid[:, None] - phi_samples[None, :]
    diff = phi_grid[:, np.newaxis] - phi_samples[np.newaxis, :]
    
    # Numerically stable: exp(kappa*(cos(diff)-1))
    # Since cos(diff) ≤ 1, exponent is always ≤ 0 (no overflow)
    kernel = norm * np.exp(kappa * (np.cos(diff) - 1.0))
    
    # Weight and sum across samples
    pdf = np.sum(kernel * weights[np.newaxis, :], axis=1)
    
    return pdf


def estimate_kappa_simple(phi_samples: np.ndarray, weights: np.ndarray,
                         grid_spacing: float) -> float:
    """
    Simple estimate of concentration parameter κ for von Mises KDE.
    
    Uses method of moments: relates κ to weighted mean resultant length.
    Includes grid resolution constraint to prevent under-resolved kernels.
    
    Args:
        phi_samples: Circular samples in [-π, π]
        weights: Normalized weights
        grid_spacing: Spacing between grid points (2π/N_grid)
        
    Returns:
        kappa: Concentration parameter (clamped to grid resolution)
    """
    # Compute weighted mean direction
    cos_mean = np.average(np.cos(phi_samples), weights=weights)
    sin_mean = np.average(np.sin(phi_samples), weights=weights)
    
    # Mean resultant length
    R = np.sqrt(cos_mean**2 + sin_mean**2)
    
    # Clip R away from 1 to prevent infinite kappa
    R = np.clip(R, 0.0, 0.999)
    
    # Approximate κ from R (for R not too close to 1)
    if R < 0.53:
        kappa = 2.0 * R + R**3 + 5.0 * R**5 / 6.0
    elif R < 0.85:
        kappa = -0.4 + 1.39 * R + 0.43 / (1.0 - R)
    else:
        kappa = 1.0 / (R**3 - 4.0 * R**2 + 3.0 * R)
    
    # Apply bandwidth scaling similar to Scott's rule
    # For circular KDE, scale κ by n^(-2/5)
    # Protect against overflow when weights are very small
    with np.errstate(over='ignore', invalid='ignore'):
        weight_sum_sq = np.sum(weights**2)
        weight_sum_sq = np.clip(weight_sum_sq, 1e-10, 1.0)  # Prevent division by zero/overflow
        n_eff = 1.0 / weight_sum_sq  # Effective sample size
        n_eff = np.clip(n_eff, 1.0, 1e6)  # Reasonable bounds
        scale_factor = n_eff**(-2.0/5.0)
    kappa = kappa / scale_factor
    
    # Grid resolution constraint: kernel width ~ 1/sqrt(kappa)
    # Require at least 5 grid points across kernel width
    # dphi < width/5 → dphi < 1/(5*sqrt(kappa)) → kappa < (5/dphi)^2
    kappa_max = (5.0 / grid_spacing) ** 2
    
    if kappa > kappa_max:
        print(f"Warning: Clamping kappa from {kappa:.1f} to {kappa_max:.1f} (grid resolution limit)")
        print(f"  Grid spacing: {grid_spacing:.4f} rad, kernel width would be {1.0/np.sqrt(kappa):.4f} rad")
        kappa = kappa_max
    
    return float(kappa)


def wrap_angle_to_range(angle: float, low: float, high: float) -> float:
    """
    Wrap angle to specified range [low, high).
    
    Args:
        angle: Angle value
        low: Lower bound
        high: Upper bound
        
    Returns:
        Wrapped angle in [low, high)
    """
    span = high - low
    wrapped = low + np.fmod(angle - low, span)
    if wrapped < low:
        wrapped += span
    return float(wrapped)


def wrap_samples_to_mode(psi_raw: np.ndarray, mode_canonical: float) -> np.ndarray:
    """
    Wrap ψ₀ samples to be continuous around mode via π shifts.
    
    For each sample ψ₀_raw, find integer m that minimizes:
        |ψ₀_raw + m·π - mode_canonical|
    
    This creates a continuous distribution centered on mode_canonical.
    Samples may extend beyond [-π/2, π/2] if mode is near boundary.
    
    Args:
        psi_raw: Raw ψ₀ samples (may be scattered due to π-periodicity)
        mode_canonical: Mode position in [-π/2, π/2]
        
    Returns:
        psi_wrapped: Samples shifted to be near mode_canonical
    """
    psi_wrapped = np.zeros_like(psi_raw)
    
    for i, psi in enumerate(psi_raw):
        # Find best integer shift
        m = np.round((mode_canonical - psi) / np.pi)
        psi_wrapped[i] = psi + m * np.pi
    
    return psi_wrapped


def assign_samples_to_modes_axial(psi_raw: np.ndarray, 
                                  mode_values: List[float]) -> np.ndarray:
    """
    Assign each sample to nearest mode using π-periodic distance.
    
    For axial parameters with π periodicity, distance between sample and mode
    is computed using optimal integer shift.
    
    Args:
        psi_raw: Raw ψ₀ samples in radians
        mode_values: List of mode positions in radians
        
    Returns:
        mode_assignment: Array of mode indices (0 to N-1) for each sample
    """
    n_samples = len(psi_raw)
    n_modes = len(mode_values)
    mode_assignment = np.zeros(n_samples, dtype=int)
    
    for i, psi in enumerate(psi_raw):
        distances = np.zeros(n_modes)
        for j, mode_psi in enumerate(mode_values):
            # Find integer number of π shifts that puts sample closest to mode
            m = np.round((mode_psi - psi) / np.pi)
            shifted_psi = psi + m * np.pi
            distances[j] = abs(shifted_psi - mode_psi)
        
        mode_assignment[i] = np.argmin(distances)
    
    return mode_assignment


def assign_samples_to_modes_linear(x: np.ndarray, 
                                   mode_values: List[float]) -> np.ndarray:
    """
    Assign each sample to nearest mode using linear distance.
    
    Args:
        x: Sample values
        mode_values: List of mode positions
        
    Returns:
        mode_assignment: Array of mode indices (0 to N-1) for each sample
    """
    n_samples = len(x)
    n_modes = len(mode_values)
    mode_assignment = np.zeros(n_samples, dtype=int)
    
    for i, sample in enumerate(x):
        distances = np.abs(np.array(mode_values) - sample)
        mode_assignment[i] = np.argmin(distances)
    
    return mode_assignment


def unwrap_samples_by_assignment(psi_raw: np.ndarray,
                                mode_assignment: np.ndarray,
                                mode_values: List[float]) -> np.ndarray:
    """
    Unwrap axial samples according to their mode assignments.
    
    Each sample is shifted by integer multiples of π to be near its assigned mode.
    
    Args:
        psi_raw: Raw ψ₀ samples in radians
        mode_assignment: Mode index for each sample
        mode_values: List of mode positions in radians
        
    Returns:
        psi_unwrapped: Unwrapped samples in radians
    """
    n_samples = len(psi_raw)
    psi_unwrapped = np.zeros_like(psi_raw)
    
    for i in range(n_samples):
        assigned_mode_idx = mode_assignment[i]
        mode_psi = mode_values[assigned_mode_idx]
        
        # Find integer number of π shifts that puts sample closest to mode
        m = np.round((mode_psi - psi_raw[i]) / np.pi)
        psi_unwrapped[i] = psi_raw[i] + m * np.pi
    
    return psi_unwrapped


def process_axial_parameter(a: np.ndarray, b: np.ndarray, weights: np.ndarray,
                           config: ProcessingConfig,
                           logl: Optional[np.ndarray] = None,
                           param_name: str = 'psi_0') -> List[Dict]:
    """
    Process axial angle ψ₀ with multi-modal detection via circular φ space.
    
    Workflow:
    1. PHASE 1: Detect candidate modes in φ-space
    2. PHASE 2: Merge nearby modes and filter by significance
    3. PHASE 3: Assign ALL samples to FINAL modes only
    4. PHASE 4: Unwrap samples according to assignments
    5. PHASE 5: Compute statistics per mode
    6. PHASE 6: Store results consistently
    
    Args:
        a: Cartesian Q component samples
        b: Cartesian U component samples  
        weights: Normalized weights
        config: ProcessingConfig with multi-modal settings
        logl: Optional log-likelihood values per sample
        
    Returns:
        List of mode dictionaries with complete assignments and statistics
    """
    # =======================================================================
    # PHASE 1: MODE DETECTION
    # =======================================================================
    phi = np.arctan2(b, a)
    
    # Compute circular mean and recenter
    sin_mean = np.average(np.sin(phi), weights=weights)
    cos_mean = np.average(np.cos(phi), weights=weights)
    phi_center = np.arctan2(sin_mean, cos_mean)
    phi_recentered = np.angle(np.exp(1j * (phi - phi_center)))
    
    # Compute KDE in φ-space
    phi_grid_recentered = np.linspace(-np.pi, np.pi, config.circular_grid_n, endpoint=False)
    grid_spacing = 2.0 * np.pi / config.circular_grid_n
    kappa = estimate_kappa_simple(phi_recentered, weights, grid_spacing)
    phi_pdf_recentered = von_mises_kde(phi_recentered, weights, 
                                       phi_grid_recentered, kappa)
    
    # Store grid for plotting
    phi_grid_original = np.angle(np.exp(1j * (phi_grid_recentered + phi_center)))
    sort_idx = np.argsort(phi_grid_original)
    phi_grid_plot = phi_grid_original[sort_idx]
    phi_pdf_plot = phi_pdf_recentered[sort_idx]
    
    psi_raw = phi / 2.0
    
    # Find peaks in φ-space (WITHOUT valley checking - will check in ψ-space)
    if not config.detect_multiple_modes:
        # Force single mode
        peak_indices = [np.argmax(phi_pdf_recentered)]
        peak_phi_recentered = [phi_grid_recentered[peak_indices[0]]]
        peak_heights_phi = [phi_pdf_recentered[peak_indices[0]]]
    else:
        peak_indices, peak_phi_recentered, peak_heights_phi = find_kde_peaks(
            phi_grid_recentered, phi_pdf_recentered, config.mode_height_threshold,
            valley_ratio=config.mode_valley_ratio,
            skip_valley_check=True,  # Skip valley check in φ-space
            verbose=False  # Suppress φ-space peak listing
        )
        
        if len(peak_indices) == 0:
            # No peaks found, use global maximum
            peak_indices = [np.argmax(phi_pdf_recentered)]
            peak_phi_recentered = [phi_grid_recentered[peak_indices[0]]]
            peak_heights_phi = [phi_pdf_recentered[peak_indices[0]]]
    
    # Convert peaks to ψ-space and canonicalize
    candidate_modes = []
    for peak_phi_rec, peak_height_phi in zip(peak_phi_recentered, peak_heights_phi):
        mode_phi = np.angle(np.exp(1j * (peak_phi_rec + phi_center)))
        mode_psi = mode_phi / 2.0
        mode_psi_canonical = wrap_angle_to_range(mode_psi, -np.pi/2, np.pi/2)
        
        candidate_modes.append({
            'value': mode_psi_canonical,
            'value_degrees': np.degrees(mode_psi_canonical),
            'peak_height': 2.0 * peak_height_phi,  # Jacobian for ψ = φ/2
            'mode_phi': mode_phi
        })
    
    if config.verbose:
        # Format parameter name: S1_psi_0 -> S1_ψ₀
        display_name = param_name.replace('psi_0', 'ψ₀')
        print(f"\n{'~'*70}")
        print(f"MULTI-MODAL PROCESSING: {display_name} (AXIAL PARAMETER)")
        print(f"{'~'*70}")
        print(f"\nInitial Peak Detection in φ-space")
        print(f"  Detected {len(candidate_modes)} candidate peak(s) before valley check:")
        for i, m in enumerate(candidate_modes):
            print(f"    Peak {i+1}: ψ={m['value']:.4f} rad ({m['value_degrees']:.2f}°), φ={m['mode_phi']:.4f} rad")
    
    # Check valley isolation in ψ-space (physical space)
    if config.detect_multiple_modes and len(candidate_modes) > 1:
        candidate_modes = check_psi_valley_isolation(
            candidate_modes, psi_raw, weights, config.mode_valley_ratio,
            grid_n=config.circular_grid_n, verbose=config.verbose
        )
    
    # =======================================================================
    # PHASE 2: MERGE AND FILTER
    # =======================================================================
    if len(candidate_modes) > 1:
        # Compute merge distance as fraction of 99% sample range
        # For axial parameters, use unwrapped samples to avoid periodicity issues
        q_low = dyfunc.quantile(psi_raw, [0.005], weights=weights)[0]
        q_high = dyfunc.quantile(psi_raw, [0.995], weights=weights)[0]
        range_99_rad = q_high - q_low
        
        # Handle wrap-around: if range spans more than π, data wraps around boundaries
        # In that case, use full range π as the scale
        if range_99_rad > np.pi:
            range_99_rad = np.pi
        
        range_99_deg = np.degrees(range_99_rad)
        merge_dist = config.mode_merge_fraction * range_99_deg
        
        if config.verbose:
            print(f"\nMode Merging")
            print(f"  Sample range (99%): [{np.degrees(q_low):.2f}°, {np.degrees(q_high):.2f}°] (span = {range_99_deg:.2f}°)")
            print(f"  Merge threshold: {config.mode_merge_fraction:.1%} × span = {merge_dist:.2f}°")
        
        # Merge nearby modes (using axial distance)
        modes_for_merge = []
        for m in candidate_modes:
            modes_for_merge.append({
                'value': m['value_degrees'],
                'posterior_mass': 1.0 / len(candidate_modes),  # Placeholder
                'peak_height': m['peak_height'],
                'original': m
            })
        
        merged = merge_nearby_modes(modes_for_merge, merge_dist, 
                                   config.mode_merge_strategy, param_type='axial',
                                   verbose=config.verbose)
        candidate_modes = [m['original'] for m in merged]
        
        if config.verbose:
            print(f"  Result: {len(candidate_modes)} mode(s) after merging")
            for i, m in enumerate(candidate_modes):
                print(f"    Mode {i+1}: ψ={m['value']:.4f} rad ({m['value_degrees']:.2f}°)")
    
    # If filtering would leave us with nothing, keep at least one mode
    if len(candidate_modes) == 0:
        mode_phi_recentered_idx = np.argmax(phi_pdf_recentered)
        mode_phi = np.angle(np.exp(1j * (phi_grid_recentered[mode_phi_recentered_idx] + phi_center)))
        mode_psi = mode_phi / 2.0
        mode_psi_canonical = wrap_angle_to_range(mode_psi, -np.pi/2, np.pi/2)
        
        candidate_modes = [{
            'value': mode_psi_canonical,
            'value_degrees': np.degrees(mode_psi_canonical),
            'peak_height': 2.0 * phi_pdf_recentered[mode_phi_recentered_idx],
            'mode_phi': mode_phi
        }]
    
    # Extract final mode values
    final_mode_values = [m['value'] for m in candidate_modes]
    
    # =======================================================================
    # PHASE 3: ASSIGN ALL SAMPLES TO FINAL MODES
    # =======================================================================
    mode_assignment = assign_samples_to_modes_axial(psi_raw, final_mode_values)
    
    if config.verbose:
        print(f"[psi_0] PHASE 3: Assigned {len(psi_raw)} samples to {len(final_mode_values)} mode(s)")
        for mode_idx in range(len(final_mode_values)):
            n_assigned = np.sum(mode_assignment == mode_idx)
            print(f"  Mode {mode_idx+1}: {n_assigned} samples assigned")
    
    # =======================================================================
    # PHASE 4: UNWRAP SAMPLES ACCORDING TO ASSIGNMENTS
    # =======================================================================
    psi_unwrapped = unwrap_samples_by_assignment(psi_raw, mode_assignment, final_mode_values)
    psi_unwrapped_degrees = np.degrees(psi_unwrapped)
    
    # =======================================================================
    # PHASE 5: COMPUTE STATISTICS PER MODE
    # =======================================================================
    final_modes = []
    
    for mode_idx, mode_data in enumerate(candidate_modes):
        mask = (mode_assignment == mode_idx)
        n_assigned = np.sum(mask)
        
        if n_assigned == 0:
            # This shouldn't happen with proper assignment, but handle gracefully
            continue
        
        # Get assigned samples
        psi_mode_degrees = psi_unwrapped_degrees[mask]
        weights_mode = weights[mask]
        logl_mode = logl[mask] if logl is not None else None
        
        # Posterior mass = sum of assigned weights
        posterior_mass = float(np.sum(weights_mode))
        
        # Compute statistics
        stats = compute_mode_statistics(
            mode_data['value_degrees'],
            psi_mode_degrees,
            weights_mode,
            logl_mode,
            config,
            None
        )
        
        # =======================================================================
        # PHASE 6: STORE RESULTS
        # =======================================================================
        mode_dict = {
            'value': mode_data['value'],
            'value_degrees': mode_data['value_degrees'],
            'peak_height': mode_data['peak_height'],
            'posterior_mass': posterior_mass,
            'hdi_68': (np.radians(stats['hdi_68'][0]), np.radians(stats['hdi_68'][1])),
            'hdi_99': (np.radians(stats['hdi_99'][0]), np.radians(stats['hdi_99'][1])),
            'hdi_68_degrees': stats['hdi_68'],
            'hdi_99_degrees': stats['hdi_99'],
            'likelihood_max': stats['likelihood_max'],
            'likelihood_std': stats['likelihood_std'],
            'n_samples': stats['n_samples'],
            'n_in_hdi_68': stats.get('n_in_hdi_68', stats['n_samples']),
            # Store full arrays for all modes (needed for plotting)
            'psi_raw': psi_raw,
            'psi_unwrapped': psi_unwrapped,
            'mode_assignment': mode_assignment,
            'unwrapped_samples': psi_unwrapped,  # In RADIANS
            # Plotting data
            'phi_grid': phi_grid_plot,
            'phi_pdf': phi_pdf_plot,
            'mode_phi': mode_data['mode_phi'],
            'kappa': kappa,
            'type': 'axial',
            'param_name': 'psi_0'
        }
        
        final_modes.append(mode_dict)
    
    # =======================================================================
    # FINAL FILTER BY MASS THRESHOLD
    # =======================================================================
    if len(final_modes) > 1:
        modes_for_filter = []
        for m in final_modes:
            modes_for_filter.append({
                'value': m['value_degrees'],
                'posterior_mass': m['posterior_mass'],
                'likelihood_max': m['likelihood_max'],
                'peak_height': m.get('peak_height', 0.0),  # Include for ranking
                'original': m
            })
        
        filtered = rank_and_filter_modes(
            modes_for_filter,
            config.mode_mass_threshold,
            config.mode_ranking,
            config.max_modes,
            verbose=config.verbose
        )
        
        if len(filtered) > 0:
            if config.verbose:
                print(f"\nMass Filtering & Reassignment")
                print(f"  Mass threshold: {config.mode_mass_threshold:.1%}")
                print(f"  Kept {len(filtered)}/{len(final_modes)} mode(s):")
                for i, m in enumerate(filtered):
                    print(f"    Mode {i+1}: ψ={m['original']['value']:.4f} rad ({m['original']['value_degrees']:.2f}°), mass={m['posterior_mass']:.3f}")
            
            # CRITICAL: Reassign samples if modes were filtered out
            filtered_modes = [m['original'] for m in filtered]
            filtered_mode_values = [m['value'] for m in filtered_modes]
            
            # Reassign all samples to remaining modes
            mode_assignment = assign_samples_to_modes_axial(psi_raw, filtered_mode_values)
            psi_unwrapped = unwrap_samples_by_assignment(psi_raw, mode_assignment, filtered_mode_values)
            psi_unwrapped_degrees = np.degrees(psi_unwrapped)
            
            # Recompute statistics with new assignments
            for mode_idx, mode_dict in enumerate(filtered_modes):
                mask = (mode_assignment == mode_idx)
                psi_mode_degrees = psi_unwrapped_degrees[mask]
                weights_mode = weights[mask]
                logl_mode = logl[mask] if logl is not None else None
                
                posterior_mass = float(np.sum(weights_mode))
                
                stats = compute_mode_statistics(
                    mode_dict['value_degrees'],
                    psi_mode_degrees,
                    weights_mode,
                    logl_mode,
                    config,
                    None
                )
                
                # Update mode dictionary
                mode_dict['posterior_mass'] = posterior_mass
                mode_dict['hdi_68'] = (np.radians(stats['hdi_68'][0]), np.radians(stats['hdi_68'][1]))
                mode_dict['hdi_99'] = (np.radians(stats['hdi_99'][0]), np.radians(stats['hdi_99'][1]))
                mode_dict['hdi_68_degrees'] = stats['hdi_68']
                mode_dict['hdi_99_degrees'] = stats['hdi_99']
                mode_dict['likelihood_max'] = stats['likelihood_max']
                mode_dict['likelihood_std'] = stats['likelihood_std']
                mode_dict['n_samples'] = stats['n_samples']
                mode_dict['n_in_hdi_68'] = stats.get('n_in_hdi_68', stats['n_samples'])
                # Update arrays
                mode_dict['psi_unwrapped'] = psi_unwrapped
                mode_dict['mode_assignment'] = mode_assignment
                mode_dict['unwrapped_samples'] = psi_unwrapped
            
            return filtered_modes
        else:
            # All modes filtered out - return single mode at global peak
            mode_phi_recentered_idx = np.argmax(phi_pdf_recentered)
            mode_phi = np.angle(np.exp(1j * (phi_grid_recentered[mode_phi_recentered_idx] + phi_center)))
            mode_psi = mode_phi / 2.0
            mode_psi_canonical = wrap_angle_to_range(mode_psi, -np.pi/2, np.pi/2)
            mode_psi_degrees = np.degrees(mode_psi_canonical)
            
            psi_wrapped = wrap_samples_to_mode(psi_raw, mode_psi_canonical)
            psi_wrapped_degrees = np.degrees(psi_wrapped)
            mode_assignment = np.zeros(len(psi_raw), dtype=int)
            
            stats = compute_mode_statistics(mode_psi_degrees, psi_wrapped_degrees,
                                           weights, logl, config, None)
            
            mode_dict = {
                'value': mode_psi_canonical,
                'value_degrees': mode_psi_degrees,
                'peak_height': 2.0 * phi_pdf_recentered[mode_phi_recentered_idx],
                'posterior_mass': 1.0,
                'hdi_68': (np.radians(stats['hdi_68'][0]), np.radians(stats['hdi_68'][1])),
                'hdi_99': (np.radians(stats['hdi_99'][0]), np.radians(stats['hdi_99'][1])),
                'hdi_68_degrees': stats['hdi_68'],
                'hdi_99_degrees': stats['hdi_99'],
                'likelihood_max': stats['likelihood_max'],
                'likelihood_std': stats['likelihood_std'],
                'n_samples': stats['n_samples'],
                'n_in_hdi_68': stats.get('n_in_hdi_68', stats['n_samples']),
                'psi_raw': psi_raw,
                'psi_unwrapped': psi_wrapped,
                'mode_assignment': mode_assignment,
                'unwrapped_samples': psi_wrapped,
                'phi_grid': phi_grid_plot,
                'phi_pdf': phi_pdf_plot,
                'mode_phi': mode_phi,
                'kappa': kappa,
                'type': 'axial',
                'param_name': 'psi_0'
            }
            
            return [mode_dict]
    
    return final_modes


# =============================================================================
# MAIN PROCESSING FUNCTION
# =============================================================================

def convert_ref_samples_to_intrinsic(samples: np.ndarray,
                                     param_names: list,
                                     model,
                                     lambda_sq_ref: float) -> np.ndarray:
    """
    Convert reference-band (a_ref, b_ref) columns in posterior samples to
    intrinsic (a_int, b_int) using the per-component transfer factor.

    For each sample and each component, the shape parameters are used to
    reconstruct a temporary component, compute K = compute_transfer_factor(λ²_ref),
    and then A_int = (a_ref + i·b_ref) / K is stored back into the samples array.

    Args:
        samples: Posterior samples array (n_samples, n_params)
        param_names: Parameter names corresponding to columns of samples
        model: CompositeModel carrying component types and fixed attributes
        lambda_sq_ref: Reference wavelength squared in m²

    Returns:
        New samples array with (a, b) columns replaced by intrinsic values
    """
    samples_int = samples.copy()
    name_to_idx = {name: i for i, name in enumerate(param_names)}

    tied_phi_rm_idx = name_to_idx.get('tied_phi_rm')

    for ci, comp in enumerate(model.components):
        comp_name = comp.name if comp.name else f'component_{ci}'

        a_idx = name_to_idx.get(f'{comp_name}_a')
        b_idx = name_to_idx.get(f'{comp_name}_b')
        if a_idx is None or b_idx is None:
            continue

        a_ref = samples[:, a_idx]
        b_ref = samples[:, b_idx]

        if isinstance(comp, ThinPowerLawComponent):
            phi_rm_idx = name_to_idx.get(f'{comp_name}_phi_rm', tied_phi_rm_idx)
            beta_idx = name_to_idx.get(f'{comp_name}_beta')
            for s in range(len(samples)):
                tmp = ThinPowerLawComponent(
                    a=1.0, b=0.0,
                    phi_rm=samples[s, phi_rm_idx],
                    beta=samples[s, beta_idx],
                    lambda_sq_0=comp.lambda_sq_0
                )
                K = tmp.compute_transfer_factor(lambda_sq_ref)
                A_int = (a_ref[s] + 1j * b_ref[s]) / K
                samples_int[s, a_idx] = A_int.real
                samples_int[s, b_idx] = A_int.imag

        elif isinstance(comp, ThinComponent):
            phi_rm_idx = name_to_idx.get(f'{comp_name}_phi_rm', tied_phi_rm_idx)
            for s in range(len(samples)):
                tmp = ThinComponent(a=1.0, b=0.0, phi_rm=samples[s, phi_rm_idx])
                K = tmp.compute_transfer_factor(lambda_sq_ref)
                A_int = (a_ref[s] + 1j * b_ref[s]) / K
                samples_int[s, a_idx] = A_int.real
                samples_int[s, b_idx] = A_int.imag

        elif isinstance(comp, ThickComponent):
            phi_peak_idx = name_to_idx.get(f'{comp_name}_phi_peak')
            sigma_phi_idx = name_to_idx.get(f'{comp_name}_sigma_phi')
            N_idx = name_to_idx.get(f'{comp_name}_N')
            for s in range(len(samples)):
                tmp = ThickComponent(
                    a=1.0, b=0.0,
                    phi_peak=samples[s, phi_peak_idx],
                    sigma_phi=samples[s, sigma_phi_idx],
                    N=samples[s, N_idx]
                )
                K = tmp.compute_transfer_factor(lambda_sq_ref)
                A_int = (a_ref[s] + 1j * b_ref[s]) / K
                samples_int[s, a_idx] = A_int.real
                samples_int[s, b_idx] = A_int.imag

    return samples_int


def process_posterior_modes(results: 'FaradayFitter',
                           config: ProcessingConfig) -> 'FaradayFitter':
    """
    Process posterior samples with multi-modal detection.
    
    Handles both linear parameters (p₀, φ_rm, σ_φ, N) and axial angles (ψ₀).
    Detects multiple modes per parameter when config.detect_multiple_modes=True.
    
    Augments the results object with:
    - results.mode_results: List of modes per parameter
    - results.processed_parameters: Stats with mode_1, mode_2, ... mode_N
    - results.samples_processed: Unwrapped samples array (from mode_1)
    - results.param_summary[param]['mode']: Primary mode value (mode_1)
    - results.param_summary[param]['hdi_68']: HDI around mode_1
    
    Args:
        results: FaradayFitter object from faraday_model.py
        config: ProcessingConfig with multi-modal settings
        
    Returns:
        Augmented FaradayFitter object with mode processing results
    """
    print("\n" + "="*70)
    print("POSTERIOR MODE DETECTION & HDI COMPUTATION")
    print("="*70)
    
    # Extract data from results object
    samples = results.samples
    weights = results.weights
    param_names = results.setup.param_names
    logl = results.logl

    # Convert reference-band samples to intrinsic parameterisation if needed
    lambda_sq_ref = getattr(results.setup, 'lambda_sq_ref', None)
    if lambda_sq_ref is not None:
        print(f"Converting posterior samples to intrinsic parameterisation "
              f"(\u03bb\u00b2_ref = {lambda_sq_ref:.4f} m\u00b2)...")
        samples = convert_ref_samples_to_intrinsic(
            samples, param_names, results.setup.model, lambda_sq_ref
        )

    # Normalize weights
    weights = normalize_weights(weights)
    
    # Create output directory if plotting
    if config.plot_diagnostics:
        config.plot_dir.mkdir(parents=True, exist_ok=True)
        print(f"Diagnostic plots will be saved to: {config.plot_dir}")
    
    mode_results = {}  # Store mode detection results here

    # Compute effective noise on polarized amplitude via inverse-variance combination
    # across channels: σ_eff = 1/sqrt(Σ 1/σ_i²)  →  gives sqrt(N) improvement for uniform noise
    _q_err = results.setup.data.Q_err
    _u_err = results.setup.data.U_err
    _q_valid = _q_err[_q_err > 0]
    _u_valid = _u_err[_u_err > 0]
    if len(_q_valid) > 0 and len(_u_valid) > 0:
        sigma_q_eff = 1.0 / np.sqrt(np.sum(1.0 / _q_valid**2))
        sigma_u_eff = 1.0 / np.sqrt(np.sum(1.0 / _u_valid**2))
        sigma_p_sq = 0.8 * max(sigma_q_eff, sigma_u_eff)**2 + 0.2 * min(sigma_q_eff, sigma_u_eff)**2
        sigma_p = np.sqrt(sigma_p_sq)
        print(f"\nMAS debiasing: σ_q_eff={sigma_q_eff:.4e}  σ_u_eff={sigma_u_eff:.4e}  "
              f"σ_p={sigma_p:.4e}  (using {len(_q_valid)} Q / {len(_u_valid)} U channels)")
    else:
        sigma_p = None
        sigma_p_sq = None
        print("\nMAS debiasing: Q_err/U_err not available — skipping bias correction")

    print(f"\nProcessing {len(param_names)} sampled parameters...")
    
    # Identify parameter types
    # Look for (a, b) pairs to compute ψ₀
    # Use endswith() to avoid confusion with parameters like *_beta
    a_indices = [i for i, name in enumerate(param_names) if name.endswith('_a')]
    b_indices = [i for i, name in enumerate(param_names) if name.endswith('_b')]
    
    # Process each parameter
    for i, param_name in enumerate(param_names):
        x = samples[:, i]
        
        # Skip fixed parameters entirely - no processing needed
        if i in results.setup.fixed_params:
            print(f"\n{param_name}: FIXED at {results.setup.fixed_params[i]:.6g} (excluded from processing)")
            continue
        
        # Skip if this is a 'b' parameter (will be processed with matching 'a')
        # Use endswith() to avoid confusion with *_beta
        if param_name.endswith('_b'):
            continue
        
        # Check if this is an 'a' parameter → compute ψ₀
        if param_name.endswith('_a') and i in a_indices:
            # Find matching 'b' parameter
            comp_name = param_name.replace('_a', '')
            b_name = comp_name + '_b'
            
            if b_name in param_names:
                b_idx = param_names.index(b_name)
                
                # Skip derived parameters if a or b are fixed
                if i in results.setup.fixed_params or b_idx in results.setup.fixed_params:
                    fixed_val_a = results.setup.fixed_params.get(i)
                    fixed_val_b = results.setup.fixed_params.get(b_idx)
                    print(f"\n{comp_name}_p0 and {comp_name}_psi_0: derived from fixed (a={fixed_val_a}, b={fixed_val_b}) - excluded from processing")
                    continue
                
                a_samples = samples[:, i]
                b_samples = samples[:, b_idx]
                
                # Process ψ₀ (axial angle)
                psi_name = comp_name + '_psi_0'
                print(f"\n{psi_name} (polarization angle):")
                result_list = process_axial_parameter(a_samples, b_samples, 
                                                weights, config, logl, param_name=psi_name)
                mode_results[psi_name] = result_list
                
                for mode_idx, mode_result in enumerate(result_list, 1):
                    print(f"  Mode {mode_idx}: {mode_result['value_degrees']:.4f}° | HDI 68%: [{mode_result['hdi_68_degrees'][0]:.2f}°, {mode_result['hdi_68_degrees'][1]:.2f}°]")
                    print(f"    Posterior mass: {mode_result['posterior_mass']:.3f}")
                    print(f"    Max log(L): {mode_result['likelihood_max']:.2f}")
                    print(f"    Std log(L): {mode_result['likelihood_std']:.2f}")
                    print(f"    N assigned: {mode_result['n_samples']}, N in HDI_68: {mode_result.get('n_in_hdi_68', mode_result['n_samples'])}")
                
                # Process p₀ (polarization fraction) with MAS bias correction
                p_raw = np.sqrt(a_samples**2 + b_samples**2)
                if sigma_p is not None and sigma_p > 0 and np.median(p_raw) / sigma_p >= 3.0:
                    p0_samples = np.sqrt(np.maximum(p_raw**2 - sigma_p_sq, 0.0))
                    print(f"\n  MAS correction applied to all samples "
                          f"(median SNR={np.median(p_raw)/sigma_p:.2f} >= 3)")
                else:
                    p0_samples = p_raw
                    if sigma_p is not None and sigma_p > 0:
                        print(f"\n  MAS correction skipped "
                              f"(median SNR={np.median(p_raw)/sigma_p:.2f} < 3)")
                p0_name = comp_name + '_p0'
                print(f"\n\n{p0_name} (polarization fraction):")
                p0_result_list = process_linear_parameter(p0_samples, weights,
                                                    p0_name, config, logl)
                mode_results[p0_name] = p0_result_list
                
                for mode_idx, mode_result in enumerate(p0_result_list, 1):
                    print(f"  Mode {mode_idx}: {mode_result['value']:.6f} | HDI 68%: [{mode_result['hdi_68'][0]:.6f}, {mode_result['hdi_68'][1]:.6f}]")
                    print(f"    Posterior mass: {mode_result['posterior_mass']:.3f}")
                    print(f"    Max log(L): {mode_result['likelihood_max']:.2f}")
                    print(f"    Std log(L): {mode_result['likelihood_std']:.2f}")
                    print(f"    N assigned: {mode_result['n_samples']}, N in HDI_68: {mode_result.get('n_in_hdi_68', mode_result['n_samples'])}")
        
        # Process other linear parameters (φ_rm, σ_φ, N, β)
        elif not param_name.endswith('_a') and not param_name.endswith('_b'):
            # Add descriptive label for parameter
            if 'phi_rm' in param_name or 'phi_peak' in param_name:
                description = "(rotation measure)"
            elif 'sigma_phi' in param_name:
                description = "(Faraday dispersion)"
            elif param_name.endswith('_N'):
                description = "(shape parameter)"
            elif 'beta' in param_name:
                description = "(spectral index)"
            else:
                description = ""
            
            print(f"\n\n{param_name} {description}:".replace(" :", ":"))
            result_list = process_linear_parameter(x, weights, param_name, config, logl)
            mode_results[param_name] = result_list
            
            for mode_idx, mode_result in enumerate(result_list, 1):
                print(f"  Mode {mode_idx}: {mode_result['value']:.4f} | HDI 68%: [{mode_result['hdi_68'][0]:.4f}, {mode_result['hdi_68'][1]:.4f}]")
                print(f"    Posterior mass: {mode_result['posterior_mass']:.3f}")
                print(f"    Max log(L): {mode_result['likelihood_max']:.2f}")
                print(f"    Std log(L): {mode_result['likelihood_std']:.2f}")
                print(f"    N assigned: {mode_result['n_samples']}, N in HDI_68: {mode_result.get('n_in_hdi_68', mode_result['n_samples'])}")
    
    print(f"\nProcessed {len(mode_results)} parameters")
    
    # Generate diagnostic plots if requested
    if config.plot_diagnostics:
        print(f"\nGenerating diagnostic plots...")
        _save_diagnostic_plots(mode_results, samples, weights, param_names, config)
    
    # Store mode_results in results object
    results.mode_results = mode_results
    
    # Build samples_processed (unwrapped) and compute comprehensive statistics
    print(f"\nComputing comprehensive statistics from unwrapped samples...")
    samples_processed_list = []
    param_names_processed_list = []
    processed_parameters = {}
    
    # Filter samples with negligible weights to prevent overflow in variance computation
    weight_threshold = np.max(weights) * 1e-10
    weight_mask = weights > weight_threshold
    
    if np.sum(weight_mask) < 100:
        # If too few samples pass threshold, use top N samples instead
        n_keep = min(10000, len(weights))
        top_indices = np.argsort(weights)[-n_keep:]
        weight_mask = np.zeros(len(weights), dtype=bool)
        weight_mask[top_indices] = True
    
    weights_filtered = weights[weight_mask]
    weight_sum_filtered = np.sum(weights_filtered)
    if weight_sum_filtered <= 0:
        raise ValueError("No valid weights after filtering - this should not happen")
    weights_filtered = weights_filtered / weight_sum_filtered
    
    for param_name, mode_result_list in mode_results.items():
        if len(mode_result_list) > 0 and 'unwrapped_samples' in mode_result_list[0]:
            unwrapped = mode_result_list[0]['unwrapped_samples']
            samples_processed_list.append(unwrapped)
            param_names_processed_list.append(param_name)
            
            unwrapped_filtered = unwrapped[weight_mask]
            
            quantiles = [0.005, 0.16, 0.5, 0.84, 0.995]
            with np.errstate(over='ignore', invalid='ignore'):
                q_values = dyfunc.quantile(unwrapped_filtered, quantiles, weights=weights_filtered)
            
            with np.errstate(over='ignore', invalid='ignore'):
                mean, cov = dyfunc.mean_and_cov(unwrapped_filtered.reshape(-1, 1), weights=weights_filtered)
            mean_val = float(mean[0])
            cov_val = np.clip(cov[0, 0], 0, 1e10)
            std_val = float(np.sqrt(cov_val))
            
            processed_parameters[param_name] = {
                'mean': mean_val,
                'median': float(q_values[2]),
                'std': std_val,
                'eti_68': [float(q_values[1]), float(q_values[3])],
                'eti_99': [float(q_values[0]), float(q_values[4])]
            }
            
            for mode_idx, mode_result in enumerate(mode_result_list, 1):
                mode_key = f'mode_{mode_idx}'
                
                mode_val = mode_result.get('value_degrees', mode_result.get('value'))
                hdi_68_val = mode_result.get('hdi_68_degrees', mode_result.get('hdi_68'))
                hdi_99_val = mode_result.get('hdi_99_degrees', mode_result.get('hdi_99'))
                peak_height = mode_result.get('peak_height', np.nan)
                
                processed_parameters[param_name][mode_key] = {
                    'mode': float(mode_val),
                    'peak_height': float(peak_height),
                    'hdi_68': [float(hdi_68_val[0]), float(hdi_68_val[1])],
                    'hdi_99': [float(hdi_99_val[0]), float(hdi_99_val[1])],
                    'posterior_mass': float(mode_result['posterior_mass']),
                    'likelihood_max': float(mode_result['likelihood_max']),
                    'likelihood_std': float(mode_result['likelihood_std']),
                    'n_samples_in_hdi': int(mode_result['n_samples'])
                }
    
    # Stack samples_processed if we have any
    if samples_processed_list:
        # Reorder parameters for better corner plot layout:
        # For each component: p0, psi_0, phi_rm/phi_peak, (sigma_phi, N)
        # Also filter out fixed parameters
        
        reordered_samples = []
        reordered_names = []
        
        # Get fixed parameters from setup
        fixed_indices_in_cartesian = set(results.setup.fixed_params.keys())
        fixed_params_derived = set()
        
        # Identify which derived parameters correspond to fixed Cartesian ones
        for i, cart_name in enumerate(results.setup.param_names):
            if i in fixed_indices_in_cartesian:
                # Map to derived names
                # Use endswith() to avoid confusion with *_beta
                if cart_name.endswith('_a') or cart_name.endswith('_b'):
                    # a or b fixed -> p0 and psi_0 are fixed
                    comp_prefix = cart_name.rsplit('_', 1)[0]  # e.g., 'S1' from 'S1_a'
                    fixed_params_derived.add(f'{comp_prefix}_p0')
                    fixed_params_derived.add(f'{comp_prefix}_psi_0')
                else:
                    # Direct parameter like phi_rm, beta, etc. - map to itself
                    fixed_params_derived.add(cart_name)
        
        # Group parameters by component using suffix-based parsing (consistent with model)
        # Maintain same suffix list as model: ['_a','_b','_phi_rm','_phi_peak','_sigma_phi','_N','_beta']
        component_params = {}  # {comp_name: [(param_type, idx), ...]}
        global_params = {}  # {param_name: (param_type, idx)} for tied_phi_rm, etc.
        
        # Define parameter suffixes and global parameters (matching model logic)
        component_suffixes = [
            ('_p0', 'p0'),
            ('_psi_0', 'psi_0'),
            ('_phi_rm', 'phi_rm'),
            ('_phi_peak', 'phi_peak'),
            ('_sigma_phi', 'sigma_phi'),
            ('_N', 'N'),
            ('_beta', 'beta'),
        ]
        global_param_names = ['tied_phi_rm']
        
        for i, name in enumerate(param_names_processed_list):
            # Skip fixed parameters (already excluded but double-check)
            if name in fixed_params_derived:
                continue
            
            # Check if this is a global parameter (not component-specific)
            if name in global_param_names:
                # Store in global bucket
                global_params[name] = ('global', i)
                continue
            
            # Extract component name and parameter type using suffix matching
            comp_name = None
            param_type = None
            
            for suffix, ptype in component_suffixes:
                if name.endswith(suffix):
                    comp_name = name[:-len(suffix)]
                    param_type = ptype
                    break
            
            if comp_name is None:
                # Unknown parameter, keep as-is
                comp_name = name
                param_type = 'unknown'
            
            if comp_name not in component_params:
                component_params[comp_name] = []
            
            component_params[comp_name].append((param_type, i))
        
        # Define parameter order priority
        # For each component: p0, psi_0, then global params (on first component only), then other params
        param_order = ['p0', 'psi_0', 'phi_rm', 'phi_peak', 'sigma_phi', 'N', 'beta']
        
        # Add parameters in correct order for each component
        sorted_components = sorted(component_params.keys())
        for comp_idx, comp_name in enumerate(sorted_components):
            # First, add p0 and psi_0
            for param_type in ['p0', 'psi_0']:
                for stored_type, idx in component_params[comp_name]:
                    if stored_type == param_type:
                        reordered_samples.append(samples_processed_list[idx])
                        reordered_names.append(param_names_processed_list[idx])
            
            # After first component's p0/psi_0, add global parameters (e.g., tied_phi_rm)
            if comp_idx == 0 and global_params:
                for param_name in sorted(global_params.keys()):
                    param_type, idx = global_params[param_name]
                    reordered_samples.append(samples_processed_list[idx])
                    reordered_names.append(param_names_processed_list[idx])
            
            # Then add remaining component-specific parameters
            for param_type in ['phi_rm', 'phi_peak', 'sigma_phi', 'N', 'beta']:
                for stored_type, idx in component_params[comp_name]:
                    if stored_type == param_type:
                        reordered_samples.append(samples_processed_list[idx])
                        reordered_names.append(param_names_processed_list[idx])
        
        if reordered_samples:
            results.samples_processed = np.column_stack(reordered_samples)
            results.param_names_processed = reordered_names
            
            n_fixed = len(fixed_params_derived)
            if n_fixed > 0:
                print(f"\nFiltered {n_fixed} fixed parameters from processed samples")
                print(f"Fixed parameters: {sorted(fixed_params_derived)}")
        else:
            results.samples_processed = None
            results.param_names_processed = None
    else:
        results.samples_processed = None
        results.param_names_processed = None
    
    # Store processed_parameters
    results.processed_parameters = processed_parameters
    
    # Also add mode/hdi to param_summary for backward compatibility
    for param_name, mode_result_list in mode_results.items():
        if len(mode_result_list) > 0:
            mode_result = mode_result_list[0]
            
            if param_name not in results.param_summary:
                results.param_summary[param_name] = {}
            
            # CRITICAL: For psi_0, use degrees; for other parameters, use value as-is
            if 'psi_0' in param_name:
                mode_val = mode_result.get('value_degrees', mode_result.get('value'))
            else:
                mode_val = mode_result.get('value')
            results.param_summary[param_name]['mode'] = mode_val
            
            # CRITICAL: For psi_0, use degrees; for other parameters, use value as-is
            if 'psi_0' in param_name:
                if 'hdi_68_degrees' in mode_result:
                    results.param_summary[param_name]['hdi_68'] = mode_result['hdi_68_degrees']
                elif 'hdi_68' in mode_result:
                    results.param_summary[param_name]['hdi_68'] = mode_result['hdi_68']
                
                if 'hdi_99_degrees' in mode_result:
                    results.param_summary[param_name]['hdi_99'] = mode_result['hdi_99_degrees']
                elif 'hdi_99' in mode_result:
                    results.param_summary[param_name]['hdi_99'] = mode_result['hdi_99']
            else:
                if 'hdi_68' in mode_result:
                    results.param_summary[param_name]['hdi_68'] = mode_result['hdi_68']
                
                if 'hdi_99' in mode_result:
                    results.param_summary[param_name]['hdi_99'] = mode_result['hdi_99']
    
    # =======================================================================
    # FINAL STEP: Convert psi_0 from radians to degrees for ALL outputs
    # All processing done in radians (native units), convert only at the very end
    # =======================================================================
    
    # Convert samples_processed
    if results.samples_processed is not None and results.param_names_processed is not None:
        for i, name in enumerate(results.param_names_processed):
            if 'psi_0' in name:
                results.samples_processed[:, i] = np.degrees(results.samples_processed[:, i])
    
    # Convert processed_parameters stats
    for param_name in processed_parameters:
        if 'psi_0' in param_name:
            processed_parameters[param_name]['mean'] = np.degrees(processed_parameters[param_name]['mean'])
            processed_parameters[param_name]['std'] = np.degrees(processed_parameters[param_name]['std'])
            processed_parameters[param_name]['median'] = np.degrees(processed_parameters[param_name]['median'])
            processed_parameters[param_name]['eti_68'] = [
                np.degrees(processed_parameters[param_name]['eti_68'][0]),
                np.degrees(processed_parameters[param_name]['eti_68'][1])
            ]
            processed_parameters[param_name]['eti_99'] = [
                np.degrees(processed_parameters[param_name]['eti_99'][0]),
                np.degrees(processed_parameters[param_name]['eti_99'][1])
            ]
    
    print("\n" + "="*70)
    print("MODE PROCESSING COMPLETE")
    print("="*70)
    print(f"Mode detection: results.mode_results")
    print(f"Processed stats: results.processed_parameters")
    print(f"Unwrapped samples: results.samples_processed")
    print(f"Updated summaries: results.param_summary (includes mode/HDI)")
    
    return results


# =============================================================================
# DIAGNOSTIC PLOTTING
# =============================================================================

def plot_linear_diagnostic(x: np.ndarray, weights: np.ndarray,
                          result_list: List[Dict], config: ProcessingConfig) -> plt.Figure:
    """
    Create TWO-PANEL diagnostic plot for linear parameter.
    
    TOP PANEL: Raw input samples (no color coding)
    BOTTOM PANEL: Samples COLOR-CODED by mode assignment
    
    Args:
        x: Sample values
        weights: Normalized weights
        result_list: List of mode dictionaries from process_linear_parameter
        config: ProcessingConfig
        
    Returns:
        Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    
    result = result_list[0]
    param_name = result['param_name']
    param_label = format_param_label(param_name)
    
    # Resample for visualization
    x_resampled = resample_equal(x[:, np.newaxis], weights)[:, 0]
    
    # Adaptive binning
    q75, q25 = np.percentile(x_resampled, [75, 25])
    iqr = q75 - q25
    if iqr > 0:
        bin_width = 2 * iqr / (len(x_resampled) ** (1/3))
        data_range = x_resampled.max() - x_resampled.min()
        n_bins = max(50, int(np.ceil(data_range / bin_width)))
        n_bins = min(n_bins, 100)
    else:
        n_bins = 50
    
    # =======================================================================
    # TOP PANEL: RAW INPUT SAMPLES (NO COLOR)
    # =======================================================================
    ax = axes[0]
    ax.hist(x_resampled, bins=n_bins, density=True, alpha=0.6, 
           color='steelblue', label='Resampled')
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'{param_label} - Raw Samples', fontsize=13)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # =======================================================================
    # BOTTOM PANEL: SAMPLES COLORED BY MODE ASSIGNMENT
    # =======================================================================
    ax = axes[1]
    
    # Check if multi-modal
    has_mode_assignment = 'mode_assignment' in result and len(result_list) > 1
    
    if has_mode_assignment:
        # Multi-modal: use original samples with importance weights
        mode_assignment = result['mode_assignment']
        
        # Plot each mode with its own color and adaptive binning
        mode_colors = ['green', 'orange', 'purple', 'cyan', 'magenta']
        
        # Shared bins across all modes for consistent scaling
        x_min, x_max = result['grid'][0], result['grid'][-1]
        n_bins_global = 60  # Fixed bin count for consistency
        bins = np.linspace(x_min, x_max, n_bins_global + 1)
        bin_width = bins[1] - bins[0]
        
        for mode_idx, mode_result in enumerate(result_list):
            mask = (mode_assignment == mode_idx)
            n_mode = np.sum(mask)
            
            if n_mode > 0:
                color = mode_colors[mode_idx % len(mode_colors)]
                mass = mode_result['posterior_mass']
                x_mode = x[mask]
                weights_mode = weights[mask]
                
                # Histogram with importance weights, density=False (mass per bin)
                ax.hist(x_mode, bins=bins, density=False, weights=weights_mode, 
                       alpha=0.6, color=color, label=f'Mode {mode_idx + 1} (mass={mass:.2f})')
    else:
        # Single mode: just one color
        ax.hist(x_resampled, bins=n_bins, density=True, alpha=0.6,
               color='green', label=f'Mode 1 (mass=1.00)')
    
    # Overlay KDE (scale to match histogram if using weighted counts)
    if has_mode_assignment:
        # Scale KDE to mass per bin: pdf * total_weight * bin_width
        kde_scaled = result['pdf'] * np.sum(weights) * bin_width
        ax.plot(result['grid'], kde_scaled, 'k-', linewidth=2, label='Weighted KDE')
        ylabel = 'Posterior mass per bin'
    else:
        # Single mode uses density, keep KDE as-is
        ax.plot(result['grid'], result['pdf'], 'k-', linewidth=2, label='Weighted KDE')
        ylabel = 'Density'
    
    # Mark mode positions
    mode_line_colors = ['red', 'blue', 'green', 'orange', 'purple']
    for mode_idx, mode_result in enumerate(result_list):
        mode = mode_result['value']
        mass = mode_result['posterior_mass']
        color = mode_line_colors[mode_idx % len(mode_line_colors)]
        
        ax.axvline(mode, color=color, linestyle='--', linewidth=2, 
                  label=f"Mode {mode_idx+1} = {mode:.4f} (mass={mass:.2f})")
    
    # HDI for primary mode
    if 'hdi_68' in result:
        low, high = result['hdi_68']
        ax.axvspan(low, high, alpha=0.3, color='lightblue', label='68% HDI (mode 1)')
    
    if 'hdi_99' in result:
        low, high = result['hdi_99']
        ax.axvspan(low, high, alpha=0.2, color='lightblue', label='99% HDI (mode 1)')
    
    ax.set_xlabel(param_label, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f'{param_label} - Colored by Mode', fontsize=13)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(result['grid'][0], result['grid'][-1])
    
    fig.tight_layout()
    
    return fig


def plot_axial_diagnostic(result_list: List[Dict], weights: np.ndarray,
                         config: ProcessingConfig) -> plt.Figure:
    """
    Create TWO-PANEL diagnostic plot for axial parameter ψ₀.
    
    TOP PANEL: Raw ψ₀ samples (no color coding, may be split at ±π/2 boundaries)
    BOTTOM PANEL: UNWRAPPED ψ₀ samples COLOR-CODED by mode assignment
    
    Args:
        result_list: List of mode dictionaries from process_axial_parameter
        weights: Normalized weights
        config: ProcessingConfig
        
    Returns:
        Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=False)
    
    result = result_list[0]
    
    psi_raw = result['psi_raw']
    psi_unwrapped = result['unwrapped_samples']
    mode = result['value']
    
    # Check if multi-modal (has mode_assignment)
    has_mode_assignment = 'mode_assignment' in result and len(result_list) > 1
    
    if has_mode_assignment:
        mode_assignment = result['mode_assignment']
        
        # Resample indices to preserve pairing between psi_unwrapped and mode_assignment
        n_resample = 10000
        weights_norm = weights / np.sum(weights)
        resampled_indices = np.random.choice(len(psi_unwrapped), size=n_resample, 
                                            replace=True, p=weights_norm)
        
        psi_unwrapped_resampled = psi_unwrapped[resampled_indices]
        mode_assignment_resampled = mode_assignment[resampled_indices]
    else:
        # Single mode case
        psi_unwrapped_resampled = resample_equal(psi_unwrapped[:, np.newaxis], weights)[:, 0]
        mode_assignment_resampled = None
    
    psi_raw_resampled = resample_equal(psi_raw[:, np.newaxis], weights)[:, 0]
    
    # Adaptive binning for raw
    q75, q25 = np.percentile(psi_raw_resampled, [75, 25])
    iqr = q75 - q25
    if iqr > 0:
        bin_width = 2 * iqr / (len(psi_raw_resampled) ** (1/3))
        data_range = psi_raw_resampled.max() - psi_raw_resampled.min()
        n_bins_raw = max(50, int(np.ceil(data_range / bin_width)))
        n_bins_raw = min(n_bins_raw, 100)
    else:
        n_bins_raw = 50
    
    # Adaptive binning for unwrapped
    q75, q25 = np.percentile(psi_unwrapped_resampled, [75, 25])
    iqr = q75 - q25
    if iqr > 0:
        bin_width = 2 * iqr / (len(psi_unwrapped_resampled) ** (1/3))
        data_range = psi_unwrapped_resampled.max() - psi_unwrapped_resampled.min()
        n_bins_unwrapped = max(50, int(np.ceil(data_range / bin_width)))
        n_bins_unwrapped = min(n_bins_unwrapped, 100)
    else:
        n_bins_unwrapped = 50
    
    # =======================================================================
    # TOP PANEL: RAW ψ₀ SAMPLES (NO COLOR)
    # =======================================================================
    ax = axes[0]
    ax.hist(psi_raw_resampled, bins=n_bins_raw, density=True, alpha=0.6,
           color='steelblue', label='Raw ψ₀')
    ax.axvline(-np.pi/2, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax.axvline(np.pi/2, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Before Wrapping (raw ψ₀)', fontsize=13)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # =======================================================================
    # BOTTOM PANEL: UNWRAPPED ψ₀ COLORED BY MODE ASSIGNMENT
    # =======================================================================
    ax = axes[1]
    
    # Get range for consistent binning from original unwrapped samples
    psi_min = np.min(psi_unwrapped)
    psi_max = np.max(psi_unwrapped)
    
    # Color by mode assignment if multi-modal
    if has_mode_assignment:
        mode_colors = ['green', 'orange', 'purple', 'cyan', 'magenta']
        
        # Shared bins across all modes for consistent scaling
        n_bins_global = 60
        bins = np.linspace(psi_min, psi_max, n_bins_global + 1)
        bin_width = bins[1] - bins[0]
        
        for mode_idx, mode_result in enumerate(result_list):
            mask = (mode_assignment == mode_idx)
            n_mode = np.sum(mask)
            
            if n_mode > 0:
                color = mode_colors[mode_idx % len(mode_colors)]
                mass = mode_result['posterior_mass']
                psi_mode = psi_unwrapped[mask]
                weights_mode = weights[mask]
                
                # Histogram with importance weights, density=False (mass per bin)
                ax.hist(psi_mode, bins=bins, density=False, weights=weights_mode, 
                       alpha=0.6, color=color, label=f'Mode {mode_idx + 1} (mass={mass:.2f})')
    else:
        ax.hist(psi_unwrapped_resampled, bins=n_bins_unwrapped, density=True, alpha=0.6,
               color='green', label='Wrapped ψ₀')
    
    # Plot KDE (scale to match histogram if using weighted counts)
    if 'phi_grid' in result and 'phi_pdf' in result:
        psi_grid = result['phi_grid'] / 2.0
        psi_pdf_base = 2.0 * result['phi_pdf']
        
        if has_mode_assignment:
            # Scale KDE to mass per bin: pdf * total_weight * bin_width
            psi_pdf = psi_pdf_base * np.sum(weights) * bin_width
            ylabel_unwrapped = 'Posterior mass per bin'
        else:
            # Single mode uses density
            psi_pdf = psi_pdf_base
            ylabel_unwrapped = 'Density'
        
        ax.plot(psi_grid, psi_pdf, 'k-', linewidth=2, label='von Mises KDE')
        
        psi_min, psi_max = -np.pi/2, np.pi/2
        if psi_grid.min() > psi_min:
            pad_left = np.linspace(psi_min, psi_grid.min(), 10)
            ax.plot(pad_left, np.zeros_like(pad_left), 'k-', linewidth=2, alpha=0.3)
        if psi_grid.max() < psi_max:
            pad_right = np.linspace(psi_grid.max(), psi_max, 10)
            ax.plot(pad_right, np.zeros_like(pad_right), 'k-', linewidth=2, alpha=0.3)
    
    mode_colors = ['red', 'blue', 'green', 'orange', 'purple']
    for mode_idx, mode_result in enumerate(result_list):
        mode_val = mode_result['value']
        mass = mode_result['posterior_mass']
        color = mode_colors[mode_idx] if mode_idx < len(mode_colors) else 'gray'
        
        ax.axvline(mode_val, color=color, linestyle='--', linewidth=2, 
                  label=f"Mode {mode_idx+1} = {mode_val:.4f} rad (mass={mass:.2f})")
    
    if 'hdi_68' in result:
        low, high = result['hdi_68']
        ax.axvspan(low, high, alpha=0.3, color='lightblue', label='68% HDI (mode 1)')
    
    if 'hdi_99' in result:
        low, high = result['hdi_99']
        ax.axvspan(low, high, alpha=0.2, color='lightblue', label='99% HDI (mode 1)')
    
    ax.axvline(-np.pi/2, color='red', linestyle=':', linewidth=2, alpha=0.7, label='Boundaries')
    ax.axvline(np.pi/2, color='red', linestyle=':', linewidth=2, alpha=0.7)
    ax.set_xlabel('ψ₀ (rad)', fontsize=12)
    ax.set_ylabel(ylabel_unwrapped if 'ylabel_unwrapped' in locals() else 'Density', fontsize=12)
    ax.set_title('After Wrapping (continuous around mode)', fontsize=13)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()
    
    return fig


def _save_diagnostic_plots(mode_results: Dict, samples: np.ndarray,
                          weights: np.ndarray, param_names: List[str],
                          config: ProcessingConfig):
    """
    Generate and save diagnostic plots showing all detected modes.
    
    Args:
        mode_results: Dictionary mapping param names to lists of mode dicts
        samples: Original samples array
        weights: Normalized weights
        param_names: Parameter names
        config: ProcessingConfig
    """
    if not config.plot_diagnostics:
        return
    
    for param_name, result_list in mode_results.items():
        if len(result_list) == 0:
            continue
        
        result = result_list[0]
        
        if result['type'] == 'linear':
            if param_name.endswith('_p0'):
                comp_name = param_name[:-3]
                a_name = comp_name + '_a'
                b_name = comp_name + '_b'
                if a_name in param_names and b_name in param_names:
                    a_idx = param_names.index(a_name)
                    b_idx = param_names.index(b_name)
                    x = np.sqrt(samples[:, a_idx]**2 + samples[:, b_idx]**2)
                else:
                    continue
            elif param_name in param_names:
                idx = param_names.index(param_name)
                x = samples[:, idx]
            else:
                continue
            
            fig = plot_linear_diagnostic(x, weights, result_list, config)
            
            if config.filename_prefix:
                filename = f"{config.filename_prefix}_{param_name}_mode.png"
            else:
                filename = f"{param_name}_mode.png"
            
            filepath = config.plot_dir / filename
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved diagnostic: {filepath}")
        
        elif result['type'] == 'axial':
            fig = plot_axial_diagnostic(result_list, weights, config)
            
            if config.filename_prefix:
                filename = f"{config.filename_prefix}_{param_name}_mode.png"
            else:
                filename = f"{param_name}_mode.png"
            
            filepath = config.plot_dir / filename
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved diagnostic: {filepath}")
