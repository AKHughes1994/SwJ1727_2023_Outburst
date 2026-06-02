# faraday_data.py
"""
Data loading, handling, and I/O for Faraday rotation analysis.

This module provides:
- PolarizationData container
- Multiple data loading modes (physical/fractional)
- Stokes I fitting (power-law, polynomial)
- Data validation and error propagation
- Data saving utilities
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, Literal, Union, List, Dict
from pathlib import Path
from scipy.optimize import curve_fit
from scipy import stats
from faraday_utils import C_LIGHT

@dataclass
class FlaggingConfig:
    """
    Configuration for channel flagging during data loading.

    Flagging is applied in this order:
    0. Flag non-finite / invalid values (deterministic)
    1. Flag gross outliers in all Stokes I/Q/U/V flux and errors (deterministic)
    2. Clip extreme Stokes I values (deterministic)
    3. Flag frequency ranges (deterministic)
    4. Sigma-clip error outliers (statistical, iterative)
    5. Spectral baseline clipping on specified Stokes (statistical, iterative)
    6. Sigma-clip Stokes V outliers (significance-based, iterative)

    Attributes
    ----------
    sigma_clip_extreme : float or None
        Sigma threshold for single-pass gross outlier rejection applied to
        every available Stokes flux density (I, Q, U, V) and every error
        array (dI, dQ, dU, dV) independently.  Flags channels where
        |value - median| > sigma_clip_extreme * MAD_sigma.
        Intended to catch corrupt/RFI channels that would otherwise distort
        the statistics of later iterative steps.
        - 50.0 : Flag values more than 50 MAD-sigma from the median (default)
        - None : Disable this step
        Default: 50.0

    clip_stokes_i : int, list of int, or None
        Clip extreme Stokes I values. Examples:
        - 2: Clip 2 highest I values
        - -3: Clip 3 lowest I values
        - [-1, 2]: Clip 1 lowest AND 2 highest
        - None: No clipping (default)

    freq_ranges : str, list of str, or None
        Frequency ranges to flag in GHz. Examples:
        - '1.7~1.8': Single range
        - ['1.7~1.8', '0.5~0.63']: Multiple ranges
        - None: No range flagging (default)

    sigma_clip_errors : float or None
        Sigma threshold for error outliers (iterative). Examples:
        - 2.0: Flag errors > median + 2σ (linear) OR log-errors > median + 2σ (log)
        - None: No error clipping (default)

    sigma_clip_errors_space : str
        Error clipping space: 'linear' or 'log' (also accepts 'logspace'/'ln').
        - 'linear': clip on error values directly (original behaviour)
        - 'log': clip on log(errors) (often better for lognormal-like errors)
        Default: 'linear'

    sigma_clip_spectrum : float or None
        Sigma threshold for spectral baseline clipping (iterative).
        Residuals are evaluated as (see spectrum_whiten_residuals):
            z = (X - X_hat) / MAD_sigma(X - X_hat)          [whitened, default]
            z = (X - X_hat) / sqrt(dX^2 + spectrum_jitter^2) [formal errors]
        Examples:
        - 5.0: Flag |z| > 5σ (two-sided, unless spectrum_clip_sided changes)
        - None: No spectral baseline clipping (default)

    spectrum_stokes : str, list of str, or None
        Which Stokes parameter(s) to apply spectral baseline clipping to.
        Examples:
        - 'I': Only Stokes I
        - ['I']: Only Stokes I (list form)
        - ['Q', 'U']: Both Q and U
        - ['i', 'q', 'u']: Case-insensitive
        - 'QU': Fit a polynomial to sqrt(Q²+U²) (linear polarization amplitude).
                Useful when high-RM wrapping makes per-parameter polynomial fitting
                unreliable. Cannot be combined with 'Q' or 'U' individually.
        - ['I', 'QU']: Clip I with polynomial baseline and sqrt(Q²+U²) amplitude
        - None: Defaults to ['I']
        Default: ['I']

    spectrum_deg : int
        Degree of the log-polynomial baseline used for spectral clipping:
            log(X) = poly( log(nu) )
        Typical: 2--4
        Default: 3

    spectrum_jitter : float
        Absolute jitter floor (same units as Stokes parameter) added in quadrature
        in the spectral residual denominator.
        Default: 0.0

    spectrum_clip_sided : str
        Which side(s) to clip for spectral residuals:
        - 'two'  : |z| > threshold
        - 'low'  : z < -threshold (dips-only)
        - 'high' : z > threshold (spikes-only)
        Default: 'two'

    spectrum_whiten_residuals : bool
        If True (default), whiten the spectral residuals before clipping by
        normalising them with the MAD-estimated scatter of the residuals rather
        than the formal per-channel errors:
            z = (X - X_hat) / MAD_sigma(X - X_hat)
        This makes the sigma threshold directly interpretable as a Gaussian
        significance regardless of how well-calibrated the formal errors are.
        If False, the formal per-channel errors are used instead:
            z = (X - X_hat) / sqrt(dX^2 + spectrum_jitter^2)
        Default: True

    sigma_clip_v : float, 'auto', or None
        Sigma threshold for Stokes V significance clipping (iterative).
        Tests: |(V - reference) / σ_V| > threshold
        Examples:
        - 3.0: Flag significance > 3σ (manual threshold)
        - 'auto': Automatic threshold for 1 expected false positive
        - None: No V clipping (default)

    clip_v_against : str
        Reference for V clipping: 'median' or 'zero'
        Default: 'median'

    max_iter : int
        Maximum iterations for statistical clipping (errors, spectrum, and V)
        Default: 10

    verbose : bool
        Print flagging statistics
        Default: True
    """
    sigma_clip_extreme: Optional[float] = 50.0

    clip_stokes_i: Optional[Union[int, List[int]]] = None
    freq_ranges: Optional[Union[str, List[str]]] = None

    sigma_clip_errors: Optional[float] = None
    sigma_clip_errors_space: str = 'linear'

    sigma_clip_spectrum: Optional[float] = None
    spectrum_stokes: Optional[Union[str, List[str]]] = None
    spectrum_deg: int = 3
    spectrum_jitter: float = 0.0
    spectrum_clip_sided: str = 'two'
    spectrum_whiten_residuals: bool = True

    sigma_clip_v: Optional[Union[float, str]] = None
    clip_v_against: str = 'median'

    max_iter: int = 10
    verbose: bool = True


@dataclass
class PolarizationData:
    """
    Container for polarization observations or simulations.
    
    Attributes:
        lambda_sq (np.ndarray): Wavelength squared in m². Shape: (n_chan,)
        Q (np.ndarray): Stokes Q (fractional if I normalized). Shape: (n_chan,)
        U (np.ndarray): Stokes U (fractional if I normalized). Shape: (n_chan,)
        Q_err (np.ndarray): Uncertainty on Q. Shape: (n_chan,)
        U_err (np.ndarray): Uncertainty on U. Shape: (n_chan,)
        I (np.ndarray): Optional Stokes I (for physical units). Shape: (n_chan,)
        I_err (np.ndarray): Optional uncertainty on I. Shape: (n_chan,)
    """
    lambda_sq: np.ndarray
    Q: np.ndarray
    U: np.ndarray
    Q_err: Optional[np.ndarray] = None
    U_err: Optional[np.ndarray] = None
    I: Optional[np.ndarray] = None
    I_err: Optional[np.ndarray] = None
    
    def __post_init__(self):
        """Validate shapes and set default errors."""
        n = len(self.lambda_sq)
        assert self.Q.shape == (n,), "Q must have same length as lambda_sq"
        assert self.U.shape == (n,), "U must have same length as lambda_sq"
        
        if self.Q_err is None:
            self.Q_err = np.zeros(n)
        if self.U_err is None:
            self.U_err = np.zeros(n)
        
        if self.I is not None:
            assert self.I.shape == (n,), "I must have same length as lambda_sq"
            if self.I_err is None:
                self.I_err = np.zeros(n)
    
    @property
    def q(self) -> np.ndarray:
        """Fractional Stokes q (assuming I=1 normalization or Q already fractional)."""
        return self.Q
    
    @property
    def u(self) -> np.ndarray:
        """Fractional Stokes u (assuming I=1 normalization or U already fractional)."""
        return self.U
    
    @property
    def p(self) -> np.ndarray:
        """Fractional polarization p = sqrt(q^2 + u^2)."""
        return np.sqrt(self.Q**2 + self.U**2)
    
    @property
    def psi(self) -> np.ndarray:
        """Polarization angle in radians, psi = 0.5 * arctan2(u, q)."""
        return 0.5 * np.arctan2(self.U, self.Q)


def power_law(lambda_sq: np.ndarray, I0: float, lambda_sq_0: float, alpha: float) -> np.ndarray:
    """
    Power-law model for Stokes I as function of λ².
    
    I(λ²) = I₀ (λ²/λ²₀)^α
    
    This corresponds to spectral index in frequency space:
    I(ν) ∝ ν^β where β = -2α
    
    Args:
        lambda_sq: Wavelength squared in m²
        I0: Normalization at λ²₀
        lambda_sq_0: Reference wavelength squared in m²
        alpha: Power-law index
        
    Returns:
        Stokes I values
    """
    return I0 * (lambda_sq / lambda_sq_0)**alpha


def polynomial(lambda_sq: np.ndarray, *coeffs) -> np.ndarray:
    """
    Polynomial model for Stokes I as function of λ².
    
    I(λ²) = c₀ + c₁λ² + c₂(λ²)² + ... + cₙ(λ²)ⁿ
    
    Args:
        lambda_sq: Wavelength squared in m²
        *coeffs: Polynomial coefficients [c₀, c₁, c₂, ..., cₙ]
        
    Returns:
        Stokes I values
    """
    return np.polyval(coeffs[::-1], lambda_sq)


class IFitter:
    """
    Fit Stokes I spectrum with smooth function.
    
    Supports power-law and polynomial models with weighted fitting.
    """
    
    def __init__(self, model: Literal['power-law', 'polynomial'] = 'power-law',
                 poly_degree: int = 3,
                 central_freq: Optional[float] = None):
        """
        Initialize I fitter.

        Args:
            model: 'power-law' or 'polynomial'
            poly_degree: Degree of polynomial (if model='polynomial')
            central_freq: Reference frequency in GHz for power-law normalization.
                         If None, λ²₀ is fixed at the median of the data (default).
                         If a float, λ²₀ = (c / central_freq_Hz)² is used instead.
        """
        self.model = model
        self.poly_degree = poly_degree
        self.central_freq = central_freq
        self.fit_params = None
        self.fit_cov = None
        self.I_fit = None
        self.I_fit_err = None
        self.lambda_sq = None  # Store for plotting
        
    def fit(self, lambda_sq: np.ndarray, I: np.ndarray, I_err: np.ndarray,
            compute_err: bool = True, alpha_guess: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
        """time domain
        Fit Stokes I and return fitted values with uncertainties.
        
        Args:
            lambda_sq: Wavelength squared in m²
            I: Stokes I data
            I_err: Uncertainties on I
            compute_err: If True, compute I_fit uncertainties (default: True)
            alpha_guess: Initial guess for spectral index α (power-law only)
                        If None, computed from log-log fit at band edges
                        Ignored for polynomial model
            
        Returns:
            Tuple of (I_fit, I_fit_err)
            If compute_err=False, I_fit_err is all zeros
        """
        # Weights for fitting (inverse variance)
        weights = 1.0 / I_err**2
        
        # =====================================================================
        # STEP 1: Fit model and get parameters + covariance
        # =====================================================================
        if self.model == 'power-law':
            # Fix reference wavelength at central_freq (if set) or data median
            if self.central_freq is not None:
                lambda_sq_0 = (C_LIGHT / (self.central_freq * 1e9))**2
            else:
                lambda_sq_0 = np.median(lambda_sq)
            
            # Initial guess for I0
            I0_guess = np.median(I)
            
            # Smart alpha guess from log-log fit at band edges (if not provided)
            if alpha_guess is None:
                # Use top and bottom 20% of wavelength range (band edges)
                n = len(lambda_sq)
                edge_fraction = 0.2
                n_edge = max(3, int(n * edge_fraction))  # At least 3 points
                
                # Bottom and top edges
                bottom_idx = np.arange(min(n_edge, n // 2))
                top_idx = np.arange(max(n // 2, n - n_edge), n)
                edge_idx = np.concatenate([bottom_idx, top_idx])
                
                # Fit line in log-log space: ln(I) = const + α · ln(λ²)
                # Avoid log(0) or log(negative)
                valid_mask = (lambda_sq[edge_idx] > 0) & (I[edge_idx] > 0)
                if np.sum(valid_mask) >= 3:
                    ln_lsq = np.log(lambda_sq[edge_idx][valid_mask])
                    ln_I = np.log(I[edge_idx][valid_mask])
                    # Weight by SNR
                    snr = I[edge_idx][valid_mask] / I_err[edge_idx][valid_mask]
                    snr_weight = np.clip(snr, 0.1, 100)  # Avoid extreme weights
                    
                    # Fit: slope is α
                    try:
                        poly_coeffs = np.polyfit(ln_lsq, ln_I, deg=1, w=snr_weight)
                        alpha_auto = poly_coeffs[0]
                    except:
                        # Fallback if polyfit fails
                        alpha_auto = -0.7  # Typical radio spectral index in λ² space
                else:
                    alpha_auto = -0.7  # Fallback
                
                alpha_init = alpha_auto
            else:
                alpha_init = alpha_guess
            
            # 2-parameter fit: [I0, alpha] (lsq0 is fixed)
            p0 = [I0_guess, alpha_init]
            
            # Define 2-parameter power-law with fixed reference
            def power_law_2param(lsq, I0, alpha):
                return I0 * (lsq / lambda_sq_0)**alpha
            
            # Fit power-law with bounds
            bounds = ([0, -5], [np.inf, 5])
            try:
                self.fit_params, self.fit_cov = curve_fit(
                    power_law_2param, lambda_sq, I, p0=p0, sigma=I_err, 
                    absolute_sigma=False, bounds=bounds
                )
            except RuntimeError:
                print("Warning: Power-law fit failed. Using polynomial instead.")
                self.model = 'polynomial'
                return self.fit(lambda_sq, I, I_err, compute_err=compute_err)
            
            # Store reference wavelength for later use
            self.lambda_sq_0 = lambda_sq_0
            
            # Compute fitted I
            I_fit = power_law_2param(lambda_sq, *self.fit_params)
            
        else:  # polynomial
            # Fit polynomial with optional covariance matrix
            if compute_err:
                self.fit_params, self.fit_cov = np.polyfit(
                    lambda_sq, I, self.poly_degree, w=np.sqrt(weights), cov=True
                )
            else:
                self.fit_params = np.polyfit(lambda_sq, I, self.poly_degree, w=np.sqrt(weights))
                self.fit_cov = None
            
            I_fit = np.polyval(self.fit_params, lambda_sq)
        
        # =====================================================================
        # STEP 2: Unified error propagation via Jacobian
        # =====================================================================
        if compute_err and self.fit_cov is not None:
            # Compute model-specific Jacobian
            if self.model == 'power-law':
                # Jacobian of power_law w.r.t. [I0, alpha]
                I0, alpha = self.fit_params
                ratio = lambda_sq / self.lambda_sq_0
                
                dI_dI0 = ratio**alpha
                dI_dalpha = I0 * ratio**alpha * np.log(ratio)
                
                jac = np.column_stack([dI_dI0, dI_dalpha])
                
            else:  # polynomial
                # Jacobian of polynomial w.r.t. [cₙ, cₙ₋₁, ..., c₁, c₀]
                # ∂I/∂cᵢ = (λ²)ⁱ (in descending order for np.polyval convention)
                jac_list = []
                for coeff_idx in range(self.poly_degree, -1, -1):
                    jac_list.append(lambda_sq**coeff_idx)
                
                jac = np.column_stack(jac_list)
            
            # Universal error propagation: Var(I) = J · Cov(θ) · J^T
            I_fit_var = np.sum(jac @ self.fit_cov * jac, axis=1)
            I_fit_err = np.sqrt(np.abs(I_fit_var))
        else:
            I_fit_err = np.zeros_like(I_fit)
        
        # Store fit results for later access (e.g., plotting)
        self.lambda_sq = np.asarray(lambda_sq)
        self.I_fit = I_fit
        self.I_fit_err = I_fit_err
        
        return I_fit, I_fit_err
    
    def __repr__(self):
        if self.fit_params is None:
            return f"IFitter(model='{self.model}', not fitted)"
        
        if self.model == 'power-law':
            I0, alpha = self.fit_params
            return f"IFitter(power-law: I₀={I0:.3e}, λ²₀={self.lambda_sq_0:.4f} (fixed), α_λ²={alpha:.3f}, α_ν={-2*alpha:.3f})"
        else:
            return f"IFitter(polynomial: degree={self.poly_degree})"


# =============================================================================
# Channel Flagging Functions
# =============================================================================

def _mad_sigma(x: np.ndarray) -> float:
    """
    Compute robust sigma estimate using Median Absolute Deviation (MAD).

    MAD is resistant to outliers unlike standard deviation.

    Parameters
    ----------
    x : np.ndarray
        Input array

    Returns
    -------
    float
        Robust sigma estimate (1.4826 * MAD)
    """
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return 0.0
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return 1.4826 * mad  # Convert MAD to sigma for normal distribution





def flag_channels(data_dict: Dict, config: 'FlaggingConfig') -> Dict:
    """
    Flag bad channels using sequential statistical criteria.

    Flagging order:
    0. Flag non-finite / invalid values (deterministic)
    1. Flag gross outliers in all Stokes I/Q/U/V flux and errors (deterministic)
    2. Clip extreme Stokes I values (deterministic)
    3. Flag frequency ranges (deterministic)
    4. Sigma-clip error outliers (statistical, iterative)
    5. Spectral baseline clipping on specified Stokes (statistical, iterative)
    6. Sigma-clip Stokes V outliers (significance-based, iterative)

    See FlaggingConfig for detailed parameter documentation.

    Parameters
    ----------
    data_dict : dict
        Must contain: 'freq', 'I', 'Q', 'U', 'dI', 'dQ', 'dU'
        Optional: 'V', 'dV'
    config : FlaggingConfig
        Flagging configuration (see FlaggingConfig docstring for all options)

    Returns
    -------
    flagging_info : dict
        {
            'mask': bool array (True = flagged),
            'kept_indices': int array (indices of kept channels),
            'flag_counts': dict (flagging statistics per criterion),
            'raw_data': dict (copy of input data),
            'n_total': int,
            'n_flagged': int,
            'n_kept': int,
            'flagged_regions': list or None (explicit frequency regions for plot highlighting)
        }
    """
    # Validate input
    required_keys = ['freq', 'I', 'Q', 'U', 'dI', 'dQ', 'dU']
    for key in required_keys:
        if key not in data_dict:
            raise KeyError(f"data_dict missing required key '{key}'")

    n = len(data_dict['freq'])

    # Basic length validation
    for k, v in data_dict.items():
        if len(v) != n:
            raise ValueError(f"data_dict key '{k}' has length {len(v)} != {n}")

    mask = np.zeros(n, dtype=bool)  # False = good, True = flagged
    flag_counts = {}

    if config.verbose:
        print(f"\nChannel Flagging:")
        print(f"  Total channels: {n}")

    # =========================================================================
    # Step 0: Flag non-finite / invalid values (deterministic)
    # =========================================================================
    flags_before = mask.sum()

    for k, v in data_dict.items():
        arr = np.asarray(v, dtype=float)

        # Flag non-finite values for any provided array
        nonfinite = ~np.isfinite(arr)
        if nonfinite.any():
            mask |= nonfinite

        # For uncertainty arrays, flag non-positive values too
        if k.startswith('d'):
            bad_err = np.isfinite(arr) & (arr <= 0)
            if bad_err.any():
                mask |= bad_err

    flags_total = mask.sum() - flags_before
    if flags_total > 0:
        flag_counts['Non-finite / invalid'] = {
            'count': flags_total,
            'notes': 'Flagged NaN/Inf values (and dX <= 0 for uncertainties)'
        }
        if config.verbose:
            print(f"  - Non-finite / invalid values: {flags_total} channels")
    else:
        if config.verbose:
            print(f"  - Non-finite / invalid values: 0 channels")

    # =========================================================================
    # Step 1: Gross outlier rejection across all Stokes flux and errors
    # =========================================================================
    if config.sigma_clip_extreme is not None:
        flags_before = mask.sum()
        threshold = float(config.sigma_clip_extreme)

        stokes_keys = [k for k in ['I', 'Q', 'U', 'V', 'dI', 'dQ', 'dU', 'dV']
                       if k in data_dict]

        for key in stokes_keys:
            arr = np.asarray(data_dict[key], dtype=float)
            valid = (~mask) & np.isfinite(arr)
            if valid.sum() < 5:
                continue

            med = np.median(arr[valid])
            sig = _mad_sigma(arr[valid])
            if sig == 0:
                sig = np.std(arr[valid])
            if sig == 0:
                continue

            outliers = np.isfinite(arr) & (~mask) & (np.abs(arr - med) > threshold * sig)
            mask |= outliers

        flags_total = mask.sum() - flags_before
        if flags_total > 0:
            flag_counts['Gross outliers'] = {
                'count': flags_total,
                'threshold': f'{threshold:.0f}σ (MAD)',
            }
            if config.verbose:
                print(f"  - Gross outliers ({threshold:.0f}σ MAD, all IQUV flux+errors): "
                      f"{flags_total} channels")
        else:
            if config.verbose:
                print(f"  - Gross outliers ({threshold:.0f}σ MAD): 0 channels")
    else:
        if config.verbose:
            print(f"  - Gross outliers: Not applied (None)")

    # =========================================================================
    # Step 2: Clip extreme Stokes I values (deterministic)
    # =========================================================================
    if config.clip_stokes_i is not None:
        flags_before = mask.sum()

        I = np.asarray(data_dict['I']).astype(float)
        unflagged_idx = np.where(~mask)[0]

        if unflagged_idx.size > 0:
            # Parse clip specification
            if isinstance(config.clip_stokes_i, int):
                clip_spec = [config.clip_stokes_i]
            else:
                clip_spec = list(config.clip_stokes_i)

            low_n = int(abs(sum(x for x in clip_spec if x < 0)))
            high_n = int(sum(x for x in clip_spec if x > 0))

            # Flag lowest values
            if low_n > 0:
                unflagged_idx_now = np.where(~mask)[0]
                if unflagged_idx_now.size > 0:
                    unflagged_I = I[unflagged_idx_now]
                    n_to_flag = min(low_n, unflagged_idx_now.size)
                    if n_to_flag > 0:
                        order = np.argsort(unflagged_I)
                        mask[unflagged_idx_now[order[:n_to_flag]]] = True

            # Flag highest values
            if high_n > 0:
                unflagged_idx_now = np.where(~mask)[0]
                if unflagged_idx_now.size > 0:
                    unflagged_I_now = I[unflagged_idx_now]
                    n_to_flag = min(high_n, unflagged_idx_now.size)
                    if n_to_flag > 0:
                        order = np.argsort(unflagged_I_now)
                        mask[unflagged_idx_now[order[-n_to_flag:]]] = True

        flags_total = mask.sum() - flags_before
        if flags_total > 0:
            flag_counts['Stokes I clipping'] = {
                'count': flags_total,
                'spec': config.clip_stokes_i
            }
            if config.verbose:
                print(f"  - Stokes I clipping (spec: {config.clip_stokes_i}): {flags_total} channels")
        else:
            if config.verbose:
                print(f"  - Stokes I clipping: 0 channels")
    else:
        if config.verbose:
            print(f"  - Stokes I clipping: Not applied (None)")

    # =========================================================================
    # Step 3: Flag specific frequency ranges (deterministic)
    # =========================================================================
    if config.freq_ranges is not None:
        flags_before = mask.sum()

        freq_hz = np.asarray(data_dict['freq']).astype(float)
        freq_ghz = freq_hz / 1e9

        # Handle single string or list of strings
        if isinstance(config.freq_ranges, str):
            ranges_list = [config.freq_ranges]
        else:
            ranges_list = list(config.freq_ranges)

        for freq_range in ranges_list:
            if '~' in freq_range:
                try:
                    freq_min_str, freq_max_str = freq_range.split('~')
                    freq_min = float(freq_min_str.strip())
                    freq_max = float(freq_max_str.strip())

                    range_flags = (freq_ghz >= freq_min) & (freq_ghz <= freq_max)
                    mask |= range_flags

                except ValueError:
                    if config.verbose:
                        print(f"  Warning: Could not parse frequency range '{freq_range}'")
                    continue
            else:
                if config.verbose:
                    print(f"  Warning: Frequency range '{freq_range}' missing '~' separator")

        flags_total = mask.sum() - flags_before
        if flags_total > 0:
            flag_counts['Frequency ranges'] = {
                'count': flags_total,
                'ranges': config.freq_ranges
            }
            if config.verbose:
                print(f"  - Frequency ranges ({config.freq_ranges}): {flags_total} channels")
        else:
            if config.verbose:
                print(f"  - Frequency ranges: 0 channels")
    else:
        if config.verbose:
            print(f"  - Frequency ranges: Not applied (None)")

    # =========================================================================
    # Step 4: Sigma-clip error outliers (statistical, iterative)
    # =========================================================================
    if config.sigma_clip_errors is not None:
        flags_before = mask.sum()

        iterations_used = 0
        errors_clip_space = str(getattr(config, 'sigma_clip_errors_space', 'linear')).lower()

        # Iterative clipping to refine statistics
        for iteration in range(int(getattr(config, 'max_iter', 0))):
            iterations_used = iteration + 1
            newly_flagged_this_iter = False

            for stokes in ['I', 'Q', 'U', 'V']:
                err_key = f'd{stokes}'
                if err_key not in data_dict:
                    continue

                errs = np.asarray(data_dict[err_key]).astype(float)

                # Only consider finite, unflagged, and positive errors (<=0 already flagged in Step 0)
                finite = np.isfinite(errs) & (~mask) & (errs > 0)
                if not finite.any():
                    continue

                if errors_clip_space in ('log', 'logspace', 'ln'):
                    # Log-space clipping (lognormal-like errors)
                    log_err = np.log(errs[finite])
                    med_log = np.median(log_err)
                    sigma_log = _mad_sigma(log_err)

                    if sigma_log == 0:
                        sigma_log = np.std(log_err) if finite.sum() > 1 else 0.0
                    if sigma_log == 0:
                        continue

                    with np.errstate(divide='ignore', invalid='ignore'):
                        err_flags = np.log(errs) > (med_log + config.sigma_clip_errors * sigma_log)

                else:
                    # Linear-space clipping (original behaviour)
                    med_err = np.median(errs[finite])
                    sigma_err_est = _mad_sigma(errs[finite])

                    if sigma_err_est == 0:
                        sigma_err_est = np.std(errs[finite]) if finite.sum() > 1 else 0.0
                    if sigma_err_est == 0:
                        continue

                    err_flags = errs > (med_err + config.sigma_clip_errors * sigma_err_est)

                newly_flagged = err_flags & (~mask) & np.isfinite(errs) & (errs > 0)

                if newly_flagged.any():
                    mask |= newly_flagged
                    newly_flagged_this_iter = True

            # Stop if no new flags in this iteration
            if not newly_flagged_this_iter:
                break

        flags_total = mask.sum() - flags_before
        if flags_total > 0:
            threshold_desc = (f'{config.sigma_clip_errors:.1f}σ above median in log-space'
                              if errors_clip_space in ('log', 'logspace', 'ln')
                              else f'{config.sigma_clip_errors:.1f}σ above median')
            flag_counts['Error outliers'] = {
                'count': flags_total,
                'threshold': threshold_desc,
                'space': errors_clip_space,
                'iterations': iterations_used
            }
            if config.verbose:
                print(f"  - Error outliers (threshold: {threshold_desc}, "
                      f"{iterations_used} iter): {flags_total} channels")
        else:
            if config.verbose:
                print(f"  - Error outliers: 0 channels")
    else:
        if config.verbose:
            print(f"  - Error outliers: Not applied (None)")

    # =========================================================================
    # Step 5: Spectral baseline clipping on specified Stokes (statistical, iterative)
    # =========================================================================
    sigma_clip_spectrum = getattr(config, 'sigma_clip_spectrum', None)

    if sigma_clip_spectrum is not None:
        flags_before = mask.sum()

        freq_hz = np.asarray(data_dict['freq'], dtype=float)
        spectrum_deg = int(getattr(config, 'spectrum_deg', 3))
        spectrum_jitter = float(getattr(config, 'spectrum_jitter', 0.0))
        spectrum_clip_sided = str(getattr(config, 'spectrum_clip_sided', 'two')).lower()
        
        # Parse which Stokes parameters to clip (case-insensitive)
        spectrum_stokes = getattr(config, 'spectrum_stokes', None)
        if spectrum_stokes is None:
            spectrum_stokes_list = ['I']  # Default to Stokes I only
        elif isinstance(spectrum_stokes, str):
            spectrum_stokes_list = [spectrum_stokes.upper()]
        else:
            spectrum_stokes_list = [s.upper() for s in spectrum_stokes]

        # Validate: QU cannot be combined with Q or U individually
        if 'QU' in spectrum_stokes_list and ('Q' in spectrum_stokes_list or 'U' in spectrum_stokes_list):
            raise ValueError(
                "spectrum_stokes: 'QU' (linear polarization amplitude) cannot be combined "
                "with 'Q' or 'U' individually — use one or the other."
            )

        do_qu_clip = 'QU' in spectrum_stokes_list
        regular_stokes = [s for s in spectrum_stokes_list if s != 'QU']

        # Ensure 'I' is always first in regular stokes (for cumulative flagging)
        if 'I' in regular_stokes:
            regular_stokes = ['I'] + [s for s in regular_stokes if s != 'I']

        whiten = bool(getattr(config, 'spectrum_whiten_residuals', True))
        iterations_used = 0
        final_fits = {}  # Store final polynomial fits for plotting (regular stokes only)
        total_flags_this_step = 0

        for iteration in range(int(getattr(config, 'max_iter', 0))):
            iterations_used = iteration + 1
            newly_flagged_this_iter = False

            # --- Regular per-Stokes polynomial baseline clipping ---
            for stokes in regular_stokes:
                if stokes not in data_dict or f'd{stokes}' not in data_dict:
                    if config.verbose and iteration == 0:
                        print(f"  Warning: Stokes {stokes} or d{stokes} not in data, skipping spectral clipping for {stokes}")
                    continue

                X = np.asarray(data_dict[stokes], dtype=float)
                dX = np.asarray(data_dict[f'd{stokes}'], dtype=float)

                fit_ok = (~mask) & np.isfinite(freq_hz) & np.isfinite(X) & np.isfinite(dX) & (dX > 0)

                if fit_ok.sum() < max(spectrum_deg + 2, 10):
                    continue

                current_jitter = spectrum_jitter if stokes == 'I' else 0.0

                try:
                    freq_norm = (freq_hz - np.mean(freq_hz[fit_ok])) / np.std(freq_hz[fit_ok])
                    weights = 1.0 / np.sqrt(dX[fit_ok]**2 + current_jitter**2)
                    coef = np.polyfit(freq_norm[fit_ok], X[fit_ok], deg=spectrum_deg, w=weights)
                    final_fits[stokes] = {
                        'coef': coef,
                        'freq_mean': np.mean(freq_hz[fit_ok]),
                        'freq_std': np.std(freq_hz[fit_ok])
                    }
                except Exception:
                    continue

                freq_norm_all = (freq_hz - final_fits[stokes]['freq_mean']) / final_fits[stokes]['freq_std']
                X_hat = np.polyval(coef, freq_norm_all)

                r = X - X_hat
                with np.errstate(divide='ignore', invalid='ignore'):
                    if whiten:
                        r_finite = r[fit_ok & np.isfinite(r)]
                        r_mad = _mad_sigma(r_finite) if r_finite.size > 1 else 0.0
                        denom = r_mad if r_mad > 0 else np.sqrt(dX**2 + current_jitter**2)
                    else:
                        denom = np.sqrt(dX**2 + current_jitter**2)
                    z = r / denom

                if spectrum_clip_sided == 'low':
                    new_flags = (z < -float(sigma_clip_spectrum))
                elif spectrum_clip_sided == 'high':
                    new_flags = (z > float(sigma_clip_spectrum))
                else:
                    new_flags = (np.abs(z) > float(sigma_clip_spectrum))

                valid_test = (~mask) & np.isfinite(z) & np.isfinite(denom) & (denom > 0)
                newly_flagged = new_flags & valid_test

                if newly_flagged.any():
                    mask[newly_flagged] = True
                    newly_flagged_this_iter = True

            # --- QU: polynomial baseline clipping on linear polarization amplitude P ---
            # Fits a polynomial to P = sqrt(Q²+U²) rather than Q and U individually,
            # which is more robust when high-RM wrapping makes per-parameter fitting unreliable.
            # The fit is not stored in final_fits (no QU panel in the diagnostic plot).
            if do_qu_clip and all(k in data_dict for k in ('Q', 'U', 'dQ', 'dU')):
                Q_arr = np.asarray(data_dict['Q'], dtype=float)
                U_arr = np.asarray(data_dict['U'], dtype=float)
                dQ_arr = np.asarray(data_dict['dQ'], dtype=float)
                dU_arr = np.asarray(data_dict['dU'], dtype=float)

                P = np.sqrt(Q_arr**2 + U_arr**2)
                with np.errstate(divide='ignore', invalid='ignore'):
                    dP = np.where(P > 0,
                                  np.sqrt((Q_arr * dQ_arr)**2 + (U_arr * dU_arr)**2) / P,
                                  np.inf)

                fit_ok = (~mask) & np.isfinite(freq_hz) & np.isfinite(P) & np.isfinite(dP) & (dP > 0)

                if fit_ok.sum() >= max(spectrum_deg + 2, 10):
                    try:
                        freq_norm = (freq_hz - np.mean(freq_hz[fit_ok])) / np.std(freq_hz[fit_ok])
                        weights = 1.0 / dP[fit_ok]
                        coef = np.polyfit(freq_norm[fit_ok], P[fit_ok], deg=spectrum_deg, w=weights)
                        freq_mean_qu = np.mean(freq_hz[fit_ok])
                        freq_std_qu = np.std(freq_hz[fit_ok])
                    except Exception:
                        coef = None

                    if coef is not None:
                        freq_norm_all = (freq_hz - freq_mean_qu) / freq_std_qu
                        P_hat = np.polyval(coef, freq_norm_all)

                        r = P - P_hat
                        with np.errstate(divide='ignore', invalid='ignore'):
                            if whiten:
                                r_finite = r[fit_ok & np.isfinite(r)]
                                r_mad = _mad_sigma(r_finite) if r_finite.size > 1 else 0.0
                                denom = r_mad if r_mad > 0 else dP
                            else:
                                denom = dP
                            z = r / denom

                        if spectrum_clip_sided == 'low':
                            new_flags = (z < -float(sigma_clip_spectrum))
                        elif spectrum_clip_sided == 'high':
                            new_flags = (z > float(sigma_clip_spectrum))
                        else:
                            new_flags = (np.abs(z) > float(sigma_clip_spectrum))

                        valid_test = (~mask) & np.isfinite(z) & np.isfinite(denom) & (denom > 0)
                        newly_flagged = new_flags & valid_test

                        if newly_flagged.any():
                            mask[newly_flagged] = True
                            newly_flagged_this_iter = True

            # Stop if no new flags in this iteration
            if not newly_flagged_this_iter:
                break

        flags_total = mask.sum() - flags_before
        if flags_total > 0:
            stokes_str = ', '.join(spectrum_stokes_list)
            qu_note = ' (QU: poly on amplitude)' if do_qu_clip else ''
            flag_counts['Spectral baseline clipping'] = {
                'count': flags_total,
                'threshold': f'{float(sigma_clip_spectrum):.1f}σ',
                'stokes': stokes_str,
                'model': f'poly(deg={spectrum_deg}){qu_note}',
                'jitter': spectrum_jitter,
                'sided': spectrum_clip_sided,
                'iterations': iterations_used,
                'fit': final_fits  # Polynomial fits for regular stokes only (no QU panel)
            }
            if config.verbose:
                whiten_str = 'whitened' if whiten else 'absolute'
                print(f"  - Spectral baseline clipping on {stokes_str} (threshold: {float(sigma_clip_spectrum):.1f}σ, "
                      f"residuals: {whiten_str}, deg: {spectrum_deg}, sided: {spectrum_clip_sided}, {iterations_used} iter): {flags_total} channels")
        else:
            if config.verbose:
                print(f"  - Spectral baseline clipping: 0 channels")
    else:
        if config.verbose:
            print(f"  - Spectral baseline clipping: Not applied (None)")

    # =========================================================================
    # Step 6: Sigma-clip Stokes V outliers (significance-based, iterative)
    # =========================================================================
    if config.sigma_clip_v is not None:
        if 'V' not in data_dict:
            if config.verbose:
                print(f"  - Stokes V clipping: Skipped (V not in data)")
        elif 'dV' not in data_dict:
            if config.verbose:
                print(f"  - Stokes V clipping: Skipped (dV errors not available)")
        else:
            flags_before = mask.sum()

            V = np.asarray(data_dict['V']).astype(float)
            dV = np.asarray(data_dict['dV']).astype(float)

            valid = (~mask) & np.isfinite(V) & np.isfinite(dV) & (dV > 0)

            if not valid.any():
                if config.verbose:
                    print(f"  - Stokes V clipping: 0 channels (no valid V/dV)")
            else:
                # Calculate threshold: auto or manual
                if config.sigma_clip_v == 'auto':
                    N_initial = int(valid.sum())
                    p_per_test = 1.0 / max(N_initial, 1)  # 1 expected false positive
                    threshold = stats.norm.ppf(1 - p_per_test / 2)

                    if config.verbose:
                        print(f"  - Stokes V auto threshold:")
                        print(f"    → N valid channels: {N_initial}")
                        print(f"    → Per-channel p-value: {p_per_test:.4f} (1/{N_initial})")
                        print(f"    → Two-tailed threshold: {threshold:.2f}σ")
                        print(f"    → Expected false positives: 1.0")
                else:
                    threshold = float(config.sigma_clip_v)

                iterations_used = 0

                for iteration in range(int(getattr(config, 'max_iter', 0))):
                    iterations_used = iteration + 1

                    good_idx = (~mask) & np.isfinite(V) & np.isfinite(dV) & (dV > 0)
                    if good_idx.sum() == 0:
                        break

                    V_good = V[good_idx]

                    # Choose reference (iteratively updated median)
                    if config.clip_v_against == 'zero':
                        ref_v = 0.0
                    else:  # 'median'
                        ref_v = np.median(V_good)

                    # Significance test: |(V - ref) / σ_V| > threshold
                    with np.errstate(divide='ignore', invalid='ignore'):
                        significance = np.abs(V - ref_v) / dV

                    new_flags = (significance > threshold)
                    newly_flagged = new_flags & good_idx & np.isfinite(significance)

                    if not newly_flagged.any():
                        break

                    mask[newly_flagged] = True

                flags_total = mask.sum() - flags_before
                if flags_total > 0:
                    flag_counts['Stokes V clipping'] = {
                        'count': flags_total,
                        'threshold': f'{threshold:.2f}σ' if config.sigma_clip_v == 'auto' else f'{threshold:.1f}σ',
                        'reference': config.clip_v_against,
                        'iterations': iterations_used
                    }
                    if config.verbose:
                        print(f"  - Stokes V clipping (threshold: {threshold:.2f}σ, "
                              f"ref: {config.clip_v_against}, {iterations_used} iter): {flags_total} channels")
                else:
                    if config.verbose:
                        print(f"  - Stokes V clipping: 0 channels")
    else:
        if config.verbose:
            print(f"  - Stokes V clipping: Not applied (None)")

    # =========================================================================
    # Compile results
    # =========================================================================
    kept_indices = np.where(~mask)[0]
    n_flagged = mask.sum()
    n_kept = len(kept_indices)

    if config.verbose:
        print(f"  Total flagged: {n_flagged}/{n} ({100*n_flagged/n:.1f}%)")
        print(f"  Total kept: {n_kept}/{n} ({100*n_kept/n:.1f}%)")

    # Build explicit flagged regions list for plot highlighting
    # Only include deterministic frequency ranges, not individual outliers
    flagged_regions = None
    if config.freq_ranges is not None:
        flagged_regions = []
        
        # Handle single string or list of strings
        if isinstance(config.freq_ranges, str):
            ranges_list = [config.freq_ranges]
        else:
            ranges_list = list(config.freq_ranges)
        
        for freq_range in ranges_list:
            if '~' in freq_range:
                try:
                    freq_min_str, freq_max_str = freq_range.split('~')
                    freq_min_ghz = float(freq_min_str.strip())
                    freq_max_ghz = float(freq_max_str.strip())
                    
                    # Convert to Hz for plot_flagging_diagnostic
                    freq_min_hz = freq_min_ghz * 1e9
                    freq_max_hz = freq_max_ghz * 1e9
                    
                    flagged_regions.append({
                        'freq_range': (freq_min_hz, freq_max_hz),
                        'reason': f'RFI band {freq_min_ghz}-{freq_max_ghz} GHz'
                    })
                except ValueError:
                    # Skip unparseable ranges
                    continue

    # Extract spectral fit if it was performed
    spectral_fit = None
    if 'Spectral baseline clipping' in flag_counts:
        spectral_fit = flag_counts['Spectral baseline clipping'].get('fit', None)
    
    flagging_info = {
        'mask': mask,
        'kept_indices': kept_indices,
        'flag_counts': flag_counts,
        'raw_data': {k: np.asarray(v).copy() for k, v in data_dict.items()},
        'n_total': n,
        'n_flagged': n_flagged,
        'n_kept': n_kept,
        'flagged_regions': flagged_regions,  # For plot_flagging_diagnostic highlighting
        'spectral_fit': spectral_fit  # Polynomial fit for plotting on diagnostic
    }

    return flagging_info


class DataLoader:
    """
    Load and prepare polarization data for fitting.
    
    Supports three loading modes:
    A) Physical units → fit I → fractional
    B) Physical units → direct fractional (I provided)
    C) Already fractional (I=1 assumed)
    """
    
    @staticmethod
    def load_physical_with_fit(
        data_source,
        I_model: Literal['power-law', 'polynomial'] = 'power-law',
        poly_degree: int = 3,
        central_freq: Optional[float] = None,
        columns: Optional[Tuple[int, ...]] = None,
        flagging: Optional[FlaggingConfig] = None,
        load_stokes_v: bool = False
    ) -> Union[Tuple[PolarizationData, IFitter], 
               Tuple[PolarizationData, IFitter, Dict]]:
        """
        Load data in physical units and fit Stokes I with optional flagging.
        
        Mode A: Load I,Q,U,dI,dQ,dU → (optionally flag) → fit I(λ²) → compute fractional q,u
        
        Args:
            data_source: Either:
                - str/Path: Path to data file to load
                - tuple: Direct data as (freq_hz, I, Q, U, dI, dQ, dU) or
                         (lambda_sq, I, Q, U, dI, dQ, dU) - auto-detected
            I_model: 'power-law' or 'polynomial'
            poly_degree: Polynomial degree (if using polynomial)
            columns: Optional tuple of column indices (only used if loading from file)
                    Default: (0, 1, 2, 3, 4, 5, 6) for I,Q,U,dI,dQ,dU
                    With V: (0, 1, 2, 3, 4, 5, 6, 7, 8, 9) for I,Q,U,V,dI,dQ,dU,dV
            flagging: Optional FlaggingConfig for channel flagging
            load_stokes_v: If True, load Stokes V from file (column 4 if V present)
                    
        Returns:
            Without flagging: (PolarizationData, IFitter)
            With flagging: (PolarizationData, IFitter, flagging_info)
            
            - PolarizationData contains fractional q,u with propagated errors
            - IFitter object with fit parameters
            - flagging_info (if flagging applied): dict with mask, raw_data, etc.
        """
        # Determine if data_source is file or direct data
        if isinstance(data_source, (str, Path)):
            # Load from file
            if columns is None:
                if load_stokes_v:
                    columns = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)  # freq, I, Q, U, V, dI, dQ, dU, dV
                else:
                    columns = (0, 1, 2, 3, 4, 5, 6)  # freq, I, Q, U, dI, dQ, dU
            
            data = np.loadtxt(data_source)
            freq_or_lambda = data[:, columns[0]]
            I = data[:, columns[1]]
            Q = data[:, columns[2]]
            U = data[:, columns[3]]
            
            if load_stokes_v and len(columns) >= 9:
                V = data[:, columns[4]]
                dI = data[:, columns[5]]
                dQ = data[:, columns[6]]
                dU = data[:, columns[7]]
                dV = data[:, columns[8]] if len(columns) >= 9 else None
            else:
                V = None
                dV = None
                dI = data[:, columns[4]]
                dQ = data[:, columns[5]]
                dU = data[:, columns[6]]
                
        elif isinstance(data_source, tuple):
            # Direct data arrays
            if len(data_source) == 7:
                freq_or_lambda, I, Q, U, dI, dQ, dU = data_source
                V = None
                dV = None
            elif len(data_source) == 9:
                freq_or_lambda, I, Q, U, V, dI, dQ, dU, dV = data_source
            else:
                raise ValueError("Direct data must be tuple of 7 or 9 arrays")
        else:
            raise TypeError("data_source must be filename (str/Path) or tuple of arrays")
        
        # Auto-detect if freq_or_lambda is frequency (Hz) or lambda_sq (m²)
        # Typical ranges: freq ~1e8-1e10 Hz, lambda_sq ~0.01-0.1 m²
        if np.median(freq_or_lambda) > 1e6:
            # Looks like frequency in Hz
            freq_hz = freq_or_lambda
            lambda_sq = (C_LIGHT / freq_or_lambda)**2
        else:
            # Already lambda_sq in m²
            lambda_sq = freq_or_lambda
            freq_hz = C_LIGHT / np.sqrt(lambda_sq)
        
        # =================================================================
        # Apply flagging if requested
        # =================================================================
        flagging_info = None
        if flagging is not None:
            # Build data dict for flagging
            data_dict = {
                'freq': freq_hz,
                'I': I,
                'Q': Q,
                'U': U,
                'dI': dI,
                'dQ': dQ,
                'dU': dU
            }
            
            if V is not None:
                data_dict['V'] = V
                if dV is not None:
                    data_dict['dV'] = dV
            
            # Apply flagging
            flagging_info = flag_channels(data_dict, flagging)
            
            # Extract kept channels
            kept_idx = flagging_info['kept_indices']
            freq_hz = freq_hz[kept_idx]
            lambda_sq = lambda_sq[kept_idx]
            I = I[kept_idx]
            Q = Q[kept_idx]
            U = U[kept_idx]
            dI = dI[kept_idx]
            dQ = dQ[kept_idx]
            dU = dU[kept_idx]
            
            if V is not None:
                V = V[kept_idx]
                if dV is not None:
                    dV = dV[kept_idx]
        
        # Fit Stokes I (on flagged data if flagging applied)
        fitter = IFitter(model=I_model, poly_degree=poly_degree, central_freq=central_freq)
        I_fit, I_fit_err = fitter.fit(lambda_sq, I, dI)
        
        # Convert to fractional with full error propagation
        q = Q / I_fit
        u = U / I_fit
        
        # Error propagation: σ²_q = (σ_Q/I)² + (Q·σ_I/I²)²
        dq = np.sqrt((dQ / I_fit)**2 + (Q * I_fit_err / I_fit**2)**2)
        du = np.sqrt((dU / I_fit)**2 + (U * I_fit_err / I_fit**2)**2)
        
        pol_data = PolarizationData(
            lambda_sq=lambda_sq,
            Q=q,
            U=u,
            Q_err=dq,
            U_err=du,
            I=I,
            I_err=dI
        )
        
        # Print summary
        if isinstance(data_source, (str, Path)):
            print(f"\nLoaded data from {data_source}")
        else:
            print(f"\nLoaded data from direct arrays")
        print(f"  Mode: Physical → Fitted I → Fractional")
        print(f"  I model: {fitter}")
        print(f"  Channels: {len(lambda_sq)}")
        print(f"  λ² range: [{lambda_sq.min():.6f}, {lambda_sq.max():.6f}] m²")
        
        if flagging_info is not None:
            return pol_data, fitter, flagging_info
        else:
            return pol_data, fitter
    
    @staticmethod
    def load_physical_direct(
        filename: str,
        columns: Optional[Tuple[int, ...]] = None,
        flagging: Optional[FlaggingConfig] = None,
        load_stokes_v: bool = False
    ) -> Union[PolarizationData, Tuple[PolarizationData, Dict]]:
        """
        Load data in physical units and compute fractional directly with optional flagging.
        
        Mode B: Load I,Q,U,dI,dQ,dU → (optionally flag) → compute q=Q/I, u=U/I directly
        
        Args:
            filename: Path to data file
            columns: Optional tuple of column indices
                    Default: (0, 1, 2, 3, 4, 5, 6) for I,Q,U,dI,dQ,dU
                    With V: (0, 1, 2, 3, 4, 5, 6, 7, 8, 9) for I,Q,U,V,dI,dQ,dU,dV
            flagging: Optional FlaggingConfig for channel flagging
            load_stokes_v: If True, load Stokes V from file
                    
        Returns:
            Without flagging: PolarizationData
            With flagging: (PolarizationData, flagging_info)
        """
        if columns is None:
            if load_stokes_v:
                columns = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)  # freq, I, Q, U, V, dI, dQ, dU, dV
            else:
                columns = (0, 1, 2, 3, 4, 5, 6)  # freq, I, Q, U, dI, dQ, dU
        
        # Load data
        data = np.loadtxt(filename)
        freq_hz = data[:, columns[0]]
        I = data[:, columns[1]]
        Q = data[:, columns[2]]
        U = data[:, columns[3]]
        
        if load_stokes_v and len(columns) >= 9:
            V = data[:, columns[4]]
            dI = data[:, columns[5]]
            dQ = data[:, columns[6]]
            dU = data[:, columns[7]]
            dV = data[:, columns[8]] if len(columns) >= 9 else None
        else:
            V = None
            dV = None
            dI = data[:, columns[4]]
            dQ = data[:, columns[5]]
            dU = data[:, columns[6]]
        
        # Convert frequency to lambda_sq
        lambda_sq = (C_LIGHT / freq_hz)**2

        # =================================================================
        # Apply flagging if requested
        # =================================================================
        flagging_info = None
        if flagging is not None:
            # Build data dict for flagging
            data_dict = {
                'freq': freq_hz,
                'I': I,
                'Q': Q,
                'U': U,
                'dI': dI,
                'dQ': dQ,
                'dU': dU
            }
            
            if V is not None:
                data_dict['V'] = V
                if dV is not None:
                    data_dict['dV'] = dV
            
            # Apply flagging
            flagging_info = flag_channels(data_dict, flagging)
            
            # Extract kept channels
            kept_idx = flagging_info['kept_indices']
            freq_hz = freq_hz[kept_idx]
            lambda_sq = lambda_sq[kept_idx]
            I = I[kept_idx]
            Q = Q[kept_idx]
            U = U[kept_idx]
            dI = dI[kept_idx]
            dQ = dQ[kept_idx]
            dU = dU[kept_idx]
            
            if V is not None:
                V = V[kept_idx]
                if dV is not None:
                    dV = dV[kept_idx]
        
        # Convert to fractional with error propagation
        q = Q / I
        u = U / I
        
        # Error propagation: σ²_q = (σ_Q/I)² + (Q·σ_I/I²)²
        dq = np.sqrt((dQ / I)**2 + (Q * dI / I**2)**2)
        du = np.sqrt((dU / I)**2 + (U * dI / I**2)**2)
        
        pol_data = PolarizationData(
            lambda_sq=lambda_sq,
            Q=q,
            U=u,
            Q_err=dq,
            U_err=du,
            I=I,
            I_err=dI
        )
        
        print(f"\nLoaded data from {filename}")
        print(f"  Mode: Physical → Direct Fractional")
        print(f"  Channels: {len(lambda_sq)}")
        print(f"  λ² range: [{lambda_sq.min():.6f}, {lambda_sq.max():.6f}] m²")
        
        if flagging_info is not None:
            return pol_data, flagging_info
        else:
            return pol_data
    
    @staticmethod
    def load_fractional(
        filename: str,
        columns: Optional[Tuple[int, ...]] = None,
        flagging: Optional[FlaggingConfig] = None,
        load_stokes_v: bool = False
    ) -> Union[PolarizationData, Tuple[PolarizationData, Dict]]:
        """
        Load data already in fractional form (I=1 assumed) with optional flagging.
        
        Mode C: Load q,u,dq,du → (optionally flag) → assume I=1
        
        Args:
            filename: Path to data file
            columns: Optional tuple of column indices
                    Default: (0, 1, 2, 3, 4) for q,u,dq,du
                    With V: (0, 1, 2, 3, 4, 5, 6) for q,u,v,dq,du,dv
            flagging: Optional FlaggingConfig for channel flagging
            load_stokes_v: If True, load Stokes V from file
                    
        Returns:
            Without flagging: PolarizationData
            With flagging: (PolarizationData, flagging_info)
        """
        if columns is None:
            if load_stokes_v:
                columns = (0, 1, 2, 3, 4, 5, 6)  # freq, q, u, v, dq, du, dv
            else:
                columns = (0, 1, 2, 3, 4)  # freq, q, u, dq, du
        
        # Load data
        data = np.loadtxt(filename)
        freq_hz = data[:, columns[0]]
        q = data[:, columns[1]]
        u = data[:, columns[2]]
        
        if load_stokes_v and len(columns) >= 7:
            v = data[:, columns[3]]
            dq = data[:, columns[4]]
            du = data[:, columns[5]]
            dv = data[:, columns[6]]
        else:
            v = None
            dv = None
            dq = data[:, columns[3]]
            du = data[:, columns[4]]
        
        # Convert frequency to lambda_sq
        lambda_sq = (C_LIGHT / freq_hz)**2

        # =================================================================
        # Apply flagging if requested
        # =================================================================
        flagging_info = None
        if flagging is not None:
            # For fractional data, I=1 so we construct I and error arrays
            I = np.ones_like(q)
            dI = np.zeros_like(q)  # No I error for fractional data (I=1 exactly)
            
            # Convert fractional Q,U to absolute for flagging consistency
            Q = q * I  # = q since I=1
            U = u * I  # = u since I=1
            dQ = dq * I  # = dq since I=1
            dU = du * I  # = du since I=1
            
            # Build data dict for flagging
            data_dict = {
                'freq': freq_hz,
                'I': I,
                'Q': Q,
                'U': U,
                'dI': dI,
                'dQ': dQ,
                'dU': dU
            }
            
            if v is not None:
                V = v * I  # = v since I=1
                data_dict['V'] = V
                if dv is not None:
                    dV = dv * I  # = dv since I=1
                    data_dict['dV'] = dV
            
            # Apply flagging
            flagging_info = flag_channels(data_dict, flagging)
            
            # Extract kept channels
            kept_idx = flagging_info['kept_indices']
            freq_hz = freq_hz[kept_idx]
            lambda_sq = lambda_sq[kept_idx]
            q = q[kept_idx]
            u = u[kept_idx]
            dq = dq[kept_idx]
            du = du[kept_idx]
            
            if v is not None:
                v = v[kept_idx]
                if dv is not None:
                    dv = dv[kept_idx]
        
        pol_data = PolarizationData(
            lambda_sq=lambda_sq,
            Q=q,
            U=u,
            Q_err=dq,
            U_err=du
        )
        
        print(f"\nLoaded data from {filename}")
        print(f"  Mode: Already Fractional (I=1)")
        print(f"  Channels: {len(lambda_sq)}")
        print(f"  λ² range: [{lambda_sq.min():.6f}, {lambda_sq.max():.6f}] m²")
        
        if flagging_info is not None:
            return pol_data, flagging_info
        else:
            return pol_data


def save_polarization_data(data: PolarizationData, filename: str, 
                           output_dir: Path = None, include_I: bool = True):
    """
    Save polarization data to ASCII file.
    
    Format depends on include_I flag:
    - If True: freq(Hz), I, Q, U, dI, dQ, dU (physical or fractional)
    - If False: freq(Hz), Q, U, dQ, dU (fractional only)
    
    Args:
        data: PolarizationData object
        filename: Name of output file (with .txt extension)
        output_dir: Directory to save file (default: current directory)
        include_I: If True and I is available, save I column
    """
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert lambda_sq to frequency in Hz
    freq_hz = C_LIGHT / np.sqrt(data.lambda_sq)
    
    if include_I and data.I is not None:
        # Save with I column
        output = np.column_stack([
            freq_hz,
            data.I,
            data.Q,
            data.U,
            data.I_err,
            data.Q_err,
            data.U_err
        ])
    else:
        # Save without I column (fractional only)
        output = np.column_stack([
            freq_hz,
            data.Q,
            data.U,
            data.Q_err,
            data.U_err
        ])
    
    filepath = output_dir / filename
    np.savetxt(filepath, output, fmt='%.10e')
    
    print(f"Saved data: {filepath}")


def save_flagged_data(flagging_info: Dict, filename: str,
                      output_dir: Path = None):
    """
    Save flagged data to ASCII file with flagging information in header.
    
    This creates a new data file containing only the kept (unflagged) channels,
    with comprehensive flagging statistics in commented header lines. The output
    file can be loaded directly in subsequent runs to skip re-flagging.
    
    Args:
        flagging_info: Dictionary from flag_channels() containing:
            - 'raw_data': Original data before flagging
            - 'kept_indices': Indices of kept channels
            - 'flag_counts': Statistics per flagging criterion
            - 'n_total', 'n_flagged', 'n_kept'
        filename: Name of output file (with .txt extension)
        output_dir: Directory to save file (default: current directory)
    
    Output file format:
        Header lines (all start with #):
            - Flagging summary statistics
            - Flagging criteria applied
            - Column descriptions
        Data columns:
            - Depends on what was in raw_data
            - Typical: freq(Hz), I, Q, U, dI, dQ, dU
            - With V: freq(Hz), I, Q, U, V, dI, dQ, dU, dV
    """
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract data
    raw_data = flagging_info['raw_data']
    kept_idx = flagging_info['kept_indices']
    n_total = flagging_info['n_total']
    n_flagged = flagging_info['n_flagged']
    n_kept = flagging_info['n_kept']
    flag_counts = flagging_info['flag_counts']
    
    # Build header with flagging information
    header_lines = []
    header_lines.append("="*70)
    header_lines.append("FLAGGED DATA - KEPT CHANNELS ONLY")
    header_lines.append("="*70)
    header_lines.append("")
    header_lines.append("FLAGGING SUMMARY:")
    header_lines.append(f"  Total channels: {n_total}")
    header_lines.append(f"  Flagged: {n_flagged} ({100*n_flagged/n_total:.1f}%)")
    header_lines.append(f"  Kept: {n_kept} ({100*n_kept/n_total:.1f}%)")
    header_lines.append("")
    
    if flag_counts:
        header_lines.append("FLAGGING CRITERIA APPLIED:")
        for criterion, info in flag_counts.items():
            count = info['count']
            header_lines.append(f"  - {criterion}: {count} channels")
            
            # Add detailed info for each criterion
            if criterion == 'Stokes I clipping':
                header_lines.append(f"      Specification: {info['spec']}")
            elif criterion == 'Frequency ranges':
                header_lines.append(f"      Ranges: {info['ranges']}")
            elif criterion == 'Error outliers':
                header_lines.append(f"      Threshold: {info['threshold']}")
                header_lines.append(f"      Iterations: {info['iterations']}")
            elif criterion == 'Stokes V clipping':
                header_lines.append(f"      Threshold: {info['threshold']}")
                header_lines.append(f"      Reference: {info['reference']}")
                header_lines.append(f"      Iterations: {info['iterations']}")
        header_lines.append("")
    else:
        header_lines.append("FLAGGING CRITERIA APPLIED: None")
        header_lines.append("")
    
    # Determine column format from raw_data keys
    has_V = 'V' in raw_data
    if has_V:
        header_lines.append("COLUMN FORMAT: freq(Hz) I Q U V dI dQ dU dV")
        columns_to_save = ['freq', 'I', 'Q', 'U', 'V', 'dI', 'dQ', 'dU', 'dV']
    else:
        header_lines.append("COLUMN FORMAT: freq(Hz) I Q U dI dQ dU")
        columns_to_save = ['freq', 'I', 'Q', 'U', 'dI', 'dQ', 'dU']
    
    header_lines.append("="*70)
    
    # Build header string
    header = "\n".join(header_lines)
    
    # Extract kept channels from raw data
    data_arrays = []
    for col in columns_to_save:
        if col in raw_data:
            data_arrays.append(np.asarray(raw_data[col])[kept_idx])
        else:
            # Handle missing columns gracefully
            if col == 'V' or col == 'dV':
                continue  # Skip optional V columns if not present
            else:
                raise KeyError(f"Required column '{col}' not in raw_data")
    
    # Stack into output array
    output = np.column_stack(data_arrays)
    
    # Save to file
    filepath = output_dir / filename
    np.savetxt(filepath, output, header=header, fmt='%.10e')
    
    print(f"Saved flagged data: {filepath}")
    print(f"  Kept {n_kept}/{n_total} channels ({100*n_kept/n_total:.1f}%)")


def validate_data(data: PolarizationData) -> bool:
    """
    Validate polarization data for fitting.
    
    Checks:
    - No NaN or Inf values
    - Positive errors
    - Reasonable fractional polarization (|p| <= 1)
    
    Args:
        data: PolarizationData to validate
        
    Returns:
        True if valid, raises ValueError if invalid
    """
    # Check for NaN/Inf
    for name, arr in [('lambda_sq', data.lambda_sq), ('Q', data.Q), ('U', data.U),
                      ('Q_err', data.Q_err), ('U_err', data.U_err)]:
        if np.any(~np.isfinite(arr)):
            raise ValueError(f"{name} contains NaN or Inf values")
    
    # Check positive errors
    if np.any(data.Q_err <= 0):
        raise ValueError("Q_err must be positive")
    if np.any(data.U_err <= 0):
        raise ValueError("U_err must be positive")
    
    # Check reasonable polarization fraction
    p = data.p
    if np.any(p > 1.5):  # Allow some margin for noise
        n_bad = np.sum(p > 1.5)
        print(f"Warning: {n_bad} channels have p > 1.5 (possibly noisy)")
    
    return True
