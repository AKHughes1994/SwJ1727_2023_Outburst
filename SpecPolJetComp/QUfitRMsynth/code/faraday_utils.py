# faraday_utils.py
"""
Utilities for Faraday rotation modeling using phenomenological components.

This module implements the super-Gaussian approach from Anderson et al. (2016)
for modeling broadband polarization spectra.

Physical conventions:
- Complex polarization: P = Q + iU = p*I*exp(2i*psi)
- Faraday depth φ in rad m^-2
- Wavelength squared λ² in m²
- Polarization angle ψ in radians
"""

from __future__ import annotations

import numpy as np
from typing import List, Tuple, Optional, TYPE_CHECKING
from scipy import signal
from scipy.special import gamma

C_LIGHT = 299792458.0  # m/s

try:
    import numba
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    import warnings
    warnings.warn("numba not available - computation will be slower", RuntimeWarning)

if TYPE_CHECKING:
    from faraday_data import PolarizationData


def cartesian_to_polar(a: float, b: float) -> Tuple[float, float]:
    """
    Convert Cartesian (a, b) to polar (p0, psi_0).
    
    Args:
        a: Q component of intrinsic polarization
        b: U component of intrinsic polarization
    
    Returns:
        p0: Total polarization amplitude
        psi_0: Intrinsic polarization angle in radians
    """
    p0 = np.sqrt(a**2 + b**2)
    psi_0 = 0.5 * np.arctan2(b, a)
    return p0, psi_0


def polar_to_cartesian(p0: float, psi_0: float) -> Tuple[float, float]:
    """
    Convert polar (p0, psi_0) to Cartesian (a, b).
    
    Args:
        p0: Total polarization amplitude
        psi_0: Intrinsic polarization angle in radians
    
    Returns:
        a: Q component of intrinsic polarization
        b: U component of intrinsic polarization
    """
    a = p0 * np.cos(2 * psi_0)
    b = p0 * np.sin(2 * psi_0)
    return a, b


class ThinComponent:
    """
    Faraday-thin polarized emission component (delta function in Faraday space).
    
    Represents a single Faraday depth with no internal depolarization.
    Models idealized Faraday rotation: P(λ²) = (a + ib) * exp(i*2*φ_RM*λ²)
    
    Parameters:
        a (float): Q component of intrinsic polarization
        b (float): U component of intrinsic polarization
        phi_rm (float): Rotation measure in rad m^-2
        name (str): Optional name for tracking (e.g., 'ISM', 'Internal')
    """
    
    def __init__(self, a: float, b: float, phi_rm: float, name: Optional[str] = None):
        self.a = a
        self.b = b
        self.phi_rm = phi_rm
        self.name = name
    
    def compute_fdf(self, phi_grid: np.ndarray) -> np.ndarray:
        """
        Compute Faraday dispersion function on a Faraday depth grid.
        
        For a thin component, F(φ) is a delta function at φ = phi_rm.
        On a discrete grid, this is a single non-zero point at the nearest
        grid location.
        
        Args:
            phi_grid (np.ndarray): Faraday depth grid in rad m^-2
            
        Returns:
            np.ndarray: Complex FDF values
        """
        # Delta function on discrete grid
        fdf = np.zeros(len(phi_grid), dtype=complex)
        
        # Find nearest grid point to phi_rm
        idx = np.argmin(np.abs(phi_grid - self.phi_rm))
        
        # Place amplitude at nearest grid point (no scaling)
        fdf[idx] = self.a + 1j * self.b
        
        return fdf
    
    def compute_polarization(self, lambda_sq: np.ndarray) -> np.ndarray:
        """
        Compute complex polarization P(λ²) = Q + iU.
        
        For a thin component: P(λ²) = (a + ib) * exp(i*2*φ_rm*λ²)
        
        Args:
            lambda_sq (np.ndarray): Wavelength squared in m²
            
        Returns:
            np.ndarray: Complex polarization P(λ²)
        """
        return (self.a + 1j * self.b) * np.exp(2j * self.phi_rm * lambda_sq)

    def compute_transfer_factor(self, lambda_sq_ref: float) -> complex:
        """
        Compute the transfer factor K at a reference wavelength squared.

        K is the model response at λ²_ref for unit intrinsic amplitude, so that
        the intrinsic prefactor A_int = A_ref / K, where A_ref = q_ref + i*u_ref
        is the complex amplitude sampled at the reference wavelength.

        For a thin component K is a pure phase rotation with |K| = 1, so p₀ is
        unchanged by the conversion and only ψ₀ is shifted.

        Args:
            lambda_sq_ref (float): Reference wavelength squared in m²

        Returns:
            complex: Transfer factor K(λ²_ref)
        """
        return np.exp(2j * self.phi_rm * lambda_sq_ref)

    def __repr__(self):
        name_str = f"'{self.name}'" if self.name else "unnamed"
        p0, psi_0 = cartesian_to_polar(self.a, self.b)
        return f"ThinComponent({name_str}: a={self.a:.4f}, b={self.b:.4f}, RM={self.phi_rm:.1f} rad/m² | p₀={p0:.4f}, ψ₀={psi_0:.3f} rad)"


class ThinPowerLawComponent:
    """
    Faraday-thin polarized emission with power-law fractional polarization spectrum.
    
    Represents a single Faraday depth with power-law frequency-dependent fractional polarization.
    
    In frequency space:   P(ν) = (a + ib) * (ν/ν₀)^β * exp(i*2*φ_RM*(c/ν)²)
    In λ² space:          P(λ²) = (a + ib) * (λ²/λ²₀)^(-β/2) * exp(i*2*φ_RM*λ²)
    
    where β is the spectral index of fractional polarization: p(ν) ∝ ν^β
    
    This is a phenomenological model useful for fitting sources with frequency-dependent
    fractional polarization across the observing band. The FDF is computed numerically
    since no analytical form exists.
    
    Parameters:
        a (float): Q component of intrinsic polarization at reference λ²₀ (or ν₀)
        b (float): U component of intrinsic polarization at reference λ²₀ (or ν₀)
        phi_rm (float): Rotation measure in rad m^-2
        beta (float): Spectral index of fractional polarization (p(ν) ∝ ν^β)
                     Typical range: -2 to +2 (radio sources typically β ~ -0.5 to -1.0)
        lambda_sq_0 (float, optional): Reference wavelength squared in m²
        name (str): Optional name for tracking (e.g., 'Core', 'Jet')
    
    Notes:
        IMPORTANT: Unlike ThinComponent, the parameters (a, b) and derived p0
        represent fractional polarization at the REFERENCE WAVELENGTH λ²₀ (or ν₀), NOT at λ² = 0.
        This is required because the power-law spectrum cannot be evaluated at λ² = 0:
        - If β < 0: P(0) = (a+ib) * 0^(-β/2) → ∞
        - If β > 0: P(0) = (a+ib) * 0^(β/2) → 0
        
        At reference: P(λ²₀) = (a + ib) * exp(2i*φ_rm*λ²₀)
                      p0 = sqrt(a² + b²) is fractional polarization at λ²₀ (or ν₀)
        
        If lambda_sq_0 is not provided during initialization, it must be set
        before calling compute_polarization() or compute_fdf(). Typically
        set to the variance-weighted mean of the data: λ²₀ = Σ(w_i*λ²_i)/Σ(w_i)
        where w_i are inverse variance weights.
    """
    
    def __init__(self, a: float, b: float, phi_rm: float, beta: float,
                 lambda_sq_0: Optional[float] = None, name: Optional[str] = None):
        self.a = a
        self.b = b
        self.phi_rm = phi_rm
        self.beta = beta
        self.lambda_sq_0 = lambda_sq_0
        self.name = name
    
    def compute_fdf(self, phi_grid: np.ndarray,
                   lambda_sq_min: Optional[float] = None,
                   lambda_sq_max: Optional[float] = None,
                   n_lambda_sq: int = 2000) -> np.ndarray:
        """
        Compute Faraday dispersion function on a Faraday depth grid.
        
        For a power-law component, the FDF must be computed numerically via
        inverse Fourier transform:
        
        F(φ) = ∫ P(λ²) exp(-i*2*φ*λ²) dλ²
        
        where P(λ²) = (a + ib) * (λ²/λ²₀)^(-β/2) * exp(i*2*φ_rm*λ²)
        
        This is computed using trapezoidal integration over a dense λ² grid.
        The integration range defaults to 0.1 GHz - 50 GHz if not specified.
        
        Args:
            phi_grid (np.ndarray): Faraday depth grid in rad m^-2
            lambda_sq_min (float, optional): Minimum λ² for integration (m²)
            lambda_sq_max (float, optional): Maximum λ² for integration (m²)
            n_lambda_sq (int): Number of λ² grid points for integration (default: 2000)
            
        Returns:
            np.ndarray: Complex FDF values
            
        Raises:
            ValueError: If lambda_sq_0 has not been set
        """
        if self.lambda_sq_0 is None:
            raise ValueError(
                "lambda_sq_0 must be set before computing FDF. "
                "Set via component.lambda_sq_0 = value or provide during initialization."
            )
        
        # Set default integration range if not provided
        # From ~0.1 GHz to ~3 GHz (covers most radio observations)

        if lambda_sq_min is None:
            freq_max = 3e9  # Hz
            lambda_sq_min = (C_LIGHT / freq_max)**2
        if lambda_sq_max is None:
            freq_min = 0.1e9  # Hz
            lambda_sq_max = (C_LIGHT / freq_min)**2
        
        # Create dense λ² grid for numerical integration
        lambda_sq_grid = np.linspace(lambda_sq_min, lambda_sq_max, n_lambda_sq)
        
        # Compute polarization on this grid
        P_lsq = self.compute_polarization(lambda_sq_grid)
        
        # Inverse Fourier transform to get FDF
        # F(φ) = ∫ P(λ²) exp(-i*2*φ*λ²) dλ²
        # Use numba-optimized version for ~100x speedup with parallelization
        fdf = inverse_transform_fdf_numba(P_lsq, lambda_sq_grid, phi_grid)
        
        return fdf
    
    def compute_polarization(self, lambda_sq: np.ndarray) -> np.ndarray:
        """
        Compute complex polarization P(λ²) = Q + iU.
        
        In frequency space:   P(ν) = (a + ib) * (ν/ν₀)^β * exp(i*2*φ_rm*(c/ν)²)
        In λ² space:          P(λ²) = (a + ib) * (λ²/λ²₀)^(-β/2) * exp(i*2*φ_rm*λ²)
        
        where β is the spectral index of fractional polarization: p(ν) ∝ ν^β
        
        The spectral dependence is (λ²/λ²₀)^(-β/2) because:
        - ν ∝ 1/λ, so λ² ∝ 1/ν²
        - p(ν) ∝ ν^β implies p(λ²) ∝ (λ²)^(-β/2)
        
        Args:
            lambda_sq (np.ndarray): Wavelength squared in m²
            
        Returns:
            np.ndarray: Complex polarization P(λ²)
            
        Raises:
            ValueError: If lambda_sq_0 has not been set
        """
        if self.lambda_sq_0 is None:
            raise ValueError(
                "lambda_sq_0 must be set before computing polarization. "
                "Set via component.lambda_sq_0 = value or provide during initialization."
            )
        
        # Power-law spectral dependence: (λ²/λ²₀)^(-β/2)
        spectral_factor = (lambda_sq / self.lambda_sq_0)**(-self.beta / 2.0)
        
        # Full polarization with Faraday rotation
        return (self.a + 1j * self.b) * spectral_factor * np.exp(2j * self.phi_rm * lambda_sq)

    def compute_transfer_factor(self, lambda_sq_ref: float) -> complex:
        """
        Compute the transfer factor K at a reference wavelength squared.

        K is the model response at λ²_ref for unit intrinsic amplitude, so that
        the intrinsic prefactor A_int = A_ref / K, where A_ref = q_ref + i*u_ref
        is the complex amplitude sampled at the reference wavelength.

        For a power-law component |K| ≠ 1 in general (the spectral factor shifts
        the amplitude), so both p₀ and ψ₀ are affected by the conversion.

        Requires lambda_sq_0 to be set.

        Args:
            lambda_sq_ref (float): Reference wavelength squared in m²

        Returns:
            complex: Transfer factor K(λ²_ref)

        Raises:
            ValueError: If lambda_sq_0 has not been set
        """
        if self.lambda_sq_0 is None:
            raise ValueError(
                "lambda_sq_0 must be set before computing transfer factor. "
                "Set via component.lambda_sq_0 = value or provide during initialization."
            )
        spectral_factor = (lambda_sq_ref / self.lambda_sq_0) ** (-self.beta / 2.0)
        return spectral_factor * np.exp(2j * self.phi_rm * lambda_sq_ref)

    def __repr__(self):
        name_str = f"'{self.name}'" if self.name else "unnamed"
        p0, psi_0 = cartesian_to_polar(self.a, self.b)
        if self.lambda_sq_0 is not None:
            lambda_sq_0_str = f"{self.lambda_sq_0:.6f} m²"
            # Convert to frequency: ν₀ = c/√(λ²₀)
    
            nu_0_hz = C_LIGHT / np.sqrt(self.lambda_sq_0)
            nu_0_ghz = nu_0_hz / 1e9
            nu_0_str = f"{nu_0_ghz:.3f} GHz"
        else:
            lambda_sq_0_str = "not set"
            nu_0_str = "not set"
        return (f"ThinPowerLawComponent({name_str}: a={self.a:.4f}, b={self.b:.4f}, "
                f"RM={self.phi_rm:.1f} rad/m², β={self.beta:.2f}, λ²₀={lambda_sq_0_str}, ν₀={nu_0_str} | "
                f"p₀={p0:.4f}, ψ₀={psi_0:.3f} rad)")


class ThickComponent:
    """
    Faraday-thick polarized emission component (super-Gaussian in Faraday space).
    
    Represents emission with internal Faraday dispersion. The super-Gaussian
    form allows modeling various depolarization scenarios (Anderson et al. 2016):
    - N=2: Burn foreground screen (turbulent cells)
    - N>>2: Burn slab (mixed emission/rotation)
    - Intermediate N: RM gradients, complex structures
    
    Parameters:
        a (float): Q component of intrinsic polarization
        b (float): U component of intrinsic polarization
        phi_peak (float): Mean Faraday depth in rad m^-2
        sigma_phi (float): Faraday dispersion width in rad m^-2
        N (float): Shape parameter (N=2 is Gaussian)
        name (str): Optional name for tracking (e.g., 'Foreground', 'Halo')
    """
    
    def __init__(self, a: float, b: float, phi_peak: float, sigma_phi: float, 
                 N: float, name: Optional[str] = None):
        # Validate N parameter
        if N <= 0:
            raise ValueError(f"Shape parameter N must be positive, got N={N}")
        
        # Warn for unusually large N values
        if N > 50:
            import warnings
            warnings.warn(
                f"ThickComponent initialized with N={N:.1f} (>50). "
                f"Such large N values produce very sharp box-like profiles and approach "
                f"numerical overflow regions. While our log-space implementation handles "
                f"this stably, N>50 is uneeded. Anderson et al. (2016) typically "
                f"uses N=2 (Gaussian) to N~30 (sharp slab).",
                UserWarning,
                stacklevel=2
            )
        
        self.a = a
        self.b = b
        self.phi_peak = phi_peak
        self.sigma_phi = sigma_phi
        self.N = N
        self.name = name
    
    def compute_fdf(self, phi_grid: np.ndarray) -> np.ndarray:
        """
        Normalized super-Gaussian FDF with numerically stable log-space computation.
        
        F(φ) = (1/C_norm) * exp(-|φ - φ_peak|^N / (2*σ_φ^N)) * (a + ib)
        
        where C_norm = 2^(1+1/N) * σ_φ * Γ(1/N) / N ensures ∫|F(φ)|dφ = √(a² + b²)
        
        For large N (>50), computing (δ/σ)^N directly causes overflow. Instead,
        we compute it in log-space: (δ/σ)^N = exp(N * log(δ/σ)), which remains
        numerically stable for arbitrary N values.
        
        Args:
            phi_grid (np.ndarray): Faraday depth grid in rad m^-2
            
        Returns:
            np.ndarray: Complex FDF values
        """
        # Normalization: C_norm ensures integral of |F(φ)| equals √(a² + b²)
        C_norm = 2**(1 + 1/self.N) * self.sigma_phi * gamma(1/self.N) / self.N
        
        # Distance from peak
        delta = np.abs(phi_grid - self.phi_peak)
        
        # Compute (delta/sigma)^N in log-space to avoid overflow for large N
        ratio = delta / self.sigma_phi
        
        # Suppress warnings from log(0) at peak - handled explicitly below
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio_N = self.N * np.log(ratio)
        
        # Clip extreme values to prevent overflow in final exp()
        log_ratio_N_clipped = np.clip(log_ratio_N, -np.inf, 700.0)
        
        # Compute exponent: -0.5 * (delta/sigma)^N
        exponent = -0.5 * np.exp(log_ratio_N_clipped)
        
        # At peak (delta=0), exponent should be 0 to give profile=1
        exponent = np.where(delta > 0, exponent, 0.0)
        
        # Compute profile (safely underflows to 0 for very negative exponents)
        profile = np.exp(exponent)
        
        # Apply normalization and complex phase using Cartesian form
        # Normalization ensures ∫|F(φ)|dφ = p0 = √(a² + b²)
        fdf = (profile / C_norm) * (self.a + 1j * self.b)
        return fdf
    
    def compute_polarization(self, lambda_sq: np.ndarray) -> np.ndarray:
        """
        Compute complex polarization P(λ²) = Q + iU via Fourier transform.
        
        P(λ²) = ∫ F(φ) exp(i*2*φ*λ²) dφ
        
        For N ≈ 2 (Gaussian case), uses analytical formula for speed and accuracy.
        Otherwise, computes numerically by creating a Faraday depth grid.
        
        Args:
            lambda_sq (np.ndarray): Wavelength squared in m²
            
        Returns:
            np.ndarray: Complex polarization P(λ²)
        """
        # Use analytical Gaussian formula if N is very close to 2
        if np.abs(self.N - 2.0) < 1e-5:
            # Analytical solution for Gaussian (N=2):
            # P(λ²) = (a + ib) * exp(2i*φ_peak*λ²) * exp(-2*σ_φ²*λ⁴)
            phase = np.exp(2j * self.phi_peak * lambda_sq)
            depol = np.exp(-2 * self.sigma_phi**2 * lambda_sq**2)
            pol = (self.a + 1j * self.b) * phase * depol
            return pol
        
        # Numerical computation for non-Gaussian cases.
        # Grid parameters from swept calibration on σ_φ=50 rad/m², band 0.85–1.8 GHz.
        # Fit model: f(N) = a * N^b + c  (pinned through N=2, N_PHI inflated 20%).
        # To recalibrate, rerun optimize_num_int.py and update the six constants below.
        _nsig_a, _nsig_b, _nsig_c = 19.556, -2.190, 1.316
        _nphi_a, _nphi_b, _nphi_c = 7252.291, -4.310, 55.407

        n_sigma   = _nsig_a * self.N ** _nsig_b + _nsig_c
        phi_range = n_sigma * self.sigma_phi

        _nphi_raw = _nphi_a * self.N ** _nphi_b + _nphi_c
        _nphi_int = max(int(round(_nphi_raw)), 3)
        n_phi     = _nphi_int if _nphi_int % 2 == 1 else _nphi_int + 1

        phi_grid = np.linspace(self.phi_peak - phi_range,
                               self.phi_peak + phi_range, n_phi)
        dphi = phi_grid[1] - phi_grid[0]
        
        # Compute FDF on grid
        fdf = self.compute_fdf(phi_grid)
        
        # Fourier transform: P(λ²) = ∫ F(φ) exp(i*2*φ*λ²) dφ
        # Use numba-optimized version for ~100x speedup
        pol = fourier_transform_fdf_numba(fdf, phi_grid, lambda_sq)
        
        # Old pure Python version (commented for easy rollback):
        # pol = np.zeros(len(lambda_sq), dtype=complex)
        # for i, lsq in enumerate(lambda_sq):
        #     integrand = fdf * np.exp(2j * phi_grid * lsq)
        #     pol[i] = np.trapz(integrand, dx=dphi)
        
        return pol

    def compute_transfer_factor(self, lambda_sq_ref: float) -> complex:
        """
        Compute the transfer factor K at a reference wavelength squared.

        K is the model response at λ²_ref for unit intrinsic amplitude, so that
        the intrinsic prefactor A_int = A_ref / K, where A_ref = q_ref + i*u_ref
        is the complex amplitude sampled at the reference wavelength.

        For a thick component |K| < 1 in general due to depolarisation, so both
        p₀ and ψ₀ are affected by the conversion. K is computed numerically
        using the same integral as compute_polarization(), exploiting the fact
        that the FDF is linear in the intrinsic amplitude (a + ib).

        Args:
            lambda_sq_ref (float): Reference wavelength squared in m²

        Returns:
            complex: Transfer factor K(λ²_ref)
        """
        A_int = self.a + 1j * self.b
        pol = self.compute_polarization(np.array([lambda_sq_ref]))[0]
        return pol / A_int

    def __repr__(self):
        name_str = f"'{self.name}'" if self.name else "unnamed"
        p0, psi_0 = cartesian_to_polar(self.a, self.b)
        return (f"ThickComponent({name_str}: a={self.a:.4f}, b={self.b:.4f}, φ_peak={self.phi_peak:.1f} rad/m², "
                f"σ_φ={self.sigma_phi:.1f} rad/m², N={self.N:.2f} | p₀={p0:.4f}, ψ₀={psi_0:.3f} rad)")


class CompositeModel:
    """
    Composite model combining multiple thin, thick, and/or power-law components.
    
    Represents the total polarized emission as a sum of individual components:
    F(φ) = F_1(φ) + F_2(φ) + ... + F_n(φ)
    P(λ²) = P_1(λ²) + P_2(λ²) + ... + P_n(λ²)
    
    Example:
        model = CompositeModel()
        model.add_component(ThinComponent(a=0.04, b=0.03, phi_rm=50))
        model.add_component(ThickComponent(a=0.02, b=-0.02, phi_peak=120, 
                                          sigma_phi=80, N=2))
        model.add_component(ThinPowerLawComponent(a=0.03, b=0.02, phi_rm=80, beta=-0.7))
    """
    
    def __init__(self):
        self.components: List[ThinComponent | ThinPowerLawComponent | ThickComponent] = []
    
    def add_component(self, component: ThinComponent | ThinPowerLawComponent | ThickComponent):
        """Add a thin, power-law, or thick component to the model."""
        self.components.append(component)
    
    @property
    def n_components(self) -> int:
        """Number of components in the model."""
        return len(self.components)
    
    @property
    def model_type(self) -> str:
        """String representation of model type (e.g., 'ST', 'TPT', 'SSP')."""
        type_str = ""
        for comp in self.components:
            if isinstance(comp, ThinPowerLawComponent):
                type_str += "P"
            elif isinstance(comp, ThinComponent):
                type_str += "S"
            elif isinstance(comp, ThickComponent):
                type_str += "T"
        return type_str
    
    def compute_fdf(self, phi_grid: np.ndarray) -> np.ndarray:
        """
        Compute total Faraday dispersion function.
        
        F(φ) = Σ_i F_i(φ)
        
        Args:
            phi_grid (np.ndarray): Faraday depth grid in rad m^-2
            
        Returns:
            np.ndarray: Complex FDF values
        """
        fdf_total = np.zeros(len(phi_grid), dtype=complex)
        for component in self.components:
            fdf_total += component.compute_fdf(phi_grid)
        return fdf_total
    
    def compute_polarization(self, lambda_sq: np.ndarray) -> np.ndarray:
        """
        Compute total complex polarization P(λ²).
        
        P(λ²) = Σ_i P_i(λ²)
        
        Args:
            lambda_sq (np.ndarray): Wavelength squared in m²
            
        Returns:
            np.ndarray: Complex polarization P(λ²)
        """
        pol_total = np.zeros(len(lambda_sq), dtype=complex)
        for component in self.components:
            pol_total += component.compute_polarization(lambda_sq)
        return pol_total
    
    def simulate_data(self, lambda_sq: np.ndarray, 
                     significance: float = 10.0,
                     I_spectrum: Optional[np.ndarray] = None,
                     seed: Optional[int] = None) -> PolarizationData:
        """
        Simulate polarization observations with receiver noise.
        
        Adds equal additive noise to Stokes I, Q, U (realistic receiver behavior).
        Noise level set by peak polarization SNR: p_peak / σ = significance.
        
        Args:
            lambda_sq: Wavelength squared array (m²)
            significance: SNR at peak polarization channel (default: 10)
            I_spectrum: Noiseless Stokes I spectrum. If None, assumes I=1
            seed: Random seed for reproducibility
            
        Returns:
            PolarizationData with noisy I, Q, U in fractional space
        """
        if seed is not None:
            np.random.seed(seed)
        
        # Noiseless I spectrum
        if I_spectrum is None:
            I_true = np.ones_like(lambda_sq)
        else:
            I_true = np.asarray(I_spectrum, dtype=float)
        
        # Noiseless fractional polarization
        P_frac = self.compute_polarization(lambda_sq)
        
        # Convert to absolute space
        P_abs = P_frac * I_true
        
        # Noise level from peak polarization
        P_peak = np.max(np.abs(P_abs))
        sigma_abs = P_peak / significance
        
        # Add same noise to I, Q, U (receiver noise is independent of Stokes parameter)
        I_obs = I_true + np.random.normal(0, sigma_abs, len(lambda_sq))
        Q_abs = np.real(P_abs) + np.random.normal(0, sigma_abs, len(lambda_sq))
        U_abs = np.imag(P_abs) + np.random.normal(0, sigma_abs, len(lambda_sq))
        
        # Convert to fractional
        Q_frac = Q_abs / I_obs
        U_frac = U_abs / I_obs
        
        # Fractional errors (same absolute noise divided by observed I)
        Q_err_frac = sigma_abs / I_obs
        U_err_frac = sigma_abs / I_obs
        I_err = np.full_like(lambda_sq, sigma_abs)

        # Local import avoids circular import at module load time.
        from faraday_data import PolarizationData
        
        return PolarizationData(lambda_sq, Q_frac, U_frac, Q_err_frac, U_err_frac, 
                               I=I_obs, I_err=I_err)

    def __repr__(self):
        if len(self.components) == 0:
            return f"CompositeModel (type: {self.model_type}, no components)"
        
        comp_reprs = []
        for i, comp in enumerate(self.components, 1):
            comp_reprs.append(f"  [{i}] {repr(comp)}")
        
        comp_list = "\n".join(comp_reprs)
        return f"CompositeModel (type: {self.model_type}, {len(self.components)} components):\n{comp_list}"


# Helper functions

def create_lambda_sq_grid(freq_min: float, freq_max: float, 
                          n_channels: int,
                          gap: tuple = None) -> np.ndarray:
    """
    Create wavelength squared grid from frequency range.
    
    Args:
        freq_min (float): Minimum frequency in GHz
        freq_max (float): Maximum frequency in GHz
        n_channels (int): Number of spectral channels
        gap (tuple): Optional (gap_start_GHz, gap_end_GHz) to create
                     discontinuous frequency coverage
        
    Returns:
        np.ndarray: Wavelength squared in m²
    """
    if gap is None:
        # Contiguous frequency coverage
        freq_hz = np.linspace(freq_min * 1e9, freq_max * 1e9, n_channels)
    else:
        # Discontinuous frequency coverage with gap
        gap_start, gap_end = gap
        
        # Calculate bandwidth in each band
        low_bandwidth = gap_start - freq_min
        high_bandwidth = freq_max - gap_end
        total_usable = low_bandwidth + high_bandwidth
        
        # Distribute channels proportionally to bandwidth
        n_low = int(np.round(n_channels * low_bandwidth / total_usable))
        n_high = n_channels - n_low
        
        # Create two frequency ranges
        freq_low = np.linspace(freq_min * 1e9, gap_start * 1e9, n_low)
        freq_high = np.linspace(gap_end * 1e9, freq_max * 1e9, n_high)
        
        # Concatenate
        freq_hz = np.concatenate([freq_low, freq_high])
    
    lambda_m = C_LIGHT / freq_hz
    return lambda_m**2


def freq_to_lambda_sq(freq_ghz: float) -> float:
    """Convert frequency in GHz to wavelength squared in m²."""
    lambda_m = C_LIGHT / (freq_ghz * 1e9)
    return lambda_m**2


def compute_rmtf(lambda_sq: np.ndarray, phi_grid: np.ndarray,
                 weights: Optional[np.ndarray] = None,
                 weight_type: str = 'uniform',
                 Q_err: Optional[np.ndarray] = None,
                 U_err: Optional[np.ndarray] = None) -> Tuple[np.ndarray, dict]:
    """
    Compute the Rotation Measure Transfer Function (RMTF).
    
    Uses λ₀² = 0 for full Faraday depth resolution following Rudnick & Cotton (2023).
    This provides ~2-3× better resolution than the traditional Brentjens & de Bruyn 
    (2005) approach using weighted mean λ₀².
    
    RMTF: R(φ) = K Σᵢ wᵢ exp(-2iφλᵢ²)
    
    Args:
        lambda_sq (np.ndarray): Wavelength squared sampling in m²
        phi_grid (np.ndarray): Faraday depth grid in rad m^-2
        weights (Optional[np.ndarray]): Custom weights for each channel
        weight_type (str): 'uniform' or 'variance' (inverse-variance weighting)
            Only used if weights=None
        Q_err (Optional[np.ndarray]): Q errors for variance weighting
        U_err (Optional[np.ndarray]): U errors for variance weighting
    
    Returns:
        Tuple[np.ndarray, dict]: 
            - Complex RMTF values on phi_grid
            - Dictionary with RMTF properties (fwhm, etc.)
    """
    if weights is None:
        if weight_type == 'uniform':
            weights = np.ones(len(lambda_sq))
        elif weight_type == 'variance':
            if Q_err is None or U_err is None:
                raise ValueError("weight_type='variance' requires Q_err and U_err")
            var_Q = Q_err**2
            var_U = U_err**2
            weights = 1.0 / (var_Q + var_U)
        else:
            raise ValueError(f"Unknown weight_type: {weight_type}")
    
    # Normalization constant K (Eq. 38 from Brentjens & de Bruyn 2005)
    K = 1.0 / np.sum(weights)
    
    # Use λ₀² = 0 for full resolution (Rudnick & Cotton 2023)
    # This gives narrower RMTF and better discrimination of Faraday components
    lambda_0_sq = 0.0
    
    # Traditional Brentjens & de Bruyn (2005) approach uses weighted mean:
    # This broadens the real/imaginary RMTF by ~factor 2-3 but doesn't affect amplitude spectrum
    lambda_0_sq = np.sum(weights * lambda_sq) / np.sum(weights)
    
    # Compute RMTF - vectorized
    # Note: -2j is correct (negative sign for inverse FT)
    a = -2.0j * phi_grid  # shape: (n_phi,)
    b = lambda_sq - lambda_0_sq  # shape: (n_lambda,)
    rmtf = K * np.sum(weights * np.exp(np.outer(a, b)), axis=1)
    
    # Calculate RMTF properties
    lambda_max_sq = np.max(lambda_sq)
    lambda_min_sq = np.min(lambda_sq)
    delta_lambda_sq = lambda_max_sq - lambda_min_sq
    
    # FWHM formulas from Rudnick & Cotton (2023)
    # Nominal (Brentjens & de Bruyn approach with weighted mean λ₀²)
    expected_fwhm_nominal = 3.8 / delta_lambda_sq
    
    # Full resolution (λ₀² = 0 approach)
    expected_fwhm_full = 2.0 / (lambda_max_sq + lambda_min_sq)
    
    max_scale = np.pi / lambda_min_sq
    
    rmtf_props = {
        'lambda_0_sq': lambda_0_sq,
        'lambda_min_sq': lambda_min_sq,
        'lambda_max_sq': lambda_max_sq,
        'delta_lambda_sq': delta_lambda_sq,
        'expected_fwhm_nominal': expected_fwhm_nominal,
        'expected_fwhm_full': expected_fwhm_full,
        'max_scale': max_scale,
        'K': K
    }
    
    return rmtf, rmtf_props



def rm_synthesis(data: PolarizationData, phi_grid: np.ndarray,
                 weights: Optional[np.ndarray] = None,
                 weight_type: str = 'uniform') -> np.ndarray:
    """
    Perform RM synthesis on observed data to get FDF.
    
    Uses λ₀² = 0 for full Faraday depth resolution following Rudnick & Cotton (2023).
    
    F̃(φ) = K Σᵢ P̃(λᵢ²) exp(-2iφλᵢ²)
    
    This directly transforms observed polarization data (including noise)
    from λ² space to Faraday depth space.
    
    Args:
        data (PolarizationData): Observed polarization data with Q, U
        phi_grid (np.ndarray): Faraday depth grid in rad m^-2
        weights (Optional[np.ndarray]): Custom weights for each channel
        weight_type (str): 'uniform' or 'variance' for weighting
    
    Returns:
        np.ndarray: Observed FDF F̃(φ) (complex)
    """
    lambda_sq = data.lambda_sq
    
    # Complex polarization P = Q + iU (with noise!)
    P_obs = data.Q + 1j * data.U
    
    if weights is None:
        if weight_type == 'uniform':
            weights = np.ones(len(lambda_sq))
        elif weight_type == 'variance':
            # Inverse variance weighting
            if data.Q_err is not None and data.U_err is not None:
                var_Q = data.Q_err**2
                var_U = data.U_err**2
                weights = 1.0 / (var_Q + var_U)
            else:
                weights = np.ones(len(lambda_sq))
        else:
            raise ValueError(f"Unknown weight_type: {weight_type}")
    
    # Normalization constant K (Eq. 38 from Brentjens & de Bruyn 2005)
    K = 1.0 / np.sum(weights)
    
    # Use λ₀² = 0 for full resolution (Rudnick & Cotton 2023)
    lambda_0_sq = 0.0
    
    # Traditional Brentjens & de Bruyn (2005) uses weighted mean:
    lambda_0_sq = np.sum(weights * lambda_sq) / np.sum(weights)
    
    # RM synthesis - vectorized
    a = -2.0j * phi_grid  # shape: (n_phi,)
    b = lambda_sq - lambda_0_sq  # shape: (n_lambda,)
    fdf_obs = K * np.sum(weights * P_obs * np.exp(np.outer(a, b)), axis=1)
    
    return fdf_obs


def compute_polarization_quantities(Q: np.ndarray, U: np.ndarray, 
                                    dQ: np.ndarray, dU: np.ndarray,
                                    method: str = 'mc', 
                                    n_samples: int = 100,
                                    debias: bool = True) -> Tuple:
    """
    Compute polarization fraction (P) and EVPA (χ) with errors.
    
    Uses either Monte Carlo sampling or analytic error propagation.
    MC method properly handles nonlinear transformations and provides
    asymmetric errors, especially important for low S/N where Rice bias matters.
    
    Args:
        Q: Stokes Q parameter (array)
        U: Stokes U parameter (array)
        dQ: Error on Q (array)
        dU: Error on U (array)
        method: 'mc' (Monte Carlo, default) or 'analytic' (Gaussian propagation)
        n_samples: Number of MC samples (default: 100, only used if method='mc')
        debias: Apply Ricean debiasing to polarization intensity (default: True).
                Uses the Hales (2012) estimator:
                  σ_QU² = 0.8·σ_max² + 0.2·σ_min²
                where σ_max = max(σ_Q, σ_U) per channel, then
                  P0 = sqrt(max(P² − σ_QU², 0))
                Does not affect EVPA.
        
    Returns:
        P: Polarization fraction (array).  If debias=True this is the
           de-biased P0; otherwise the raw Rice-biased |P|.
        P_err: Error on P
               - If MC: tuple (P_err_lower, P_err_upper) for asymmetric errors
               - If analytic: array (symmetric errors)
        chi: EVPA in radians (array)
        chi_err: Error on chi
               - If MC: tuple (chi_err_lower, chi_err_upper) for asymmetric errors
               - If analytic: array (symmetric errors)
    
    Notes:
        - MC method is recommended for accurate error representation
        - Analytic method is faster but assumes Gaussian errors
        - EVPA (chi) is in range [-π/2, π/2]
        - Polarization fraction P is always positive (zero where debiased value would be negative)
    """
    # Hales (2012) noise estimator: weight towards the larger noise axis
    # σ_QU² = 0.8·σ_max² + 0.2·σ_min²
    if debias:
        sigma_max2 = np.maximum(dQ**2, dU**2)
        sigma_min2 = np.minimum(dQ**2, dU**2)
        sigma_QU2 = 0.8 * sigma_max2 + 0.2 * sigma_min2  # (n_channels,)

    if method == 'mc':
        # Monte Carlo: Sample Q, U and compute derived quantities
        # Shape: (n_channels, n_samples)
        Q_samples = np.random.normal(Q[:, np.newaxis], dQ[:, np.newaxis], 
                                     (len(Q), n_samples))
        U_samples = np.random.normal(U[:, np.newaxis], dU[:, np.newaxis], 
                                     (len(U), n_samples))
        
        # Compute P and chi for all samples (vectorized)
        P_samples = np.sqrt(Q_samples**2 + U_samples**2)
        chi_samples = 0.5 * np.arctan2(U_samples, Q_samples)
        
        # Ricean debiasing inside the sample distribution (Hales 2012)
        # Subtract the known noise bias from each sample before computing statistics.
        # Channels where P² < σ_QU² are consistent with zero polarisation → clamp to 0.
        if debias:
            P_sq_debiased = np.maximum(P_samples**2 - sigma_QU2[:, np.newaxis], 0.0)
            P_samples = np.sqrt(P_sq_debiased)

        # Extract percentiles (16%, 50%, 84% for 1σ equivalent)
        P = np.percentile(P_samples, 50, axis=1)
        P_lower = P - np.percentile(P_samples, 16, axis=1)
        P_upper = np.percentile(P_samples, 84, axis=1) - P
        
        chi = np.percentile(chi_samples, 50, axis=1)
        chi_lower = chi - np.percentile(chi_samples, 16, axis=1)
        chi_upper = np.percentile(chi_samples, 84, axis=1) - chi
        
        return P, (P_lower, P_upper), chi, (chi_lower, chi_upper)
    
    elif method == 'analytic':
        # Analytic error propagation (Gaussian approximation)
        P_biased = np.sqrt(Q**2 + U**2)
        
        # Avoid division by zero
        P_safe = np.where(P_biased > 0, P_biased, 1.0)
        dP = np.sqrt((Q*dQ)**2 + (U*dU)**2) / P_safe
        dP = np.where(P_biased > 0, dP, 0.0)
        
        if debias:
            # P0 = sqrt(max(P² − σ_QU², 0)); propagate: dP0 ≈ dP · P_biased / P0
            P_sq_debiased = np.maximum(P_biased**2 - sigma_QU2, 0.0)
            P = np.sqrt(P_sq_debiased)
            P0_safe = np.where(P > 0, P, 1.0)
            dP = np.where(P > 0, dP * P_biased / P0_safe, dP)
        else:
            P = P_biased
        
        chi = 0.5 * np.arctan2(U, Q)
        
        # Error on chi (handle P=0 case)
        P_sq = Q**2 + U**2
        P_sq_safe = np.where(P_sq > 0, P_sq, 1.0)
        dchi = 0.5 * np.sqrt((Q*dU)**2 + (U*dQ)**2) / P_sq_safe
        dchi = np.where(P_sq > 0, dchi, 0.0)
        
        return P, dP, chi, dchi
    
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'mc' or 'analytic'.")


# =============================================================================
# NUMBA SPEED UPS
# =============================================================================

if NUMBA_AVAILABLE:
    @numba.jit(nopython=True, parallel=True)
    def inverse_transform_fdf_numba(P_lsq, lambda_sq_grid, phi_grid):
        """
        Compute F(φ) = ∫ P(λ²) exp(-i*2*φ*λ²) dλ² using numba for speed.
        
        Performs inverse Fourier transform of polarization spectrum to 
        Faraday depth space using trapezoidal integration. Parallelizes
        over Faraday depths for significant speedup.
        
        Used by ThinPowerLawComponent for FDF computation.
        
        Args:
            P_lsq: Complex polarization P(λ²) on lambda_sq_grid (1D array)
            lambda_sq_grid: Wavelength squared grid in m² (1D array)
            phi_grid: Faraday depth values in rad m^-2 (1D array)
            
        Returns:
            Complex FDF F(φ) for each Faraday depth (1D array)
            
        Performance:
            ~100x faster than pure Python loop with parallel=True
            (500 phi points × 2000 lambda_sq grid points)
        """
        n_phi = len(phi_grid)
        n_lambda_sq = len(lambda_sq_grid)
        d_lambda_sq = lambda_sq_grid[1] - lambda_sq_grid[0]
        
        fdf = np.zeros(n_phi, dtype=np.complex128)
        
        # Parallel over Faraday depths
        for i in numba.prange(n_phi):
            phi = phi_grid[i]
            
            # Compute integrand: P(λ²) * exp(-i*2*φ*λ²)
            phase = np.exp(-2j * phi * lambda_sq_grid)
            integrand = P_lsq * phase
            
            # Trapezoidal integration (manual implementation for numba)
            integral = integrand[0] + integrand[-1]  # End points
            integral += 2.0 * np.sum(integrand[1:-1])  # Interior points
            integral *= d_lambda_sq / 2.0
            
            fdf[i] = integral
        
        return fdf
    
    @numba.jit(nopython=True, parallel=False)
    def fourier_transform_fdf_numba(fdf, phi_grid, lambda_sq):
        """
        Compute P(λ²) = ∫ F(φ) exp(i*2*φ*λ²) dφ using numba for speed.
        
        Performs Fourier transform of Faraday dispersion function to 
        wavelength space using trapezoidal integration. Parallelizes
        over wavelengths for additional speedup.
        
        Args:
            fdf: Complex FDF values on phi_grid (1D array)
            phi_grid: Faraday depth grid in rad m^-2 (1D array)
            lambda_sq: Wavelength squared values in m² (1D array)
            
        Returns:
            Complex polarization P(λ²) for each wavelength (1D array)
            
        Performance:
            ~100x faster than pure Python loop for typical cases
            (100 wavelengths × 101 phi grid points)
        """
        n_lambda = len(lambda_sq)
        n_phi = len(phi_grid)
        dphi = phi_grid[1] - phi_grid[0]
        
        pol = np.zeros(n_lambda, dtype=np.complex128)
        
        # Parallel over wavelengths
        for i in numba.prange(n_lambda):
            lsq = lambda_sq[i]
            
            # Compute integrand: F(φ) * exp(i*2*φ*λ²)
            phase = np.exp(2j * phi_grid * lsq)
            integrand = fdf * phase
            
            # Trapezoidal integration (manual implementation for numba)
            integral = integrand[0] + integrand[-1]  # End points
            integral += 2.0 * np.sum(integrand[1:-1])  # Interior points
            integral *= dphi / 2.0
            
            pol[i] = integral
        
        return pol
else:
    # Fallback: pure numpy implementations (slower)
    def inverse_transform_fdf_numba(P_lsq, lambda_sq_grid, phi_grid):
        """Fallback implementation when numba is not available."""
        d_lambda_sq = lambda_sq_grid[1] - lambda_sq_grid[0]
        fdf = np.zeros(len(phi_grid), dtype=complex)
        for i, phi in enumerate(phi_grid):
            integrand = P_lsq * np.exp(-2j * phi * lambda_sq_grid)
            fdf[i] = np.trapz(integrand, dx=d_lambda_sq)
        return fdf
    
    def fourier_transform_fdf_numba(fdf, phi_grid, lambda_sq):
        """Fallback implementation when numba is not available."""
        dphi = phi_grid[1] - phi_grid[0]
        pol = np.zeros(len(lambda_sq), dtype=complex)
        for i, lsq in enumerate(lambda_sq):
            integrand = fdf * np.exp(2j * phi_grid * lsq)
            pol[i] = np.trapz(integrand, dx=dphi)
        return pol