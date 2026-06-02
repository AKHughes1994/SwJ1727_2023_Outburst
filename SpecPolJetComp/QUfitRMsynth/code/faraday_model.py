# faraday_model.py
"""
Bayesian fitting for Faraday rotation analysis using dynesty.

Fits models in Cartesian (a, b) parameterization and provides preliminary
parameter extraction using dynesty's built-in statistical functions.
"""

import numpy as np
import pickle
import copy
import time
import json
import multiprocessing
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Union, Literal
from scipy import stats
from pathlib import Path
from multiprocessing import Pool
from faraday_utils import (ThinComponent, ThinPowerLawComponent, ThickComponent,
                           CompositeModel, C_LIGHT)
from faraday_data import PolarizationData, validate_data

import dynesty
from dynesty import utils as dyfunc

try:
    import numba
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    import warnings
    warnings.warn("numba not available - chi-squared computation will be slower", RuntimeWarning)


# =============================================================================
# PRIOR SETUP & VALIDATION
# =============================================================================

SUPPORTED_PRIOR_TYPES = ['uniform', 'log_uniform', 'gaussian', 'truncated_gaussian', 'fixed']
PSI_0_PEAKED_K = 1.0


def validate_prior_spec(prior_spec: Tuple) -> Tuple[str, Tuple]:
    """
    Validate a single prior specification.
    
    Args:
        prior_spec: Tuple defining prior (type, param1, param2, ...)
        
    Returns:
        Tuple of (prior_type, parameters)
        
    Raises:
        ValueError: If prior specification is invalid
    """
    if not isinstance(prior_spec, tuple) or len(prior_spec) < 2:
        raise ValueError(f"Prior must be tuple (type, params...), got: {prior_spec}")
    
    prior_type = prior_spec[0]
    params = prior_spec[1:]
    
    if prior_type not in SUPPORTED_PRIOR_TYPES:
        raise ValueError(
            f"Unknown prior type '{prior_type}'. "
            f"Supported types: {SUPPORTED_PRIOR_TYPES}"
        )
    
    # Validate parameter counts and values
    if prior_type == 'uniform':
        if len(params) != 2:
            raise ValueError(f"uniform prior needs (min, max), got {len(params)} params")
        if params[0] >= params[1]:
            raise ValueError(f"uniform prior: min={params[0]} >= max={params[1]}")
    
    elif prior_type == 'log_uniform':
        if len(params) != 2:
            raise ValueError(f"log_uniform prior needs (min, max), got {len(params)} params")
        if params[0] <= 0 or params[1] <= 0:
            raise ValueError(f"log_uniform prior: both min and max must be > 0")
        if params[0] >= params[1]:
            raise ValueError(f"log_uniform prior: min={params[0]} >= max={params[1]}")
    
    elif prior_type == 'gaussian':
        if len(params) != 2:
            raise ValueError(f"gaussian prior needs (mean, std), got {len(params)} params")
        if params[1] <= 0:
            raise ValueError(f"gaussian prior: std={params[1]} must be > 0")
    
    elif prior_type == 'truncated_gaussian':
        if len(params) != 4:
            raise ValueError(
                f"truncated_gaussian prior needs (mean, std, min, max), got {len(params)} params"
            )
        if params[1] <= 0:
            raise ValueError(f"truncated_gaussian prior: std={params[1]} must be > 0")
        if params[2] >= params[3]:
            raise ValueError(
                f"truncated_gaussian prior: min={params[2]} >= max={params[3]}"
            )
    
    elif prior_type == 'fixed':
        if len(params) != 1:
            raise ValueError(f"fixed prior needs (value,), got {len(params)} params")
    
    return prior_type, params


def generate_prior_template(model: CompositeModel) -> Dict:
    """
    Generate a prior template dictionary for a given model.
    
    This provides suggested priors that the user can customize.
    
    Args:
        model: CompositeModel to generate template for
        
    Returns:
        Dictionary with suggested prior structure
    """
    priors = {}
    
    for i, comp in enumerate(model.components):
        comp_key = f'component_{i}'
        
        if isinstance(comp, ThinPowerLawComponent):
            priors[comp_key] = {
                'type': 'ThinPowerLawComponent',
                'polar_conversion': {
                    'p0_bounds': (0.0, 0.7, 'radial_uniform'),
                    'psi_0_bounds': (-np.pi/2, np.pi/2, 'uniform')
                },
                'phi_rm': ('uniform', -500.0, 500.0),
                'beta': ('uniform', -2.0, 2.0)
            }
        
        elif isinstance(comp, ThinComponent):
            priors[comp_key] = {
                'type': 'ThinComponent',
                'polar_conversion': {
                    'p0_bounds': (0.0, 0.7, 'radial_uniform'),
                    'psi_0_bounds': (-np.pi/2, np.pi/2, 'uniform')
                },
                'phi_rm': ('uniform', -500.0, 500.0)
            }
        
        elif isinstance(comp, ThickComponent):
            priors[comp_key] = {
                'type': 'ThickComponent',
                'polar_conversion': {
                    'p0_bounds': (0.0, 0.7, 'radial_uniform'),
                    'psi_0_bounds': (-np.pi/2, np.pi/2, 'uniform')
                },
                'phi_peak': ('uniform', -500.0, 500.0),
                'sigma_phi': ('log_uniform', 1.0, 200.0),
                'N': ('uniform', 2.0, 30.0)
            }
    
    return priors


def validate_priors(model: CompositeModel, priors: Dict, p0_physical_max: float = 0.7,
                    tie_phi_rm: bool = False) -> bool:
    """
    Validate that priors match the model structure.
    
    Supports two styles:
    1. Cartesian: Direct priors on (a, b)
    2. Polar conversion: Bounds on (p0, psi_0) converted internally to (a, b)
    
    For polar conversion, p0_bounds accepts three prior distributions:
    
    **'radial_uniform'**: p(p0) = const
        - Uniform in radius; any p0 is equally likely
        
    **'log_uniform'**: p(p0) ∝ 1/p0
        - Log uniform in radius; smaller p0 values are more likely
        
    **'area_uniform'**: p(p0) ∝ p0
        - Uniform in area; any (Q, U) is equally likely
        - Natural uninformative prior for Q-U plane sampling
    
    Example:
        priors = {
            'ISM': {
                'type': 'ThickComponent',
                'polar_conversion': {
                    'p0_bounds': (0.01, 0.5, 'area_uniform'),
                    'psi_0_bounds': (-np.pi/2, np.pi/2)
                },
                ...
            }
        }
    
    Checks:
    - Correct number of components
    - Component types match
    - All required parameters have priors
    - Prior specifications are valid
    - Physical constraints (p0 ≤ p0_physical_max)
    - Mutual exclusivity (cannot specify both styles)
    
    Args:
        model: CompositeModel to fit
        priors: Prior dictionary to validate
        p0_physical_max: Maximum allowed polarization fraction (default: 0.7)
        
    Returns:
        True if valid
        
    Raises:
        ValueError: If validation fails (with helpful error message)
    """
    import warnings
    
    n_components = len(model.components)
    
    special_keys = {'tied_phi_rm'}
    component_keys = [k for k in priors.keys() if k not in special_keys]
    n_priors = len(component_keys)
    
    # Print validation header
    print("\nValidating priors...")
    
    # Print all priors for debugging
    print("\n" + "="*70)
    print("PRIOR SPECIFICATION")
    print("="*70)
    for key, value in priors.items():
        if key in special_keys:
            print(f"\n{key}: {value}")
        else:
            print(f"\n{key}:")
            for param, prior_spec in value.items():
                if param == 'polar_conversion':
                    print(f"  {param}:")
                    for sub_param, sub_spec in prior_spec.items():
                        print(f"    {sub_param}: {sub_spec}")
                else:
                    print(f"  {param}: {prior_spec}")
    print("="*70)
    
    # Check number of components
    if n_priors != n_components:
        template = generate_prior_template(model)
        raise ValueError(
            f"Prior dictionary has {n_priors} components but model has {n_components}.\n\n"
            f"Expected prior template:\n{_format_prior_dict(template)}"
        )
    
    # Validate each component
    print("\n" + "="*70)
    print("VALIDATING COMPONENT PRIORS")
    print("="*70)
    
    for i, comp in enumerate(model.components):
        # Check for component name first, fall back to indexed name
        comp_name = comp.name if comp.name else f'component_{i}'
        indexed_name = f'component_{i}'
        
        # Try both possible keys
        if comp_name in priors:
            comp_key = comp_name
        elif indexed_name in priors:
            comp_key = indexed_name
        else:
            raise ValueError(
                f"Missing prior for component {i} (name='{comp.name}'). "
                f"Expected key '{comp_name}' or '{indexed_name}'"
            )
        
        comp_priors = priors[comp_key]
        
        print(f"\n• Component: {comp_key} (index {i})")
        
        # Check component type specified
        if 'type' not in comp_priors:
            raise ValueError(f"Prior for '{comp_key}' must specify 'type'")
        
        # Determine expected type
        if isinstance(comp, ThinPowerLawComponent):
            expected_type = 'ThinPowerLawComponent'
        elif isinstance(comp, ThinComponent):
            expected_type = 'ThinComponent'
        elif isinstance(comp, ThickComponent):
            expected_type = 'ThickComponent'
        else:
            raise ValueError(f"Unknown component type: {type(comp)}")
        
        if comp_priors['type'] != expected_type:
            raise ValueError(
                f"Prior '{comp_key}' has type '{comp_priors['type']}' but "
                f"model component is {expected_type}"
            )
        
        # Detect prior style
        has_cartesian = ('a' in comp_priors and 'b' in comp_priors)
        has_polar = 'polar_conversion' in comp_priors
        
        # Check mutual exclusivity
        if has_cartesian and has_polar:
            raise ValueError(
                f"Component '{comp_key}': specify EITHER (a, b) OR polar_conversion, not both"
            )
        
        if not has_cartesian and not has_polar:
            raise ValueError(
                f"Component '{comp_key}': must specify EITHER (a, b) OR polar_conversion"
            )
        
        # Validate style-specific requirements
        if has_cartesian:
            # Cartesian style validation
            if 'a' not in comp_priors or 'b' not in comp_priors:
                raise ValueError(
                    f"Component '{comp_key}': must specify both 'a' AND 'b' for Cartesian priors"
                )
            
            # Validate prior specifications
            validate_prior_spec(comp_priors['a'])
            validate_prior_spec(comp_priors['b'])
            
            # Check physical constraint
            a_bounds = _get_prior_bounds(comp_priors['a'])
            b_bounds = _get_prior_bounds(comp_priors['b'])
            
            if a_bounds is not None and b_bounds is not None:
                # Maximum p0 is achieved at the bounds with largest absolute values
                # For asymmetric bounds like a ∈ [-0.9, 0.2], need to check both ends
                max_abs_a = max(abs(a_bounds[0]), abs(a_bounds[1]))
                max_abs_b = max(abs(b_bounds[0]), abs(b_bounds[1]))
                max_possible_p0 = np.sqrt(max_abs_a**2 + max_abs_b**2)
                
                if max_possible_p0 > p0_physical_max:
                    warnings.warn(
                        f"Component '{comp_key}': Cartesian bounds allow p0 up to {max_possible_p0:.2f} > {p0_physical_max:.2f}. "
                        f"Samples with p0 > {p0_physical_max:.2f} will be rejected during sampling. "
                        f"This can cause boundary pile-ups that masquerade as modes.",
                        UserWarning
                    )
                
                # Print confirmation of Cartesian style
                print(f"   {comp_key}: Using Cartesian priors")
                print(f"      a: [{a_bounds[0]:.4f}, {a_bounds[1]:.4f}]")
                print(f"      b: [{b_bounds[0]:.4f}, {b_bounds[1]:.4f}]")
                print(f"      (max p0: {max_possible_p0:.3f})")
                print(f"      Note: p(p0) ∝ p0 → favors higher polarizations")
        
        elif has_polar:
            # Polar conversion style validation
            polar_spec = comp_priors['polar_conversion']
            
            if 'p0_bounds' not in polar_spec:
                raise ValueError(
                    f"Component '{comp_key}': polar_conversion must include 'p0_bounds'"
                )
            
            if 'psi_0_bounds' not in polar_spec:
                raise ValueError(
                    f"Component '{comp_key}': polar_conversion must include 'psi_0_bounds'"
                )
            
            p0_bounds = polar_spec['p0_bounds']
            psi_0_bounds = polar_spec['psi_0_bounds']

            # Validate p0_bounds format
            if not isinstance(p0_bounds, tuple) or len(p0_bounds) < 2:
                raise ValueError(
                    f"Component '{comp_key}': p0_bounds must be tuple of (min, max) or (min, max, distribution)"
                )
            
            p0_min, p0_max = p0_bounds[0], p0_bounds[1]
            distribution = p0_bounds[2] if len(p0_bounds) == 3 else 'radial_uniform'
            
            # Validate distribution type
            valid_distributions = ['radial_uniform', 'log_uniform', 'area_uniform']
            if distribution not in valid_distributions:
                raise ValueError(
                    f"Component '{comp_key}': distribution must be one of {valid_distributions}, got '{distribution}'"
                )
            
            # Check p0 bounds
            if p0_min < 0:
                raise ValueError(
                    f"Component '{comp_key}': p0_min must be >= 0, got {p0_min}"
                )
            
            if p0_min >= p0_max:
                raise ValueError(
                    f"Component '{comp_key}': p0_min ({p0_min}) must be < p0_max ({p0_max})"
                )
            
            # Log-uniform requires p0_min > 0
            if distribution == 'log_uniform' and p0_min <= 0:
                raise ValueError(
                    f"Component '{comp_key}': log_uniform distribution requires p0_min > 0, got {p0_min}"
                )
            
            # Check physical constraint
            if p0_max > p0_physical_max:
                warnings.warn(
                    f"Component '{comp_key}': p0_max ({p0_max:.2f}) > physical maximum ({p0_physical_max:.2f}). "
                    f"Samples with p0 > {p0_physical_max:.2f} will be rejected during sampling.",
                    UserWarning
                )

            # Validate psi_0_bounds format
            if not isinstance(psi_0_bounds, tuple) or len(psi_0_bounds) not in (2, 3):
                raise ValueError(
                    f"Component '{comp_key}': psi_0_bounds must be tuple of (min, max) "
                    f"or (min, max, angle_prior)"
                )

            # Validate optional angle prior type
            _valid_angle_priors = ('uniform', 'peaked')
            angle_prior = psi_0_bounds[2] if len(psi_0_bounds) == 3 else 'uniform'
            if angle_prior not in _valid_angle_priors:
                raise ValueError(
                    f"Component '{comp_key}': psi_0_bounds angle_prior must be one of "
                    f"{_valid_angle_priors}, got '{angle_prior}'"
                )
            if angle_prior == 'peaked':
                _pk  = polar_spec.get('peaked_k', PSI_0_PEAKED_K)
                _pa  = polar_spec.get('peaked_angle_deg', 0.0)
                print(f"      psi_0 angle prior: peaked  (k={_pk}, peak_angle={_pa}°, also {_pa + 90.0}°)")
            else:
                print(f"      psi_0 angle prior: {angle_prior}")
            
            # Print confirmation of polar conversion style
            print(f"   {comp_key}: Using polar conversion priors → converting to Cartesian (a, b)")
            print(f"      p0: [{p0_min:.3f}, {p0_max:.3f}] ({distribution})")
            
            # Add note about distribution properties
            if distribution == 'radial_uniform':
                print(f"      Note: Uniform in radius; any p0 is equally likely")
            elif distribution == 'log_uniform':
                print(f"      Note: Log uniform in radius; smaller p0 values are more likely")
            elif distribution == 'area_uniform':
                print(f"      Note: Uniform in area; any (Q, U) is equally likely")
            
            psi_range = psi_0_bounds[1] - psi_0_bounds[0]
            unbounded_str = " (unbounded)" if psi_range >= 0.9 * np.pi else ""
            print(f"      ψ0: [{psi_0_bounds[0]:.2f}, {psi_0_bounds[1]:.2f}] rad{unbounded_str}")
            
            # Compute resulting a, b ranges
            # a = p0 * cos(2*psi_0), b = p0 * sin(2*psi_0)
            if psi_range >= 0.9 * np.pi:
                # Unbounded angle → a, b can reach full ±p0_max
                a_min, a_max = -p0_max, p0_max
                b_min, b_max = -p0_max, p0_max
                print(f"      → a: [{a_min:.4f}, {a_max:.4f}]")
                print(f"      → b: [{b_min:.4f}, {b_max:.4f}]")
            else:
                # Bounded angle → compute actual ranges
                psi_vals = np.linspace(psi_0_bounds[0], psi_0_bounds[1], 1000)
                a_vals_min = p0_min * np.cos(2 * psi_vals)
                a_vals_max = p0_max * np.cos(2 * psi_vals)
                b_vals_min = p0_min * np.sin(2 * psi_vals)
                b_vals_max = p0_max * np.sin(2 * psi_vals)
                
                a_min = min(a_vals_min.min(), a_vals_max.min())
                a_max = max(a_vals_min.max(), a_vals_max.max())
                b_min = min(b_vals_min.min(), b_vals_max.min())
                b_max = max(b_vals_min.max(), b_vals_max.max())
                
                print(f"      → a: [{a_min:.4f}, {a_max:.4f}]")
                print(f"      → b: [{b_min:.4f}, {b_max:.4f}]")
        
        # Validate other component-specific parameters
        if isinstance(comp, ThinPowerLawComponent):
            # beta always required
            if 'beta' not in comp_priors:
                raise ValueError(f"Component '{comp_key}': missing required parameter 'beta'")
            validate_prior_spec(comp_priors['beta'])
            
            # phi_rm: only validate if not tying (if tying, we use 'tied_phi_rm' instead)
            if not tie_phi_rm:
                if 'phi_rm' not in comp_priors:
                    raise ValueError(f"Component '{comp_key}': missing required parameter 'phi_rm'")
                validate_prior_spec(comp_priors['phi_rm'])
            elif 'phi_rm' in comp_priors:
                print(f"WARNING: Component '{comp_key}' phi_rm prior will be IGNORED (using 'tied_phi_rm' instead)")
        
        elif isinstance(comp, ThinComponent):
            # phi_rm: only validate if not tying (if tying, we use 'tied_phi_rm' instead)
            if not tie_phi_rm:
                if 'phi_rm' not in comp_priors:
                    raise ValueError(f"Component '{comp_key}': missing required parameter 'phi_rm'")
                validate_prior_spec(comp_priors['phi_rm'])
            elif 'phi_rm' in comp_priors:
                print(f"WARNING: Component '{comp_key}' phi_rm prior will be IGNORED (using 'tied_phi_rm' instead)")
        
        else:  # ThickComponent
            # phi_peak, sigma_phi, N required
            required = ['phi_peak', 'sigma_phi', 'N']
            for param in required:
                if param not in comp_priors:
                    raise ValueError(
                        f"Component '{comp_key}': missing required parameter '{param}'"
                    )
                validate_prior_spec(comp_priors[param])
            
            # Check for effective_width_mode flag (optional)
            effective_width_mode = comp_priors.get('effective_width_mode', False)
            
            if effective_width_mode:
                # Print confirmation
                sigma_phi_bounds = _get_prior_bounds(comp_priors['sigma_phi'])
                N_bounds = _get_prior_bounds(comp_priors['N'])
                
                print(f"   {comp_key}: Using EFFECTIVE WIDTH MODE")
                print(f"      Input 'sigma_phi': [{sigma_phi_bounds[0]:.1f}, {sigma_phi_bounds[1]:.1f}] rad/m²")
                print(f"      → Interpreted as φ_eff (99% flux containment)")
                print(f"      N: [{N_bounds[0]:.1f}, {N_bounds[1]:.1f}]")
                print(f"      Conversion: σ_φ = φ_eff / (2 ln(100))^(1/N)")
                
                # Show example conversion
                if N_bounds is not None and sigma_phi_bounds is not None:
                    N_example = N_bounds[0]
                    phi_eff_example = sigma_phi_bounds[1]
                    sigma_phi_example = phi_eff_example / (2 * np.log(100))**(1.0/N_example)
                    print(f"      Example: N={N_example:.1f}, φ_eff={phi_eff_example:.1f} → σ_φ={sigma_phi_example:.1f} rad/m²")
    
    # Validate tied_phi_rm if RM tying is enabled
    if tie_phi_rm:
        if 'tied_phi_rm' not in priors:
            raise ValueError(
                "tie_phi_rm=True requires 'tied_phi_rm' prior specification. "
                "Add to priors: 'tied_phi_rm': ('uniform', min, max)"
            )
        
        print("\n" + "="*70)
        print("TIED PHI_RM VALIDATION")
        print("="*70)
        tied_spec = priors['tied_phi_rm']
        validate_prior_spec(tied_spec)
        print("   Tied φ_RM: All thin/power-law components share this prior")
        bounds = _get_prior_bounds(tied_spec)
        if bounds:
            print(f"      φ_RM: [{bounds[0]:.1f}, {bounds[1]:.1f}] rad/m²")
    else:
        # Warn if tied_phi_rm specified but not being used
        if 'tied_phi_rm' in priors:
            print("WARNING: 'tied_phi_rm' prior specified but will be IGNORED "
                  "(tie_phi_rm=False, using per-component phi_rm priors instead)")
    
    print("="*70)
    print("✓ Priors validated successfully")
    print("="*70 + "\n")
    return True


def _get_prior_bounds(prior_spec: Tuple) -> Optional[Tuple[float, float]]:
    """
    Extract bounds from prior specification.
    
    Returns (min, max) for bounded priors, None for unbounded.
    """
    prior_type = prior_spec[0]
    
    if prior_type == 'uniform':
        return (prior_spec[1], prior_spec[2])
    elif prior_type == 'log_uniform':
        return (prior_spec[1], prior_spec[2])
    elif prior_type == 'truncated_gaussian':
        return (prior_spec[3], prior_spec[4])
    else:
        # Gaussian, fixed - no definite bounds
        return None


def _format_prior_dict(priors: Dict) -> str:
    """Format prior dictionary for display."""
    lines = ["priors = {"]
    for comp_key, comp_priors in priors.items():
        lines.append(f"    '{comp_key}': {{")
        for param_key, param_prior in comp_priors.items():
            lines.append(f"        '{param_key}': {param_prior},")
        lines.append("    },")
    lines.append("}")
    return "\n".join(lines)


# =============================================================================
# LIKELIHOOD FUNCTION
# =============================================================================

def log_likelihood(params: np.ndarray,
                   model: CompositeModel,
                   data: PolarizationData,
                   param_names: List[str],
                   p0_physical_max: float = 0.7,
                   lambda_sq_ref: Optional[float] = None) -> float:
    """
    Compute log-likelihood for given parameters.

    ln L = -0.5 * Σ[chi² + log(2π σ²)]

    Physical constraint: p0 = sqrt(a² + b²) must be <= p0_physical_max

    Args:
        params: Parameter values (flattened array, in Cartesian space a, b)
        model: CompositeModel (structure only, params will be updated)
        data: PolarizationData with observations
        param_names: List of parameter names (for debugging)
        p0_physical_max: Maximum allowed polarization fraction (default: 0.7)
        lambda_sq_ref: If set, sampled (a, b) are interpreted as the complex
            amplitude at this reference λ² (m²) and converted to the intrinsic
            amplitude via the per-component transfer factor before evaluating
            the forward model.

    Returns:
        Log-likelihood value (or -inf if invalid parameters)
    """
    try:
        try:
            update_model_from_params(model, params, param_names,
                                     lambda_sq_ref=lambda_sq_ref)
        except Exception:
            return -np.inf

        # Physical constraint: check p0 = sqrt(a² + b²) <= p0_physical_max
        for comp in model.components:
            p0 = np.sqrt(comp.a**2 + comp.b**2)
            if p0 > p0_physical_max:
                return -np.inf

        # Compute model prediction
        P_model = model.compute_polarization(data.lambda_sq)
        Q_model = P_model.real
        U_model = P_model.imag

        if not np.all(np.isfinite(Q_model)) or not np.all(np.isfinite(U_model)):
            return -np.inf

        sigma_Q = data.Q_err
        sigma_U = data.U_err

        chi2_Q = np.sum((data.Q - Q_model)**2 / sigma_Q**2)
        chi2_U = np.sum((data.U - U_model)**2 / sigma_U**2)

        log_norm_Q = np.sum(np.log(2 * np.pi * sigma_Q**2))
        log_norm_U = np.sum(np.log(2 * np.pi * sigma_U**2))

        return -0.5 * (chi2_Q + chi2_U + log_norm_Q + log_norm_U)

    except Exception:
        return -np.inf


def update_model_from_params(model: CompositeModel,
                             params: np.ndarray,
                             param_names: List[str],
                             lambda_sq_ref: Optional[float] = None):
    """
    Update model components with parameter values.
    
    Handles both indexed names ('component_0_a') and named components ('ISM_a').
    Parameters are always in Cartesian space (a, b).
    
    Args:
        model: CompositeModel to update (modified in place)
        params: Parameter values (Cartesian space)
        param_names: Parameter names (e.g., 'ISM_a' or 'component_0_a')
        lambda_sq_ref: If set, sampled (a, b) for each component are treated as
            the complex amplitude at this reference λ² (m²) and converted to the
            intrinsic amplitude via the component's transfer factor.
    """
    # Build a mapping from component name/index to component object
    comp_map = {}
    comp_name_list = []
    for i, comp in enumerate(model.components):
        comp_name = comp.name if comp.name else f'component_{i}'
        comp_map[comp_name] = comp
        comp_map[f'component_{i}'] = comp  # Also support index-based lookup
        comp_name_list.append(comp_name)
    
    # Parse parameter names and group by component
    param_dict = {}
    
    # Extract tied_phi_rm if present
    tied_phi_rm_value = None
    if 'tied_phi_rm' in param_names:
        idx = param_names.index('tied_phi_rm')
        tied_phi_rm_value = params[idx]
    
    # Component parameter suffixes (Cartesian form)
    param_suffixes = ['_a', '_b', '_phi_rm', '_phi_peak', '_sigma_phi', '_N', '_beta']
    
    non_component_params = ['tied_phi_rm']

    for name, value in zip(param_names, params):
        if name in non_component_params:
            continue
        
        # Find which parameter suffix matches
        comp_identifier = None
        param_name = None
        
        for suffix in param_suffixes:
            if name.endswith(suffix):
                comp_identifier = name[:-len(suffix)]
                param_name = suffix[1:]  # Remove leading underscore
                break
        
        if comp_identifier is None or param_name is None:
            raise ValueError(f"Could not parse parameter name: {name}")
        
        if comp_identifier not in param_dict:
            param_dict[comp_identifier] = {}
        param_dict[comp_identifier][param_name] = value
    
    # Update each component
    for comp_identifier, comp_params in param_dict.items():
        if comp_identifier not in comp_map:
            raise ValueError(
                f"Component '{comp_identifier}' not found in model. "
                f"Available components: {comp_name_list}"
            )
        
        comp = comp_map[comp_identifier]
        
        if isinstance(comp, ThinPowerLawComponent):
            comp.beta = comp_params['beta']
            if tied_phi_rm_value is not None:
                comp.phi_rm = tied_phi_rm_value
            else:
                comp.phi_rm = comp_params['phi_rm']
            comp.a = comp_params['a']
            comp.b = comp_params['b']
            if lambda_sq_ref is not None:
                K = comp.compute_transfer_factor(lambda_sq_ref)
                A_int = (comp_params['a'] + 1j * comp_params['b']) / K
                comp.a = A_int.real
                comp.b = A_int.imag

        elif isinstance(comp, ThinComponent):
            if tied_phi_rm_value is not None:
                comp.phi_rm = tied_phi_rm_value
            else:
                comp.phi_rm = comp_params['phi_rm']
            comp.a = comp_params['a']
            comp.b = comp_params['b']
            if lambda_sq_ref is not None:
                K = comp.compute_transfer_factor(lambda_sq_ref)
                A_int = (comp_params['a'] + 1j * comp_params['b']) / K
                comp.a = A_int.real
                comp.b = A_int.imag

        elif isinstance(comp, ThickComponent):
            comp.phi_peak = comp_params['phi_peak']
            comp.sigma_phi = comp_params['sigma_phi']
            comp.N = comp_params['N']
            comp.a = comp_params['a']
            comp.b = comp_params['b']
            if lambda_sq_ref is not None:
                K = comp.compute_transfer_factor(lambda_sq_ref)
                A_int = (comp_params['a'] + 1j * comp_params['b']) / K
                comp.a = A_int.real
                comp.b = A_int.imag


# =============================================================================
# PRIOR TRANSFORM WITH ORDERING
# =============================================================================

def ordered_phi_uniform(u_block: np.ndarray, phi_min: float, phi_max: float) -> np.ndarray:
    """
    Generate ordered parameters from uniform prior using sorting.
    
    Transforms K i.i.d. uniform samples to K ordered values preserving
    the uniform prior marginals (order-statistics distribution).
    
    Args:
        u_block: Array of K values in [0,1]
        phi_min: Lower bound of uniform prior
        phi_max: Upper bound of uniform prior
        
    Returns:
        Array of K ordered values: phi[0] <= phi[1] <= ... <= phi[K-1]
    """
    v = np.sort(u_block)
    return phi_min + (phi_max - phi_min) * v


def ordered_phi_gaussian(u_block: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """
    Generate ordered parameters from Gaussian prior using sorting + inverse CDF.
    
    Transforms K i.i.d. uniform samples to K ordered values preserving
    the Gaussian prior marginals (order-statistics distribution).
    
    Uses stable erfinv implementation for normal PPF.
    
    Args:
        u_block: Array of K values in [0,1]
        mu: Mean of Gaussian prior
        sigma: Standard deviation of Gaussian prior
        
    Returns:
        Array of K ordered values: phi[0] <= phi[1] <= ... <= phi[K-1]
    """
    from scipy.special import erfinv
    
    # Clip to avoid +/- inf at exactly 0 or 1
    eps = 1e-15
    v = np.sort(np.clip(u_block, eps, 1.0 - eps))
    
    # Stable normal PPF via erfinv
    # norm_ppf(p) = sqrt(2) * erfinv(2*p - 1)
    z = np.sqrt(2.0) * erfinv(2.0 * v - 1.0)
    
    return mu + sigma * z


def ordered_phi_log_uniform(u_block: np.ndarray, pmin: float, pmax: float) -> np.ndarray:
    """
    Generate ordered parameters from log-uniform prior using sorting + inverse CDF.
    
    Transforms K i.i.d. uniform samples to K ordered values preserving
    the log-uniform prior marginals (order-statistics distribution).
    
    Args:
        u_block: Array of K values in [0,1]
        pmin: Lower bound (must be > 0)
        pmax: Upper bound (must be > pmin)
        
    Returns:
        Array of K ordered values: phi[0] <= phi[1] <= ... <= phi[K-1]
    """
    v = np.sort(u_block)
    log_pmin = np.log(pmin)
    log_pmax = np.log(pmax)
    return np.exp(log_pmin + v * (log_pmax - log_pmin))


def ordered_phi_truncated_gaussian(u_block: np.ndarray, mu: float, sigma: float, 
                                   pmin: float, pmax: float) -> np.ndarray:
    """
    Generate ordered parameters from truncated Gaussian prior using sorting + inverse CDF.
    
    Transforms K i.i.d. uniform samples to K ordered values preserving
    the truncated Gaussian prior marginals (order-statistics distribution).
    
    Args:
        u_block: Array of K values in [0,1]
        mu: Mean of underlying Gaussian
        sigma: Std dev of underlying Gaussian
        pmin: Lower truncation bound
        pmax: Upper truncation bound
        
    Returns:
        Array of K ordered values: phi[0] <= phi[1] <= ... <= phi[K-1]
    """
    # Clip to avoid edge case issues at exactly 0 or 1
    eps = 1e-15
    u_clipped = np.clip(u_block, eps, 1.0 - eps)
    v = np.sort(u_clipped)
    a = (pmin - mu) / sigma
    b = (pmax - mu) / sigma
    return stats.truncnorm.ppf(v, a, b, loc=mu, scale=sigma)


def _analyze_ordering_requirements(model: CompositeModel, priors: Dict, 
                                   tie_phi_rm: bool = False) -> Dict:
    """
    Analyze which parameter groups need ordering and how to implement it.
    
    Groups parameters by (component_type, param_name) to ensure ordering only
    applies within exchangeable families. Returns single source of truth for
    all ordering logic.
    
    Args:
        model: CompositeModel
        priors: Prior dictionary
        tie_phi_rm: If True, phi_rm is tied (use beta for ThinPowerLaw ordering)
        
    Returns:
        Dictionary with ordering configuration keyed by (component_type, param_name):
        {
            (ThinComponent, 'phi_rm'): {'needs_ordering': True, 'method': 'uniform', ...},
            (ThinPowerLawComponent, 'beta'): {'needs_ordering': True, 'method': 'gaussian', ...},
            ...
        }
    """
    # Collect priors grouped by (component_type, param_name)
    param_groups = {}  # {(comp_type, param_name): [prior1, prior2, ...]}
    
    for i, comp in enumerate(model.components):
        comp_key = comp.name if comp.name else f'component_{i}'
        comp_type = type(comp)
        
        if comp_key in priors:
            comp_priors = priors[comp_key]
        elif f'component_{i}' in priors:
            comp_priors = priors[f'component_{i}']
        else:
            continue
        
        # Determine orderable parameter for this component type
        # B) Beta ordering conditional on tie_phi_rm
        if comp_type == ThinComponent:
            orderable = ['phi_rm'] if not tie_phi_rm else []
        elif comp_type == ThinPowerLawComponent:
            orderable = ['phi_rm'] if not tie_phi_rm else ['beta']
        elif comp_type == ThickComponent:
            orderable = ['phi_peak']
        else:
            orderable = []
        
        for param_name in orderable:
            if param_name in comp_priors:
                group_key = (comp_type, param_name)
                if group_key not in param_groups:
                    param_groups[group_key] = []
                param_groups[group_key].append(comp_priors[param_name])
    
    # Analyze each group
    config = {}
    
    for (comp_type, param_name), prior_list in param_groups.items():
        if len(prior_list) <= 1:
            config[(comp_type, param_name)] = {
                'needs_ordering': False, 'method': None, 'params': None
            }
        else:
            config[(comp_type, param_name)] = {
                'needs_ordering': True, 'method': None, 'params': None
            }
            
            # Check if all identical
            first = prior_list[0]
            if not all(p == first for p in prior_list):
                config[(comp_type, param_name)]['method'] = 'likelihood'
                continue
            
            prior_type, params = validate_prior_spec(first)
            
            # C) Extended order-statistics support
            if prior_type == 'uniform':
                config[(comp_type, param_name)]['method'] = 'uniform'
                config[(comp_type, param_name)]['params'] = {
                    'phi_min': params[0], 'phi_max': params[1]
                }
            elif prior_type == 'gaussian':
                config[(comp_type, param_name)]['method'] = 'gaussian'
                config[(comp_type, param_name)]['params'] = {
                    'mu': params[0], 'sigma': params[1]
                }
            elif prior_type == 'log_uniform':
                config[(comp_type, param_name)]['method'] = 'log_uniform'
                config[(comp_type, param_name)]['params'] = {
                    'pmin': params[0], 'pmax': params[1]
                }
            elif prior_type == 'truncated_gaussian':
                config[(comp_type, param_name)]['method'] = 'truncated_gaussian'
                config[(comp_type, param_name)]['params'] = {
                    'mu': params[0], 'sigma': params[1],
                    'pmin': params[2], 'pmax': params[3]
                }
            else:
                config[(comp_type, param_name)]['method'] = 'likelihood'
    
    return config


def _build_peaked_angle_ppf(k, peak_angle=0.0, n_grid=100_000):
    """
    Build numerical inverse CDF for the peaked polarization angle prior.

    Prior form: f(theta) proportional to exp(k * cos(4 * (theta - peak_angle)))
      - Peaks at peak_angle and peak_angle ± 90 deg
      - Troughs at peak_angle ± 45 deg
      - k=0 recovers uniform; larger k concentrates more mass at peaks
      - peak_angle=0 (default) reproduces the original behaviour

    Grid is built once at initialisation; each sample is a fast interpolation lookup.

    Args:
        k:           Concentration parameter (see PSI_0_PEAKED_K)
        peak_angle:  Angle in radians where the prior peaks (default 0.0).
                     The prior also peaks at peak_angle ± pi/2.
        n_grid:      Number of grid points for numerical integration (default 100,000)

    Returns:
        Callable mapping u in [0, 1] -> theta in radians over [-pi/2, pi/2)
    """
    from scipy.interpolate import interp1d

    theta = np.linspace(-np.pi / 2, np.pi / 2, n_grid, endpoint=False)
    log_p = k * np.cos(4 * (theta - peak_angle))
    p = np.exp(log_p - log_p.max())
    cdf = np.cumsum(p)
    cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    return interp1d(cdf, theta, bounds_error=False,
                    fill_value=(theta[0], theta[-1]))


class PriorTransform:
    """
    Transform unit cube to physical parameters.
    
    Handles fixed parameters by sampling only free parameters from dynesty,
    then inserting fixed values at correct positions before returning.
    Picklable for multiprocessing.
    """
    
    def __init__(self, model: CompositeModel, priors: Dict, 
                 enforce_ordering: bool = True, tie_phi_rm: bool = False):
        """
        Initialize prior transform.
        
        Args:
            model: CompositeModel structure
            priors: Prior dictionary (may include fixed parameters)
            enforce_ordering: If True, enforce component ordering
            tie_phi_rm: If True, tie phi_rm across all thin/power-law components
        """
        # Build full parameter info
        self.param_info_full = _build_param_info(model, priors, tie_phi_rm=tie_phi_rm)
        
        # Separate free and fixed parameters
        self.param_info_free = [p for p in self.param_info_full if not p['is_fixed']]
        
        # Mapping from sampled space (free params) to full space (all params)
        self.free_to_full = [i for i, p in enumerate(self.param_info_full) 
                             if not p['is_fixed']]
        
        # Store fixed parameter values
        self.fixed_params = {i: p['fixed_value'] 
                            for i, p in enumerate(self.param_info_full) 
                            if p['is_fixed']}
        
        # Detect periodic boundaries for dynesty
        self.periodic = []
        for i, info in enumerate(self.param_info_free):
            if info.get('is_polar_conversion', False) and info['param_name'] == 'b':
                polar_spec = info['prior_spec']
                psi_0_bounds = polar_spec['psi_0_bounds']
                psi_0_min, psi_0_max = psi_0_bounds[0], psi_0_bounds[1]
                if abs(psi_0_max - psi_0_min) >= 0.9 * np.pi:
                    self.periodic.append(i)

        self._psi_ppf = {}
        for i, info in enumerate(self.param_info_free):
            if info.get('is_polar_conversion', False) and info['param_name'] == 'b':
                polar_spec = info['prior_spec']
                psi_0_bounds = polar_spec['psi_0_bounds']
                angle_prior = psi_0_bounds[2] if len(psi_0_bounds) == 3 else 'uniform'
                if angle_prior == 'peaked':
                    peaked_k     = polar_spec.get('peaked_k', PSI_0_PEAKED_K)
                    peaked_angle = np.deg2rad(polar_spec.get('peaked_angle_deg', 0.0))
                    self._psi_ppf[i] = _build_peaked_angle_ppf(peaked_k, peak_angle=peaked_angle)

        self.enforce_ordering = enforce_ordering
        self.model = model
        self.tie_phi_rm = tie_phi_rm
        
        # Analyze ordering requirements
        if enforce_ordering:
            self.ordering_config = _analyze_ordering_requirements(model, priors, tie_phi_rm)
            
            # Check if ANY parameter group needs likelihood enforcement
            self.needs_likelihood_ordering = any(
                cfg['needs_ordering'] and cfg['method'] == 'likelihood'
                for cfg in self.ordering_config.values()
            )
        else:
            self.ordering_config = {}
            self.needs_likelihood_ordering = False
        
        # Build index maps for ordered parameters (for efficient lookup in __call__)
        self._build_ordering_indices()
    
    def _build_ordering_indices(self):
        """
        Build index maps for parameters that need ordering.
        
        Groups by (component_type, param_name) to ensure ordering only within
        exchangeable families. Creates self.ordering_indices = {(type, param): [indices]}
        """
        self.ordering_indices = {}
        
        for group_key in self.ordering_config.keys():
            self.ordering_indices[group_key] = []
        
        # Scan through free parameters and collect indices
        for i, info in enumerate(self.param_info_free):
            comp_idx = info['component_idx']
            param_name = info['param_name']
            
            # Skip non-component parameters (like tied_phi_rm)
            if comp_idx is None:
                continue
            
            # Determine component type
            comp = self.model.components[comp_idx]
            comp_type = type(comp)
            
            # Check if this (comp_type, param_name) pair is in our ordering config
            group_key = (comp_type, param_name)
            if group_key in self.ordering_indices:
                self.ordering_indices[group_key].append(i)
    
    def __call__(self, u: np.ndarray) -> np.ndarray:
        """
        Transform from unit cube to physical parameters.
        
        Handles:
        - Cartesian (a, b) priors
        - Polar conversion (p0, psi_0) → (a, b)
        - Order-statistics transforms for phi_rm and phi_peak (when applicable)
        
        Args:
            u: Array in [0,1]^n_free (free parameters only)
            
        Returns:
            Free physical parameters (n_free dimensions)
            Full parameter vector is constructed internally for ordering check.
        """
        params_free = np.zeros(len(self.param_info_free))
        
        # Track which parameters have been transformed
        transformed = np.zeros(len(self.param_info_free), dtype=bool)
        
        # STEP -1: Apply order-statistics transforms if applicable
        # Loop over all parameter groups that use order-statistics
        for group_key, indices in self.ordering_indices.items():
            config = self.ordering_config.get(group_key, {})
            
            if (config.get('needs_ordering', False) and 
                config.get('method') in ['uniform', 'gaussian', 'log_uniform', 'truncated_gaussian'] and
                len(indices) > 0):
                
                u_block = u[indices]
                method = config['method']
                params = config['params']
                
                if method == 'uniform':
                    ordered_values = ordered_phi_uniform(u_block, params['phi_min'], params['phi_max'])
                elif method == 'gaussian':
                    ordered_values = ordered_phi_gaussian(u_block, params['mu'], params['sigma'])
                elif method == 'log_uniform':
                    ordered_values = ordered_phi_log_uniform(u_block, params['pmin'], params['pmax'])
                elif method == 'truncated_gaussian':
                    ordered_values = ordered_phi_truncated_gaussian(
                        u_block, params['mu'], params['sigma'], params['pmin'], params['pmax'])
                
                # Insert ordered values
                for i, idx in enumerate(indices):
                    params_free[idx] = ordered_values[i]
                    transformed[idx] = True
        
        # STEP 0: Handle polar_conversion pairs (a, b)
        # Must process these together before other transformations
        i = 0
        while i < len(self.param_info_free):
            info = self.param_info_free[i]
            
            if info.get('is_polar_conversion', False) and info['param_name'] == 'a':
                # Found 'a' with polar_conversion, next should be 'b'
                if i + 1 >= len(self.param_info_free):
                    raise ValueError(f"Parameter '{info['name']}' has polar_conversion but no 'b' follows")
                
                info_b = self.param_info_free[i + 1]
                if not (info_b.get('is_polar_conversion', False) and info_b['param_name'] == 'b'):
                    raise ValueError(f"Parameter '{info['name']}' (a) must be followed by 'b'")
                
                # Both 'a' and 'b' share the same polar_conversion spec
                polar_spec = info['prior_spec']
                p0_bounds = polar_spec['p0_bounds']
                psi_0_bounds = polar_spec['psi_0_bounds']
                
                p0_min, p0_max = p0_bounds[0], p0_bounds[1]
                distribution = p0_bounds[2] if len(p0_bounds) == 3 else 'radial_uniform'
                psi_0_min, psi_0_max = psi_0_bounds[0], psi_0_bounds[1]
                
                # Sample p0 based on distribution type
                if distribution == 'radial_uniform':
                    # Uniform in radius; any p0 is equally likely
                    p0 = p0_min + (p0_max - p0_min) * u[i]
                    
                elif distribution == 'log_uniform':
                    # Log uniform in radius; smaller p0 values are more likely
                    log_p0 = np.log(p0_min) + u[i] * (np.log(p0_max) - np.log(p0_min))
                    p0 = np.exp(log_p0)
                    
                elif distribution == 'area_uniform':
                    # Uniform in area; any (Q, U) is equally likely
                    # Jacobian: |∂(Q,U)/∂(p0,ψ)| = 2p0
                    # CDF: F(p0) = (p0² - p_min²) / (p_max² - p_min²)
                    # Inverse CDF: p0 = sqrt(u * (p_max² - p_min²) + p_min²)
                    p0_sq = u[i] * (p0_max**2 - p0_min**2) + p0_min**2
                    p0 = np.sqrt(p0_sq)
                    
                else:
                    raise ValueError(f"Unknown distribution: {distribution}")
                
                # Sample psi_0 (check for unbounded angle)
                if abs(psi_0_max - psi_0_min) >= 0.9 * np.pi:
                    angle_prior = psi_0_bounds[2] if len(psi_0_bounds) == 3 else 'uniform'
                    if angle_prior == 'peaked' and (i + 1) in self._psi_ppf:
                        psi_0 = float(self._psi_ppf[i + 1](u[i + 1]))
                    else:
                        psi_0 = np.pi * (u[i + 1] - 0.5)
                else:
                    # Bounded wedge
                    psi_0 = psi_0_min + (psi_0_max - psi_0_min) * u[i + 1]
                
                # Convert to Cartesian
                a = p0 * np.cos(2 * psi_0)
                b = p0 * np.sin(2 * psi_0)
                
                params_free[i] = a
                params_free[i + 1] = b
                transformed[i] = True
                transformed[i + 1] = True
                
                i += 2  # Skip the 'b' we just processed
            else:
                i += 1
        
        # STEP 1: Transform sigma_phi parameters (with effective width conversion if needed)
        # If effective_width_mode is enabled, sigma_phi is interpreted as φ_eff and
        # converted to σ_phi after N is known.
        sigma_phi_deferred = []
        for i, info in enumerate(self.param_info_free):
            if not transformed[i] and info['name'].endswith('_sigma_phi'):
                if info.get('effective_width_mode', False):
                    # Using effective width mode - defer conversion until N is sampled
                    sigma_phi_deferred.append(i)
                else:
                    # Standard sigma_phi
                    params_free[i] = _transform_parameter(u[i], info['prior_spec'])
                    transformed[i] = True
        
        # STEP 2: Transform all remaining parameters (including N)
        for i, info in enumerate(self.param_info_free):
            if not transformed[i]:
                params_free[i] = _transform_parameter(u[i], info['prior_spec'])
                transformed[i] = True
        
        # STEP 3: Convert φ_eff → σ_phi after N has been set
        for sigma_idx in sigma_phi_deferred:
            sigma_info = self.param_info_free[sigma_idx]
            
            # Find corresponding N parameter for this component in FULL space
            comp_idx = sigma_info['component_idx']
            N_value = None
            
            # First check if N is in free parameters
            for j, info in enumerate(self.param_info_free):
                if info['component_idx'] == comp_idx and info['param_name'] == 'N':
                    N_value = params_free[j]
                    break
            
            # If not found in free params, check if it's fixed
            if N_value is None:
                for full_idx, info in enumerate(self.param_info_full):
                    if info['component_idx'] == comp_idx and info['param_name'] == 'N':
                        if info['is_fixed']:
                            N_value = self.fixed_params[full_idx]
                        break
            
            if N_value is None:
                raise ValueError(
                    f"Could not find N parameter for component {comp_idx}. "
                    f"effective_width_mode requires N to be either free or fixed."
                )
            
            # Sample phi_eff from prior (using sigma_phi prior specification)
            phi_eff = _transform_parameter(u[sigma_idx], sigma_info['prior_spec'])
            
            # Convert phi_eff to sigma_phi
            # phi_eff is defined at 99% flux containment (1% of peak)
            # f(phi_eff) = exp(-|phi_eff|^N / (2*sigma_phi^N)) = 0.01
            # Solving: sigma_phi = phi_eff / (2 * ln(100))^(1/N)
            sigma_phi = phi_eff / (2.0 * np.log(100.0))**(1.0 / N_value)
            
            params_free[sigma_idx] = sigma_phi
            transformed[sigma_idx] = True
        
        # Expand to full parameter space for ordering check
        params_full = np.zeros(len(self.param_info_full))
        
        for i, full_idx in enumerate(self.free_to_full):
            params_full[full_idx] = params_free[i]
        
        for full_idx, value in self.fixed_params.items():
            params_full[full_idx] = value
        
        # NOTE: Ordering constraints using order-statistics are already satisfied by construction.
        # Only groups with 'likelihood' enforcement method need checking in log_likelihood.
        # prior_transform should always return valid finite values.
        # Returning NaNs can corrupt dynesty's internal geometry.
        
        # Return only free parameters (dynesty expects n_free dimensions)
        return params_free


def create_prior_transform(model: CompositeModel, 
                          priors: Dict,
                          enforce_ordering: bool = True,
                          tie_phi_rm: bool = False):
    """
    Create prior transform function for dynesty.
    
    Transforms from unit cube [0,1]^N to physical parameters,
    with optional ordering constraints and RM tying.
    
    Args:
        model: CompositeModel structure
        priors: Prior dictionary
        enforce_ordering: If True, enforce component ordering
        tie_phi_rm: If True, tie phi_rm across all thin/power-law components
        
    Returns:
        PriorTransform callable (picklable for multiprocessing)
    """
    return PriorTransform(model, priors, enforce_ordering, tie_phi_rm)


# =============================================================================
# FITSETUP CLASS
# =============================================================================

class FitSetup:
    """
    Container for fitting setup: model + priors + data.
    
    This class packages everything needed for Bayesian fitting:
    - CompositeModel (defines structure)
    - Prior dictionary (defines parameter search space)
    - PolarizationData (observations to fit)
    
    It validates that everything is consistent and provides
    convenient access to fitting components.
    """
    
    def __init__(self,
                 model: CompositeModel,
                 priors: Dict,
                 data: PolarizationData,
                 enforce_ordering: bool = True,
                 p0_physical_max: float = 0.7,
                 tie_phi_rm: bool = False,
                 lambda_sq_ref: Optional[float] = None):
        """
        Initialize fit setup.
        
        Args:
            model: CompositeModel to fit
            priors: Prior dictionary (can use component names or indices)
            data: PolarizationData to fit
            enforce_ordering: If True, enforce component ordering constraints
            p0_physical_max: Maximum allowed polarization fraction (default: 0.7)
            tie_phi_rm: If True, tie phi_rm across all thin/power-law components
                       (all share same Rotation Measure). Default: False
            
        Raises:
            ValueError: If validation fails
        """
        # Validate p0_physical_max
        if p0_physical_max <= 0:
            raise ValueError(f"p0_physical_max must be > 0, got {p0_physical_max}")
        
        self.model = model
        self.data = data
        self.enforce_ordering = enforce_ordering
        self.p0_physical_max = p0_physical_max
        self.tie_phi_rm = tie_phi_rm
        self.lambda_sq_ref = lambda_sq_ref
        if lambda_sq_ref is not None:
            print(f"Reference-band parameterisation active: λ²_ref = {lambda_sq_ref:.4f} m²")

        # Set lambda_sq_0 for ThinPowerLawComponent instances if not already set
        # Use variance-weighted mean of the data
        self._setup_power_law_reference()
        
        # Validate RM tying if enabled
        if self.tie_phi_rm:
            self._validate_rm_tying()
        
        # Normalize prior keys to use component names if available
        self.priors = self._normalize_priors(priors)
        
        # Validate everything
        self.validate()
        
        # Pre-compute useful quantities
        self.ndim = get_ndim(model, tie_phi_rm=tie_phi_rm)  # Will be overridden if there are fixed params
        self.component_names = self._get_component_names()
        
        # Create prior transform first (needed to build param_names)
        self.prior_transform = create_prior_transform(
            model, self.priors, enforce_ordering=enforce_ordering, tie_phi_rm=tie_phi_rm
        )
        
        self.param_names = [p['name'] for p in self.prior_transform.param_info_full]
        
        # Parse fixed vs free parameters (updates self.ndim, needs param_names)
        self._parse_fixed_params()
        
        # Report fixed parameters if any
        if len(self.fixed_params) > 0:
            print("\n" + "="*70)
            print("FIXED PARAMETERS")
            print("="*70)
            for idx, value in self.fixed_params.items():
                param_name = self.param_names[idx]
                print(f"  {param_name:20s} = {value:.6g}")
            print(f"\nTotal: {len(self.fixed_params)} fixed, {self.ndim} free parameters")
            print("="*70)
        
        # Report periodic boundaries if any
        if len(self.prior_transform.periodic) > 0:
            print("\n" + "="*70)
            print("PERIODIC BOUNDARIES")
            print("="*70)
            for idx in self.prior_transform.periodic:
                # Get the parameter info to show which component this belongs to
                param_info = self.prior_transform.param_info_free[idx]
                comp_idx = param_info['component_idx']
                # Get component name from model
                comp = self.model.components[comp_idx]
                comp_name = comp.name if comp.name else f'component_{comp_idx}'
                print(f"  Unit cube index {idx} (samples psi_0 for {comp_name})")
            print(f"\nTotal: {len(self.prior_transform.periodic)} periodic dimensions")
            print("="*70)
        
        # Report ordering strategy
        if enforce_ordering:
            ordering_config = self.prior_transform.ordering_config
            
            # Check if any ordering is needed
            any_ordering = any(cfg['needs_ordering'] for cfg in ordering_config.values())
            
            if any_ordering:
                print("\n" + "="*70)
                print("ORDERING STRATEGY")
                print("="*70)
                
                for group_key, config in ordering_config.items():
                    if config['needs_ordering']:
                        comp_type, param_name = group_key
                        method = config['method']
                        
                        # Format display name
                        type_name = comp_type.__name__
                        display = f"  {type_name}.{param_name}"
                        
                        if method in ['uniform', 'gaussian', 'log_uniform', 'truncated_gaussian']:
                            print(f"{display}: Order-statistics prior applied")
                            print(f"    Method: {method} (mathematically exact)")
                            print(f"    Location: prior transform (efficient)")
                        else:
                            print(f"{display}: Likelihood enforcement")
                            print(f"    Ordering will be enforced via likelihood because priors differ across")
                            print(f"    components or are not supported for order-statistics")
                            print(f"    Location: likelihood (less efficient, more complex geometry)")
                
                print("="*70)
    
    def _setup_power_law_reference(self):
        """
        Set lambda_sq_0 for ThinPowerLawComponent instances if not already set.

        Uses lambda_sq_ref when provided, ensuring the spectral normalisation
        point matches the amplitude reference point so the transfer factor
        reduces to a pure phase for power-law components. Falls back to the
        variance-weighted mean of the data otherwise.
        """
        if self.lambda_sq_ref is not None:
            lambda_sq_0 = self.lambda_sq_ref
        else:
            lambda_sq_0 = 0.05

        for comp in self.model.components:
            if isinstance(comp, ThinPowerLawComponent):
                if comp.lambda_sq_0 is None:
                    comp.lambda_sq_0 = lambda_sq_0
                    comp_name = comp.name if comp.name else "unnamed"
                    print(f"Set λ²₀ = {lambda_sq_0:.6f} m² for ThinPowerLawComponent '{comp_name}'")

    def _validate_rm_tying(self):
        """
        Validate that RM tying is appropriate and count tieable components.
        
        Checks:
        - At least one thin/power-law component exists
        - If enforce_ordering is also True, ALL components must be ThinPowerLawComponent
        - Warns if thick components are present (tying may not be physically meaningful)
        - Warns if per-source phi_rm priors are specified (they'll be ignored)
        """
        # Count component types
        thin_components = []
        powerlaw_components = []
        thick_components = []
        
        for i, comp in enumerate(self.model.components):
            comp_name = comp.name if comp.name else f'component_{i}'
            if isinstance(comp, ThinPowerLawComponent):
                powerlaw_components.append((i, comp_name, comp))
            elif isinstance(comp, ThinComponent):
                thin_components.append((i, comp_name, comp))
            elif isinstance(comp, ThickComponent):
                thick_components.append((i, comp_name, comp))
        
        n_tieable = len(thin_components) + len(powerlaw_components)
        n_powerlaw = len(powerlaw_components)
        n_thick = len(thick_components)
        
        # Check that we have components to tie
        if n_tieable == 0:
            raise ValueError(
                "tie_phi_rm=True requires at least one ThinComponent or ThinPowerLawComponent. "
                f"Model has {n_thick} ThickComponent(s) only."
            )
        
        # CRITICAL: If enforce_ordering AND tie_phi_rm, all components must be ThinPowerLawComponent
        # This allows beta ordering (spectral index) which only exists for power-law components
        if self.enforce_ordering:
            n_total = len(self.model.components)
            if n_powerlaw != n_total:
                raise ValueError(
                    "tie_phi_rm=True AND enforce_ordering=True requires ALL components "
                    f"to be ThinPowerLawComponent for beta ordering. "
                    f"Found: {n_powerlaw} ThinPowerLawComponent, {len(thin_components)} ThinComponent, "
                    f"{n_thick} ThickComponent. Convert all to ThinPowerLawComponent."
                )
            print(f"\nBeta ordering enabled: {n_powerlaw} ThinPowerLawComponent(s) with tied phi_rm")
        
        if n_tieable == 1:
            print(f"WARNING: tie_phi_rm=True but only 1 thin/power-law component exists. "
                  f"Tying has no effect (nothing to tie to). Consider setting tie_phi_rm=False.")
        
        # Warn about thick components (only relevant when not enforcing ordering)
        if n_thick > 0 and not self.enforce_ordering:
            print(f"WARNING: tie_phi_rm=True with {n_thick} ThickComponent(s) present: "
                  f"Tying phi_rm to phi_peak may not be physically meaningful. "
                  f"Thick components represent extended Faraday depth distributions, "
                  f"while thin components represent discrete Faraday depths. "
                  f"Consider using separate RMs for thick components.")
        
        # Print confirmation
        if not self.enforce_ordering:
            print(f"\nRM Tying enabled: phi_rm will be tied across {n_tieable} thin/power-law component(s)")
            for i, comp_name, comp in (thin_components + powerlaw_components):
                comp_type = "ThinPowerLawComponent" if isinstance(comp, ThinPowerLawComponent) else "ThinComponent"
                print(f"  - {comp_name} ({comp_type})")
    
    def _normalize_priors(self, priors: Dict) -> Dict:
        """
        Normalize prior keys to use component names.
        
        Supports both:
        - Named keys: {'ISM': {...}, 'Halo': {...}}
        - Indexed keys: {'component_0': {...}, 'component_1': {...}}
        - Special keys: {'tied_phi_rm': (...)}
        
        Args:
            priors: Input prior dictionary
            
        Returns:
            Normalized prior dictionary with component names/indices as keys
        """
        normalized = {}
        
        for i, comp in enumerate(self.model.components):
            # Try component name first
            comp_name = comp.name if comp.name else f'component_{i}'
            
            # Check if prior exists with this name
            if comp_name in priors:
                normalized[comp_name] = priors[comp_name]
            # Fall back to indexed name
            elif f'component_{i}' in priors:
                normalized[comp_name] = priors[f'component_{i}']
            else:
                raise ValueError(
                    f"No prior found for component {i} (name='{comp.name}'). "
                    f"Expected key '{comp_name}' or 'component_{i}'"
                )
        
        # Copy over non-component keys (like 'tied_phi_rm')
        component_keys = set()
        for i, comp in enumerate(self.model.components):
            comp_name = comp.name if comp.name else f'component_{i}'
            component_keys.add(comp_name)
            component_keys.add(f'component_{i}')
        
        for key, value in priors.items():
            if key not in component_keys:
                normalized[key] = value
        
        return normalized
    
    def _get_component_names(self) -> List[str]:
        """Get list of component names (or indices if unnamed)."""
        names = []
        for i, comp in enumerate(self.model.components):
            names.append(comp.name if comp.name else f'component_{i}')
        return names
    
    def _parse_fixed_params(self):
        """
        Separate fixed and free parameters.
        
        Builds mappings between sampled space (free params only) and full 
        parameter space. Updates self.ndim to reflect only free parameters.
        
        Fixed parameters are specified as ('fixed', value) in the prior dict.
        """
        # Build parameter info to identify fixed vs free
        param_info = _build_param_info(self.model, self.priors, tie_phi_rm=self.tie_phi_rm)
        
        n_params = len(self.param_names)
        self.fixed_params = {}  # {param_idx: fixed_value}
        self.free_param_indices = []
        self.is_fixed = np.zeros(n_params, dtype=bool)
        
        for i, pinfo in enumerate(param_info):
            if pinfo['is_fixed']:
                self.fixed_params[i] = pinfo['fixed_value']
                self.is_fixed[i] = True
            else:
                self.free_param_indices.append(i)
        
        # Ensure at least one parameter is free
        if len(self.free_param_indices) == 0:
            raise ValueError("All parameters are fixed - nothing to fit!")
        
        # Validate fixed RM ordering if enforce_ordering is True
        if self.enforce_ordering:
            self._validate_fixed_rm_ordering(param_info)
        
        # Convert fixed sigma_phi from phi_eff to sigma_phi in effective_width_mode
        self._convert_fixed_sigma_phi_in_effective_width_mode(param_info)
        
        # Update dimensionality (dynesty samples only free parameters)
        self.ndim = len(self.free_param_indices)
        self.free_param_names = [self.param_names[i] for i in self.free_param_indices]
    
    def _validate_fixed_rm_ordering(self, param_info):
        """
        Check that fixed parameter values respect ordering constraints per group.
        
        Only validates groups that actually enforce ordering. Uses non-decreasing
        constraint (<=) to match runtime _check_ordering behavior. Ignores tied_phi_rm
        (component_idx=None) which has no component index.
        
        Args:
            param_info: List of parameter info dicts from _build_param_info
            
        Raises:
            ValueError: If fixed values violate ordering constraint
        """
        # Get ordering groups from prior transform config
        if not hasattr(self, 'prior_transform') or not hasattr(self.prior_transform, 'ordering_config'):
            return  # Can't validate without ordering config
        
        ordering_config = self.prior_transform.ordering_config
        
        # For each group that needs ordering, collect fixed values
        for group_key, config in ordering_config.items():
            if not config.get('needs_ordering', False):
                continue
            
            comp_type, param_name = group_key
            fixed_values = []
            
            # Collect fixed values for this group
            for i, pinfo in enumerate(param_info):
                if (pinfo['is_fixed'] and 
                    pinfo['param_name'] == param_name and
                    pinfo['component_idx'] is not None):  # Skip tied_phi_rm
                    
                    # Check component type matches
                    comp_idx = pinfo['component_idx']
                    if type(self.model.components[comp_idx]) == comp_type:
                        fixed_values.append((i, pinfo['fixed_value'], comp_idx))
            
            # Check ordering (non-decreasing, by component index)
            if len(fixed_values) > 1:
                fixed_values.sort(key=lambda x: x[2])  # Sort by component index
                
                for j in range(len(fixed_values) - 1):
                    idx1, val1, comp1 = fixed_values[j]
                    idx2, val2, comp2 = fixed_values[j + 1]
                    
                    if val1 > val2:  # Non-decreasing: allow equality
                        raise ValueError(
                            f"Fixed {comp_type.__name__}.{param_name} values violate ordering!\n"
                            f"  {self.param_names[idx1]} = {val1:.3g}\n"
                            f"  {self.param_names[idx2]} = {val2:.3g}\n"
                            f"  Required: {val1:.3g} <= {val2:.3g}"
                        )
    
    def _convert_fixed_sigma_phi_in_effective_width_mode(self, param_info):
        """
        Convert fixed sigma_phi from phi_eff to sigma_phi when effective_width_mode is True.
        
        When a user fixes sigma_phi in effective_width_mode, they're specifying phi_eff,
        not sigma_phi. We need to convert it using the component's N value.
        
        Args:
            param_info: List of parameter info dicts from _build_param_info
        """
        for i, pinfo in enumerate(param_info):
            # Check if this is a fixed sigma_phi in effective_width_mode
            if (pinfo['is_fixed'] and 
                pinfo['param_name'] == 'sigma_phi' and 
                pinfo.get('effective_width_mode', False)):
                
                phi_eff_fixed = pinfo['fixed_value']
                comp_idx = pinfo['component_idx']
                
                # Find N for this component (must be either free or fixed)
                N_value = None
                for j, pinfo_N in enumerate(param_info):
                    if (pinfo_N['component_idx'] == comp_idx and 
                        pinfo_N['param_name'] == 'N'):
                        if pinfo_N['is_fixed']:
                            N_value = pinfo_N['fixed_value']
                        else:
                            # N is free - can't convert fixed phi_eff without knowing N
                            raise ValueError(
                                f"Component {comp_idx}: Cannot fix sigma_phi in "
                                f"effective_width_mode when N is free. Either:\n"
                                f"  1. Fix both sigma_phi and N, or\n"
                                f"  2. Make sigma_phi free (let it be fitted), or\n"
                                f"  3. Turn off effective_width_mode"
                            )
                        break
                
                if N_value is None:
                    raise ValueError(
                        f"Component {comp_idx}: Could not find N parameter. "
                        f"effective_width_mode requires N to be defined."
                    )
                
                # Convert phi_eff to sigma_phi
                # phi_eff is defined at 99% flux containment (1% of peak)
                # Solving: sigma_phi = phi_eff / (2 * ln(100))^(1/N)
                sigma_phi_converted = phi_eff_fixed / (2.0 * np.log(100.0))**(1.0 / N_value)
                
                # Update the fixed value
                self.fixed_params[i] = sigma_phi_converted
    
    def validate(self) -> bool:
        """
        Validate the fit setup.
        
        Checks:
        - Model and priors are consistent
        - Data is valid
        - All required components present
        
        Returns:
            True if valid
            
        Raises:
            ValueError: If validation fails
        """
        # Validate priors match model (using normalized priors)
        validate_priors(self.model, self.priors, p0_physical_max=self.p0_physical_max,
                       tie_phi_rm=self.tie_phi_rm)
        
        # Validate data
        validate_data(self.data)
        
        return True
    
    def log_likelihood(self, params: np.ndarray) -> float:
        """
        Compute log-likelihood for given parameters.
        
        Expands free parameters to full parameter space if needed,
        enforces ordering constraints (returns -inf if violated),
        then calls the module-level log_likelihood function.
        
        Args:
            params: Parameter values (n_free dimensions from dynesty)
            
        Returns:
            Log-likelihood value (-inf if ordering violated or model invalid)
        """
        # If there are fixed parameters, expand to full space
        if len(self.fixed_params) > 0:
            params_full = np.zeros(len(self.param_names))
            
            # Insert free parameters at correct positions
            for i, full_idx in enumerate(self.free_param_indices):
                params_full[full_idx] = params[i]
            
            # Insert fixed values
            for full_idx, value in self.fixed_params.items():
                params_full[full_idx] = value
            
            params = params_full
        
        # Enforce ordering constraints ONLY if likelihood enforcement is needed
        # (Parameters using order-statistics are already guaranteed to be ordered)
        if self.enforce_ordering and self.prior_transform.needs_likelihood_ordering:
            if not _check_ordering(params, self.prior_transform.param_info_full, 
                                  self.model, self.prior_transform.ordering_config):
                return -np.inf
        
        return log_likelihood(params, self.model, self.data, self.param_names,
                            p0_physical_max=self.p0_physical_max,
                            lambda_sq_ref=self.lambda_sq_ref)
    
    def __repr__(self):
        return (f"FitSetup(\n"
                f"  model: {self.model.model_type} with {len(self.model.components)} components\n"
                f"  components: {self.component_names}\n"
                f"  ndim: {self.ndim}\n"
                f"  data: {len(self.data.lambda_sq)} channels\n"
                f"  ordering: {self.enforce_ordering}\n"
                f")")


def _build_param_info(model: CompositeModel, priors: Dict, tie_phi_rm: bool = False) -> List[Dict]:
    """
    Build parameter information list.
    
    Supports both Cartesian (a, b) and polar_conversion styles.
    Always builds parameter list in Cartesian space (a, b) for sampling.
    
    When tie_phi_rm=True, only adds ONE phi_rm parameter that applies to all
    thin/power-law components (reduces dimensionality by N-1).
    
    Returns list of dicts with:
    - name: parameter name (e.g., 'ISM_a' or 'component_0_a')
    - component_idx: which component (or None for tied_phi_rm)
    - component_type: 'ThinComponent' or 'ThickComponent'
    - param_name: parameter within component (e.g., 'a', 'b', 'phi_rm')
    - prior_spec: prior specification tuple (or polar_conversion dict)
    - is_fixed: whether parameter is fixed (not sampled)
    - fixed_value: value if fixed, None otherwise
    - is_polar_conversion: True if using polar_conversion style
    """
    param_info = []
    tied_phi_rm_added = False  # Track if we've added the tied parameter
    
    for i, comp in enumerate(model.components):
        # Use component name if available, otherwise use index
        comp_key = comp.name if comp.name else f'component_{i}'
        
        # Find the prior for this component
        if comp_key in priors:
            comp_priors = priors[comp_key]
        elif f'component_{i}' in priors:
            comp_priors = priors[f'component_{i}']
        else:
            raise ValueError(f"No prior found for component {i} (name='{comp.name}')")
        
        # Determine component type string
        if isinstance(comp, ThinPowerLawComponent):
            comp_type = 'ThinPowerLawComponent'
        elif isinstance(comp, ThinComponent):
            comp_type = 'ThinComponent'
        elif isinstance(comp, ThickComponent):
            comp_type = 'ThickComponent'
        else:
            raise ValueError(f"Unknown component type: {type(comp)}")
        
        # Detect style
        has_polar = 'polar_conversion' in comp_priors
        
        # Always sample in Cartesian space (a, b)
        if isinstance(comp, ThinPowerLawComponent):
            cartesian_params = ['a', 'b']
            other_params = ['phi_rm', 'beta']
        elif isinstance(comp, ThinComponent):
            cartesian_params = ['a', 'b']
            other_params = ['phi_rm']
        else:  # ThickComponent
            cartesian_params = ['a', 'b']
            other_params = ['phi_peak', 'sigma_phi', 'N']
        
        # Add a, b parameters (with polar_conversion marker if applicable)
        for param_name in cartesian_params:
            if has_polar:
                # Store polar_conversion dict as "prior_spec"
                prior_spec = comp_priors['polar_conversion']
                is_fixed = False
                fixed_value = None
            else:
                # Direct Cartesian prior
                prior_spec = comp_priors[param_name]
                is_fixed = isinstance(prior_spec, tuple) and prior_spec[0] == 'fixed'
                fixed_value = prior_spec[1] if is_fixed else None
            
            param_info.append({
                'name': f'{comp_key}_{param_name}',
                'component_idx': i,
                'component_type': comp_type,
                'param_name': param_name,
                'prior_spec': prior_spec,
                'is_fixed': is_fixed,
                'fixed_value': fixed_value,
                'is_polar_conversion': has_polar,
                'effective_width_mode': False  # Only applies to sigma_phi
            })
        
        # Add other parameters
        for param_name in other_params:
            # Special handling for phi_rm when tying is enabled
            if param_name == 'phi_rm' and tie_phi_rm:
                # Only add tied_phi_rm once for all thin/power-law components
                if not tied_phi_rm_added:
                    # Get prior from special 'tied_phi_rm' key
                    prior_spec = priors['tied_phi_rm']
                    is_fixed = isinstance(prior_spec, tuple) and prior_spec[0] == 'fixed'
                    fixed_value = prior_spec[1] if is_fixed else None
                    
                    param_info.append({
                        'name': 'tied_phi_rm',
                        'component_idx': None,  # Applies to all thin/power-law
                        'component_type': 'tied_parameter',
                        'param_name': 'phi_rm',
                        'prior_spec': prior_spec,
                        'is_fixed': is_fixed,
                        'fixed_value': fixed_value,
                        'is_polar_conversion': False,
                        'effective_width_mode': False
                    })
                    tied_phi_rm_added = True
                # Skip adding per-component phi_rm
                continue
            
            prior_spec = comp_priors[param_name]
            is_fixed = isinstance(prior_spec, tuple) and prior_spec[0] == 'fixed'
            fixed_value = prior_spec[1] if is_fixed else None
            
            # Check for effective_width_mode flag (only for sigma_phi in ThickComponents)
            effective_width_mode = False
            if param_name == 'sigma_phi' and comp_type == 'ThickComponent':
                effective_width_mode = comp_priors.get('effective_width_mode', False)
            
            param_info.append({
                'name': f'{comp_key}_{param_name}',
                'component_idx': i,
                'component_type': comp_type,
                'param_name': param_name,
                'prior_spec': prior_spec,
                'is_fixed': is_fixed,
                'fixed_value': fixed_value,
                'is_polar_conversion': False,
                'effective_width_mode': effective_width_mode  # Flag for phi_eff conversion
            })
    
    return param_info


def _transform_parameter(u: float, prior_spec: Tuple) -> float:
    """
    Transform single parameter from unit interval to physical value.
    
    Args:
        u: Value in [0, 1]
        prior_spec: Prior specification tuple
        
    Returns:
        Physical parameter value
        
    Note:
        'fixed' priors should never reach here (they're excluded from sampled space),
        but we handle them defensively to avoid crashes.
    """
    prior_type, params = validate_prior_spec(prior_spec)
    
    if prior_type == 'fixed':
        # Fixed parameters should not be transformed (they're not sampled)
        # But handle defensively to avoid crashes
        return params[0]  # Return the fixed value
    
    elif prior_type == 'uniform':
        pmin, pmax = params
        return pmin + u * (pmax - pmin)
    
    elif prior_type == 'log_uniform':
        pmin, pmax = params
        log_min = np.log(pmin)
        log_max = np.log(pmax)
        return np.exp(log_min + u * (log_max - log_min))
    
    elif prior_type == 'gaussian':
        mean, std = params
        # Clip u to avoid ±inf at exactly 0 or 1
        eps = 1e-15
        u_clipped = np.clip(u, eps, 1.0 - eps)
        return stats.norm.ppf(u_clipped, loc=mean, scale=std)
    
    elif prior_type == 'truncated_gaussian':
        mean, std, pmin, pmax = params
        # Clip u to avoid edge case issues at exactly 0 or 1
        eps = 1e-15
        u_clipped = np.clip(u, eps, 1.0 - eps)
        a = (pmin - mean) / std
        b = (pmax - mean) / std
        return stats.truncnorm.ppf(u_clipped, a, b, loc=mean, scale=std)
    
    else:
        raise ValueError(f"Unknown prior type: {prior_type}")


def _check_ordering(params: np.ndarray, 
                   param_info: List[Dict],
                   model: CompositeModel,
                   ordering_config: Dict = None) -> bool:
    """
    Check if component ordering constraints are satisfied.
    
    Only enforces ordering for parameter groups where ordering_config indicates
    'likelihood' method. Groups using order-statistics in prior transform are
    already ordered by construction.
    
    Uses ordering_config as single source of truth for what needs checking.
    
    Args:
        params: Physical parameter values
        param_info: Parameter information
        model: CompositeModel
        ordering_config: Ordering configuration from PriorTransform (required)
        
    Returns:
        True if ordering satisfied, False otherwise
    """
    if ordering_config is None:
        return True
    
    # Determine which groups need likelihood enforcement
    groups_to_check = {}  # {(comp_type, param_name): should_enforce}
    
    for group_key, config in ordering_config.items():
        needs_ordering = config.get('needs_ordering', False)
        method = config.get('method')
        groups_to_check[group_key] = (needs_ordering and method == 'likelihood')
    
    # Early return if no enforcement needed
    if not any(groups_to_check.values()):
        return True
    
    # Collect parameter values for each group that needs checking
    group_values = {gk: [] for gk, check in groups_to_check.items() if check}
    
    # Scan components and collect values
    for i, comp in enumerate(model.components):
        comp_type = type(comp)
        
        # Check each group that needs enforcement
        for (check_type, param_name), values_list in group_values.items():
            # Only collect if component type matches
            if comp_type == check_type:
                # Find the parameter in param_info
                for j, info in enumerate(param_info):
                    if info['component_idx'] == i and info['param_name'] == param_name:
                        values_list.append(params[j])
                        break
    
    # Check ordering for each group (non-decreasing: allow equality)
    for group_key, values in group_values.items():
        if len(values) > 1:
            if not all(values[i] <= values[i+1] for i in range(len(values)-1)):
                return False
    
    return True



# =============================================================================
# FARADAYFITTER CLASS
# =============================================================================

class FaradayFitter:
    """
    Interface to dynesty nested sampling for Faraday rotation fitting.
    
    This class handles:
    - Running dynesty with sensible defaults
    - Automatic parallelization
    - Progress tracking
    - Returning results
    """
    
    def __init__(self, setup: FitSetup, 
                 sampler: str = 'static',
                 dynesty_kwargs: Optional[Dict] = None,
                 use_pool: bool = True, n_cores: Optional[int] = None):
        """
        Initialize fitter.
        
        Args:
            setup: FitSetup object with model, priors, data
            sampler: Sampler type - 'static' (NestedSampler) or 'dynamic' (DynamicNestedSampler)
            dynesty_kwargs: Optional dictionary with 'sampler' and/or 'run' keys.
                          'sampler' contains kwargs for NestedSampler initialization.
                          'run' contains kwargs for run_nested() method.
                          Example: {'sampler': {'nlive': 1000, 'slices': 5},
                                   'run': {'dlogz': 0.01, 'maxiter': 10000}}
            use_pool: If True, use multiprocessing
            n_cores: Number of cores (None = auto-detect and use N-2)
        """
        # Validate sampler choice
        if sampler not in ['static', 'dynamic']:
            raise ValueError(
                f"sampler must be 'static' or 'dynamic', got '{sampler}'"
            )
        
        self.setup = setup
        self.sampler_type = sampler  # Store as sampler_type, not sampler (which will be the object)
        self.use_pool = use_pool
        
        # Determine number of cores for parallelization
        if use_pool:
            total_cores = multiprocessing.cpu_count()
            if n_cores is None:
                # Auto: use N-2 cores, minimum 1
                self.n_cores = max(1, total_cores - 2)
            else:
                self.n_cores = min(n_cores, total_cores)
            
            print(f"Parallelization: Using {self.n_cores} of {total_cores} available cores")
        else:
            self.n_cores = 1
            print("Parallelization: Disabled")
        
        # Require explicit dynesty_kwargs
        if dynesty_kwargs is None:
            raise ValueError(
                "dynesty_kwargs is REQUIRED. Must specify explicit sampler and run parameters.\n"
                "Example: dynesty_kwargs = {\n"
                "    'sampler': {'nlive': 1000, 'bound': 'multi', 'sample': 'rslice'},\n"
                "    'run': {'dlogz': 0.01, 'maxiter': None, 'maxcall': None}\n"
                "}"
            )
        
        if 'sampler' not in dynesty_kwargs or 'run' not in dynesty_kwargs:
            raise ValueError(
                "dynesty_kwargs must contain BOTH 'sampler' AND 'run' keys.\n"
                "Example: dynesty_kwargs = {\n"
                "    'sampler': {'nlive': 1000, 'bound': 'multi', 'sample': 'rslice'},\n"
                "    'run': {'dlogz': 0.01, 'maxiter': None, 'maxcall': None}\n"
                "}"
            )
        
        # Direct assignment - no defaults, no merging
        self.dynesty_kwargs = dynesty_kwargs
        
        print(f"Dynesty settings:")
        print(f"  Sampler:")
        for key, value in self.dynesty_kwargs['sampler'].items():
            print(f"    {key}: {value}")
        print(f"  Run:")
        for key, value in self.dynesty_kwargs['run'].items():
            print(f"    {key}: {value}")
        
        # Sampler storage (populated during fit)
        self.sampler = None
        
        # Summary statistics (computed after fit)
        self.param_summary = None
        self.chi_squared_best = None
        self.chi_squared_mean = None
        self.chi_squared_std = None
        self.chi_squared_history = None
        self.chi_squared_mean_history = None
        self.chi_squared_std_history = None
        self.bic = None
        self.aic = None
        # Derived samples (populated by compute_summaries in model module)
        # Replaces (a, b) → (p0, psi_0) via direct math transformation
        self.samples_derived = None
        self.param_names_derived = None
        
        # Processed samples (populated by process_posterior_modes in processing module)
        # Unwrapped/continuous around modes for coherent statistics
        self.samples_processed = None
        self.param_names_processed = None
        
        # Mode detection results (populated by process_posterior_modes in processing module)
        self.mode_results = None
        self.processed_parameters = None  # Comprehensive statistics from unwrapped samples
    
    def fit(self, verbose: bool = True, progress_bar: bool = True) -> 'FaradayFitter':
        """
        Run nested sampling fit.
        
        Args:
            verbose: If True, print progress information
            progress_bar: If True, show dynesty progress bar (requires matplotlib)
            
        Returns:
            self (FaradayFitter) with populated sampler and computed summaries
        """
        
        print("\n" + "="*70)
        print("STARTING NESTED SAMPLING FIT")
        print("="*70)
        print(f"Model: {self.setup.model.model_type}")
        print(f"Components: {self.setup.component_names}")
        print(f"Parameters: {self.setup.ndim}")
        print(f"Data channels: {len(self.setup.data.lambda_sq)}")
        print(f"Ordering constraints: {self.setup.enforce_ordering}")
        print(f"Sampler type: {self.sampler_type}")
        print("="*70)
        
        # Prepare sampler kwargs
        sampler_kwargs = self.dynesty_kwargs['sampler'].copy()
        
        # Prepare run kwargs
        run_kwargs = self.dynesty_kwargs['run'].copy()
        if verbose and progress_bar:
            run_kwargs['print_progress'] = True
        else:
            run_kwargs['print_progress'] = False
        
        # Track time
        start_time = time.time()
        
        # Run with or without multiprocessing
        if self.use_pool and self.n_cores > 1:
            
            print(f"\nUsing {self.n_cores} CPUs for parallel processing")
            print("\nRunning nested sampling...")
            
            # Use context manager for proper cleanup
            with Pool(self.n_cores) as pool:
                # Select sampler class based on type
                if self.sampler_type == 'dynamic':
                    SamplerClass = dynesty.DynamicNestedSampler
                else:  # 'static'
                    SamplerClass = dynesty.NestedSampler
                
                # Get periodic boundaries from prior transform
                periodic = self.setup.prior_transform.periodic if len(self.setup.prior_transform.periodic) > 0 else None
                
                sampler = SamplerClass(
                    loglikelihood=self.setup.log_likelihood,
                    prior_transform=self.setup.prior_transform,
                    ndim=self.setup.ndim,
                    periodic=periodic,
                    # periodic=None, # hard code as a check
                    pool=pool,
                    queue_size=self.n_cores,
                    **sampler_kwargs
                )
                
                self.sampler = sampler
                sampler.run_nested(**run_kwargs)
        else:
            print("\nRunning nested sampling (single-core)...")
            
            # Select sampler class based on type
            if self.sampler_type == 'dynamic':
                SamplerClass = dynesty.DynamicNestedSampler
            else:  # 'static'
                SamplerClass = dynesty.NestedSampler
            
            # Get periodic boundaries from prior transform
            periodic = self.setup.prior_transform.periodic if len(self.setup.prior_transform.periodic) > 0 else None
            
            sampler = SamplerClass(
                loglikelihood=self.setup.log_likelihood,
                prior_transform=self.setup.prior_transform,
                ndim=self.setup.ndim,
                periodic=periodic,
                **sampler_kwargs
            )
            
            self.sampler = sampler
            sampler.run_nested(**run_kwargs)
        
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*70)
        print("SAMPLING COMPLETE")
        print("="*70)
        print(f"Total iterations: {self.sampler.results.niter}")
        print(f"Log(Z): {self.sampler.results.logz[-1]:.2f} ± {self.sampler.results.logzerr[-1]:.2f}")
        print(f"Elapsed time: {elapsed_time:.1f} seconds ({elapsed_time/60:.1f} minutes)")
        print("="*70)
        
        # Compute summaries
        print("\nProcessing results...")
        self.compute_summaries()
        
        return self
    
    # =========================================================================
    # Properties for sampler access
    # =========================================================================
    
    @property
    def samples(self) -> np.ndarray:
        """
        Get posterior samples with fixed parameters expanded.
        
        Returns base samples from dynesty, expanded to include fixed parameters
        at their correct positions if any were fixed during fitting.
        """
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        
        base_samples = self.sampler.results.samples
        
        # Expand for fixed params if needed
        if len(self.setup.fixed_params) > 0:
            n_samples = len(base_samples)
            n_total = len(self.setup.param_names)
            expanded = np.zeros((n_samples, n_total))
            
            # Insert free parameter samples at correct positions
            for i, full_idx in enumerate(self.setup.free_param_indices):
                expanded[:, full_idx] = base_samples[:, i]
            
            # Insert fixed values (constant across all samples)
            for full_idx, value in self.setup.fixed_params.items():
                expanded[:, full_idx] = value
            
            return expanded
        
        return base_samples
    
    @property
    def weights(self) -> np.ndarray:
        """Get importance weights from dynesty."""
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        return self.sampler.results.importance_weights()
    
    @property
    def n_eff(self) -> float:
        """Get effective sample size from dynesty."""
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        return self.sampler.n_effective
    
    @property
    def logz(self) -> float:
        """Get log evidence."""
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        return self.sampler.results.logz[-1]
    
    @property
    def logzerr(self) -> float:
        """Get log evidence uncertainty."""
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        return self.sampler.results.logzerr[-1]
    
    @property
    def n_iter(self) -> int:
        """Get number of iterations."""
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        results = self.sampler.results
        return results.niter if hasattr(results, 'niter') else len(results.samples)
    
    @property
    def logl(self) -> np.ndarray:
        """Get log-likelihood values for each sample."""
        if self.sampler is None:
            raise ValueError("No sampler available. Run fit() first.")
        return self.sampler.results.logl
    
    @property
    def n_data(self) -> int:
        """Get number of data points (Q and U for each channel)."""
        return len(self.setup.data.lambda_sq) * 2
    
    @property
    def n_params(self) -> int:
        """Get number of free parameters."""
        return self.setup.ndim
    
    def compute_summaries(self):
        """
        Compute summary statistics from posterior samples.
        
        Computes and stores:
        - Parameter summaries (mean, median, quantiles)
        - Chi-squared statistics (best, mean, history)
        - Model comparison metrics (BIC, AIC)

        Can be called multiple times after continuing sampling with add_batch().
        """
        print("Computing parameter summaries...")
        t_start = time.time()
        self.param_summary = self._compute_param_summary()
        print(f"   Parameter summaries: {time.time()-t_start:.1f}s")
        
        # Compute and add derived parameters (p0, psi_0) from (a, b)
        t_start = time.time()
        self._add_derived_parameters()
        print(f"   Derived parameters (p0, psi_0): {time.time()-t_start:.1f}s")
        
        # Build samples_derived array (replace a/b with p0/psi_0)
        t_start = time.time()
        self._build_samples_derived()
        print(f"   Derived samples array: {time.time()-t_start:.1f}s")
        
        print("Computing model comparison metrics...")
        t_start = time.time()
        self.bic = self._compute_bic()
        self.aic = self._compute_aic()
        print(f"   BIC/AIC: {time.time()-t_start:.1f}s")
        
        print("Computing chi-squared statistics...")
        t0 = time.time()
        
        # Best sample chi-squared
        t_start = time.time()
        chi2_best = self.compute_chi_squared(method='best')
        self.chi_squared_best = chi2_best['total']
        self.chi_squared_Q_best = chi2_best['Q']
        self.chi_squared_U_best = chi2_best['U']
        self.reduced_chi_squared_best = self.chi_squared_best / (self.n_data - self.n_params)
        print(f"   Best-fit χ²: {time.time()-t_start:.1f}s")
        
        # History (thinned samples for convergence diagnostics)
        t_start = time.time()
        chi2_history = self.compute_chi_squared(method='history')
        self.chi_squared_history = chi2_history['chi_squared']
        self.chi_squared_mean_history = chi2_history['mean_history']
        self.chi_squared_std_history = chi2_history['std_history']
        print(f"   χ² convergence history: {time.time()-t_start:.1f}s")
        
        # Weighted mean from top samples
        t_start = time.time()
        n_final = min(1000, len(self.samples))
        chi2_final = self.compute_chi_squared(method='mean', n_samples=n_final)
        self.chi_squared_mean = chi2_final['mean']
        self.chi_squared_std = chi2_final['std']
        self.reduced_chi_squared_mean = self.chi_squared_mean / (self.n_data - self.n_params)
        print(f"   Weighted mean χ²: {time.time()-t_start:.1f}s")
        
        print(f"Total χ² computation: {time.time()-t0:.1f}s")
        print("Summary computation complete!")
    
    def _build_samples_derived(self):
        """
        Build samples_derived array by replacing (a, b) with (p0, psi_0).
        
        Creates interpretable parameter space for plotting where:
        - (a, b) pairs are replaced with (p0, psi_0) 
        - All other parameters remain unchanged
        - Values are NOT unwrapped (simple math transformation only)
        
        Sets:
        - self.samples_derived: Array with replaced parameters
        - self.param_names_derived: Names for derived columns
        """
        param_names = self.setup.param_names
        samples = self.samples
        
        processed_samples_list = []
        processed_names_list = []
        skip_indices = set()
        ab_replacements = {}
        
        # First pass: identify a/b pairs and mark for replacement
        for i, param_name in enumerate(param_names):
            if param_name.endswith('_a'):
                comp_name = param_name[:-2]
                b_name = comp_name + '_b'
                
                if b_name in param_names:
                    b_idx = param_names.index(b_name)
                    a_samples = samples[:, i]
                    b_samples = samples[:, b_idx]
                    
                    # Compute derived parameters
                    p0_samples = np.sqrt(a_samples**2 + b_samples**2)
                    psi_samples = 0.5 * np.arctan2(b_samples, a_samples) * 180 / np.pi
                    
                    p0_name = comp_name + '_p0'
                    psi_name = comp_name + '_psi_0'
                    
                    # Store replacement info
                    ab_replacements[i] = (p0_samples, psi_samples, p0_name, psi_name)
                    
                    # Mark both a and b for skipping
                    skip_indices.add(i)
                    skip_indices.add(b_idx)
        
        # Second pass: build processed array
        for i, param_name in enumerate(param_names):
            if i in skip_indices:
                # Replace a/b pair with p0/psi_0
                if i in ab_replacements:
                    p0_samples, psi_samples, p0_name, psi_name = ab_replacements[i]
                    processed_samples_list.append(p0_samples)
                    processed_names_list.append(p0_name)
                    processed_samples_list.append(psi_samples)
                    processed_names_list.append(psi_name)
            else:
                # Keep parameter as-is
                processed_samples_list.append(samples[:, i])
                processed_names_list.append(param_name)
        
        # Stack into array and store
        self.samples_derived = np.column_stack(processed_samples_list)
        self.param_names_derived = processed_names_list
    
    def _add_derived_parameters(self):
        """
        Compute and add derived parameters (p0, psi_0) to param_summary.
        
        For each (a, b) pair in sampled parameters, computes:
        - p0 = sqrt(a^2 + b^2) - fractional polarization amplitude
        - psi_0 = 0.5 * arctan2(b, a) - polarization angle in degrees
        
        Adds them to param_summary with full statistics.
        """
        param_names = self.setup.param_names
        
        for i, name in enumerate(param_names):
            if name.endswith('_a'):
                comp_name = name[:-2]
                b_name = comp_name + '_b'
                
                if b_name in param_names:
                    a_idx = i
                    b_idx = param_names.index(b_name)
                    
                    # Extract samples
                    a_samples = self.samples[:, a_idx]
                    b_samples = self.samples[:, b_idx]
                    
                    # Compute p0 = sqrt(a^2 + b^2)
                    p0_samples = np.sqrt(a_samples**2 + b_samples**2)
                    p0_stats = self._compute_derived_param_stats(p0_samples)
                    p0_name = comp_name + '_p0'
                    self.param_summary[p0_name] = p0_stats
                    
                    # Compute psi_0 = 0.5 * arctan2(b, a) in degrees
                    psi_samples = 0.5 * np.arctan2(b_samples, a_samples) * 180 / np.pi
                    psi_stats = self._compute_derived_param_stats(psi_samples)
                    psi_name = comp_name + '_psi_0'
                    self.param_summary[psi_name] = psi_stats
    
    def _compute_param_summary(self) -> Dict:
        """
        Compute parameter statistics using dynesty's built-in functions.
        
        Returns Cartesian parameters only with:
        - Mean and covariance from dyfunc.mean_and_cov() (full array)
        - Median, 68% CI, and 99% range from dyfunc.quantile() (per parameter)
        
        Returns:
            Dictionary with parameter names as keys and summary dicts as values
        """
        # Filter samples
        good_mask = np.all(np.isfinite(self.samples), axis=1) & np.isfinite(self.weights) & (self.weights > 0)
        samples_clean = self.samples[good_mask]
        weights_clean = self.weights[good_mask]
        
        # Compute mean and covariance for all parameters at once
        mean_all, cov_all = dyfunc.mean_and_cov(samples_clean, weights=weights_clean)
        std_all = np.sqrt(np.diag(cov_all))
        
        # Compute quantiles per parameter (dyfunc.quantile requires 1D input)
        quantiles_all = np.array([
            dyfunc.quantile(samples_clean[:, i], [0.005, 0.16, 0.5, 0.84, 0.995], weights=weights_clean)
            for i in range(samples_clean.shape[1])
        ])
        
        summary = {}
        
        for i, param_name in enumerate(self.setup.param_names):
            if self.setup.is_fixed[i]:
                fixed_value = self.setup.fixed_params[i]
                summary[param_name] = {
                    'mean': fixed_value,
                    'std': 0.0,
                    'median': fixed_value,
                    'eti_68': [fixed_value, fixed_value],
                    'eti_99': [fixed_value, fixed_value],
                    'fixed': True
                }
            else:
                summary[param_name] = {
                    'mean': float(mean_all[i]),
                    'std': float(std_all[i]),
                    'median': float(quantiles_all[i, 2]),
                    'eti_68': [float(quantiles_all[i, 1]), float(quantiles_all[i, 3])],
                    'eti_99': [float(quantiles_all[i, 0]), float(quantiles_all[i, 4])],
                    'fixed': False
                }
                
                if param_name.endswith('_N'):
                    summary[param_name]['min'] = float(np.min(samples_clean[:, i]))
        
        # If RM tying is enabled, populate per-component phi_rm from tied_phi_rm
        if self.setup.tie_phi_rm and 'tied_phi_rm' in summary:
            tied_stats = summary['tied_phi_rm']
            for i, comp in enumerate(self.setup.model.components):
                if isinstance(comp, (ThinComponent, ThinPowerLawComponent)):
                    comp_name = comp.name if comp.name else f'component_{i}'
                    comp_phi_rm_name = f'{comp_name}_phi_rm'
                    # Duplicate the tied_phi_rm stats for this component
                    summary[comp_phi_rm_name] = tied_stats.copy()
        
        return summary
    
    def _compute_bic(self) -> float:
        """
        Compute Bayesian Information Criterion.
        
        BIC = -2 * ln(L_max) + k * ln(n)
        where k = number of parameters, n = number of data points
        
        We approximate ln(L_max) using the maximum log-likelihood from samples.
        """
        # Get maximum log-likelihood
        max_loglike = np.max(self.sampler.results.logl)
        
        bic = -2 * max_loglike + self.n_params * np.log(self.n_data)
        return float(bic)
    
    def _compute_aic(self) -> float:
        """
        Compute Akaike Information Criterion.
        
        AIC = -2 * ln(L_max) + 2 * k
        where k = number of parameters
        """
        max_loglike = np.max(self.sampler.results.logl)
        aic = -2 * max_loglike + 2 * self.n_params
        return float(aic)
    
    def compute_chi_squared(self, method: Literal['best', 'mean', 'history'] = 'best', 
                           n_samples: Optional[int] = None) -> Dict:
        """
        Compute chi-squared goodness of fit.

        Args:
            method: Computation method
                'best': Chi² for MAP sample (highest likelihood)
                'mean': Importance-weighted mean ± std across posterior
                'history': Running chi² through posterior samples
            n_samples: For method='mean', use only last n_samples (default: all)
        
        Returns:
            Dict with method-specific keys:
                method='best': {'total', 'Q', 'U'}
                method='mean': {'mean', 'std'}
                method='history': {'chi_squared', 'mean_history', 'std_history'}
        """
        Q_obs = self.setup.data.Q
        U_obs = self.setup.data.U
        sigma_Q_data = self.setup.data.Q_err
        sigma_U_data = self.setup.data.U_err
        lambda_sq = self.setup.data.lambda_sq
        
        if method == 'best':
            # Find best sample from valid likelihood values
            valid_mask = np.isfinite(self.logl)
            valid_indices = np.where(valid_mask)[0]
            
            if len(valid_indices) == 0:
                raise ValueError("No valid samples found in posterior")
            
            best_idx_in_valid = np.argmax(self.logl[valid_mask])
            best_sample_idx = valid_indices[best_idx_in_valid]
            best_sample_params = self.samples[best_sample_idx]
            
            sigma_Q_total = sigma_Q_data
            sigma_U_total = sigma_U_data

            # Compute chi-squared from best sample
            best_model = copy.deepcopy(self.setup.model)
            update_model_from_params(best_model, best_sample_params, self.setup.param_names, lambda_sq_ref=self.setup.lambda_sq_ref)
            
            P_best = best_model.compute_polarization(lambda_sq)
            chi2_Q = float(np.sum(((Q_obs - P_best.real) / sigma_Q_total)**2))
            chi2_U = float(np.sum(((U_obs - P_best.imag) / sigma_U_total)**2))
            chi2_total = chi2_Q + chi2_U
            
            return {
                'total': chi2_total,
                'Q': chi2_Q,
                'U': chi2_U
            }
        
        elif method == 'mean':
            import time
            t_start = time.time()
            
            # Use top n_samples by importance weight
            total_samples = len(self.samples)
            if n_samples is None:
                n_use = total_samples
                selected_indices = np.arange(total_samples)
            else:
                n_use = min(n_samples, total_samples)
                selected_indices = np.argsort(self.weights)[-n_use:]
            
            # Create single model object to reuse sample
            model_reusable = copy.deepcopy(self.setup.model)
            
            # Pre-allocate arrays for batch computation
            n_wavelengths = len(lambda_sq)
            Q_models = np.zeros((n_use, n_wavelengths))
            U_models = np.zeros((n_use, n_wavelengths))
            sigma_Q_array = np.zeros((n_use, n_wavelengths))
            sigma_U_array = np.zeros((n_use, n_wavelengths))
            
            # Loop to compute model predictions (still need Python for object manipulation)
            for i, sample_idx in enumerate(selected_indices):
                params = self.samples[sample_idx]

                sigma_Q_array[i] = sigma_Q_data
                sigma_U_array[i] = sigma_U_data

                # Update model parameters
                update_model_from_params(model_reusable, params, self.setup.param_names, lambda_sq_ref=self.setup.lambda_sq_ref)

                # Compute predictions
                P_model = model_reusable.compute_polarization(lambda_sq)
                Q_models[i] = P_model.real
                U_models[i] = P_model.imag
            
            print(f"      Computed {n_use} model predictions: {time.time()-t_start:.1f}s")
            
            # Batch compute chi-squared with numba
            t_chi2 = time.time()
            chi_sq_samples = compute_chi2_batch_numba(Q_obs, U_obs, Q_models, U_models,
                                                      sigma_Q_array, sigma_U_array)
            print(f"      Computed chi-squared: {time.time()-t_chi2:.3f}s")
            
            # Get corresponding weights
            weights_selected = self.weights[selected_indices]
            
            # Compute importance-weighted mean and std
            finite_mask = np.isfinite(chi_sq_samples)
            
            if np.any(finite_mask):
                weights_finite = weights_selected[finite_mask]
                chi_sq_finite = chi_sq_samples[finite_mask]
                
                # Normalize weights
                weight_sum = np.sum(weights_finite)
                if weight_sum > 1e-300:
                    weights_norm = weights_finite / weight_sum
                    
                    chi_sq_mean = np.sum(weights_norm * chi_sq_finite)
                    chi_sq_var = np.sum(weights_norm * (chi_sq_finite - chi_sq_mean)**2)
                    chi_sq_std = np.sqrt(chi_sq_var)
                else:
                    chi_sq_mean = np.nan
                    chi_sq_std = np.nan
            else:
                chi_sq_mean = np.nan
                chi_sq_std = np.nan
            
            return {
                'mean': float(chi_sq_mean),
                'std': float(chi_sq_std)
            }
        
        elif method == 'history':
            # Apply adaptive thinning (max 100 samples for speed)
            n_samples = len(self.samples)
            thin_factor = max(1, n_samples // 100)
            thinned_indices = np.arange(0, n_samples, thin_factor)
            
            n_thinned = len(thinned_indices)
            chi2_history = np.zeros(n_thinned)
            mean_history = np.zeros(n_thinned)
            std_history = np.zeros(n_thinned)
            
            import time
            
            # Create single model object to reuse
            model_reusable = copy.deepcopy(self.setup.model)
            
            # Pre-allocate arrays for batch computation
            n_wavelengths = len(lambda_sq)
            Q_models = np.zeros((n_thinned, n_wavelengths))
            U_models = np.zeros((n_thinned, n_wavelengths))
            sigma_Q_array = np.zeros((n_thinned, n_wavelengths))
            sigma_U_array = np.zeros((n_thinned, n_wavelengths))
            
            # Compute chi-square at thinned indices
            t_start = time.time()
            for idx, sample_idx in enumerate(thinned_indices):
                params = self.samples[sample_idx]

                sigma_Q_array[idx] = sigma_Q_data
                sigma_U_array[idx] = sigma_U_data

                # Update model parameters
                update_model_from_params(model_reusable, params, self.setup.param_names, lambda_sq_ref=self.setup.lambda_sq_ref)

                # Compute predictions
                P_model = model_reusable.compute_polarization(lambda_sq)
                Q_models[idx] = P_model.real
                U_models[idx] = P_model.imag
            
            print(f"      Computed {n_thinned} model predictions: {time.time()-t_start:.1f}s")
            
            # Batch compute chi-squared with numba
            t_chi2 = time.time()
            chi2_at_thinned = compute_chi2_batch_numba(Q_obs, U_obs, Q_models, U_models,
                                                       sigma_Q_array, sigma_U_array)
            chi2_history = chi2_at_thinned.copy()
            print(f"      Computed chi-squared: {time.time()-t_chi2:.3f}s")
            
            # Compute rolling window statistics (last 100 samples)
            t_start = time.time()
            window_size = 100
            for idx in range(n_thinned):
                start_idx = max(0, idx + 1 - window_size)
                end_idx = idx + 1
                
                chi2_window = chi2_at_thinned[start_idx:end_idx]
                weights_window = self.weights[thinned_indices[start_idx:end_idx]]
                
                # Filter out NaN/inf values
                finite_mask = np.isfinite(chi2_window)
                if np.any(finite_mask):
                    weights_finite = weights_window[finite_mask]
                    chi2_finite = chi2_window[finite_mask]
                    
                    # Normalize weights
                    weight_sum = np.sum(weights_finite)
                    if weight_sum > 1e-300:
                        weights_norm = weights_finite / weight_sum
                        mean_history[idx] = np.sum(weights_norm * chi2_finite)
                        var = np.sum(weights_norm * (chi2_finite - mean_history[idx])**2)
                        std_history[idx] = np.sqrt(var)
                    else:
                        mean_history[idx] = np.nan
                        std_history[idx] = np.nan
                else:
                    mean_history[idx] = np.nan
                    std_history[idx] = np.nan
            
            print(f"      Computed rolling statistics: {time.time()-t_start:.1f}s")
            
            return {
                'chi_squared': chi2_history,
                'mean_history': mean_history,
                'std_history': std_history
            }
        
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'best', 'mean', or 'history'.")
    
    def save(self, filename: str) -> None:
        """
        Save fitter with sampler and all results.
        
        Saves entire FaradayFitter object including sampler, setup, and 
        computed summaries. Can be loaded with load_results() for continued
        sampling or analysis.
        
        Args:
            filename: Path to output file (e.g., 'fits/results.pkl')
            
        Example:
            fitter.save('fits/results.pkl')
            fitter2 = load_results('fits/results.pkl')
            fitter2.sampler.add_batch(nlive=100)
            fitter2.compute_summaries()
        """
        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
        
        print(f"Saved results: {filepath}")
    
    def get_best_fit_params(self, method: Literal['median', 'mean', 'map', 'mode'] = 'median') -> np.ndarray:
        """
        Get best-fit parameter values.
        
        Args:
            method: How to define "best"
                - 'median': Use weighted median (default, robust)
                - 'mean': Use weighted mean
                - 'map': Maximum a posteriori (highest likelihood sample)
                - 'mode': Use histogram-based mode (falls back to median if not available)
                
        Returns:
            Array of best-fit parameter values
        """
        if method == 'mode':
            # Check if mode exists, otherwise fall back to median
            if 'mode' in self.param_summary[self.setup.param_names[0]]:
                return np.array([self.param_summary[name]['mode'] 
                               for name in self.setup.param_names])
            else:
                return np.array([self.param_summary[name]['median'] 
                               for name in self.setup.param_names])
        
        elif method == 'median':
            return np.array([self.param_summary[name]['median'] 
                           for name in self.setup.param_names])
        
        elif method == 'mean':
            return np.array([self.param_summary[name]['mean'] 
                           for name in self.setup.param_names])
        
        elif method == 'map':
            # Find sample with highest likelihood
            max_idx = np.argmax(self.sampler.results.logl)
            return self.samples[max_idx]
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def get_best_fit_model(self, method: Literal['median', 'mean', 'map', 'mode'] = 'median') -> CompositeModel:
        """
        Get CompositeModel with best-fit parameters.
        
        Args:
            method: How to define "best" (see get_best_fit_params)
                - 'median': Use weighted median (default)
                - 'mean': Use weighted mean
                - 'map': Maximum a posteriori (highest likelihood sample)
                - 'mode': Use histogram-based mode (falls back to median if not available)
            
        Returns:
            CompositeModel with updated parameters
        """
        # Make a deep copy of the model
        best_model = copy.deepcopy(self.setup.model)
        
        # Get best-fit parameters
        best_params = self.get_best_fit_params(method=method)
        
        # Update model
        update_model_from_params(best_model, best_params, self.setup.param_names, lambda_sq_ref=self.setup.lambda_sq_ref)
        
        return best_model
    
    def get_representative_sample(self, percentile: float = 1.0, 
                                  method: Literal['median', 'mean', 'mode'] = 'median') -> np.ndarray:
        """
        Get representative sample from high-likelihood region near best-fit.
        
        This addresses the issue that individual parameter best-fits may not correspond
        to a high-likelihood sample due to parameter covariances. Instead, find
        an actual posterior sample with both high likelihood and proximity to best-fit.
        
        Algorithm:
        1. Identify top percentile of samples by likelihood
        2. Calculate normalized distance from each to best-fit
        3. Return sample with minimum distance to best-fit
        
        Args:
            percentile: Likelihood percentile threshold (default: 1.0 = top 1%)
            method: Method to define best-fit reference point
                - 'median': Use weighted median (default)
                - 'mean': Use weighted mean
                - 'mode': Use histogram-based mode (falls back to median if not available)
            
        Returns:
            Array of parameter values from actual posterior sample
        """
        # Get best-fit for each parameter
        best_fit_params = self.get_best_fit_params(method=method)
        
        # Find top percentile of samples by likelihood
        likelihood_threshold = np.percentile(self.logl, 100 - percentile)
        high_likelihood_mask = self.logl >= likelihood_threshold
        high_likelihood_samples = self.samples[high_likelihood_mask]
        
        if len(high_likelihood_samples) == 0:
            # Fallback: return MAP if no samples above threshold
            max_idx = np.argmax(self.logl)
            return self.samples[max_idx]
        
        # Compute normalized distance to best-fit for each high-likelihood sample
        # Normalize by standard deviation to account for different parameter scales
        stds = np.array([self.param_summary[name]['std'] for name in self.setup.param_names])
        
        # Avoid division by zero for fixed parameters
        stds = np.where(stds > 0, stds, 1.0)
        
        # Calculate normalized Euclidean distance
        distances = np.sqrt(np.sum(((high_likelihood_samples - best_fit_params) / stds)**2, axis=1))
        
        # Find sample with minimum distance to best-fit
        min_distance_idx = np.argmin(distances)
        representative_sample = high_likelihood_samples[min_distance_idx]
        
        return representative_sample
    
    
    def _format_sampling_config(self) -> Dict:
        """
        Format sampling configuration for JSON export.
        
        Returns:
            Dictionary with sampler_type, n_cores, sampler_kwargs, and run_kwargs
        """
        config = {
            'sampler_type': self.sampler_type,
            'n_cores': self.n_cores,
            'sampler_kwargs': {},
            'run_kwargs': {}
        }
        
        for key, value in self.dynesty_kwargs['sampler'].items():
            if isinstance(value, (int, float, str, bool, type(None))):
                config['sampler_kwargs'][key] = value
            elif isinstance(value, (list, tuple)):
                config['sampler_kwargs'][key] = list(value)
            else:
                config['sampler_kwargs'][key] = str(value)
        
        for key, value in self.dynesty_kwargs['run'].items():
            if isinstance(value, (int, float, str, bool, type(None))):
                config['run_kwargs'][key] = value
            elif isinstance(value, (list, tuple)):
                config['run_kwargs'][key] = list(value)
            else:
                config['run_kwargs'][key] = str(value)
        
        return config
    
    def _compute_derived_param_stats(self, samples: np.ndarray) -> Dict:
        """
        Compute importance-weighted statistics for derived parameters.
        
        Args:
            samples: 1D array of derived parameter samples
            
        Returns:
            Dictionary with mean, std, median, and quantiles
        """
        # Import dynesty functions
        from dynesty import utils as dyfunc
        
        # Compute importance-weighted mean and covariance
        mean, cov = dyfunc.mean_and_cov(samples.reshape(-1, 1), self.weights)
        mean_val = float(mean[0])
        std_val = float(np.sqrt(cov[0, 0]))
        
        # Compute quantiles
        quantiles = [0.005, 0.16, 0.5, 0.84, 0.995]
        q_values = dyfunc.quantile(samples, quantiles, weights=self.weights)
        
        return {
            'mean': mean_val,
            'std': std_val,
            'median': float(q_values[2]),
            'eti_68': [float(q_values[1]), float(q_values[3])],
            'eti_99': [float(q_values[0]), float(q_values[4])]
        }
    
    def save_json(self, filename: str, I_fitter=None):
        """
        Save results to JSON file.

        Args:
            filename: Path to output JSON file
            I_fitter: Optional IFitter instance; if provided and is a power-law fit,
                      alpha and its 1-sigma error are stored under 'stokes_i_fit'.
        """
        # Build output dictionary
        output = {
            'metadata': {
                'model_type': self.setup.model.model_type,
                'n_components': len(self.setup.model.components),
                'component_names': self.setup.component_names,
                'n_params': self.n_params,
                'n_data': self.n_data,
                'n_iter': int(self.n_iter),
                'n_eff': float(self.n_eff),
                'ordering_enforced': self.setup.enforce_ordering
            },
            'sampling_parameters': self._format_sampling_config(),
            'evidence': {
                'logz': float(self.logz),
                'logzerr': float(self.logzerr)
            },
            'model_comparison': {
                'BIC': float(self.bic),
                'AIC': float(self.aic)
            },
            'goodness_of_fit': {
                'best_sample': {
                    'chi_squared': float(self.chi_squared_best),
                    'chi_squared_Q': float(self.chi_squared_Q_best),
                    'chi_squared_U': float(self.chi_squared_U_best),
                    'reduced_chi_squared': float(self.reduced_chi_squared_best)
                },
                'importance_weighted': {
                    'chi_squared_mean': float(self.chi_squared_mean),
                    'chi_squared_std': float(self.chi_squared_std),
                    'reduced_chi_squared_mean': float(self.reduced_chi_squared_mean)
                },
                'dof': int(self.n_data - self.n_params)
            },
            'parameters': {}
        }
        
        # Add all parameters (sampled, derived, and fixed) from param_summary
        for param_name, summary in self.param_summary.items():
            # Determine parameter type based on name and fixed status
            if param_name.endswith('_p0') or param_name.endswith('_psi_0'):
                param_type = 'derived'
            elif summary.get('fixed', False):
                param_type = 'fixed'
            else:
                param_type = 'sampled'
            
            param_dict = {
                'parameter_type': param_type,
                'mean': float(summary['mean']),
                'std': float(summary['std']),
                'median': float(summary['median']),
                'eti_68': [float(summary['eti_68'][0]), float(summary['eti_68'][1])],
                'eti_99': [float(summary['eti_99'][0]), float(summary['eti_99'][1])]
            }
            
            if 'min' in summary:
                param_dict['min'] = float(summary['min'])
            
            output['parameters'][param_name] = param_dict
        
        # Add processed parameters if available (from process_posterior_modes)
        if hasattr(self, 'processed_parameters') and self.processed_parameters is not None:
            output['processed_parameters'] = self.processed_parameters
        
        # Add Stokes I power-law fit parameters if available
        if (I_fitter is not None and I_fitter.model == 'power-law'
                and I_fitter.fit_params is not None and I_fitter.fit_cov is not None):
            alpha_lambda2 = float(I_fitter.fit_params[1])
            alpha_lambda2_err = float(np.sqrt(I_fitter.fit_cov[1, 1]))
            # alpha_nu = -2 * alpha_lambda2  (since lambda^2 = c^2/nu^2)
            alpha_nu = -2.0 * alpha_lambda2
            alpha_nu_err = 2.0 * alpha_lambda2_err
            lambda_sq_0 = float(I_fitter.lambda_sq_0)
            freq0_GHz = float(C_LIGHT / np.sqrt(lambda_sq_0) / 1e9)
            output['stokes_i_fit'] = {
                'model': 'power-law',
                'I0': float(I_fitter.fit_params[0]),
                'I0_err': float(np.sqrt(I_fitter.fit_cov[0, 0])),
                'alpha_lambda2': alpha_lambda2,
                'alpha_lambda2_err': alpha_lambda2_err,
                'alpha_nu': alpha_nu,
                'alpha_nu_err': alpha_nu_err,
                'freq0_GHz': freq0_GHz
            }

        # Write to file
        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, 'w') as f:
            json.dump(output, f, indent=2)

        print(f"Results saved to: {filepath}")
    
    def print_summary(self):
        """Print a formatted summary of results."""
        
        print("\n" + "="*70)
        print("FIT RESULTS SUMMARY")
        print("="*70)
        
        print(f"\nModel: {self.setup.model.model_type}")
        print(f"Components: {', '.join(self.setup.component_names)}")
        print(f"Parameters: {self.n_params}")
        print(f"Data points: {self.n_data}")
        
        print(f"\nEvidence:")
        print(f"  ln(Z) = {self.logz:.2f} ± {self.logzerr:.2f}")
        
        print(f"\nModel comparison:")
        print(f"  BIC = {self.bic:.2f}")
        print(f"  AIC = {self.aic:.2f}")
        
        print(f"\nGoodness of fit (Best Sample):")
        print(f"  χ² = {self.chi_squared_best:.2f}")
        print(f"    Q: {self.chi_squared_Q_best:.2f}")
        print(f"    U: {self.chi_squared_U_best:.2f}")
        print(f"  Reduced χ² = {self.reduced_chi_squared_best:.3f}")
        
        print(f"\nGoodness of fit (Importance-Weighted):")
        print(f"  χ² = {self.chi_squared_mean:.2f} ± {self.chi_squared_std:.2f}")
        print(f"  Reduced χ² = {self.reduced_chi_squared_mean:.3f}")
        
        dof = self.n_data - self.n_params
        print(f"  Degrees of freedom = {dof}")
        
        rchi = self.reduced_chi_squared_mean
        if rchi < 0.5:
            print(f"  → χ²_red = {rchi:.3f}: likely overfitting or errors significantly overestimated")
        elif rchi < 1.5:
            print(f"  → χ²_red = {rchi:.3f}: good fit")
        elif rchi < 2.0:
            print(f"  → χ²_red = {rchi:.3f}: acceptable — marginal systematics or mild model mismatch")
        elif rchi < 3.0:
            print(f"  → χ²_red = {rchi:.3f}: questionable — likely unaccounted systematics or model inadequacy")
        else:
            print(f"  → χ²_red = {rchi:.3f}: poor fit — model is strongly inadequate or errors are severely underestimated")
        
        print(f"\nSampling:")
        print(f"  Iterations: {self.n_iter}")
        print(f"  Effective samples: {self.n_eff:.1f}")
        
        print(f"\nParameter estimates:")
        for param_name, summary in self.param_summary.items():
            mean = summary['mean']
            std = summary['std']
            median = summary['median']
            eti_68 = summary['eti_68']
            eti_99 = summary['eti_99']
            
            err_lower_1sigma = median - eti_68[0]
            err_upper_1sigma = eti_68[1] - median
            
            print(f"\n  {param_name}:")
            print(f"    Mean:   {mean:8.4f} ± {std:.4f}")
            print(f"    Median: {median:8.4f} +{err_upper_1sigma:.4f} -{err_lower_1sigma:.4f}  (68%: [{eti_68[0]:.4f}, {eti_68[1]:.4f}])")
            print(f"    99% range: [{eti_99[0]:.4f}, {eti_99[1]:.4f}]")
            
            if 'min' in summary:
                print(f"    Min: {summary['min']:.4f}")
        
        print("="*70)
    
    def __repr__(self):
        return (f"FaradayFitter(\n"
                f"  model: {self.setup.model.model_type}\n"
                f"  ln(Z): {self.logz:.2f} ± {self.logzerr:.2f}\n"
                f"  BIC: {self.bic:.2f}\n"
                f"  n_iter: {self.n_iter}\n"
                f")")


def load_results(filename: str) -> 'FaradayFitter':
    """
    Load saved fit results from pickle file.
    
    Loads FaradayFitter object with sampler, setup, and computed summaries.
    Can be used for plotting, analysis, and continuing sampling.
    
    Args:
        filename: Path to saved results file (from save())
        
    Returns:
        FaradayFitter object ready for analysis and plotting
        
    Example:
        fitter = FaradayFitter(...)
        fitter.fit()
        fitter.save('fits/results.pkl')
        
        fitter2 = load_results('fits/results.pkl')
        fitter2.sampler.add_batch(nlive=100)
        fitter2.compute_summaries()
        fitter2.save('fits/results_continued.pkl')
    """
    filepath = Path(filename)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Results file not found: {filepath}")
    
    with open(filepath, 'rb') as f:
        fitter = pickle.load(f)
    
    if not isinstance(fitter, FaradayFitter):
        raise TypeError(
            f"Loaded object is not a FaradayFitter instance. "
            f"Got {type(fitter).__name__}."
        )
    
    print(f"Loaded results: {filepath}")
    print(f"  Model: {fitter.setup.model.model_type}")
    print(f"  Components: {len(fitter.setup.model.components)}")
    print(f"  Parameters: {fitter.n_params}")
    print(f"  Iterations: {fitter.n_iter}")
    
    return fitter


def get_ndim(model: CompositeModel, tie_phi_rm: bool = False) -> int:
    """
    Get number of parameters in model.
    
    When tie_phi_rm=True, reduces count by (N-1) where N is the number
    of thin/power-law components, since they all share one phi_rm.
    
    Args:
        model: CompositeModel
        tie_phi_rm: If True, account for tied phi_rm (default: False)
        
    Returns:
        Number of free parameters
    """
    ndim = 0
    n_thin_powerlaw = 0  # Count components with phi_rm
    
    for comp in model.components:
        if isinstance(comp, ThinPowerLawComponent):
            ndim += 4  # a, b, phi_rm, beta
            n_thin_powerlaw += 1
        elif isinstance(comp, ThinComponent):
            ndim += 3  # a, b, phi_rm
            n_thin_powerlaw += 1
        elif isinstance(comp, ThickComponent):
            ndim += 5  # a, b, phi_peak, sigma_phi, N
    
    # Adjust for tied phi_rm (reduces by N-1 where N = number of thin/power-law components)
    if tie_phi_rm and n_thin_powerlaw > 1:
        ndim -= (n_thin_powerlaw - 1)
    
    return ndim


# =============================================================================
# NUMBA SPEED UPS
# =============================================================================

if NUMBA_AVAILABLE:
    @numba.jit(nopython=True, parallel=True)
    def compute_chi2_batch_numba(Q_obs, U_obs, Q_models, U_models, 
                                  sigma_Q_array, sigma_U_array):
        """
        Compute chi-squared for batch of model predictions using numba.
        
        Parallelizes chi-squared computation across samples for ~2-5x speedup.
        Most time is still spent in model.compute_polarization() (already
        optimized with numba in faraday_utils.py).
        
        Args:
            Q_obs: Observed Q values (n_wavelengths,)
            U_obs: Observed U values (n_wavelengths,)
            Q_models: Model Q predictions (n_samples, n_wavelengths)
            U_models: Model U predictions (n_samples, n_wavelengths)
            sigma_Q_array: Q uncertainties for each sample (n_samples, n_wavelengths)
            sigma_U_array: U uncertainties for each sample (n_samples, n_wavelengths)
            
        Returns:
            chi2_array: Chi-squared for each sample (n_samples,)
        """
        n_samples = Q_models.shape[0]
        chi2_array = np.zeros(n_samples)
        
        for i in numba.prange(n_samples):
            chi2_Q = np.sum(((Q_obs - Q_models[i]) / sigma_Q_array[i])**2)
            chi2_U = np.sum(((U_obs - U_models[i]) / sigma_U_array[i])**2)
            chi2_array[i] = chi2_Q + chi2_U
        
        return chi2_array
else:
    def compute_chi2_batch_numba(Q_obs, U_obs, Q_models, U_models, 
                                  sigma_Q_array, sigma_U_array):
        """Fallback implementation when numba is not available."""
        n_samples = Q_models.shape[0]
        chi2_array = np.zeros(n_samples)
        for i in range(n_samples):
            chi2_Q = np.sum(((Q_obs - Q_models[i]) / sigma_Q_array[i])**2)
            chi2_U = np.sum(((U_obs - U_models[i]) / sigma_U_array[i])**2)
            chi2_array[i] = chi2_Q + chi2_U
        return chi2_array
