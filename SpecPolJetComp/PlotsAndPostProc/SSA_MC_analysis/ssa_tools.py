#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:48:48 2026

@author: cowie
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root, newton
import scipy.special as special
import yaml
import warnings
from scipy.special import expit  # numerically stable sigmoid

from mc_tools import (mc_uncertainty, plot_mc_diagnostics, save_mc_samples,
                      draw_params, plot_mc_input_params, plot_mc_corner)


######################## Constants

# Set to a YAML file path string to override auto-detection.
# Leave as "" to search for inputs.yaml in the current working directory.
INPUT_YAML = ""

electron_mass_cgs = 9.1094e-28
electron_charge_cgs = 4.803e-10
c_cgs = 3e10
c1 = 6.27e18
kB_cgs = 1.38e-16

######################## Pseudo-constants and helper functions

def c5_func(p):
    prefactor = np.sqrt(3)/(16*np.pi) * electron_charge_cgs**3/(electron_mass_cgs * c_cgs**2)
    gammas = (p+7/3)/(p+1) * special.gamma((3*p - 1)/12) * special.gamma((3*p+7)/12)
    return prefactor*gammas

def c6_func(p):
    prefactor = np.sqrt(3) * np.pi/72 *  electron_charge_cgs * electron_mass_cgs**5 * c_cgs**10 
    gammas = (p+10/3) * special.gamma((3*p+2)/12) * special.gamma((3*p+10)/12)
    return prefactor*gammas

def k1_func(p):
    k1, _, _ = _k_funcs(p)
    return k1

def k3_e_func(p):
    _, k3_e, _ = _k_funcs(p)
    return k3_e

def k3_nu_func(p):
    _, _, k3_nu = _k_funcs(p)
    return k3_nu

def K_E_func(p, Emin, Emax):
    # np.where evaluates BOTH branches before selecting; suppress the expected
    # divide-by-zero / power warnings from the discarded branch.
    p = np.asarray(p, dtype=float)
    near_two = np.abs(p - 2.0) < 0.01
    denom = np.where(near_two, 1.0, 2.0 - p)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        return np.where(near_two,
                        np.log(Emax / Emin),
                        (Emax**(2.0 - p) - Emin**(2.0 - p)) / denom)

def K_E_dens_func(p, Emin, Emax):
    p = np.asarray(p, dtype=float)
    near_one = np.abs(p - 1.0) < 0.01
    denom = np.where(near_one, 1.0, 1.0 - p)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        return np.where(near_one,
                        np.log(Emax / Emin),
                        (Emax**(1.0 - p) - Emin**(1.0 - p)) / denom)

def K_nu_func(p, nu_min, nu_max):
    p = np.asarray(p, dtype=float)
    near_two = np.abs(p - 2.0) < 0.01
    denom = np.where(near_two, 1.0, 2.0 - p)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        return np.where(near_two,
                        0.5 * np.log(nu_max / nu_min),
                        0.5 * c1**((p - 2) / 2) * 2.0 / denom
                        * (nu_max**((2.0 - p) / 2) - nu_min**((2.0 - p) / 2)))

def K_nu_dens_func(p, nu_min, nu_max):
    p = np.asarray(p, dtype=float)
    near_one = np.abs(p - 1.0) < 0.01
    denom = np.where(near_one, 1.0, 1.0 - p)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        return np.where(near_one,
                        0.5 * np.log(nu_max / nu_min),
                        0.5 * c1**((p - 1) / 2) * 2.0 / denom
                        * (nu_max**((1.0 - p) / 2) - nu_min**((1.0 - p) / 2)))
    
def p_from_alpha(alpha_thin):
    return 1-2*alpha_thin

def gamma_to_beta(gamma):
    return np.sqrt(1-1/gamma**2)

def doppler_func(gamma, inc):
    beta = gamma_to_beta(gamma)
    return (gamma*(1-beta*np.cos(np.pi/180 * inc)))**(-1)

def optical_depth_to_be_minimised(tau, p):
    return np.exp(tau) - 1 - (p+4)/5 * tau

# --- tau_m precomputed lookup table ------------------------------------------
# tau_m depends only on p = 1 - 2*alpha_thin.  Precomputing it on a dense grid
# and using np.interp replaces one scipy root() call per MC sample with a
# trivial array lookup — the dominant cost in the original code.
print("Precomputing tau_m(p) lookup table ...", end=" ", flush=True)
_TAU_M_P    = np.linspace(0.5, 6.0, 50_000)
_TAU_M_VALS = np.array([root(optical_depth_to_be_minimised, 1.0, args=(p_,)).x[0]
                         for p_ in _TAU_M_P])
print(f"done  (p range [{_TAU_M_P[0]:.2f}, {_TAU_M_P[-1]:.2f}], "
      f"tau_m range [{_TAU_M_VALS.min():.4f}, {_TAU_M_VALS.max():.4f}]).")

def tau_m_interp(p):
    """Array-safe O(1) lookup for the synchrotron optical-depth peak tau_m(p)."""
    return np.interp(p, _TAU_M_P, _TAU_M_VALS)


def plot_tau_m_lookup(plots_dir: str = "", output_stem: str = "mc_out",
                      show_plot: bool = False, n_check: int = 25):
    """
    Diagnostic plot of the tau_m(p) precomputed lookup table.

    Panels
    ------
    Top  : tau_m(p) curve over the full grid; n_check exact root() solutions
           overlaid as scatter points to verify interpolation accuracy.
    Bottom : residual  tau_m_interp(p) − tau_m_exact(p)  at the check points.

    Saved to <plots_dir>/<stem>_tau_m_lookup.pdf.
    """
    import os

    # Sparse check points spread across the grid
    p_check  = np.linspace(_TAU_M_P[0], _TAU_M_P[-1], n_check)
    tau_exact = np.array([root(optical_depth_to_be_minimised, 1.0, args=(p_,)).x[0]
                          for p_ in p_check])
    tau_interp = tau_m_interp(p_check)
    residuals  = tau_interp - tau_exact

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(9, 7),
                                          gridspec_kw={"height_ratios": [3, 1]})
    fig.suptitle(r"$\tau_m(p)$ precomputed lookup table — diagnostic", fontsize=13)

    # --- Top panel: lookup curve + exact check points ---
    ax_top.plot(_TAU_M_P, _TAU_M_VALS, color="C0", lw=1.5, label="np.interp grid")
    ax_top.scatter(p_check, tau_exact, color="C1", zorder=5, s=30,
                   label=f"exact root()  (n={n_check})")
    ax_top.set_ylabel(r"$\tau_m$")
    ax_top.set_xlim(_TAU_M_P[0], _TAU_M_P[-1])
    ax_top.legend(fontsize=10)
    ax_top.set_title(
        f"Grid: {len(_TAU_M_P):,} points,  "
        r"$p \in$" + f"[{_TAU_M_P[0]:.2f}, {_TAU_M_P[-1]:.2f}],  "
        r"$\tau_m \in$" + f"[{_TAU_M_VALS.min():.4f}, {_TAU_M_VALS.max():.4f}]",
        fontsize=9
    )

    # --- Bottom panel: interpolation residuals ---
    ax_bot.axhline(0, color="k", lw=0.8)
    ax_bot.scatter(p_check, residuals, color="C2", s=30, zorder=5)
    ax_bot.set_xlabel(r"$p = 1 - 2\alpha_{\rm thin}$")
    ax_bot.set_ylabel(r"$\Delta\tau_m$")
    ax_bot.set_xlim(_TAU_M_P[0], _TAU_M_P[-1])
    ax_bot.set_title(
        f"Max |residual| = {np.abs(residuals).max():.2e}", fontsize=9
    )

    fig.tight_layout()

    dir_ = plots_dir.strip() if (plots_dir and plots_dir.strip()) else os.getcwd()
    os.makedirs(dir_, exist_ok=True)
    path = os.path.join(dir_, f"{output_stem}_tau_m_lookup.pdf")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved tau_m lookup diagnostic: {path}")

    if show_plot:
        plt.show()
    else:
        plt.close(fig)


# --- k-function cache --------------------------------------------------------
# k1_func and k3_e/nu_func all call c5_func / c6_func (each involving two
# special.gamma evaluations).  Within a single physics-function call the same p
# is passed to several k-functions; a 1-slot cache means c5/c6 are evaluated
# exactly once per unique p value.  Scalars are keyed by value (helps inside
# the gamma_min root-solver loop where p is fixed); arrays are keyed by
# object identity (same variable reused across k-function calls in one
# vectorised physics-function evaluation).
_KFUNC_LAST_P    = None
_KFUNC_LAST_VALS = None

def _k_funcs(p):
    global _KFUNC_LAST_P, _KFUNC_LAST_VALS
    p_key = float(p) if np.ndim(p) == 0 else id(p)
    if _KFUNC_LAST_P == p_key:
        return _KFUNC_LAST_VALS
    c5 = c5_func(p)
    c6 = c6_func(p)
    k1    = (np.pi * c5/c6)**2 * (2*c1)**(-5) * (1 - np.exp(-1))**2
    k3_e  = (11/(2*(p+1)) * 1/(8*np.pi) * 4/3 * c6)**(-1/(1+2*(p+6))) * (2*c1)**(-(p+4)/(2+4*(p+6)))
    k3_nu = (2*c1)**(-(p+4)/34) * (1/(8*np.pi) * 4/3 * 11/6 * c6)**(-1/17)
    vals  = (k1, k3_e, k3_nu)
    _KFUNC_LAST_P, _KFUNC_LAST_VALS = p_key, vals
    return vals

def nu1_from_numax(numax, p, tau_max):
    return numax * tau_max**(2/(p+4))

def fnu1_from_fmax(fmax, p, nu_max, nu1):
    factor = (1-np.exp(-1))/((nu_max/nu1)**(5/2) * (1-np.exp(-(nu_max/nu1)**(-(p+4)/2))))
    return fmax * factor

def gamma_to_energy(gamma):
    return electron_mass_cgs * c_cgs**2 * gamma

def kpc_to_cm(dist_kpc):
    return dist_kpc * 3.086e21

def nu_crit(electron_energy,B):
    return c1 * B * electron_energy**2


######################## Base quantities - energy form - checked against previous equations - to be checked against simulated spectra

#Total energy in the magnetic field and non-thermal particles (including protons) of the emitting region
def E_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    
    deviation_factor = (2*(1+p))/(1+2*(p+6)) * equip_deviation**(11/(1+2*(p+6))) + 11/(1+2*(p+6)) * equip_deviation**(-(2*(1+p))/(1+2*(p+6)))
    
    energy = ( (1+2*(p+6))/(2*(p+1)) * 4/3 * 1/8 * k3_e_func(p)**11 * k1_func(p)**(-(3*p+14)/(2*(1+2*(p+6)))) * K_E_func(p, Emin, Emax)**(11/(1+2*(p+6))) * 
            nu1**(-1) * dist_cm**((4+6*(p+4))/(1+2*(p+6))) * fnu1**((2+3*(p+4))/(1+2*(p+6))) * doppler**(-(1+7*(p+4))/(1+2*(p+6))) * (1+redshift)**(-(2+5*(p+5))) * 
            (1+proton_energy_ratio)**(11/(1+2*(p+6))) * deviation_factor)
    
    if log10:
        return np.log10(energy)
    else:
        return energy

#Size of (quasi-spherical) emitting region
def R_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    size = ( k3_e_func(p) * k1_func(p)**(-(2*(p+6))/(4*(1+2*(p+6)))) * K_E_func(p, Emin, Emax)**(1/(1+2*(p+6))) * nu1**(-1) * 
            dist_cm**(2*(p+6)/(1+2*(p+6))) * fnu1**((p+6)/(1+2*(p+6))) * doppler**(-(p+5)/(1+2*(p+6))) * (1+redshift)**(-(1+3*(p+6))/(1+2*(p+6))) * 
            (1+proton_energy_ratio)**(1/(1+2*(p+6))) * equip_deviation**(1/(1+2*(p+6))) ) 
    
    if log10:
        return np.log10(size)
    else:
        return size

#Magnetic field in emitting region
def B_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    B = (  k3_e_func(p)**(4) * k1_func(p)**(1/(1+2*(p+6))) * K_E_func(p, Emin, Emax)**(4/(1+2*(p+6))) * nu1 * 
         dist_cm**(-4/(1+2*(p+6))) * fnu1**(-2/(1+2*(p+6))) * doppler**(-(2*p+7)/(1+2*(p+6))) * 
         (1+redshift)**((1+2*(p+7))/(1+2*(p+6))) * (1+proton_energy_ratio)**(4/(1+2*(p+6))) * equip_deviation**(4/(1+2*(p+6))) )
    
    if log10:
        return np.log10(B)
    else:
        return B
    
######################## Base quantities - frequency form

#CHECK ALL FREQUENCY FORMS    
#Total energy in the magnetic field and non-thermal particles (including protons) of the emitting region
def E_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    
    deviation_factor = (6/17 * equip_deviation**(11/17) + 11/17 * equip_deviation**(-6/17))
    
    energy = ( 17/36 * k1_func(p)**(-10/17) * k3_nu_func(p)**(11) * K_nu_func(p, nu_min, nu_max)**(11/17) * 
              fnu1**(20/17) * dist_cm**(40/17) * nu1**((11*p - 56)/34) * doppler**(-(64+11*p)/34) * 
              (1+redshift)**((11*p + 96)/34) * (1+proton_energy_ratio)**(11/17) * deviation_factor )
    
    if log10:
        return np.log10(energy)
    else:
        return energy
    
def R_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    size = ( k1_func(p)**(-4/17) * k3_nu_func(p) * K_nu_func(p, nu_min, nu_max)**(1/17) * 
            fnu1**(8/17) * dist_cm**(16/17) * nu1**((p-36)/34) * doppler**(-(12+p)/34) * 
            (1+redshift)**((p-52)/34) * (1+proton_energy_ratio)**(1/17) * equip_deviation**(1/17) )
    
    if log10:
        return np.log10(size)
    else:
        return size
    
def B_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    B = ( k1_func(p)**(1/17) * k3_nu_func(p)**(4) * K_nu_func(p, nu_min, nu_max)**(4/17) * 
            fnu1**(-2/17) * dist_cm**(-4/17) * nu1**((2*p+13)/17) * doppler**(-(2*p+7)/17) * 
            (1+redshift)**((2*p+15)/17) * (1+proton_energy_ratio)**(4/17) * equip_deviation**(4/17) )
    
    if log10:
        return np.log10(B)
    else:
        return B

######################## Dervied quantities - energy form

#Normalisation for electron energy distribution power law
def N0_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    N0 = ( 11/(16*np.pi*(p+1)) * k3_e_func(p)**(8) * k1_func(p)**(2/(1+2*(p+6))) * K_E_func(p, Emin, Emax)**(-(5+2*p)/(1+2*(p+6))) * nu1**2 *
          dist_cm**(-8/(1+2*(p+6))) * fnu1**(-4/(1+2*(p+6))) * doppler**(-(4*p+14)/(1+2*(p+6))) * (1+redshift)**((2+4*(p+7))/(1+2*(p+6))) * 
          (1+proton_energy_ratio)**(-(5+2*p)/(1+2*(p+6))) * equip_deviation**(-(5+2*p)/(1+2*(p+6))) )
    
    if log10:
        return np.log10(N0)
    else:
        return N0
    
#NUmber density of non-thermal electrons
def ne_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    ne =  N0_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=False) * K_E_dens_func(p, Emin, Emax)

    if log10:
        return np.log10(ne)
    else:
        return ne
    
#Total number of non-thermal electrons in emitting volume
def Ne_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    Ne = (  11/(12*(p+1)) * K_E_dens_func(p, Emin, Emax) * k3_e_func(p)**11 * k1_func(p)**(-(14+3*p)/(2*(1+2*(p+6)))) * K_E_func(p, Emin, Emax)**(-(2+2*p)/(1+2*(p+6))) * 
          nu1**(-1) * dist_cm**((4+6*(p+4))/(1+2*(p+6))) * fnu1**((2+3*(p+4))/(1+2*(p+6))) * doppler**(-(1+7*(p+4))/(1+2*(p+6))) * 
          (1+redshift)**(-(2+5*(p+5))/(1+2*(p+6))) * (1+proton_energy_ratio)**(-(2+2*p)/(1+2*(p+6))) * equip_deviation**(-(2+2*p)/(1+2*(p+6)))  )
    
    if log10:
        return np.log10(Ne)
    else:
        return Ne

#Brightness temperature of emitting region
def TB_energy_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    const = c_cgs**2 / (2*np.pi*kB_cgs)
    
    TB = const * (  k3_e_func(p)**(-2) * k1_func(p)**((p+6)/(1+2*(p+6))) * K_E_func(p, Emin, Emax)**(-2/(1+2*(p+6))) * 
          dist_cm**(2/(1+2*(p+6))) * fnu1**(1/(1+2*(p+6))) * doppler**(-3/(1+2*(p+6))) * (1+redshift)**(-1/(1+2*(p+6))) * (1+proton_energy_ratio)**(-2/(1+2*(p+6))) * equip_deviation**(-2/(1+2*(p+6))) )
    
    if log10:
        return np.log10(TB)
    else:
        return TB
    
#################### Gamma min equations - energy form

def _newton_vec(func, x0, args, maxiter=150, tol=1e-9, max_step=2.0):
    """Vectorised Newton with finite-difference derivative and step clamping.

    Clamping |Δx| ≤ max_step keeps iterations away from extreme logit values
    (gamma_min → 0 or gamma_min → gamma_max) where the physics functions
    produce Inf/NaN and scipy's unclamped secant method can diverge.

    Returns (x, converged) where converged is a boolean mask indicating which
    elements reached the tolerance.  Non-converged elements should be treated
    as invalid (set to NaN) so mc_uncertainty discards them.
    """
    x = x0.copy()
    h = 1e-5 * np.maximum(np.abs(x), 1.0)
    converged = np.zeros(len(x), dtype=bool)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        for _ in range(maxiter):
            f0 = func(x,     *args)
            f1 = func(x + h, *args)
            df = (f1 - f0) / h
            safe_df = np.where(np.isfinite(df) & (df != 0.0), df, 1e-30)
            step = np.clip(f0 / safe_df, -max_step, max_step)
            x -= step
            converged |= (np.abs(step) < tol)
            if np.all(converged):
                break
    return x, converged

#Function to be minimised numerically to find lower limit of gamma_min
def gamma_min_func_to_minimise(gamma_min, flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation):
    
    p = p_from_alpha(alpha_thin)
    Emin = gamma_to_energy(gamma_min)
    Emax = gamma_to_energy(gamma_max)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    factor = (  c1**(-1) * k3_e_func(p)**(-4) * k1_func(p)**(-1/(1+2*(p+6))) * K_E_func(p, Emin, Emax)**(-4/(1+2*(p+6))) * 
              dist_cm**(4/(1+2*(p+6))) * fnu1**(2/(1+2*(p+6))) * doppler**(-6/(1+2*(p+6))) *  (1+redshift)**(-2/(1+2*(p+6))) * (1+proton_energy_ratio)**(-4/(1+2*(p+6))) * equip_deviation**(-4/(1+2*(p+6))) )
        
    
    return (electron_mass_cgs**2 * c_cgs**4 * gamma_min**2 - factor)


def gamma_min_func_to_minimise_x(x, flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation):
    
    gamma_min = gamma_max * expit(x)  # always 0 < gamma_min < gamma_max
    
    return gamma_min_func_to_minimise(
        gamma_min, flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin,
        gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio,
        equip_deviation
    )

# Lower limit on gamma_min — gamma_min arg is passed by the MC framework but
# not used in this calculation (the self-consistent value is solved for).
def gamma_min_constraint(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_min, gamma_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, init=10, log10=True):
    args = (flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, gamma_max,
            bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation)

    # Initial guess in logit-space: logit(gamma_min / gamma_max).
    # The sampled gamma_min is physically plausible, so this lands the
    # starting point near the true root rather than at gamma_max/2 (x=0).
    x0_scalar = float(np.log(gamma_min / (gamma_max - gamma_min))) \
        if np.ndim(gamma_min) == 0 else np.log(gamma_min / (gamma_max - gamma_min))

    if np.ndim(flux_dens_peak_mJy) == 0:
        # Scalar path — original robust solver (used in the loop fallback and
        # for standalone calls).
        x_sol = root(gamma_min_func_to_minimise_x, x0_scalar, args=args).x[0]
        gamma_min_upper_lim = gamma_max * expit(x_sol)
    else:
        # Vectorised path — step-clamped Newton so iterations can't overshoot
        # into extreme logit values (gamma_min → 0) where the physics functions
        # produce Inf/NaN and scipy's unclamped secant method diverges.
        x_sol, converged = _newton_vec(gamma_min_func_to_minimise_x, x0_scalar, args=args)
        gamma_min_upper_lim = gamma_max * expit(x_sol)
        # Non-converged samples get NaN so mc_uncertainty discards them cleanly
        # rather than returning gamma_min ≈ 0 (expit of very negative x).
        gamma_min_upper_lim = np.where(converged, gamma_min_upper_lim, np.nan)

    if log10:
        return np.log10(gamma_min_upper_lim)
    else:
        return gamma_min_upper_lim
    
######################## Derived quantities - frequency form

#Normalisation for electron energy distribution power law
def N0_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    N0 = ( 11/(48*np.pi) * (k1_func(p)**(1/17) * k3_nu_func(p)**(4) * fnu1**(-2/17) * dist_cm**(-4/17) * nu1**((2*p+13)/17) * doppler**(-(2*p+7)/17) * 
          (1+redshift)**((2*p+15)/17) )**((6-p)/2) * 
          (1+proton_energy_ratio)**(-(5+2*p)/17) * equip_deviation**(-(5+2*p)/17) * K_nu_func(p, nu_min, nu_max)**(-(5+2*p)/17) )
    
    if log10:
        return np.log10(N0)
    else:
        return N0
    
def ne_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    ne = ( 11/(48*np.pi) * K_nu_dens_func(p, nu_min, nu_max) * k1_func(p)**(5/34) * k3_nu_func(p)**(10) * K_nu_func(p, nu_min, nu_max)**(-7/17) * 
          fnu1**(-10/34) * dist_cm**(-20/34) * nu1**((10*p + 65)/34) * doppler**(-(10*p+35)/34) * (1+redshift)**((10*p+75)/34) * 
          (1+proton_energy_ratio)**(-7/17) * equip_deviation**(-7/17) )
    
    if log10:
        return np.log10(ne)
    else:
        return ne
    
def Ne_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    Ne = ( 11/36 * K_nu_dens_func(p, nu_min, nu_max) * k1_func(p)**(-19/34) * k3_nu_func(p)**(13) * K_nu_func(p, nu_min, nu_max)**(-4/17) * 
          fnu1**(19/17) * dist_cm**(38/17) * nu1**((13*p - 43)/34) * doppler**(-(13*p+71)) * (1+redshift)**((13*p-81)/34) * 
          (1+proton_energy_ratio)**(-4/17) * equip_deviation**(-4/17)  )
    
    if log10:
        return np.log10(Ne)
    else:
        return Ne
    
def TB_frequency_form(flux_dens_peak_mJy, nu_obs_Hz, dist_kpc, alpha_thin, nu_min, nu_max, bulk_gamma, cos_inclination, redshift, proton_energy_ratio, equip_deviation, log10=True):
    p = p_from_alpha(alpha_thin)
    
    tau_m = tau_m_interp(p)
    nu1 = nu1_from_numax(nu_obs_Hz, p, tau_m)
    
    flux_dens_peak_cgs = flux_dens_peak_mJy * 1e-26
    
    fnu1 = fnu1_from_fmax(flux_dens_peak_cgs, p, nu_obs_Hz, nu1)
    
    dist_cm = kpc_to_cm(dist_kpc)
    
    inclination_deg = np.arccos(cos_inclination) * 180/np.pi
    doppler = doppler_func(bulk_gamma, inclination_deg)
    
    const = c_cgs**2 / (2*np.pi*kB_cgs)
    
    TB = const * ( k1_func(p)**(8/17) * k3_nu_func(p)**(-2) * K_nu_func(p, nu_min, nu_max)**(-2/17) * fnu1**(1/17) * 
                  dist_cm**(2/17) ** nu1**((2-p)/17) * doppler**((p-5)/17) * (1+redshift)**((1-p)/17) * 
                  (1+proton_energy_ratio)**(-2/17) * equip_deviation**(-2/17) )
    
    if log10:
        return np.log10(TB)
    else:
        return TB

########################


def resolve_input_yaml_path():
    if isinstance(INPUT_YAML, str) and INPUT_YAML.strip():
        input_yaml = INPUT_YAML.strip()
        print(f"Using INPUT_YAML override: {input_yaml}")
        return input_yaml

    cwd = os.getcwd()
    pattern = os.path.join(cwd, "*inputs.yaml")
    matches = sorted(glob.glob(pattern))

    if not matches:
        raise FileNotFoundError(
            f"No inputs.yaml found in CWD: {cwd}. Set INPUT_YAML to a YAML file path."
        )

    if len(matches) > 1:
        print(
            "Multiple inputs.yaml matches found in CWD. "
            f"Using first match: {matches[0]}"
        )
    else:
        print(f"Found inputs.yaml in CWD: {matches[0]}")

    return matches[0]


input_yaml_path = resolve_input_yaml_path()
with open(input_yaml_path, "r", encoding="utf-8") as f:
    raw = yaml.safe_load(f)

# Separate run-level settings from sampling parameters
globals_cfg = raw.pop("globals", {})
params = raw

_cwd         = os.getcwd()
_plots_raw   = str(globals_cfg.get("plots_path",   "")).strip()
_samples_raw = str(globals_cfg.get("samples_path", "")).strip()
plots_dir    = _plots_raw   if _plots_raw   else os.path.join(_cwd, "plots")
samples_dir  = _samples_raw if _samples_raw else os.path.join(_cwd, "samples")
output_stem       = str(globals_cfg.get("output_stem", "mc_out"))
show_plot         = bool(globals_cfg.get("show_plot", False))
_padding_raw      = float(globals_cfg.get("gamma_min_padding", 1.0))
gamma_min_padding = max(1.0, _padding_raw)
joint_sampling    = bool(globals_cfg.get("joint_sampling", False))
joint_seed        = int(globals_cfg.get("joint_seed", 100))

os.makedirs(plots_dir,   exist_ok=True)
os.makedirs(samples_dir, exist_ok=True)

print(f"\n{'=' * 62}")
print(f"  Config loaded      : {input_yaml_path}")
print(f"  Parameters         : {len(params)}")
print(f"  plots_dir          : {plots_dir}")
print(f"  samples_dir        : {samples_dir}")
print(f"  output_stem        : {output_stem!r}")
print(f"  show_plot          : {show_plot}")
print(f"  gamma_min_padding  : {gamma_min_padding}"
      + ("  (pinned from {:.3g})".format(_padding_raw) if _padding_raw < 1.0 else ""))
print(f"  joint_sampling     : {joint_sampling}"
      + (f"  (seed = {joint_seed})" if joint_sampling else ""))
print(f"{'=' * 62}\n")

# Shared keyword arguments forwarded to every plot_mc_diagnostics call
plot_kwargs = dict(plots_dir=plots_dir, output_stem=output_stem, show_plot=show_plot)

plot_tau_m_lookup(plots_dir=plots_dir, output_stem=output_stem, show_plot=show_plot)


def _update_gamma_min(params, res_gmin, padding):
    """Extract posterior max from gamma_min_constraint and update params in-place."""
    log10_max  = float(np.max(res_gmin.samples))
    new_high   = 10.0**log10_max * padding
    print(f"\n{'─' * 62}")
    print(f"  gamma_min upper limit from posterior maximum:")
    print(f"    log10(gamma_min_max) = {log10_max:.4f}")
    print(f"    linear max (pre-padding)        = {10.0**log10_max:.4g}")
    print(f"    gamma_min_high (x{padding:.3g} padding) = {new_high:.4g}")
    print(f"  Updating params['gamma_min']['high'] → {new_high:.4g}")
    print(f"{'─' * 62}\n")
    params["gamma_min"]["high"] = new_high
    return new_high


# =============================================================================
if joint_sampling:
# =============================================================================
    N = 1_000_000

    # -- 1. Initial master draw -----------------------------------------------
    draws = draw_params(params, n=N, seed=joint_seed)

    # -- 2. gamma_min constraint on shared draw --------------------------------
    res_gmin = mc_uncertainty(gamma_min_constraint, params, ci_level=0.68,
                              fixed_kwargs={"log10": True}, fixed_draws=draws)
    plot_mc_diagnostics(res_gmin, params, **plot_kwargs)
    save_mc_samples(res_gmin, samples_dir=samples_dir, output_stem=output_stem)

    # -- 3. Per-sample gamma_min: draw uniformly up to the per-sample upper limit
    #
    # gamma_min_constraint returns the UPPER LIMIT on gamma_min for each sample.
    # We call it directly (not via mc_uncertainty) so we get the full N-length
    # array with NaN for non-converged samples, rather than the already-filtered
    # .samples array from MCResult (which would cause a shape mismatch).
    #
    # We then filter ALL draws to the valid subset so every array is consistently
    # sized — NaN gamma_min samples cannot contribute to any calculation anyway.
    #
    # OLD (incorrect — commented out for comparison):
    # _update_gamma_min(params, res_gmin, gamma_min_padding)   # global max only
    # draws["gamma_min"] = draw_params(params, n=N, seed=joint_seed + 1,
    #                                  subset=["gamma_min"])["gamma_min"]  # independent redraw

    gmin_low = params["gamma_min"]["low"]
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        gmin_upper = gamma_min_constraint(**draws, log10=False)  # shape (N,), NaN where non-converged

    valid_gmin = np.isfinite(gmin_upper)
    valid_gmin &= (gmin_upper > gmin_low)  # drop samples where SSA upper limit < prior low
    rng_gmin = np.random.default_rng(joint_seed + 1)
    gmin_draws = rng_gmin.uniform(gmin_low, np.where(valid_gmin, gmin_upper, gmin_low))
    draws["gamma_min"] = np.where(valid_gmin, gmin_draws, np.nan)

    # Filter ALL draws to the valid subset so every array has the same length.
    # This avoids shape mismatches in the sample-by-sample fallback loop inside
    # mc_uncertainty when any draw contains NaN.
    draws = {k: v[valid_gmin] for k, v in draws.items()}

    # Remove gamma_min from params so mc_uncertainty no longer prints it as a
    # prior — it is now a solved quantity passed via fixed_draws, not sampled.
    params.pop("gamma_min", None)

    n_valid = valid_gmin.sum()
    print(f"\n{'─' * 62}")
    print(f"  gamma_min: per-sample uniform draw up to self-consistent upper limit")
    print(f"    valid samples    : {n_valid:,} / {N:,} ({100*n_valid/N:.1f}%)")
    print(f"    upper limit range: [{gmin_upper[valid_gmin].min():.2f}, {gmin_upper[valid_gmin].max():.2f}]")
    print(f"    gamma_min range  : [{draws['gamma_min'].min():.2f}, {draws['gamma_min'].max():.2f}]")
    print(f"{'─' * 62}\n")

    # -- 4. Plot shared input parameters once (with refined gamma_min) --------
    plot_mc_input_params(draws, params, plots_dir=plots_dir,
                         output_stem=output_stem, show_plot=show_plot)

    # -- 4b. Save dist_kpc (D) as its own .npy --------------------------------
    _d_path = os.path.join(samples_dir, f"{output_stem}_dist_kpc_samples.npy")
    np.save(_d_path, draws["dist_kpc"])
    print(f"Saved dist_kpc samples:      {_d_path}")

    # -- 5. Evaluate E, R, n_e, B on the shared draws -------------------------
    res_E  = mc_uncertainty(E_energy_form,  params, ci_level=0.68,
                            fixed_kwargs={"log10": True}, fixed_draws=draws)
    res_R  = mc_uncertainty(R_energy_form,  params, ci_level=0.68,
                            fixed_kwargs={"log10": True}, fixed_draws=draws)
    res_ne = mc_uncertainty(ne_energy_form, params, ci_level=0.68,
                            fixed_kwargs={"log10": True}, fixed_draws=draws)
    res_B  = mc_uncertainty(B_energy_form,  params, ci_level=0.68,
                            fixed_kwargs={"log10": True}, fixed_draws=draws)

    # -- 6. Per-component posterior plots (no redundant params grid) ----------
    for _res in [res_E, res_R, res_ne, res_B]:
        plot_mc_diagnostics(_res, params, show_params=False, **plot_kwargs)
        save_mc_samples(_res, samples_dir=samples_dir, output_stem=output_stem)

    # -- 7. Corner plot of [E, R, n_e, B] -------------------------------------
    # Compute a jointly-valid chain so correlations are exact.
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        _yE  = E_energy_form( **draws, log10=True)
        _yR  = R_energy_form( **draws, log10=True)
        _yne = ne_energy_form(**draws, log10=True)
        _yB  = B_energy_form( **draws, log10=True)

    _jmask = np.isfinite(_yE) & np.isfinite(_yR) & np.isfinite(_yne) & np.isfinite(_yB)
    chain  = np.column_stack([_yE[_jmask], _yR[_jmask], _yne[_jmask], _yB[_jmask]])

    _corner_labels = [
        r"$\log_{10}(E\ /\ \mathrm{erg})$",
        r"$\log_{10}(R\ /\ \mathrm{cm})$",
        r"$\log_{10}(n_e\ /\ \mathrm{cm}^{-3})$",
        r"$\log_{10}(B\ /\ \mathrm{G})$",
    ]
    plot_mc_corner(chain, labels=_corner_labels, plots_dir=plots_dir,
                   output_stem=output_stem, show_plot=show_plot)

    # -- 8. Pickle of everything -----------------------------------------------
    import pickle
    _pkl = {
        "B":         res_B.samples,
        "E":         res_E.samples,
        "R":         res_R.samples,
        "ne":        res_ne.samples,
        "gamma_min": res_gmin.samples,
        "draws":     draws,
        "results":   {"B": res_B, "E": res_E, "R": res_R,
                      "ne": res_ne, "gamma_min": res_gmin},
    }
    _pkl_path = os.path.join(samples_dir, f"{output_stem}_all.pkl")
    with open(_pkl_path, "wb") as _f:
        pickle.dump(_pkl, _f)
    print(f"Saved full results pkl:      {_pkl_path}")

# =============================================================================
else:
# =============================================================================
    # Sammples 
    n_samp = 1_000_000
    seed = 19950325

    # -- 1. gamma_min constraint ----------------------------------------------
    res = mc_uncertainty(gamma_min_constraint, params, n=n_samp, seed=seed, ci_level=0.68,
                         fixed_kwargs={"log10": True})
    plot_mc_diagnostics(res, params, **plot_kwargs)
    save_mc_samples(res, samples_dir=samples_dir, output_stem=output_stem)

    _update_gamma_min(params, res, gamma_min_padding)

    # -- 2. Independent runs with tightened gamma_min prior -------------------
    res = mc_uncertainty(E_energy_form, params, n=n_samp, seed=seed,
                         ci_level=0.68, fixed_kwargs={"log10": True})
    plot_mc_diagnostics(res, params, **plot_kwargs)
    save_mc_samples(res, samples_dir=samples_dir, output_stem=output_stem)

    res = mc_uncertainty(R_energy_form, params, n=n_samp, seed=seed,
                         ci_level=0.68, fixed_kwargs={"log10": True})
    plot_mc_diagnostics(res, params, **plot_kwargs)
    save_mc_samples(res, samples_dir=samples_dir, output_stem=output_stem)

    res = mc_uncertainty(B_energy_form, params, n=n_samp, seed=seed,
                         ci_level=0.68, fixed_kwargs={"log10": True})
    plot_mc_diagnostics(res, params, **plot_kwargs)
    save_mc_samples(res, samples_dir=samples_dir, output_stem=output_stem)

    res = mc_uncertainty(ne_energy_form, params, n=n_samp, seed=seed,
                         ci_level=0.68, fixed_kwargs={"log10": True})
    plot_mc_diagnostics(res, params, **plot_kwargs)
    save_mc_samples(res, samples_dir=samples_dir, output_stem=output_stem)
