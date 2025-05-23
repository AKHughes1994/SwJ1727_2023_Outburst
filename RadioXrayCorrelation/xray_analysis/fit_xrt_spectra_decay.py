#!/usr/bin/env python3
"""
Script to fit Swift XRT spectra for Swift J1727.8-1613 during its decay phase.

This script processes Swift XRT spectral data, performs spectral fitting using PyXspec,
and extracts key parameters for different spectral models (power-law, disk blackbody, 
and a combination of both). It also generates diagnostic plots for each fit and saves 
the results in a JSON file.

Authors:
- Andrew Hughes (andrew.hughes@physics.ox.ac.uk)

Dependencies:
- PyXspec (XSPEC Python bindings)
- Astropy
- NumPy
- Matplotlib
- Python 3.x

Features:
- Automatically identifies and processes spectral files.
- Fits spectra using different models based on observation conditions.
- Handles both chi-squared and Cash statistics for fitting.
- Saves fit results and uncertainties in a structured JSON file.
- Generates diagnostic plots for each spectrum.

Usage:
- Ensure the required spectral files are located in the `spectra/` directory.
- Run the script directly to process and fit the spectra.
- Results are saved in `results/SwiftJ1727_swift-xrt_decay.json`.

"""

import json
import glob
import os
import sys
import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.io import fits
from xspec import *

def msg(txt):
    """
    Print a timestamped message to the console.

    Parameters:
    - txt (str): The message to print.
    """
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt)

def plot(name):
    """
    Generate diagnostic plots for the spectral fit.

    Parameters:
    - name (str): The filename to save the plot.

    Returns:
    - None
    """
    Plot.device = '/null'
    Plot.xAxis = "KeV"  # Set x-axis to energy (keV)
    Plot('data')
    chans = Plot.x()
    rates = Plot.y()
    yer = Plot.yErr()
    folded = Plot.model()
    resids = np.array(rates) - np.array(folded)
    chi = resids / yer

    # Create a two-panel plot: data + model (top) and residuals (bottom)
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax[0].errorbar(chans, rates, yerr=yer, fmt='ro', ms=1, label="Data", elinewidth=0.2)
    ax[0].plot(chans, folded, 'b', label="Model", linewidth=1)
    ax[0].set_ylabel(r'Counts/cm$^2$/sec/chan')
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")

    ax[1].plot(chans, chi, 'g', linewidth=1)
    ax[1].axhline(0.0, ls=':', color='k')
    ax[1].set_ylim(-5, 5)
    ax[1].set_ylabel(r'$\chi$')
    ax[1].set_xlabel('Energy [keV]')

    plt.savefig(name)
    plt.clf()
    plt.close()

def initialize_model(mod_index, flux_guess, plaw_frac=0.5, fix=False):
    """
    Initialize the spectral model and parameters for fitting.

    Parameters:
    - mod_index (int): Model index (0: power-law, 1: disk blackbody, 2: combined).
    - flux_guess (float): Initial flux guess to avoid local minima.
    - plaw_frac (float): Fraction of flux assigned to the power-law component (default: 0.5).
    - fix (bool): Whether to fix certain parameters (default: False).

    Returns:
    - mod (str): XSPEC model string.
    - initial_pars (dict): Dictionary of initial parameter values.
    """
    # Define energy range and initial parameter values
    Emin = 1.0
    Emax = 10.0
    nh = '0.268 -1'  # Weighted-average NH after allowing it to vary freely
    gamma = '1.7,,0.0,0.0,4.0,4.0'
    Tin = '0.5,,0.05,0.05,5.0,5.0'

    if fix:
        gamma = '1.7 -1'
        Tin = '0.5 -1'
        nh = '0.268 -1'

    # Define models based on the index
    if mod_index == 0:  # Power-law
        mod = 'tbabs * pegpwrlw'
        initial_pars = {1: nh, 2: gamma, 3: Emin, 4: Emax, 5: f'{flux_guess}'}

    elif mod_index == 1:  # Disk blackbody
        bb_flux = np.log10(flux_guess * 1e-12)
        mod = 'tbabs * cflux * diskbb'
        initial_pars = {1: nh, 2: Emin, 3: Emax, 4: f'{bb_flux}', 5: Tin, 6: '1.0 -1'}

    elif mod_index == 2:  # Combined power-law + disk blackbody
        mod = 'tbabs * (pegpwrlw + diskbb)'
        initial_pars = {
            1: nh,
            2: gamma, 3: Emin, 4: Emax, 5: f'{flux_guess * plaw_frac}',
            6: Tin, 7: 3e3
        }

    else:
        sys.exit('Error: Invalid model index. Use 0 (power-law), 1 (diskbb), or 2 (combined).')

    return mod, initial_pars

def main():
    """
    Main function to process and fit Swift XRT spectra for Swift J1727.8-1613.

    This function performs the following steps:
    1. Initializes XSPEC settings and prepares the spectral fitting environment.
    2. Reads and filters spectral files based on observation time.
    3. Sorts the spectral files by observation time for sequential processing.
    4. Iterates through each spectrum to:
        - Perform initial power-law fits to estimate flux.
        - Fit models (power-law, disk blackbody, or both) based on observation conditions.
        - Extract fit parameters, uncertainties, and flux values.
        - Generate diagnostic plots for each fit.
    5. Handles errors during fitting by appending placeholder values (-1) for failed fits.
    6. Saves the results in a structured JSON file for further analysis.

    Outputs:
    - Diagnostic plots saved in the `pngs/` directory.
    - Results saved in `results/SwiftJ1727_swift-xrt_decay.json`.

    Dependencies:
    - XSPEC Python bindings (PyXspec)
    - Astropy for time handling
    - NumPy for array operations
    - Matplotlib for plotting
    """
    # Initialize a dictionary to store spectral parameters
    xrt_dict = {}
    model_names = {0: 'powerlaw', 1: 'diskbb', 2: 'both'}  # Model names for indexing
    skiplist = ['Emin', 'Emax', 'eMin', 'eMax']  # Parameters to skip during extraction
    Emin = 1.0  # Minimum energy for fitting (keV)
    Emax = 10.0  # Maximum energy for fitting (keV)

    # Initialize XSPEC settings
    Xset.abund = "wilm"    # Elemental abundances
    Xset.xsect = "vern"    # Cross-sections
    Fit.query = "yes"      # Default response to XSPEC queries
    Xset.chatter = 0       # Suppress XSPEC output chatter
    Fit.nIterations = 100  # Maximum fit iterations
    Plot.xAxis = "KeV"     # Set x-axis to energy (keV)
    Fit.statTest = "chi"   # Default fit statistic
    Plot.device = '/null'  # Disable XSPEC plotting
    Plot.yLog = True       # Use logarithmic y-axis for plots
    Xset.parallel.error = 10  # Parallelize error calculations

    # Read in the spectral files (Grade 0 spectra)
    spectra = glob.glob('spectra/*final.pi')
    obs_isots = []  # Observation times (ISO format)
    exposure = []
    bin_counts = []  # Binned counts
    tot_counts = []  # Total counts
    spectral_files = []  # List of spectral files

    # Filter and reorder spectra by observation time
    for spectrum in spectra:
        header = fits.getheader(spectrum, ext=1)

        # Get the timing information
        mjd_start = Time(header['DATE-OBS'], format='isot').mjd
        mjd_end = Time(header['DATE-END'], format='isot').mjd
        isot_mid = Time((mjd_start + mjd_end) / 2.0, format='mjd').isot

        # Only use spectra from the hard state (after MJD 60341)
        if Time(header['DATE-OBS'], format='isot').mjd > 60341:
            obs_isots.append(isot_mid)
            exposure.append(header['ONTIME'])
            bin_counts.append(header['COUNTGRP'])
            tot_counts.append(header['COUNTTOT'])
            spectral_files.append(spectrum)

    # Sort the spectra by observation time
    sort_index = np.argsort(obs_isots)
    obs_isots = np.array(obs_isots)[sort_index]
    exposure = np.array(exposure)[sort_index]
    spectral_files = np.array(spectral_files)[sort_index]
    bin_counts = np.array(bin_counts)[sort_index]
    tot_counts = np.array(tot_counts)[sort_index]
    obs_mjd = Time(obs_isots, format='isot').mjd  # Convert ISO times to MJD

    # Iterate through each spectrum to perform fitting
    for k, spectrum_file in enumerate(spectral_files):
        # Load the spectrum and set the energy range for fitting
        AllData(spectrum_file)
        AllData.ignore('*:10.0-**')  # Ignore energies above 10 keV
        AllData.ignore('*:**-0.5')  # Ignore energies below 0.5 keV
        AllData.notice('*:0.5-10.0')  # Focus on the 0.5-10 keV range
        AllData.ignore('bad')  # Ignore bad data points

        # Determine the fit statistic (chi-squared or Cash) based on bin counts
        if bin_counts[k] == 1:
            Fit.statMethod = 'cstat'  # Use Cash statistic for low counts
        else:
            Fit.statMethod = 'chi'  # Use chi-squared for sufficient counts

        # Perform an initial power-law fit to estimate the flux
        Model('tbabs * pegpwrlw')
        AllModels(1).setPars({1: '0.268 -1', 2: '2.0 -1', 3: Emin, 4: Emax, 5: '1000.0'})
        Fit.delta = 1e-2
        Fit.perform()
        flux_guess = AllModels(1).pegpwrlw.norm.values[0]  # Initial flux estimate

        # Determine which models to fit based on fit stats
        if Fit.statMethod == 'cstat':
            models = [0]  # Use power-law for low counts
        else:
            models = [2]  # Use combined model for sufficient counts

        # Iterate through the selected models
        for mod_index in models:
            # Adjust initial power-law fraction based on observation time
            if obs_mjd[k] > 60400:
                plaw_frac = 0.95  # Hard state: high power-law fraction
            else:
                plaw_frac = 0.05  # Soft state: low power-law fraction

            # Initialize the model and systematics
            mod, initial_pars = initialize_model(mod_index, flux_guess, plaw_frac)
            Model(mod)
            AllModels(1).setPars(initial_pars)
            AllModels.systematic = 0.03  # Apply 3% systematic uncertainty

            # Perform the fit with decreasing delta values for convergence
            for delt in [1e-2, 1e-3, 1e-4, 1e-5]:
                Fit.delta = delt
                Fit.perform()

            # Initialize model dictionary if it doesn't exist
            mname = model_names[mod_index]
            if mname not in xrt_dict.keys():
                xrt_dict[mname] = {'chi2': [], 'dof': [], 'redchi2': [], 'obs_isot': [], 'obs_mjd': [], 'exposure':[], 'fitstat': []}

                # Add parameter keys for each model component
                for comp in AllModels(1).componentNames:
                    component = getattr(AllModels(1), comp)
                    for par in component.parameterNames:
                        if par not in skiplist:
                            if par == 'norm' and comp == 'diskbb':
                                pass
                            else:
                                xrt_dict[mname][par] = []
                                xrt_dict[mname][par + '_neg'] = []
                                xrt_dict[mname][par + '_pos'] = []

                if mname == 'both':
                    xrt_dict[mname]['lg10Flux'] = []
                    xrt_dict[mname]['lg10Flux_neg'] = []
                    xrt_dict[mname]['lg10Flux_pos'] = []

            # Append fit statistics and observation details
            xrt_dict[mname]['obs_isot'].append(obs_isots[k])
            xrt_dict[mname]['obs_mjd'].append(obs_mjd[k])
            xrt_dict[mname]['exposure'].append(exposure[k])
            xrt_dict[mname]['chi2'].append(Fit.testStatistic)
            xrt_dict[mname]['dof'].append(Fit.dof)
            xrt_dict[mname]['redchi2'].append(Fit.testStatistic / Fit.dof)
            xrt_dict[mname]['fitstat'].append(Fit.statMethod)

            # Generate diagnostic plots
            msg(f'Fitting {spectrum_file} ({Fit.statMethod}, {plaw_frac}) on {obs_isots[k]} (MJD {obs_mjd[k]}) with chi2(dof) = {Fit.testStatistic} ({Fit.dof})')
            # Define the directory name
            directory_name = "pngs"

            # Check if the directory exists
            if not os.path.exists(directory_name):
                # If it doesn't exist, create the directory
                os.makedirs(directory_name)
            plot_name = spectrum_file.replace('spectra', 'pngs').replace('spectra', 'pngs') + f'_mod_{mod_index}.png'
            plot(plot_name)

            # Try to extract flux uncertainties
            try:
                nParameters = AllModels(1).nParameters
                Fit.error(f'maximum 3.5 1.0 1-{nParameters}')
                pars = []
                for comp in AllModels(1).componentNames:
                    component = getattr(AllModels(1), comp)
                    for par in component.parameterNames:
                        if par not in skiplist:
                            parameter = getattr(component, par)
                            pars.append(parameter.values[0])
                            if par == 'norm' and comp == 'diskbb':
                                pass
                            else:
                                xrt_dict[mname][par].append(parameter.values[0])
                                xrt_dict[mname][par + '_neg'].append(abs(parameter.values[0] - parameter.error[0]))
                                xrt_dict[mname][par + '_pos'].append(abs(parameter.error[1] - parameter.values[0]))

                # Use CFLUX to calculate total flux for two-component models
                if mname == 'both':
                    Model('tbabs * cflux * (pegpwrlw + diskbb)')
                    initial_pars = {
                        1: f'{pars[0]} -1',
                        2: Emin, 3: Emax, 4: np.log10(flux_guess * 1e-12),
                        5: f'{pars[1]} -1', 6: Emin, 7: Emax, 8: f'{pars[2]}',
                        9: f'{pars[3]} -1', 10: f'{pars[4]} -1'
                    }
                    AllModels(1).setPars(initial_pars)
                    Fit.perform()
                    nParameters = AllModels(1).nParameters
                    Fit.error(f'maximum 3.5 1.0 1-{nParameters}')
                    xrt_dict[mname]['lg10Flux'].append(AllModels(1).cflux.lg10Flux.values[0])
                    xrt_dict[mname]['lg10Flux_neg'].append(AllModels(1).cflux.lg10Flux.values[0] - AllModels(1).cflux.lg10Flux.error[0])
                    xrt_dict[mname]['lg10Flux_pos'].append(AllModels(1).cflux.lg10Flux.error[1] - AllModels(1).cflux.lg10Flux.values[0])

            except ValueError:
                # Handle fitting failures by appending -1 for all parameters
                for comp in AllModels(1).componentNames:
                    component = getattr(AllModels(1), comp)
                    for par in component.parameterNames:
                        if par not in skiplist:
                            if par == 'norm' and comp == 'diskbb':
                                pass
                            else:
                                xrt_dict[mname][par].append(-1)
                                xrt_dict[mname][par + '_neg'].append(-1)
                                xrt_dict[mname][par + '_pos'].append(-1)

                if mname == 'both':
                    xrt_dict[mname]['lg10Flux'].append(-1)
                    xrt_dict[mname]['lg10Flux_neg'].append(-1)
                    xrt_dict[mname]['lg10Flux_pos'].append(-1)

    # Ensure the results directory exists
    results_directory = "results"
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    
    # Save the final XRT dictionary to a JSON file
    with open(f'results/SwiftJ1727_swift-xrt_decay.json', 'w') as j:
        json.dump(xrt_dict, j, indent=4)

if __name__ in "__main__":
    main()
