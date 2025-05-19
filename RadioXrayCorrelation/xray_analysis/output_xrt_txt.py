"""
Script to process Swift XRT JSON data and save it as formatted text files.

This script reads X-ray data from JSON files for the rise and decay phases of Swift J1727.8-1613, 
processes the data to compute fluxes, errors, and systematic uncertainties, and saves the results 
in text files for further analysis. The output includes both linear and logarithmic flux values 
along with their respective errors.

Features:
- Converts JSON data into NumPy arrays for efficient processing.
- Computes X-ray flux and errors (statistical and systematic) in both linear and logarithmic scales.
- Handles asymmetric errors using difference-based calculations.
- Supports optional inclusion of upper limits for flux values.
- Outputs data in a sorted text file with a clear and descriptive header.

Inputs:
- JSON files containing X-ray data for rise and decay phases:
  - 'results/SwiftJ1727_swift-xrt_rise.json'
  - 'results/SwiftJ1727_swift-xrt_decay.json'

Outputs:
- Text files with processed X-ray data:
  - '../SwiftJ1727_Files/XRT_rise.txt'
  - '../SwiftJ1727_Files/XRT_decay.txt'

Output Columns:
1. Date (MJD): Observation date in Modified Julian Date.
2. Flux (erg/s/cm^2): X-ray flux in linear scale.
3. Stat Err neg (erg/s/cm^2): Lower bound of the statistical error in linear scale.
4. Stat Err pos (erg/s/cm^2): Upper bound of the statistical error in linear scale.
5. Sys Err neg (erg/s/cm^2): Lower bound of the error (statistical + systematic) in linear scale.
6. Sys Err pos (erg/s/cm^2): Upper bound of the error (statistical + systematic) in linear scale.
7. Log10 Flux: X-ray flux in logarithmic scale.
8. Log10 Stat Err neg: Lower bound of the statistical error in logarithmic scale.
9. Log10 Stat Err pos: Upper bound of the statistical error in logarithmic scale.
10. Log10 Sys Err neg: Lower bound of the error (statistical + systematic) in logarithmic scale.
11. Log10 Sys Err pos: Upper bound of the error (statistical + systematic) in logarithmic scale.
12. Upper limit (0 - False, 1 - True): Indicates whether the flux value is an upper limit.

Dependencies:
- NumPy for numerical operations.
- JSON for reading input data.

Author:
- Andrew Hughes (andrew.hughes@physics.ox.ac.uk)

Usage:
- Ensure the input JSON files are located in the specified paths.
- Run the script to process the data and generate the output text files.
"""

import numpy as np
import json

def convert_to_numpy(data):
    """
    Recursively convert lists in a nested dictionary or list to NumPy arrays.

    Parameters:
    - data: The input data (dictionary, list, or other types).

    Returns:
    - The data with all lists converted to NumPy arrays.
    """
    if isinstance(data, dict):
        return {key: convert_to_numpy(value) for key, value in data.items()}
    elif isinstance(data, list):
        return np.array(data)
    else:
        return data


def load_json(filename):
    """
    Load a JSON file and return its content.
    Args:
        filename (str): The path to the JSON file.
    Returns:
        dict: The content of the JSON file.
    """
    with open(filename, 'r') as f:
        data = json.load(f)
    return convert_to_numpy(data)


def write_to_txt(data, filename, upper_limit=False):
    """
    Write the data to a text file.
    Args:
        data (dict): The data to write.
        filename (str): The name of the output text file.
    """

    # Load the JSON data
    json_data = load_json(data)

    # Initialize arrays to hold the data
    mjd = np.array([]).reshape(1, 0)  # (1, N)
    flux = np.array([]).reshape(1, 0)  # (1, N)
    flux_err = np.array([]).reshape(2, 0)  # (2, N)
    log_flux = np.array([]).reshape(1, 0)  # (1, N)
    log_flux_err = np.array([]).reshape(2, 0)  # (2, N)

    # Fractional systematic error for the flux calculation
    sys_frac = 0.1

    # Iterate through the JSON models and append the data to the arrays
    for model in json_data.keys():
        # For ease of use, we will define a model dictionary (m)
        m = json_data[model]

        if model == 'both':  # cflux * (diskbb + pegpwrlaw)
            # Append the data to the arrays
            mjd = np.append(mjd, m['obs_mjd'].reshape(1, -1), axis=1)
            flux = np.append(flux, (10 ** m['lg10Flux']).reshape(1, -1), axis=1)
            flux_err_neg = 10 ** m['lg10Flux'] - 10 ** (m['lg10Flux'] - m['lg10Flux_neg'])
            flux_err_pos = 10 ** (m['lg10Flux'] + m['lg10Flux_pos']) - 10 ** m['lg10Flux']
            flux_err = np.append(flux_err, np.vstack([flux_err_neg, flux_err_pos]), axis=1)
            log_flux = np.append(log_flux, m['lg10Flux'].reshape(1, -1), axis=1)
            log_flux_err = np.append(
                log_flux_err,
                np.vstack([m['lg10Flux_neg'], m['lg10Flux_pos']]),
                axis=1
            )

        if model == 'powerlaw':  # pegpwrlw
            # Append the data to the arrays
            mjd = np.append(mjd, m['obs_mjd'].reshape(1, -1), axis=1)
            flux = np.append(flux, (1e-12 * m['norm']).reshape(1, -1), axis=1)
            flux_err_neg = 1e-12 * m['norm_neg']
            flux_err_pos = 1e-12 * m['norm_pos']
            flux_err = np.append(flux_err, np.vstack([flux_err_neg, flux_err_pos]), axis=1)
            log_flux_values = np.log10(1e-12 * m['norm'])
            log_flux = np.append(log_flux, log_flux_values.reshape(1, -1), axis=1)
            log_flux_err_neg = log_flux_values - np.log10(1e-12 * (m['norm'] - m['norm_neg']))
            log_flux_err_pos = np.log10(1e-12 * (m['norm'] + m['norm_pos'])) - log_flux_values
            log_flux_err = np.append(
                log_flux_err,
                np.vstack([log_flux_err_neg, log_flux_err_pos]),
                axis=1
            )

    # Solve for systematic errors
    flux_sys_neg = (flux_err[0, :] ** 2 + (sys_frac * flux[0, :]) ** 2) ** 0.5
    flux_sys_pos = (flux_err[1, :] ** 2 + (sys_frac * flux[0, :]) ** 2) ** 0.5
    flux_sys = np.vstack([flux_sys_neg, flux_sys_pos])

    # Convert back to logarithmic scale after applying systematic in linear space
    log_flux_sys_neg = log_flux[0, :] - np.log10(flux[0, :] - flux_sys_neg) 
    log_flux_sys_pos = np.log10(flux[0, :] + flux_sys_pos) - log_flux[0, :]
    log_flux_sys = np.vstack([log_flux_sys_neg, log_flux_sys_pos])

    # Initialize the array for upper limits
    ulims = np.zeros_like(flux[0, :]).reshape(1,-1)

    # If upper_limit is True, include manually calculated upper limits (decay only)
    if upper_limit:
        mjd = np.append(mjd, np.array([60444.44583857, 60448.53271186]).reshape(1, -1), axis=1)
        flux = np.append(flux, np.array([4.83704517e-13, 1.72585874e-12]).reshape(1, -1), axis=1)
        flux_err = np.append(
            flux_err,
            np.array([[0.0, 0.0], [0.0, 0.0]]),
            axis=1
        )
        log_flux = np.append(log_flux, np.array([np.log10(4.83704517e-13), np.log10(1.72585874e-12)]).reshape(1, -1), axis=1)
        log_flux_err = np.append(
            log_flux_err,
            np.array([[0.0, 0.0], [0.0, 0.0]]),
            axis=1
        )
        flux_sys = np.append(
            flux_sys,
            np.array([[0.0, 0.0], [0.0, 0.0]]),
            axis=1
        )
        log_flux_sys = np.append(
            log_flux_sys,
            np.array([[0.0, 0.0], [0.0, 0.0]]),
            axis=1
        )
        ulims = np.append(ulims, np.array([1.0, 1.0]).reshape(1, -1), axis=1)

    # Save the data to a text file
    output_data = np.vstack([
        mjd,
        flux,
        flux_err,
        flux_sys,
        log_flux,
        log_flux_err,
        log_flux_sys,
        ulims
    ]).T

    # Sort the output data by MJD
    output_data = output_data[np.argsort(output_data[:, 0])]

    header = (
        'Date (MJD), Flux (erg s^-1 cm^-2), Stat Err neg (erg s^-1 cm^-2), '
        'Stat Err pos (erg s^-1 cm^-2), Sys Err neg (erg s^-1 cm^-2), '
        'Sys Err pos (erg s^-1 cm^-2), Log10 Flux, Log10 Stat Err neg, '
        'Log10 Stat Err pos, Log10 Sys Err neg, Log10 Sys Err pos, Upper limit (0 - False, 1 - True)'
    )

    np.savetxt(filename, output_data, header=header)
         

def main():
    """
    Main function to execute the script.
    """
    
    # Run the function with the rise JSON file    
    write_to_txt('results/SwiftJ1727_swift-xrt_rise.json', '../SWJ1717_Files/XRT_rise.txt', upper_limit=False)

    # Run the function with the decay JSON file    
    write_to_txt('results/SwiftJ1727_swift-xrt_decay.json', '../SWJ1717_Files/XRT_decay.txt', upper_limit=True)

if __name__ == "__main__":
    main()
