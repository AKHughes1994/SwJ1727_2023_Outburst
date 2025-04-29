#!/usr/bin/env python3
"""
Script to retrieve and process Swift XRT data for a given target.

This script uses the `swifttools` library to request and download XRT products 
(spectra and lightcurves) for specified targets. It also processes the downloaded 
spectral files for PyXspec spectral fitting by updating headers and grouping counts.

Authors:
- Andrew Hughes (andrew.hughes@physics.ox.ac.uk)

Dependencies:
- swifttools.xrt_prods
- astropy
- numpy
- glob
- subprocess
- Python 3.x

Usage:
- Run the script directly to process data for predefined targets.
- Ensure the required directories (`../spectra` and `../lightcurves`) exist.

"""

import time
import json
import glob
import subprocess
import os
import sys
import numpy as np
from astropy.io import fits

from swifttools.xrt_prods import XRTProductRequest

def msg(txt):
    """
    Print a timestamped message to the console.

    Parameters:
    - txt (str): The message to print.
    """
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt)

def get_xrt_prods(target_id, target_name, target_coords, segments, centroid=True, prod_type='spectrum', grade = 'all'):
    """
    Retrieve Swift XRT products (spectra or lightcurves) for a given target.

    Parameters:
    - target_id (str): The target ID for the Swift XRT request.
    - target_name (str): The name of the target.
    - target_coords (list): The [RA, Dec] coordinates of the target.
    - segments (list): List of observation segments to include.
    - centroid (bool): Whether to use centroiding for the request (default: True).
    - prod_type (str): Type of product to retrieve ('spectrum' or 'lightcurve').

    Returns:
    - None
    """
    # Validate the product type
    if prod_type not in ['spectrum', 'lightcurve']:
        sys.exit('Error: Please specify variable [prod_type] as either "spectrum" or "lightcurve".')

    # Ensure segments is at least a 1D array
    segments = np.atleast_1d(segments)

    # Initialize observation IDs
    obs_ids = ','.join([f'{target_id}{seg:03d}' for seg in segments])

    # Initialize the XRT product request
    myRequest = XRTProductRequest('hughes1@ualberta.ca', silent=False)

    # Set global parameters for the request
    myRequest.setGlobalPars(
        name=target_name,
        targ=target_id,
        SinceT0=False,
        RA=target_coords[0],
        Dec=target_coords[1],
        centroid=centroid,
        useSXPS=False,
        poserr=1
    )

    # ~~~~~~~~ #
    # Spectrum #
    # ~~~~~~~~ #
    if prod_type == 'spectrum':
        myRequest.addSpectrum(
            whichData='user',
            useObs=obs_ids,
            hasRedshift=False,
            grades=grade,
            doNotFit=True,
            srcrad=40,
            timeslice='obsid'
        )

        # Submit the XRT product request
        if myRequest.isValid()[0]:
            if not myRequest.submit():
                msg(f'I could not submit error: {myRequest.submitError}')
    
            while not myRequest.complete:
                time.sleep(10)

            # Retrieve the spectral products
            outdict = myRequest.retrieveSpectralFits(returnData=True)
        else:
            msg(f'BAD REQUEST: {myRequest.isValid()[1]}')
            exit()

        # Distribute the data files to the appropriate directories
        observations = [key for key in outdict.keys() if 'Obs' in key]
   
        # Download and extract files
        for obs in observations:
            pipeline_output = outdict[obs]['DataFile']
            subprocess.run([f'rm -rf ../spectra/*{obs}*'], shell=True)
            subprocess.run([f'wget -O ../spectra/{obs}.tar.gz {pipeline_output}'], shell=True)
            subprocess.run([f'tar -xf ../spectra/{obs}.tar.gz -C ../spectra/'], shell=True)

    # ~~~~~~~~~~ #
    # Lightcurve #
    # ~~~~~~~~~~ #
    else:
        myRequest.addLightCurve(
            whichData='user',
            useObs=obs_ids,
            binMeth='obsid',
            timeFormat='m',
            minEnergy=0.5,
            maxEnergy=10.0,
            softLo=0.5,
            softHi=2.0,
            hardLo=2.0,
            hardHi=10.0,
            allowUL='both'
        )

        # Submit the XRT product request
        if myRequest.isValid()[0]:
            if not myRequest.submit():
                msg(f'I could not submit error: {myRequest.submitError}')
    
            while not myRequest.complete:
                time.sleep(10)

            # Retrieve and extract lightcurve products
            outdict = myRequest.retrieveLightCurve(returnData=True, incbad=True)
        else:
            msg(f'BAD REQUEST: {myRequest.isValid()[1]}')
            exit()

        # Create a directory for the target and download the products
        subprocess.run([f'rm -rf ../lightcurves/{target_id}; mkdir ../lightcurves/{target_id}'], shell=True)
        myRequest.downloadProducts(f'../lightcurves/{target_id}')

        # Extract and clean up the files
        subprocess.run([f'tar -xf ../lightcurves/{target_id}/lc.tar.gz -C ../lightcurves/{target_id}/; mv ../lightcurves/{target_id}/USERPROD*/lc/* ../lightcurves/{target_id}/.'], shell=True)
        subprocess.run([f'rm -rf ../lightcurves/{target_id}/USERPROD*; rm -rf ../lightcurves/{target_id}/lc.tar.gz'], shell=True)

def main():
    """
    Main function to retrieve and process Swift XRT data for predefined targets.
    """

    # Make sure spectra and lightcurves directories exist
    if not os.path.exists('../spectra'):
        os.makedirs('../spectra')
    if not os.path.exists('../lightcurves'):
        os.makedirs('../lightcurves')

    # Define the Swift XRT parameters for Swift J1727.8-2609
    target_coords = [261.9318, -16.2066]
    target_names = ['SwiftJ1727d8m1613', 'Swift J1727.8-1613', 'GRB 230824A']
    target_ids = ['00089766', '00016584', '01186959']
    segments = [
        [2, 3, 4, 5, 6, 7, 12],
        [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17, 18, 19],
        [5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 1, 2, 20, 21, 22, 23, 24, 26, 25, 27, 28, 29, 30, 31, 32, 34, 33, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57]
    ]
    grades = ['all', 'all', '0']
    centroids = [False, False, False]

    # Iterate through target names, IDs, and segments to retrieve XRT products
    for k, target_name in enumerate(target_names):
        # get_xrt_prods(target_ids[k], target_name, target_coords, segments[k], centroid=centroids[k], prod_type='lightcurve', grade=grades[k])
        # get_xrt_prods(target_ids[k], target_name, target_coords, segments[k], centroid=centroids[k], prod_type='spectrum', grade=grades[k])
        pass

    # Process spectral files for PyXspec spectral fitting
    spectrum_directory = os.getcwd().replace('code', 'spectra')
    source_spectra = sorted(glob.glob(f'{spectrum_directory}/*source.pi'))

    for source_spectrum in source_spectra:
        # Define file names
        group_spectrum = source_spectrum.replace('source.pi', 'group.pi')
        final_spectrum = source_spectrum.replace('source.pi', 'final.pi')
        background_spectrum = source_spectrum.replace('source.pi', 'back.pi')
        rmf = source_spectrum.replace('source.pi', '.rmf')
        arf = source_spectrum.replace('source.pi', '.arf')

        print('\n', source_spectrum, fits.getheader(source_spectrum)['DATE-OBS'])

        # Update headers with spectral file references
        syscall = f'grppha {source_spectrum} !{group_spectrum} comm="chkey backfile {background_spectrum} & chkey respfile {rmf} & chkey ancrfile {arf} & exit"'
        subprocess.run([syscall], shell=True)

        # Group counts using Sivakoff GRPPHA
        syscall = f'python3 new_grppha.py -i {group_spectrum} -o {final_spectrum} -l 500 -u 10000 -c 25 -e 0.0'
        subprocess.run([syscall], shell=True)

if __name__ == "__main__":
    main()