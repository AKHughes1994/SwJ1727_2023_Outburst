"""
Simple script to convert the WAPITI text output 
to text files compatible with the QU modelling code 
"""

import glob
import numpy as np
from astropy.time import Time

def main():

    # Define input and output directories
    INPUT_DIR = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/3C286/wapiti_output'
    OUTPUT_DIR = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/3C286/QU_text'

    # Construct sorted list of files
    txt_files = sorted(glob.glob(f"{INPUT_DIR}/*.txt"))
    
    # Iterate over files
    for txt_file in txt_files[:]:

        print(f"Processing file: {txt_file}")

        # Read the header to extract the Middle MJD and source name
        with open(txt_file, 'r') as f:
            lines = f.readlines()
            middle_mjd = None
            source_name = None
            for line in lines:
                if '# Middle MJD:' in line:
                    middle_mjd = float(line.split(':')[1].strip())
                if '# Source:' in line:
                    source_name = line.split(':')[1].strip().split('_')[0]  # Extract source name before underscore
            
        if middle_mjd is None:
            print(f"Warning: Could not find Middle MJD in {txt_file}, skipping...")
            continue
        
        if source_name is None:
            source_name = "UNKNOWN"  # Fallback if source name not found
        
        # Convert MJD to ISO format
        t = Time(middle_mjd, format='mjd')
        date_isot = t.isot
        
        # Load the data (skip header lines starting with #)
        data = np.loadtxt(txt_file, comments='#')
        
        # Extract the columns: Channel Freq[GHz] I[mJy] Q[mJy] U[mJy] V[mJy] Ptot[mJy] err_I err_Q err_U err_V err_Ptot
        # We need: freq(Hz) I(Jy) Q(Jy) U(Jy) V(Jy) dI(Jy) dQ(Jy) dU(Jy) dV(Jy)
        freq_Hz = data[:, 1] * 1e9        # Freq[GHz] -> Hz
        I_Jy = data[:, 2] / 1e3           # I[mJy] -> Jy
        Q_Jy = data[:, 3] / 1e3           # Q[mJy] -> Jy
        U_Jy = data[:, 4] / 1e3           # U[mJy] -> Jy
        V_Jy = data[:, 5] / 1e3           # V[mJy] -> Jy
        dI_Jy = data[:, 7] / 1e3          # err_I -> Jy
        dQ_Jy = data[:, 8] / 1e3          # err_Q -> Jy
        dU_Jy = data[:, 9] / 1e3          # err_U -> Jy
        dV_Jy = data[:, 10] / 1e3         # err_V -> Jy
        
        # Stack the data
        text_array = np.column_stack([freq_Hz, I_Jy, Q_Jy, U_Jy, V_Jy, dI_Jy, dQ_Jy, dU_Jy, dV_Jy])

        # Save array to text file where each row is a channel, split at the T
        output_file = f"{OUTPUT_DIR}/3C286_WAPITI_{date_isot.split('T')[0].replace('-', '')}.txt"
        header = (f"# {source_name} QU data converted from WAPITI text output\n"
                  f"# Date (ISO): {date_isot}\n"
                  f"# Columns: freq(Hz) I(Jy) Q(Jy) U(Jy) V(Jy) dI(Jy) dQ(Jy) dU(Jy) dV(Jy)\n")
        np.savetxt(output_file, text_array, header=header, fmt='%.6f')

if __name__ == "__main__":
    main()