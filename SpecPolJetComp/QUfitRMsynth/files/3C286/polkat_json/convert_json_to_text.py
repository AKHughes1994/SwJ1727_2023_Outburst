"""
Simple script to convert the JSON output from polkat 
to text files compatible with the QU modelling code 
"""

import json
import glob
import numpy as np

def main():

    # Define input and output directories
    INPUT_DIR = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/3C286/polkat_json'
    OUTPUT_DIR = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/3C286/QU_text'

    # Construct sorted list of files
    json_files = sorted(glob.glob(f"{INPUT_DIR}/*.json"))
    
    # Iterate over files
    for json_file in json_files[:]:

        print(f"Processing file: {json_file}")

        # Open json file
        with open(json_file, 'r') as f:
            data = json.load(f)

        # Take only the channelised array and the time
        # The MFS section contains scan000 (not component0), with scalar values
        date_obs = data['MFS']['scan000']['date_obs']
        
        # Process scan000 channelised data
        data_array = data['CHAN']['scan000']

        # Convert lists to numpy arrays
        for key in data_array:
            data_array[key] = np.array(data_array[key])

        # Extract the data of interest - all channels (not just [0])
        text_array = np.column_stack([
            data_array['freq_GHz'] * 1e9,
            data_array['I_flux_mJy'] / 1e3,
            data_array['Q_flux_mJy'] / 1e3,
            data_array['U_flux_mJy'] / 1e3,
            data_array['V_flux_mJy'] / 1e3,
            data_array['I_rms_mJy'] / 1e3,
            data_array['Q_rms_mJy'] / 1e3,
            data_array['U_rms_mJy'] / 1e3,
            data_array['V_rms_mJy'] / 1e3
        ])

        # Save array to text file where each row is a channel, split at the T
        output_file = f"{OUTPUT_DIR}/3C286_QU_{date_obs.split('T')[0].replace('-', '')}.txt"
        header = (f"# 3C286 QU data converted from polkat JSON output\n"
                  f"# Date (ISO): {date_obs}\n"
                  f"# Columns: freq(Hz) I(Jy) Q(Jy) U(Jy) V(Jy) dI(Jy) dQ(Jy) dU(Jy) dV(Jy)\n")
        np.savetxt(output_file, text_array, header=header, fmt='%.6f')


if __name__ == "__main__":
    main()