"""
Simple script to convert the JSON output from polkat 
to text files compatible with the QU modelling code 
"""

import json
import glob
import numpy as np

def main():

    # Define input and output directories
    INPUT_DIR = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/J1727/polkat_json'
    OUTPUT_DIR = '/home/andrewhughes/Projects/X-KAT/SwiftJ1727/new_meerkat/SpectroPol/publication_work/QU_phenom/files/J1727/QU_text'

    # Construct sorted list of files
    json_files = sorted(glob.glob(f"{INPUT_DIR}/*.json"))
    
    # Iterate over files
    for json_file in json_files[:]:

        print(f"Processing file: {json_file}")

        # Open json file
        with open(json_file, 'r') as f:
            data = json.load(f)

        # Take only the channelised array and the time
        date_isot = data['MFS']['component0']['date_isot'][0]
        
        # Process component0
        data_array = data['CHAN']['component0']

        # Numpize the data array
        for key in data_array:
            data_array[key] = np.array(data_array[key])

        # Extract the data of interest
        text_array = [
            data_array['freq_GHz'][0] * 1e9,
            data_array['I_flux_mJy'][0] / 1e3,
            data_array['Q_flux_mJy'][0] / 1e3,
            data_array['U_flux_mJy'][0] / 1e3,
            data_array['V_flux_mJy'][0] / 1e3,
            data_array['I_rms_mJy'][0] / 1e3,
            data_array['Q_rms_mJy'][0] / 1e3,
            data_array['U_rms_mJy'][0] / 1e3,
            data_array['V_rms_mJy'][0] / 1e3
        ]

        # Save array to text file where each row is a channel, split at the T
        output_file = f"{OUTPUT_DIR}/SwiftJ1727_QU_{date_isot.split('T')[0].replace('-', '')}.txt"
        header = (f"# SwiftJ1727 QU data converted from polkat JSON output\n"
                  f"# Date (ISO): {date_isot}\n"
                  f"# Columns: freq(Hz) I(Jy) Q(Jy) U(Jy) V(Jy) dI(Jy) dQ(Jy) dU(Jy) dV(Jy)\n")
        np.savetxt(output_file, np.column_stack(text_array), header=header, fmt='%.6f')
        
        # Check if component1 exists and process it for ejecta
        if 'component1' in data['CHAN']:
            print(f"  Found component1, creating ejecta file...")
            data_array_ejecta = data['CHAN']['component1']
            
            # Numpize the data array
            for key in data_array_ejecta:
                data_array_ejecta[key] = np.array(data_array_ejecta[key])
            
            # Extract the data of interest
            text_array_ejecta = [
                data_array_ejecta['freq_GHz'][0] * 1e9,
                data_array_ejecta['I_flux_mJy'][0] / 1e3,
                data_array_ejecta['Q_flux_mJy'][0] / 1e3,
                data_array_ejecta['U_flux_mJy'][0] / 1e3,
                data_array_ejecta['V_flux_mJy'][0] / 1e3,
                data_array_ejecta['I_rms_mJy'][0] / 1e3,
                data_array_ejecta['Q_rms_mJy'][0] / 1e3,
                data_array_ejecta['U_rms_mJy'][0] / 1e3,
                data_array_ejecta['V_rms_mJy'][0] / 1e3
            ]
            
            # Save ejecta array to text file
            output_file_ejecta = f"{OUTPUT_DIR}/SwiftJ1727_QU_{date_isot.split('T')[0].replace('-', '')}_ejecta.txt"
            header_ejecta = (f"# SwiftJ1727 QU data (ejecta/component1) converted from polkat JSON output\n"
                            f"# Date (ISO): {date_isot}\n"
                            f"# Columns: freq(Hz) I(Jy) Q(Jy) U(Jy) V(Jy) dI(Jy) dQ(Jy) dU(Jy) dV(Jy)\n")
            np.savetxt(output_file_ejecta, np.column_stack(text_array_ejecta), header=header_ejecta, fmt='%.6f')


if __name__ == "__main__":
    main()