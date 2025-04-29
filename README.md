# Description

This repo hosts the analysis that I led (as a part of a large collaboration) for the 2023-2024 outburst of the black hole X-ray binary Swift J1727.8-1613.

## Subdirectories:

### LightCurves

This hosts the data assembled and presented in Hughes et al., submitted in 2025. The data include +200 epochs of integrated radio flux densities from 7 facilities, with radio frequencies ranging from 0.3-300 GHz. 

The main components are: 

1. The radio data (in CSV format) are found in `files/SWJ1727_Radio.csv`
2. The plotting and analysis scripts are included in the notebook
3. The MAXI/GSC X-ray data needed to create the Hardness-Intensity diagram are also included in `files`


### RadioXrayCorrelation

This hosts the data and analysis for the radio-X-ray correlation ('LrLx'), Hughes et al. submitted in 2025

The components are:

1. The directory xray_analysis can reproduce the Swift-XRT data -- this is optional as all the files to reproduce figures/analysis are in 'SWJ1717_Files' -- I include this for any curious:
   a. Run `get_xrt_from_pipeline.py` (this requires swift-tools python API and HEASoft)
   b. Run `fit_xrt_spectra_rise.py` to fit the rise of the outburst
   c. Run `fit_xrt_spectra_decay.py` to fit the decay of the outburst
   d. Run `output_xrt_txt.py` to convert from the JSON output from b/c to the .txt files used in the notebook
2. The Xray/Radio data used are found in 'SWJ1727_Files'
3. Analysis and plotting scripts are in the notebook(s)

The notebook should be self-explanatory; however, please leave a note or discussion if anything is unclear. The radio data is taken from the LightCurves, with the only caveat being that we now fit only the core emission during the decay. It's difficult to include the images for reproducibility as they are > 1 GB. However, by the time you read this, all the data and analysis scripts will likely be publicly available, allowing you to verify the results yourself. More quickly, if you add a discussion/issue I'm happy to let you know what I did, and if I still have the images (and the request is reasonable) I may be able to provide them! 
