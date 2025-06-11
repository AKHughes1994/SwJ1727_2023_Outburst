# LightCurve Subdirectory

This subdirectory contains data and scripts related to the light curve analysis of the 2023-2024 outburst of the black hole X-ray binary Swift J1727.8-1613. The paper corresponding to this data is called:

Hughes et al. 2025: "Comprehensive Radio Monitoring of the Black Hole X-ray Binary Swift J1727.8−1613 during its 2023–2024 Outburst" 

If you make use of this data please cite:

```
@ARTICLE{2025arXiv250607798H,
       author = {{Hughes}, Andrew K. and {Carotenuto}, Francesco and {Russell}, Thomas D. and {Tetarenko}, Alexandra J. and {Miller-Jones}, James C.~A. and {Bahramian}, Arash and {Bright}, Joe S. and {Cowie}, Fraser J. and {Fender}, Rob and {Gurwell}, Mark A. and {Khaulsay}, Jasvinderjit K. and {Kirby}, Anastasia and {Jones}, Serena and {Lescure}, Elodie and {McCollough}, Michael and {Plotkin}, Richard M. and {Rao}, Ramprasad and {Vrtilek}, Saeqa D. and {Williams-Baldwin}, David R.~A. and {Wood}, Callan M. and {Sivakoff}, Gregory R. and {Altamirano}, Diego and {Casella}, Piergiorgio and {Corbel}, Stephane and {DeBoer}, David R. and {Del Santo}, Melania and {Echiburu-Trujillo}, Constanza and {Farah}, Wael and {Gandhi}, Poshak and {Koljonen}, Karri I.~I. and {Maccarone}, Thomas and {Matthews}, James H. and {Markoff}, Sera B. and {Pollak}, Alexander W. and {Russell}, David M. and {Saikia}, Payaswini and {Castro Segura}, Noel and {Shaw}, Aarran W. and {Siemion}, Andrew and {Soria}, Roberto and {Tomsick}, John A. and {van den Eijnden}, Jakob},
        title = "{Comprehensive Radio Monitoring of the Black Hole X-ray Binary Swift J1727.8$-$1613 during its 2023$-$2024 Outburst}",
      journal = {arXiv e-prints},
     keywords = {High Energy Astrophysical Phenomena},
         year = 2025,
        month = jun,
          eid = {arXiv:2506.07798},
        pages = {arXiv:2506.07798},
archivePrefix = {arXiv},
       eprint = {2506.07798},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv250607798H},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

This is a "living" dataset, so as more data comes in, we will update it, below are the current additions that were not in the paper: 

**None**


Below is a breakdown of its components:

## Components

### 1. **Data Files**
   - **`SWJ1727_Radio.csv`**: This file contains over 200 epochs of integrated radio flux densities collected from seven facilities. The data is crucial for analyzing the radio emission behavior during the outburst. The columns are as follows:
   - **X-ray Data**:
     - **`Swift_J1727_2.0-6.0keV_gsclc.dat`**: Contains light curve data for the 2.0-6.0 keV energy range.
     - **`Swift_J1727_6.0-20.0keV_gsclc.dat`**: Contains light curve data for the 6.0-20.0 keV energy range.

### 2. **Notebooks**
   - **Jupyter Notebooks**: The notebooks in this directory perform the plotting and analysis of the light curves. They use the data from the provided CSV files to generate visualizations and insights.

### 3. **X-ray Data**
   - Additional MAXI/GSC X-ray data required for creating the Hardness-Intensity diagram are located here. These datasets complement the radio data and help in understanding the X-ray emission characteristics during the outburst.

---

### Usage Instructions

1. **Data Analysis**: 
   - Use the provided Jupyter notebooks to analyze the radio and X-ray data.
   - Ensure you have the necessary Python libraries installed (e.g., numpy, matplotlib, pandas).

2. **Plotting**: 
   - Run the notebooks to generate plots such as the light curves and Hardness-Intensity diagrams.

3. **Customization**:
   - Feel free to modify the notebooks to explore additional analyses or plots.

---

### Notes
If you encounter any issues or have questions about the data or analysis scripts, please leave a note or start a discussion in the repository.
