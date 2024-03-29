# surface 

This repository provides the datasets and Python scripts used in the paper "A Prescription for the Asteroseismic Surface Correction" by Yaguang Li et al. (2023). The stellar models used in this project can be accessed from [Zenodo](https://zenodo.org/record/7905521).


## Directory Structure 
The repository is organized as follows:

- `data/`: Contains intermediate data files used in this project.

- `paper_surface_plots/`: Contains scripts and data to reproduce the plots in the paper. This directory is curated for easy replication of the results.

- `sample/`: Contains the following datasets for the sample stars studied:
     - `samples.xlsx`: Stellar parameters.
     - `modes.xlsx`: Oscillation frequencies.

- `src/`: Contains the source code used in this project. This directory is not curated, so use the contents with caution.

- `mesa_work/`: Constains the MESA/GYRE inlists used to produce the stellar models.

- `fDnu.ipynb`: Contains the source code to calculate correction factors for the $\Delta\nu$ scaling relation presented in this work.

- `requirements.txt`: Lists the Python packages required for this project. Note that some code requires functions defined in two additional unreleased packages, `asteroseismology` and `grid`. These can be accessed on GitHub at the following links:
    - [asteroseismology](https://github.com/parallelpro/asteroseismology)
    - [grid](https://github.com/parallelpro/grid)


# Citation
If you use the data or routines provided in this repository for your work, please cite the following paper:

     @ARTICLE{2023MNRAS.523..916L,
          author = {{Li}, Yaguang and {Bedding}, Timothy R. and {Stello}, Dennis and {Huber}, Daniel and {Hon}, Marc and {Joyce}, Meridith and {Li}, Tanda and {Perkins}, Jean and {White}, Timothy R. and {Zinn}, Joel C. and {Howard}, Andrew W. and {Isaacson}, Howard and {Hey}, Daniel R. and {Kjeldsen}, Hans},
          title = "{A prescription for the asteroseismic surface correction}",
          journal = {\mnras},
          keywords = {stars: low-mass, stars: oscillations, stars: solar-type, Astrophysics - Solar and Stellar Astrophysics},
          year = 2023,
          month = jul,
          volume = {523},
          number = {1},
          pages = {916-927},
               doi = {10.1093/mnras/stad1445},
     archivePrefix = {arXiv},
          eprint = {2208.01176},
     primaryClass = {astro-ph.SR},
          adsurl = {https://ui.adsabs.harvard.edu/abs/2023MNRAS.523..916L},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
     }





