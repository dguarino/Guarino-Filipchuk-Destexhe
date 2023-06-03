[![DOI](https://zenodo.org/badge/498655295.svg)](https://zenodo.org/badge/latestdoi/498655295) 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dguarino/Guarino-Filipchuk-Destexhe/HEAD)

# Attractor dynamics without structural signature in the mouse cortex

Analysis code to reproduce all figures of the forthcoming paper by Guarino, Filipchuk, Destexhe      
preprint link: https://www.biorxiv.org/content/10.1101/2022.05.24.493230

The [`GuarinoFilipchukDestexhe`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/GuarinoFilipchukDestexhe.ipynb)  Jupyter notebook performs loading and selection of the MICrONS data, structural and dynamical analyses, and plots the results as in the paper panels.    

The [`cortical_electrophysiology`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/cortical_electrophysiology.ipynb) Jupyter notebook performs loading and selection of the [Stringer et al. 2019](https://www.science.org/doi/10.1126/science.aav7893) datasets, and performs supporting dynamical analyses.

The service Jupyter notebooks [`attractor_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/attractor_analysis.ipynb), [`dynamical_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/dynamical_analysis.ipynb), and [`structural_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/structural_analysis.ipynb) keep the analysis logic separated and commented.

All the code is hosted on this github repository (with a Zenodo DOI persistent identifier [here](https://zenodo.org/badge/latestdoi/498655295)) and can be interactively executed here [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dguarino/Guarino-Filipchuk-Destexhe/HEAD)

This repository also contains a copy of some required data files from the [MICrONS project phase1](https://www.microns-explorer.org/phase1), to ease the setup on Binder. The rest is downloaded from the MICrONS and Figshare repositories at binder runtime.
