[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8112490.svg)](https://doi.org/10.5281/zenodo.8112490)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dguarino/Guarino-Filipchuk-Destexhe/HEAD)
[![Open in Code Ocean](https://codeocean.com/codeocean-assets/badge/open-in-code-ocean.svg)](https://codeocean.com/capsule/9782876/tree)

# Convergent flows in modular networks generate reproducible firing patterns

Analysis code to reproduce all figures of the forthcoming paper by Guarino, Filipchuk, Destexhe      

The [`Ko`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/ko.ipynb) Jupyter notebook performs post-hoc power analysis of the data presented in [Ko et al. 2011](https://www.nature.com/articles/nature09880).    

The [`GuarinoFilipchukDestexhe`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/GuarinoFilipchukDestexhe.ipynb) notebook loads the MICrONS data, performs structural and dynamical analyses, and plots the results as in the paper figures (and more).    

Similarly, the notebooks [`stringer`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/stringer.ipynb) for [Stringer et al. 2019](https://www.science.org/doi/10.1126/science.aav7893), [`svoboda`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/svoboda.ipynb) for [Li et al. 2015](https://www.nature.com/articles/nature14178), [`goard`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/goard.ipynb) for [Franco and Goard 2021](https://www.science.org/doi/10.1126/sciadv.abf9815).

The service Jupyter notebooks [`attractor_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/attractor_analysis.ipynb), [`dynamical_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/dynamical_analysis.ipynb), [`structural_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/structural_analysis.ipynb), [`functional_analysis`](https://github.com/dguarino/Guarino-Filipchuk-Destexhe/blob/main/functional_analysis.ipynb) keep the analysis logic separated and commented.

All the code is hosted on this github repository (with a Zenodo DOI persistent identifier [here](https://zenodo.org/badge/latestdoi/498655295)) and can be interactively executed here [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dguarino/Guarino-Filipchuk-Destexhe/HEAD)

This repository also contains a copy of some required data files from the [MICrONS project phase1](https://www.microns-explorer.org/phase1), to ease the setup on Binder. The rest is downloaded from the MICrONS and Figshare repositories at binder runtime.
