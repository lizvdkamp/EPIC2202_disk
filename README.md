# EPIC 220208795
This repository contains all the code to analyse the K2 light curve and the ground-based photometry.
This is separated into several notebooks and a certain file structure which is presented at the end of this README.

The steps for data analysis are as follows.

### 1) Eclipse Modelling
<em>Goal</em>: This section attempts to find a circumplanetary disk solution to the eclipse. This is done with the help of <strong>pyPplusS</strong> a code package developed by Rein & Ofir 2019 (https://github.com/EdanRein/pyPplusS) and using an MCMC sampling method. 

#### a) Hard-edged Disk
Here we model an opaque disk (disk opacity = 1) with no edge component, based on a previously made test-run disk.
#### b) Soft-edged Disk
Here we add a translucent ring around the opaque disk in an attempt to fix the remaining residuals.

<em>Results</em>: Both models resulted in a promising solution, but the soft-edged disk did not manage to get rid of the remaining residuals


### 2) Period Folding
<em>Goal</em>: This section attempts to find a second eclipse in the ~20 years of ground-based photometric data. Lots of periods are tested, the photometry is left unbinned and binned after subtracting the stellar variation model and then inspected by performing a chi-squared ratio test.

<em>Result</em>: Two promising periods have been found.


### 3) Orbital Analysis
<em>Goal</em>: given the results from the previous sections we can determine some orbital parameters (the eccentricity, the periastron and apastron distance) and relate these to allowed masses and periods for the proposed companion around EPIC 220208795.

<em>Results</em>: these are presented in the plots in "plots 2202/parameters/"


### File Structure

  - Code : contains all the relevant functions for each of the notebooks
  - data : contains all the data files
      - photometry : light curves of K2 data
  - models : contains all the model data with mcmc backends (.h5) and best fits (.npy)
      - backends (can be requested via e-mail)
  - plots : contains all the plots
      - parameters : has the parameter maps around EPIC 220208795
  - pyPplusS : code package by Rein & Ofir 2019 (https://github.com/EdanRein/pyPplusS)


### Notebooks

  - k2_lightcurve - EPIC 2202.ipynb : contains the relevant code to extract the K2 lightcurve and save it
  - testrun_disk - EPIC 2202.ipynb : contains the MCMC run for a first iteration of an opaque disk, resulting in the parameters used as initial parameters for the hard-edged disk
  - hard_and_soft-edged_disk - EPIC 2202.ipynb : contains the MCMC run for the hard-edged and soft-edged disk
  - period_analysis - EPIC 2202.ipynb : contains the code for the period analysis as well as creating figure 4 in the paper
  - orbital_analysis - EPIC 2202.ipynb : contains the code for the orbital analysis as well as creating the parameter maps for the paper
  - figure making - EPIC 2202.ipynb : contains the relevant code used to create the rest of the figures for the paper

