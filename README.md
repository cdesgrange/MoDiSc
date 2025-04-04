# MoDiSc

MoDiSc stands for Modeling Disk in Scattered light. One can use the MoDiSc code to constrain the morphology and spectroscopy of a disk in scattered-light. It was originally used to model the total and polarized intensity data from the VLT/SPHERE instrument on the debris disk HD 120326 (HIP 67497) presented in Desgrange et al. (2025) resubm.

<br>

## Installation

First create a new environment to use MoDiSc.

> $ conda create -n modisc_env python=python3.9

Activate this environment

> $ conda activate modisc_env

Install MoDiSc via pip:

> $ pip install MoDiSc

The packages required by MoDiSc are indicated in MoDiSc/packages.py.
Most of the packages can be installed via 

> $ pip install [packagename]

**Exceptions:**

- opencv must be installed via $ pip install opencv-python

- emcee: the appropriate version of the emcee package to be used by MoDiSc is MoDiSc/emcee/

<br>

## Run MoDiSc

1. Go to the directory where you would like to run the simulations. Make sure to add the scripts run_modisc.py and plot_mcmc_results.py in this directory, and the folder "config_files/" in which are the configuration files of the simulations to be launched.

2. In the folder "config_files/", there are six example of simulations to be launched.

- job1 and job2 constrain the morphology of the debris disk HD 120326 based on SPHERE polarized intensity data at 1.6 $\mu$m (one observation), using a Nelder-Mead optimization (faster; job1) or a MCMC exploration (job2)

- job3 and job4 constrain the morphology of the debris disk HD 120326 based on SPHERE total intensity data at 1.6 $\mu$m (one observation), using a Nelder-Mead optimization (faster; job3) or a MCMC exploration (job4)

- job5 and job6 constrain the morphology of the debris disk HD 120326 based on SPHERE total and polarized intensity data at 1.6 $\mu$m (two observations in total), using a Nelder-Mead optimization (faster; job5) or a MCMC exploration (job6)

> python3 run_modisc.py config_files/job1_example_1obs_polarized_intensity_Nelder-Mead.yaml

3. Wait for the results! The results will be stored in results/*/beginning_date_of_the_simulation/final/. A logfile will be created and stored there, and a copy will also be added in the folder logs/, once the script finished to run.

4. If you run a simulation using a MCMC exploration, you can plot the cornerplot of the posteriors of the simulations by running

> python3 plot_mcmc_results.py config_files/job2_example_1obs_polarized_intensity_MCMC.yaml

the figures will be saved in results/*/beginning_date_of_the_simulation/results_mcmc/\*.pdf

<br>

## Issues

Do not hesitate to open a new issue.

<br>
<br>

*(this page and the documentation is work in progress)*

<br>
<br>

All the best, Célia Desgrange - ESO Chile



