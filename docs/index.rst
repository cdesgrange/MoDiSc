Welcome to MoDiSc's documentation!
===================================

.. image:: _figs/logo.pdf
   :align: center
   :width: 600px


``MoDiSc`` is a Python package for debris disk enthusiasts available on `GitHub <https://github.com/cdesgrange/MoDiSc>`_.
MoDiSc (Modeling Disks in Scattered light) is designed to constrain the morphology and photometry of debris disks observed in scattered light. 

In a nutshell, the tool ``MoDiSc`` takes as input a configuration file setting the paths to the datasets of the observations and defining the variables used in the simulations. By exploring the parameter space with a MCMC algorithm (``emcee`` package from `Foreman-Mackey et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_), ``MoDiSc`` searches for the disk model that best matches the observations. To initialize the values of the MCMC simulation, a first satisfying guess can be determined by running a Nelder-Mead minimization (`Nelder & Mead 1965 <https://academic.oup.com/comjnl/article-abstract/7/4/308/354237?redirectedFrom=fulltext>`_), which is implemented in MoDiSc. 


The disk model can be composed of one or several belts, each of them generated with the module ``fm.scattered_light_disk`` from `VIP_HCI <https://github.com/vortex-exoplanet/>`_ (`Gomez Gonzalez et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017AJ....154....7G/abstract>`_, `Christiaens et al. 2023a <https://ui.adsabs.harvard.edu/abs/2023JOSS....8.4774C/abstract>`_),  based on the radiative transfer code ``GRaTeR``  (`Augereau et al. 1999b <https://ui.adsabs.harvard.edu/abs/1999A&A...348..557A/abstract>`_). One or several observations in total and/or polarized intensity can be considered simultaneously, from one or several instruments.

**Origin**: I wrote the code ``MoDiSc`` based on the code ``DiskFM`` (`Mazoyer et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020SPIE11447E..59M/abstract>`_), which I adapted. I wanted a code more versatile than ``DiskFM``, easily launchable for different types of simulations (based on one or several epochs, one or several instruments, modeling one or several belts, using different number free parameters). I did not want to modify everytime the code, but simply give in input a configuration file containing all the information need for the simulation. I also wanted to have the results of the simulations automatically saved in a new folder, with a log of the simulation. I hope the code ``MoDiSc`` could be useful for other people too. Feel free to contact me (celia.desgrange@eso.org) if you would like a new feature implemented in ``MoDiSc`` for your work.

Check out the :doc:`Installation-and-dependencies` section to install MoDiSc and the tutorials in :doc:`Quick-start-to-run-MoDiSc` to run it on several examples.

.. note::

   This documentation is currently under active development (April 22, 2025).


Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   Installation-and-dependencies

.. toctree::
   :maxdepth: 3
   :caption: Tutorials

   Quick-start-to-run-MoDiSc
   Overview-of-the-configuration-file
   1A-Modeling-a-polarized-intensity-observation-Nelder-Mead
   1B-Modeling-a-polarized-intensity-observation-MCMC
   2A-Modeling-a-total-intensity-observation-Nelder-Mead
   2B-Modeling-a-total-intensity-observation-MCMC
   3A-Modeling-two-observations-Nelder-Mead
   3B-Modeling-two-observations-MCMC

.. toctree::
   :maxdepth: 2
   :caption: About

   Citation
   Contact
   
