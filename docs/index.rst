Welcome to MoDiSc's documentation!
===================================

``MoDiSc`` is a Python package for debris disk enthusiasts.
MoDiSc (Modeling Disks in Scattered light) is designed to constrain the morphology and photometry of debris disks observed in scattered light. 

In a nutshell, the tool ``MoDiSc`` takes as input a configuration file setting the paths to the datasets of the observations and defining the variables used in the simulations. By exploring the parameter space with a MCMC algorithm (``emcee`` package from `Foreman-Mackey et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_), **MoDiSc** searches for the disk model that best matches the observations. To initialize the values of the MCMC simulation, a first satisfying guess can be determined by running a Nelder-Mead minimization (`Nelder & Mead 1965 <https://academic.oup.com/comjnl/article-abstract/7/4/308/354237?redirectedFrom=fulltext>`_), which is implemented in MoDiSc. 


The disk model can be composed of one or several belts, each of them generated with the module ``fm.scattered_light_disk`` from `VIP_HCI <https://github.com/vortex-exoplanet/>`_ (`Gomez Gonzalez et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017AJ....154....7G/abstract>`_, `Christiaens et al. 2023a <https://ui.adsabs.harvard.edu/abs/2023JOSS....8.4774C/abstract>`_). The ``fm.scattered_light_disk`` module is based on the radiative transfer code GRaTeR  (`Augereau et al. 1999b <https://ui.adsabs.harvard.edu/abs/1999A&A...348..557A/abstract>`_). One or several observations in total and/or polarized intensity can be considered simultaneously, from one or several instruments.

**Origin**: I wrote the ``MoDiSc`` code by first working with the ``DiskFM`` code (`Mazoyer et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020SPIE11447E..59M/abstract>`_). I wanted to have a code more resilient for different types of simulations (based on one/several epochs, one/several instruments), without modifying everytime the code and so by simply giving in input a configuration file. I also wanted to have the results of the simulations automatically saved in a new folder, with a log of the simulation. I hope this could be useful for other people too.

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
   :maxdepth: 2
   :caption: Tutorials

   Quick-start-to-run-MoDiSc

.. toctree::
   :maxdepth: 2
   :caption: About

   Citation
   Contact
   
