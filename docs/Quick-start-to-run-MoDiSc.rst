Quick start to run MoDiSc
-------------------------

1. **Go to the directory where you would like to run the simulations**. Make sure to add the scripts run_modisc.py and plot_mcmc_results.py in this directory, and the folder "config_files/" in which are the configuration files of the simulations to be launched.

    

2. In the folder "config_files/", there are **six examples of simulations to be launched** just by running one commant.

- ``job1_*`` and ``job2_*``: constrain the morphology of the debris disk around HD 120326 based on SPHERE polarized intensity data at 1.6 $\mu$m (one observation), using a Nelder-Mead optimization (faster; ``job1_*``) 

.. code-block:: bash

  $ python3 run_modisc.py config_files/job1_example_1obs_polarized_intensity_Nelder-Mead.yaml

or a MCMC exploration (longer; ``job2_*``)

.. code-block:: bash

  $ python3 run_modisc.py config_files/job2_example_1obs_polarized_intensity_Nelder-Mead.yaml


- ``job3_*`` and ``job4_*`` constrain the morphology of the debris disk around HD 120326 based on SPHERE total intensity data at 1.6 $\mu$m (one observation), using a Nelder-Mead optimization

.. code-block:: bash

  $ python3 run_modisc.py config_files/job3_example_1obs_total_intensity_Nelder-Mead.yaml

or a MCMC exploration 

.. code-block:: bash

  $ python3 run_modisc.py config_files/job3_example_1obs_total_intensity_MCMC.yaml


- ``job5_*`` and ``job6_*`` constrain the morphology of the debris disk around HD 120326 based on SPHERE total and polarized intensity data at 1.6 $\mu$m (two observations in total), using a Nelder-Mead optimization 

.. code-block:: bash

  $ python3 run_modisc.py config_files/job5_example_2obs_polar_and_total_intensity_Nelder-Mead.yaml

or a MCMC exploration

.. code-block:: bash

  $ python3 run_modisc.py config_files/job5_example_2obs_polar_and_total_intensity_MCMC.yaml

3. **Wait for the results!** They will be stored at ``results/*/beginning_date_of_the_simulation/final/``. A logfile will be created and stored there, and a copy will also be added in the folder ``logs/``, once the Python script has finished running.



4. If you run a simulation using MCMC exploration, you can plot the **cornerplot of the posterior distributions** by running

.. code-block:: bash

  $ python3 plot_mcmc_results.py config_files/job2_example_1obs_polarized_intensity_MCMC.yaml

the figures will be saved at ``results/{name of the directory defined in the configuration file "job*.yaml" under the variable "SAVINGDIR"}/{date and time of the launch of the simulation}/results_mcmc/*.pdf``.


