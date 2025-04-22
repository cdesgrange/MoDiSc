Here, I present the information written in a configuration file to model a polarized intensity observation using a Nelder-Mead optimization.

The example configuration file corresponds to ``job1_example_1obs_polarized_intensity_Nelder-Mead.yaml`` located in the folder ``config_files/`` on `GitHub <https://github.com/cdesgrange/MoDiSc>`_,

.. note::

  - When a parameter is a list, the list is the size of the number of observations te be matched in the simulations.
  
  - Even if there is only one observation, keep the list format.

  - In principle, the modifications only have to be done in this configuration file and nothing in the scripts *.py.

  - Some keywords are not used to launch a simulation considering only one polarized intensity dataset. Don't remove unuseful keywords, because the script ``run_modisc.py`` will still load their value.


Paths
-----

In each configuration file, the user can either indicate the path to the dataset on their laptop, or on the external cluster. The type of path considered depend on the value of the keyword ``EXTERNAL_SERVER``. If ``EXTERNAL_SERVER`` is equal to 1 or True, the paths considered will be those finishing by **_EXT. Otherwise, the paths considered will be those finishing by **_INT.


.. code-block:: bash

  ###########
  ## PATHS ##
  ###########
  EXTERNAL_SERVER: 0 # 1 = True: the paths considered will be **_EXT ; 0 = False: the paths considered will be **_INT
  
  ########################
  # CASE: ONE'S COMPUTER #
  ########################
  DATADIR_INT    : ['data/example/', ] # path to the folder in which are located the files of the observations
  
  FN_SCIENCE_INT : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/science_data_Q_star_pol_subtr.fits', ] # filename(s) of the science data. Science data are pre- or post-processed, depending if they should be post-processed when doing the simulations to look for the best disk model. Example: FN_SCIENCE_ALL should indicate for SPHERE polarized intensity data the IRDAP post-processed science data, but for SPHERE pupil-stabilized observations, the pre-processed data, to take into account the self-subtraction (Milli+2012) effect.
  FN_CUBE_REF_INT: [None, None, ] # filename(s) of the reference cube data. If no reference cube, the value is [..., None, ...]. Reference cube data are pre-processed, and will be used to post-process total intensity data using the RDI(+ADI) technique(s).
  FN_PSF_INT     : [None, ]
  FN_PSF_1_INT   : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/psf_cube_flux_processed_left.fits', ]
  FN_PSF_2_INT   : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/psf_cube_flux_processed_right.fits', ]
  # (!) check that the value of TWO_PSF_FILES is consistent with the filenames for the psf given above. If for one given observation, there are two different files (= located at two different paths) for the PSF indicate them in FN_PSF_1_INT (= [..., 'fn_psf1', ...]) and FN_PSF_2_INT (= [..., 'fn_psf2', ...]). In this case, for this observation the value of FN_PSF_INT will be ignored. TWO_PSF_FILES must be set to [..., True, ...] (or [..., 1, ...]).  Otherwise, if for one given observation, there is only one file (= located at a single path) for the PSF, indicate it in FN_PSF_INT (= [..., 'fn_psf', ...]). In this case, for this observation the values of FN_PSF_1_INT and FN_PSF_2_INT will be ignored. TWO_PSF_FILES must be set to [..., False, ...] (or [..., 0, ...]). 
  
  FN_PA_INT      : [None, ]  # filename(s) of the parallactic angles file(s). Relevant for pupil-tracking stabilized observations, or observations acquired for several rolling angles. The FN_PA_ALL may be set to [..., None, ...] for other types of observations.
  FN_NOISE_INT   : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/noise_map_U_phi_star_pol_subtr_annuli.fits', ] # filename(s) of the noise map. The noise map can be provided or computed later in the script. In the latter case, FN_NOISE_ALL is set to [..., None, ...] and COMPUTE_NOISE_MAP_ALL should be set to [..., 1, ...]
  
  PATH_MASK_INT  : ['data/example/mask_inner_belt_SPHERE_IRDIS.fits', ] # Remark: The .fits file corresponding to the mask can be in a folder different than the folder where are located the data. This is why here the full path is given, and not only the filename.
  
  #########################
  # CASE: EXTERNAL SERVER #
  #########################
  # See the comments in the previous block listing the parameters "*_INT"
  DATADIR_EXT    : ['data/example/', ]
  FN_SCIENCE_EXT : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/science_data_Q_star_pol_subtr.fits', ]
  FN_CUBE_REF_EXT: [None, None, ] 
  FN_PSF_EXT     : [None, ]
  FN_PSF_1_EXT   : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/psf_cube_flux_processed_left.fits', ]
  FN_PSF_2_EXT   : ['polarized_intensity_SPHERE_IRDIS_2018-06-01_processed_with_IRDAP/psf_cube_flux_processed_right.fits', ]
  FN_NOISE_EXT   : [None, ]
  FN_PA_EXT      : [None, ] 
  PATH_MASK_EXT  : ['data/example/mask_inner_belt_SPHERE_IRDIS.fits', ] 


In every case, the results will be saved in the directory ``SAVINGDIR``. If the folder does noes exist yet, it will be created when running the simulation.
The keyword ``RESULTDIR`` is not used for simulations using Nelder-Mead optimization to determine a disk fitting the observations.


.. code-block:: bash

  #################
  # IN EVERY CASE #
  #################
  SAVINGDIR: 'results/results_BBH_polar/' # path where will be created the folder containing the outputs of the simulations. The name of the folder will correspond to the date when MoDiSc was launched and is automatically generated by MoDiSc
  RESULTDIR: 'results/results_BBH_polar/2025-2-3-16h45min21s/' # when one's wishes to plot the results of the MCMC simulations, at the end of this path must be updated the name of the folder automatically generated by MoDiSc containing the outputs of the simulations. For example: RESULTDIR: 'results/results_BBH_polar/2025-2-3-16h45min21s/'


Display Parameters
------------------

The user can decide to print more or less information regarding the simulation. It is adviced to set the keyword ``DISPLAY_GENERAL_INFO`` to 1 (or True) to print (and save in the log file) general information regarding the simulation.

The keyword ``DISPLAY_INFO_SIMU_MCMC`` should be set to 1 only when testing whether the simulation runs well. The values associated to the different walkers and the associated chisquare will be printed. This is time-consuming, so it should be set to 0 if unuseful.

The keyword ``DISPLAY_INFO_SIMU_NELDERMEAD`` should can be set to 1 to see how the values of the free parameters evolve during the optimization; as the value of the associated chisquare.

.. code-block:: bash

  ################
  ## PARAMETERS ##
  ################
  DISPLAY_GENERAL_INFO: 1   # print the information about loading the dataset(s), parameters of the simulations and general status of the MCMC or Nelder-Mead simulations
  DISPLAY_INFO_SIMU_MCMC: 0 # print the information when MCMC simulations are running (this should be set to 1 only when testing if the simulations runs well, because this is time-consuming to print the information for all the MCMC iterations and for all the walkers)
  DISPLAY_INFO_SIMU_NELDERMEAD: 1 # print the information when Nelder-Mead simulations are running
  
  
  
  
  
  











