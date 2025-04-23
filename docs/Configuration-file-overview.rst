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
  
  
System Parameter
----------------

Indicate the distance of the observed system. This would be used by the ``VIP_HCI`` function ``vip_hci.fm.scattered_light_disk.ScatteredLightDisk()`` in the function ``MoDiSc.simulations.generate_disk_model()``. This is to derive the location of the dust belt in the preprocessed (or postprocessed) science image in which the disk has to be modeled. The reference radius of the dust belt is given in au and the platescale of the image is given in arcseconds/au.
  
.. code-block:: bash

  ####################
  # SYSTEM PARAMETER #
  ####################

  DISTANCE_STAR: 113.27  # distance to the system of interest in pc


Observation Parameters
----------------------

Here are listed all the parameters of the observations, in a list. A list of one element implies that one observation will be modeled. A list of two elements implies that two observations will be modeled, and so on. 

.. code-block:: bash

  ##########################
  # OBSERVATION PARAMETERS #
  ##########################
  EPOCHS       : ['2018-06-01'] # list of epochs(s) corresponding to the observation(s). One value per observation.
  INSTRU       : ['IRDIS']      # list of instrument(s) corresponding to the observation(s). One value per observation.
  TYPE_OBS     : ['polarized_intensity'] # list of the type of the observation: total intensity ('total_intensity') or polarimetry ('polarized_intensity'). One value per observation.
  PLATE_SCALE  : [0.012255] # list of the plate scale values in arcseconds. One value per observation.
  SPECTRAL_AXIS: [0]        # if [..., 1, ...], there is a spectral axis for the psf and science data, if [..., 0, ...], there is not
  CHANNELS     : [[0,1]]    # list of list of spectral channels (one channel = one wavelength) to be considered. Example: For 3 observations, CHANNELS_ALL = [[0],[0,1],[None]] indicates that for the first observation, only the first spectral will be considered, whereas for the second observation, both the first and second channels will be considered, and for the third observation, the parameter is not relevant because there is no spectral axis (SPECTRAL_AXIS should be equal to [1, 1, 0]).
  TWO_PSF_FILES: [1]      # indicate whether there are two different files (= located at two different paths) to consider for the PSF (1 = yes, 0 = no)

.. note::

  If the observation to be modeled in a science cube with a spectral axis, ``SPECTRAL_AXIS`` should be set to [..., 1, ...], and in ``CHANNELS`` should be indicated the index of the spectral channels to be considered. If there are several indexes; the images at these wavelength will be mean summed.

.. note::

  The keyword TWO_PSF_FILES indicate whether there are two different files (= located at two different paths) to consider for the PSF (1 = yes, 0 = no). For instance, SPHERE/IRDIS polarized intensity data processed with IRDAP have two PSF, for the left and right part of the detector, stored in two different files. On the other hand; the PSF(s) of SPHERE/IRDIS total intensity data pre-processed by the High-Contrast Data Center are stored in one given file.

.. code-block:: bash
  
  CROP_PSF    : [500]   # for the PSF DATA
  CROP_SCIENCE: [412]   # list of cropping parameter for the spatial dimensions of the science cube/image. One value per observation.
  CROP_NOISE  : [412]   # list of cropping parameter for the spatial dimensions of the noise cube/image. One value per observation. 
  CROP_MASK   : [0]     # list of cropping parameter for the spatial dimensions of the mask image. One value per observation.
  CROP_REF    : [None]  # list of cropping parameter for the spatial dimensions of the ref cube. One value per observation.

  COMPUTE_NOISE_MAP: [0]   # list of booleans (or 0/1) indicating whether the noise map should be computed or is already provided. One value per observation. True = 1 means yes, compute the noise map from the science data. False = 0 means no, load it from the path DATADIR + FN_NOISE.
  NOISE_MULTIPLICATION_FACTOR: [1]  # list of floats. One value per observation. This multiplication factor can be used to artificially increase the value of the noise cube/image. See e.g. Mazoyer et al. 2020 in their SPIE paper about diskFM. Default value: 1

.. note::

   # The "CROP_*" parameters indicate the number of pixels to remove both in left-right, top-bottom directions

.. warning::

  - The spatial dimension of the science, noise, mask data should be the same. 

  - If the keyword ``COMPUTE_NOISE_MAP`` is set to 1, ``CROP_NOISE`` should be set to 0 because the noise map is computed from the cropped ``SCIENCE_DATA``.

  
.. code-block:: bash

  SPATIAL_SHIFT_PSF_DATA: [0.5]      # list of floats indicating the number of pixels to offset the psf image. One value per observation.
  SPATIAL_SHIFT_SCIENCE_DATA: [0.5]  # list of floats indicating the number of pixels to offset the science data. One value per observation. In practice, this is only use in the case of polarized intensity data and if DO_ROBUST_CONVOLUTION is set to 1, when the IM_PA image is computed. 
  
  NORM_FACTOR_SCIENCE: [1] # factor by which the science image/cube can be normalized. Default value: 1. One value per observation.

The **center of the psf and science data** is supposed to be at (n//2, n//2), where n is the size of the image in x and y directions, starting the count at 0. If this is indeed the case, the keywords ``SPATIAL_SHIFT_PSF_DATA_ALL`` and ``SPATIAL_SHIFT_SCIENCE_DATA_ALL`` should be set to 0. Otherwise, set ``SPATIAL_SHIFT_PSF_DATA_ALL`` to the number of pixels to offset the image. 

.. note:: 

  - For SPHERE/IRDIS polarized intensity data processed with IRDAP, SPATIAL_SHIFT_*_DATA = 0.5.
  - For SPHERE/IRDIS or SPHERE/IFS total intensity data preprocessed with SpeCal, SPATIAL_SHIFT_*_DATA = 0.

For example, ``SPATIAL_SHIFT_*_DATA_ALL`` = 0.5 means that the center of the image is at (n//2 + 0.5, n//2 + 0.5).

In practice, I do not renormalized the science image/cube, so I let ``NORM_FACTOR_SCIENCE`` = 1.
  










