#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script to plot the results of the MCMC simulations, which constrained the disk morphology and brightness. 

(!) Don't forget to update the keyword RESULTDIR in the configuration file job*.yaml to the path where is located the folder '/results_MCMC'. For example, if in the configuration file job*.yaml, SAVING_DIR was set to 'results/results_BBH_polar/' and the MCMC simulation was launched on the date and time '2025-2-3-16h45min21s', RESULTDIR should be set to
RESULTDIR: 'results/results_BBH_polar/2025-2-3-16h45min21s/' 
'''
 
__author__ = 'Celia Desgrange'

####### This is the MCMC plotting code for HR 4796 data #######

# Import classic packages
from MoDiSc.packages import * 

# Import MoDiSc functions 
from MoDiSc.functions.dataset_load_and_process import *
from MoDiSc.functions.noise_map_derive import *
from MoDiSc.functions.simulations import *
from MoDiSc.functions.mcmc_plot_results import *
from MoDiSc.functions.utils import *
from MoDiSc.functions.chisquare_test import *
from MoDiSc.functions.parameters_rearrange_print import *


# Define relevamt class, functions to write the info and error log files
# Info log file
class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a+")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass

# Error log file
def log_error(fn, message):
    with open("{}_error_{}.log".format(fn, date), "a") as myfile:
        myfile.write(f'ERROR : {message} \n')


if __name__ == '__main__':
    ##################################################################################################
    ### INITIALIZATION: define relevant parameters to plot and save the results of the simulations ###
    ##################################################################################################

     # Define the date (will be used to keep track of the simulation)
    L = time.localtime()
    date = "{}-{}-{}-{}h{}min{}s".format(L[0],L[1],L[2],L[3],L[4],L[5],L[6])
     
    # Uncomment the three following lines to remove warnings
    # warnings.filterwarnings("ignore", category=RuntimeWarning)
    # warnings.filterwarnings("ignore", category=UserWarning)
    # warnings.simplefilter('ignore', category=AstropyWarning)

    # Open the configuration file of the simulation
    config_filename = sys.argv[1] 
    # Note: the filename of the configuration file is given as an argument when launching the python script
    with open(config_filename, 'r') as yaml_file:
        config_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

    # Load which information are going to be printed when this script is launched
    DISPLAY_GENERAL_INFO   = config_file['DISPLAY_GENERAL_INFO']   # print the information about loading the dataset(s), parameters of the simulations and general status of the simulations
    DISPLAY_INFO_SIMU_MCMC = config_file['DISPLAY_INFO_SIMU_MCMC'] # print the information when MCMC simulations are running (this should be set to 1 only when testing if the simulations runs well, because this is time-consuming to print the information for all the MCMC iterations and walkers)
    DISPLAY_INFO_SIMU_NELDERMEAD = config_file['DISPLAY_INFO_SIMU_NELDERMEAD'] #  print the information when Nelder-Mead simulations are running

    if DISPLAY_GENERAL_INFO: print('\n=== Initialization ===')

    # Test on which machine I am
    if DISPLAY_GENERAL_INFO: print(f'Name of the machine: {socket.gethostname()}')
    EXTERNAL_SERVER = config_file['EXTERNAL_SERVER'] #  if True (or 1), the script is run on another server 
    # Advice: this script can be run on one's laptop without issues. Just make sure to have the results of the MCMC simulations (results/[date of the simulation]/results_MCMC/*backend_file_mcmc.h5) and the configuration file of the simulation (config_files/job*.yaml) locally.

    ## Initialize paths
    if EXTERNAL_SERVER: # Case: this script is run on another machine
        DATADIR_ALL = config_file['DATADIR_EXT'] # path to the folder in which are located the files of the observations
        # Regarding FN_PSF_*: Filename of the non-coronagraphic image of the star (PSF). There are three variable for the filename of the PSF, depending if there is one or two PSF files. Specifically, for SPHERE polarized intensity observations processed with IRDAP, there are two distinct PSF files, corresponding to the left and right part of the SPHERE/IRDIS detector. In this case, the filenames where to look for the PSF files are "FN_PSF_1_ALL" and "FN_PSF_2_ALL". Otherwise "FN_PSF_ALL". This is done automatically in the script "functions_load_and_process_dataset".
        FN_PSF_1_ALL   = config_file['FN_PSF_1_EXT'] # should be used in the case of SPHERE polarized intensity data processed with IRDAP
        FN_PSF_2_ALL   = config_file['FN_PSF_2_EXT'] # should be used in case of SPHERE polarized intensity data processed with IRDAP
        FN_PSF_ALL     = config_file['FN_PSF_EXT']   # probably used in all the other the cases (e.g., SPHERE data processed with SpeCal)
        FN_PA_ALL      = config_file['FN_PA_EXT']    # filename(s) of the parallactic angles file(s). Relevant for pupil-tracking stabilized observations, or observations acquired for several rolling angles. The FN_PA_ALL may be set to [..., None, ...] for other types of observations.
        FN_SCIENCE_ALL = config_file['FN_SCIENCE_EXT'] # filename(s) of the science data. Science data are pre- or post-processed, depending if they should be post-processed when doing the simulations to look for the best disk model. Example: FN_SCIENCE_ALL should indicate for SPHERE polarized intensity data the IRDAP post-processed science data, but for SPHERE pupil-stabilized observations, the pre-processed data, to take into account the self-subtraction (Milli+2012) effect.
        FN_REF_ALL     = config_file['FN_CUBE_REF_EXT'] # filename(s) of the reference cube data. If no reference cube, the value is [..., None, ...]. Reference cube data are pre-processed, and will be used to post-process total intensity data using the RDI(+ADI) technique(s).
        FN_NOISE_ALL   = config_file['FN_NOISE_EXT']   # filename(s) of the noise map. The noise map can be provided or computed later in the script. In the latter case, FN_NOISE_ALL is set to [..., None, ...] and COMPUTE_NOISE_MAP_ALL should be set to [..., 1, ...]
        PATH_MASK_ALL  = config_file['PATH_MASK_EXT'] # path(s) to the mask(s). The mask corresponds to the region where the disk model will be optimized to match the SCIENCE_DATA. 
        progress = False  # if on an external machine, do not display the MCMC progress bar

    else: # Case: this script is run on one own's machine
        # See the comments in the previous block "if EXTERNAL_SERVER"
        DATADIR_ALL = config_file['DATADIR_INT']
        FN_PSF_ALL  = config_file['FN_PSF_INT']
        FN_PSF_1_ALL     = config_file['FN_PSF_1_INT'] 
        FN_PSF_2_ALL     = config_file['FN_PSF_2_INT']
        FN_PA_ALL        = config_file['FN_PA_INT']
        FN_SCIENCE_ALL   = config_file['FN_SCIENCE_INT']
        FN_REF_ALL       = config_file['FN_CUBE_REF_EXT'] 
        FN_NOISE_ALL     = config_file['FN_NOISE_INT']
        PATH_MASK_ALL    = config_file['PATH_MASK_INT']
        progress = True  # if on my local machine, display the MCMC pro

    # Defined the saving directories: define the main and associated paths and create the corresponding folders
    RESULTDIR = config_file['RESULTDIR']
    RESULTDIR_MCMC  = os.path.join(RESULTDIR, 'results_MCMC/')
    RESULTDIR_FINAL_FILES  = os.path.join(RESULTDIR + 'final/')
    RESULTDIR_FINAL_FILES_CHI2MIN = os.path.join(RESULTDIR + 'final/chi2min/')
    RESULTDIR_FINAL_FILES_MEDIAN  = os.path.join(RESULTDIR + 'final/median/')
    os.makedirs(RESULTDIR_FINAL_FILES,exist_ok=True)
    os.makedirs(RESULTDIR_FINAL_FILES_CHI2MIN,exist_ok=True)
    os.makedirs(RESULTDIR_FINAL_FILES_MEDIAN,exist_ok=True)
   
    RESULTDIR_INPUT_FILES = os.path.join(RESULTDIR + 'inputs/')
    os.makedirs(RESULTDIR_INPUT_FILES,exist_ok=True)
    RESULTDIR_FIRSTGUESS_FILES = os.path.join(RESULTDIR + 'first_guess/')
    os.makedirs(RESULTDIR_FIRSTGUESS_FILES,exist_ok=True)
    RESULTDIR_INTERMEDIATE_FILES = os.path.join(RESULTDIR + 'intermediate/')
    os.makedirs(RESULTDIR_INTERMEDIATE_FILES,exist_ok=True)

    # Log file: Start it
    fn_log = os.path.join(RESULTDIR, f"log_plotbackend_{config_filename[len('config_files/'):-len('.yaml')]}" )
    fn_log_info = f"{fn_log}_info_{date}.log"
    sys.stdout  = Logger(fn_log_info)
    if DISPLAY_GENERAL_INFO: 
        print(f'\nWrite a logfile with all printed infos at \n{fn_log_info}')
        print(f'\nThe filename of the configuration file is \n{config_filename}')

    # Copy configuration file in the results folder to keep track of the parameters used for the simulation
    file_destination = RESULTDIR_INPUT_FILES+'/'+'job_including_path_resultdir_updated'
    os.makedirs(file_destination, exist_ok=True)
    if DISPLAY_GENERAL_INFO: print(f'\nCopy the configuration file as well at the path:\n{file_destination}')
    shutil.copy(config_filename, file_destination)

    ## Initialize variables ##
    # System
    DISTANCE_STAR = config_file['DISTANCE_STAR'] # distance to the system of interest

    ## Observation
    EPOCHS_ALL            = config_file['EPOCHS']           # list of epochs(s) corresponding to the observation(s). One value per observation.
    NB_OBS                = len(EPOCHS_ALL)                 # number of observations
    config_file['NB_OBS'] = NB_OBS
    INSTRU_ALL            = config_file['INSTRU']           # list of instrument(s) corresponding to the observation(s). One value per observation.
    if NB_OBS != len(INSTRU_ALL) : raise ValueError(f'The value of NB_OBS ({NB_OBS}) should be equal to the number of elements in INSTRU ({len(INSTRU_ALL)})...')
    SPECTRAL_AXIS_ALL = config_file['SPECTRAL_AXIS'] # if [..., 1, ...], there is a spectral axis for the PSF and SCIENCE data, if [..., 0, ...], there is not
    CHANNELS_ALL          = config_file['CHANNELS']  # list of list of spectral channels (one channel = one wavelength) to be considered. Example: For 2 observations, CHANNELS_ALL = [[0],[0,1]] indicates that for the first observation, only the first spectral will be considered, whereas for the second observation, both the first and second channels will be considered.
    PLATE_SCALE_ALL  = config_file['PLATE_SCALE']    # list of the plate scale values. One value per observation.
    TYPE_OBS_ALL      = config_file['TYPE_OBS']      # list of the type of the observation: total intensity ('total_intensity') or polarimetry ('polarized_intensity'). One value per observation.
    TWO_PSF_FILES_ALL = config_file['TWO_PSF_FILES'] # list of boolean (or 0/1) indicating whether there are two different files (= located at two different paths) to consider for the PSF (True = 1 means yes, False = 0 means no). One value per observation.

    ## Processing
    # Except the PSF image, the other images/cubes should ultimately have the same spatial shape. The "CROP_*" parameters indicate the number of pixels to remove both in left-right, top-bottom directions.
    CROP_SCIENCE_ALL = config_file['CROP_SCIENCE'] # list of cropping parameter for the spatial dimensions of the science cube/image. One value per observation.
    CROP_NOISE_ALL   = config_file['CROP_NOISE']   # list of cropping parameter for the spatial dimensions of the noise cube/image. One value per observation. 
    CROP_MASK_ALL    = config_file['CROP_MASK']    # list of cropping parameter for the spatial dimensions of the mask image. One value per observation.
    CROP_PSF_ALL     = config_file['CROP_PSF']     # list of cropping parameter for the spatial dimensions of the PSF image. One value per observation.
    CROP_REF_ALL     = config_file['CROP_REF']     # list of cropping parameter for the spatial dimensions of the ref cube. One value per observation.


    # Center of the PSF / SCIENCE DATA image: The center of the PSF and SCIENCE DATA is supposed to be at (n//2, n//2), where n is the size of the image in x and y directions, starting to count at 0. If this is the case,  SPATIAL_SHIFT_PSF_DATA_ALL and SPATIAL_SHIFT_SCIENCE_DATA_ALL should be set to 0. Otherwise, set SPATIAL_SHIFT_PSF_DATA_ALL to the number of pixels to offset the image. # Example: SPATIAL_SHIFT_*_DATA_ALL = 0.5 means that the center of the image is at (n//2 + 0.5, n//2 + 0.5).
    # Remark: for polarized intensity IRDAP post-processed data, SPATIAL_SHIFT_*_DATA = 0.5.
    # Remark: for total intensity SpeCal pre-processed data, SPATIAL_SHIFT_*_DATA = 0.
    SPATIAL_SHIFT_PSF_DATA_ALL = config_file['SPATIAL_SHIFT_PSF_DATA']          # list of floats indicating the number of pixels to offset the PSF image. One value per observation.
    SPATIAL_SHIFT_SCIENCE_DATA_ALL = config_file['SPATIAL_SHIFT_SCIENCE_DATA']  # list of floats indicating the number of pixels to offset the science data. One value per observation. In practice, this is use either 
# - in the case of polarized intensity data and if DO_ROBUST_CONVOLUTION is set to [..., 1, ...], when the IM_PA image is computed
# - in the case of total intensity data and if COMPUTE_NOISE_MAP is set to [..., 1, ...], when the NOISE_MAP is comuted
    IWA_ALL                 = config_file['IWA']     # list of the radius of the inner working angle in pixels. One value per observation.
    NORM_FACTOR_SCIENCE_ALL = config_file['NORM_FACTOR_SCIENCE'] # list of factors by which the science and noise images/cubes can be normalized. Default value: [..., 1, ...].  One value per observation.
    COMPUTE_NOISE_MAP_ALL   = config_file['COMPUTE_NOISE_MAP']   # list of booleans (or 0/1) indicating whether the noise map should be computed or is already provided. One value per observation. True = 1 means yes, compute the noise map from the SCIENCE_DATA. False = 0 means no, load it from the path DATADIR + FN_NOISE.
    NOISE_MULTIPLICATION_FACTOR_ALL = config_file['NOISE_MULTIPLICATION_FACTOR'] # list of floats. One value per observation. Thi multiplication factor can be used to artificially increase the value of the noise cube/image. See e.g. Mazoyer et al. 2020 in their SPIE paper about diskFM. Default value: 1
    RUN_POSTPROCESSING_TECHNIQUE_ALL = config_file['RUN_POSTPROCESSING_TECHNIQUE'] # algorithm used to post-process the data. In the case of total-intensity observations, acquired in pupil-stabilized mode, i.e., with various parallactic angles, the value should be set to 'PCA-ADI'. Otherwise, the value should be set to None.
    # regarding total intensity data to be postprocessed:
    NB_MODES_ALL = config_file['NB_MODES']       # number of modes/components to use when applying PCA
    # regarding polarized intensity data:
    DO_ROBUST_CONVOLUTION_ALL = config_file['DO_ROBUST_CONVOLUTION'] # list of booleans (or 0/1) indicating whether for polarized intensity data, the robust convolution should be made (see Heikamp & Keller 2019). One value per observation. True = 1 means yes, do the robust convolution. False = 0 means no, don't do it, instead it will do the classic convolution with convolve_fft() (see function chisquare() in the script functions/simulations.py)
    
    ## Disk modelling ##
    ## Names, values and bounds of the parameters
    PARAMS_NAMES_LIST_OF_DICT   = config_file['PARAMS_NAMES'] # names of the parameters of the disk to be fitted. Misc type, it is a list of dictionary.
    PARAMS_NAMES_UNFLAT         = from_param_list_of_dict_to_param_list_of_list(PARAMS_NAMES_LIST_OF_DICT) # convert the list of dictionary into a list of list. 
    NB_DISK_STRUCTURES          = len(PARAMS_NAMES_UNFLAT)    # define the number of disk structures 
    config_file['PARAMS_NAMES'] = PARAMS_NAMES_UNFLAT         # for conveniency, update in the configuration file the value of PARAMS_NAMES.
    PARAMS_NAMES_FLAT           = from_params_names_unflat_to_params_names_flat(config_file) # for conveniency (useful later in the script), also defined a flat list of PARAMS_NAMES

    ## Load the unflat list of initial parameters and their bounds. 
    try:
        PARAMS_INIT_UNFLAT   = from_params_names_to_param_init(  config_file, shape_output='unflat') # first guess of the to-be-fitted parameters
        PARAMS_BOUNDS_UNFLAT = from_params_names_to_param_bounds(config_file, shape_output='unflat') # authorized ranges of values of the to-be-fitted parameters
    
    except: # if an error is going to be raised, print some additional information before to help to understand the issue
        if DISPLAY_GENERAL_INFO: print(f'\nList of parameters: {PARAMS_NAMES_UNFLAT} \nTherefore, {NB_DISK_STRUCTURES} disk structure(s) will be modeled.')
        PARAMS_INIT_UNFLAT   = from_params_names_to_param_init(  config_file, shape_output='unflat') 
        PARAMS_BOUNDS_UNFLAT = from_params_names_to_param_bounds(config_file, shape_output='unflat') 

    # For conveniency (useful later in the script, in particular when running the MCMC- or Nelder-Mead-based simulations), also defined a flat list of PARAMS_INIT and PARAMS_BOUNDS
    PARAMS_INIT_FLAT   = from_params_names_to_param_init(  config_file, shape_output='flat') # first guess of the to-be-fitted parameters
    PARAMS_BOUNDS_FLAT = from_params_names_to_param_bounds(config_file, shape_output='flat') # authorized ranges of values of the to-be-fitted parameters
    
    if DISPLAY_GENERAL_INFO: print_params_names_and_values_and_bounds(PARAMS_NAMES_UNFLAT, PARAMS_INIT_UNFLAT, PARAMS_BOUNDS_UNFLAT, params_flat=False)

    PARAMS_LABELS_LIST_OF_DICT = config_file['PARAMS_LABELS']

    PARAMS_LABELS_UNFLAT         = from_param_list_of_dict_to_param_list_of_list(PARAMS_LABELS_LIST_OF_DICT) # convert the list of dictionary into a list of list. 
    config_file['PARAMS_LABELS'] = PARAMS_LABELS_UNFLAT         # for conveniency, update in the configuration file the value of PARAMS_NAMES.
    PARAMS_LABELS_FLAT           = from_params_names_unflat_to_params_names_flat(config_file, config_file_param='PARAMS_LABELS') # for conveniency (useful later in the script), also defined a flat list of PARAMS_NAMES

    EXPLORATION_ALGO       = config_file['EXPLORATION_ALGO']        # algorithm used to explore the parameter space: "MCMC" or "Nelder-Mead"
    NB_FREE_PARAMS_TOT     = int(config_file['NB_FREE_PARAMS_TOT']) # total number of free parameters 
    NB_FREE_PARAMS_PER_OBS = config_file['NB_FREE_PARAMS_PER_OBS']  # list of number of free parameters per observation. One value per observation.
    
    CONVENTION_UNIT      = config_file['CONVENTION_UNIT']  # two different conventions of unit are used: "MCMC" or "user-friendly". If equal CONVENTION_UNIT is set to 'MCMC', the value of the inclination corresponds to cos(inclination) and the value of the flux scaling factor to 10**(flux_scaling_factor)
    DISK_MODEL_POLAR_ALL = config_file['DISK_MODEL_POLAR'] # list of boolean (or 0/1) indicating whether Rayleigh scattering should be used in the function vip_hci.fm.scattered_light_disk.ScatteredLightDisk(). If True = 1, it means that in the function ScatteredLightDisk(), the argument spf_dico would be set to {...,'polar': 1,...}. If False = 0, no Rayleigh scattering would be considered in the function ScatteredLightDisk(), so the argument spf_dico would be set to {...,'polar': 0,...}. One value per observation.
    SAVE_SOME_RESULTS = config_file['SAVE_SOME_RESULTS']   # equal to True (or 1) if yes, otherwise equal to False (or 0).
    SAVE_FULL_RESULTS = config_file['SAVE_FULL_RESULTS']   # equal to True (or 1) if yes, otherwise equal to False (or 0).

    
    # Parameters for the first disk structure modeled
    RAD_INIT_STRUCT1     = config_file['RAD_INIT_STRUCT1']     # disk radius r0 (au)
    PA_INIT_STRUCT1      = config_file['PA_INIT_STRUCT1']      # disk position angle (deg)
    INC_INIT_STRUCT1     = config_file['INC_INIT_STRUCT1']     # disk inclination (deg or cos(inclination) value, depending of the value of CONVENTION_UNIT)
    G1_INIT_STRUCT1      = config_file['G1_INIT_STRUCT1']      # first anisotropic parameter in the Henyey-Greenstein (HG) function
    G2_INIT_STRUCT1      = config_file['G2_INIT_STRUCT1']      # second anisotropic parameter in the function
    ARGPERI_INIT_STRUCT1 = config_file['ARGPERI_INIT_STRUCT1'] # argument of periastron 
    ECC_INIT_STRUCT1     = config_file['ECC_INIT_STRUCT1']     # eccentricity of the disk
    KSI0_INIT_STRUCT1    = config_file['KSI0_INIT_STRUCT1']    # disk scale height
    GAMMA_INIT_STRUCT1   = config_file['GAMMA_INIT_STRUCT1']   # disk vertical profile
    BETA_INIT_STRUCT1    = config_file['BETA_INIT_STRUCT1']    # disk flaring
    ALPHA_INIT_STRUCT1   = config_file['ALPHA_INIT_STRUCT1']   # weight between the two anisotropic parameters in the HG function (if relevant)
    AIN_INIT_STRUCT1     = config_file['AIN_INIT_STRUCT1']     # inner slope of the disk belt
    AOUT_INIT_STRUCT1    = config_file['AOUT_INIT_STRUCT1']    # outer slope of the disk belt
    SCALING_INIT_STRUCT1 = config_file['SCALING_INIT_STRUCT1'] # flux scaling parameter of the disk belt (to match the observation)

    if NB_DISK_STRUCTURES > 1:  # Parameters for the second disk structure modeled
        RAD_INIT_STRUCT2     = config_file['RAD_INIT_STRUCT2']     # disk radius r0 (au)
        PA_INIT_STRUCT2      = config_file['PA_INIT_STRUCT2']      # disk position angle (deg)
        INC_INIT_STRUCT2     = config_file['INC_INIT_STRUCT2']     # disk inclination (deg or cos(inclination) value, depending of the value of CONVENTION_UNIT)
        G1_INIT_STRUCT2      = config_file['G1_INIT_STRUCT2']      # first anisotropic parameter in the Henyey-Greenstein (HG) function
        G2_INIT_STRUCT2      = config_file['G2_INIT_STRUCT2']      # second anisotropic parameter in the function
        ARGPERI_INIT_STRUCT2 = config_file['ARGPERI_INIT_STRUCT2'] # argument of periastron 
        ECC_INIT_STRUCT2     = config_file['ECC_INIT_STRUCT2']     # eccentricity of the disk
        KSI0_INIT_STRUCT2    = config_file['KSI0_INIT_STRUCT2']    # disk scale height
        GAMMA_INIT_STRUCT2   = config_file['GAMMA_INIT_STRUCT2']   # disk vertical profile
        BETA_INIT_STRUCT2    = config_file['BETA_INIT_STRUCT2']    # disk flaring
        ALPHA_INIT_STRUCT2   = config_file['ALPHA_INIT_STRUCT2']   # weight between the two anisotropic parameters in the HG function (if relevant)
        AIN_INIT_STRUCT2     = config_file['AIN_INIT_STRUCT2']     # inner slope of the disk belt
        AOUT_INIT_STRUCT2    = config_file['AOUT_INIT_STRUCT2']    # outer slope of the disk belt
        SCALING_INIT_STRUCT2 = config_file['SCALING_INIT_STRUCT2'] # flux scaling parameter of the disk belt (to match the observation)

    if NB_DISK_STRUCTURES > 2: # Parameters for the third disk structure modeled
        RAD_INIT_STRUCT3     = config_file['RAD_INIT_STRUCT3']     # disk radius r0 (au)
        PA_INIT_STRUCT3      = config_file['PA_INIT_STRUCT3']      # disk position angle (deg)
        INC_INIT_STRUCT3     = config_file['INC_INIT_STRUCT3']     # disk inclination (deg or cos(inclination) value, depending of the value of CONVENTION_UNIT)
        G1_INIT_STRUCT3      = config_file['G1_INIT_STRUCT3']      # first anisotropic parameter in the Henyey-Greenstein (HG) function
        G2_INIT_STRUCT3      = config_file['G2_INIT_STRUCT3']      # second anisotropic parameter in the function
        ARGPERI_INIT_STRUCT3 = config_file['ARGPERI_INIT_STRUCT3'] # argument of periastron 
        ECC_INIT_STRUCT3     = config_file['ECC_INIT_STRUCT3']     # eccentricity of the disk
        KSI0_INIT_STRUCT3    = config_file['KSI0_INIT_STRUCT3']    # disk scale height
        GAMMA_INIT_STRUCT3   = config_file['GAMMA_INIT_STRUCT3']   # disk vertical profile
        BETA_INIT_STRUCT3    = config_file['BETA_INIT_STRUCT3']    # disk flaring
        ALPHA_INIT_STRUCT3   = config_file['ALPHA_INIT_STRUCT3']   # weight between the two anisotropic parameters in the HG function (if relevant)
        AIN_INIT_STRUCT3     = config_file['AIN_INIT_STRUCT3']     # inner slope of the disk belt
        AOUT_INIT_STRUCT3    = config_file['AOUT_INIT_STRUCT3']    # outer slope of the disk belt
        SCALING_INIT_STRUCT3 = config_file['SCALING_INIT_STRUCT3'] # flux scaling parameter of the disk belt (to match the observation)

    if NB_DISK_STRUCTURES > 3: # Parameters for the fourth disk structure modeled
        RAD_INIT_STRUCT4     = config_file['RAD_INIT_STRUCT4']     # disk radius r0 (au)
        PA_INIT_STRUCT4      = config_file['PA_INIT_STRUCT4']      # disk position angle (deg)
        INC_INIT_STRUCT4     = config_file['INC_INIT_STRUCT4']     # disk inclination (deg or cos(inclination) value, depending of the value of CONVENTION_UNIT)
        G1_INIT_STRUCT4      = config_file['G1_INIT_STRUCT4']      # first anisotropic parameter in the Henyey-Greenstein (HG) function
        G2_INIT_STRUCT4      = config_file['G2_INIT_STRUCT4']      # second anisotropic parameter in the function
        ARGPERI_INIT_STRUCT4 = config_file['ARGPERI_INIT_STRUCT4'] # argument of periastron 
        ECC_INIT_STRUCT4     = config_file['ECC_INIT_STRUCT4']     # eccentricity of the disk
        KSI0_INIT_STRUCT4    = config_file['KSI0_INIT_STRUCT4']    # disk scale height
        GAMMA_INIT_STRUCT4   = config_file['GAMMA_INIT_STRUCT4']   # disk vertical profile
        BETA_INIT_STRUCT4    = config_file['BETA_INIT_STRUCT4']    # disk flaring
        ALPHA_INIT_STRUCT4   = config_file['ALPHA_INIT_STRUCT4']   # weight between the two anisotropic parameters in the HG function (if relevant)
        AIN_INIT_STRUCT4     = config_file['AIN_INIT_STRUCT4']     # inner slope of the disk belt
        AOUT_INIT_STRUCT4    = config_file['AOUT_INIT_STRUCT4']    # outer slope of the disk belt
        SCALING_INIT_STRUCT4 = config_file['SCALING_INIT_STRUCT4'] # flux scaling parameter of the disk belt (to match the observation)

    if NB_DISK_STRUCTURES > 4: # Parameters for the fifth disk structure modeled
        RAD_INIT_STRUCT5     = config_file['RAD_INIT_STRUCT5']     # disk radius r0 (au)
        PA_INIT_STRUCT5      = config_file['PA_INIT_STRUCT5']      # disk position angle (deg)
        INC_INIT_STRUCT5     = config_file['INC_INIT_STRUCT5']     # disk inclination (deg or cos(inclination) value, depending of the value of CONVENTION_UNIT)
        G1_INIT_STRUCT5      = config_file['G1_INIT_STRUCT5']      # first anisotropic parameter in the Henyey-Greenstein (HG) function
        G2_INIT_STRUCT5      = config_file['G2_INIT_STRUCT5']      # second anisotropic parameter in the function
        ARGPERI_INIT_STRUCT5 = config_file['ARGPERI_INIT_STRUCT5'] # argument of periastron 
        ECC_INIT_STRUCT5     = config_file['ECC_INIT_STRUCT5']     # eccentricity of the disk
        KSI0_INIT_STRUCT5    = config_file['KSI0_INIT_STRUCT5']    # disk scale height
        GAMMA_INIT_STRUCT5   = config_file['GAMMA_INIT_STRUCT5']   # disk vertical profile
        BETA_INIT_STRUCT5    = config_file['BETA_INIT_STRUCT5']    # disk flaring
        ALPHA_INIT_STRUCT5   = config_file['ALPHA_INIT_STRUCT5']   # weight between the two anisotropic parameters in the HG function (if relevant)
        AIN_INIT_STRUCT5     = config_file['AIN_INIT_STRUCT5']     # inner slope of the disk belt
        AOUT_INIT_STRUCT5    = config_file['AOUT_INIT_STRUCT5']    # outer slope of the disk belt
        SCALING_INIT_STRUCT5 = config_file['SCALING_INIT_STRUCT5'] # flux scaling parameter of the disk belt (to match the observation)

    if NB_DISK_STRUCTURES > 5: # Parameters for the second disk structure modeled
        raise ValueError(f'To data, MoDiSc can modeled a disk in scattered light with up to 5 disk structures. If you would like to model more than five disk structures (wow!), please reach out Célia Desgrange (celia.desgrange@eso.org). She will update MoDiSc. Note: the number of disk structures given in the configuration file {config_filename} in PARAMS_NAMES is {NB_DISK_STRUCTURES}.')
    
   
    # Stored the disk parameters in a dico
    dico_params = {'rad_init_struct1': RAD_INIT_STRUCT1,
                'pa_init_struct1'        : PA_INIT_STRUCT1,
                'inc_init_struct1'       : INC_INIT_STRUCT1,
                'g1_init_struct1'        : G1_INIT_STRUCT1,
                'g2_init_struct1'        : G2_INIT_STRUCT1,
                'argperi_init_struct1'   : ARGPERI_INIT_STRUCT1,
                'ecc_init_struct1'       : ECC_INIT_STRUCT1,
                'ksi0_init_struct1'      : KSI0_INIT_STRUCT1,
                'gamma_init_struct1'     : GAMMA_INIT_STRUCT1,
                'beta_init_struct1'      : BETA_INIT_STRUCT1,
                'alpha_init_struct1'     : ALPHA_INIT_STRUCT1,
                'ain_init_struct1'       : AIN_INIT_STRUCT1,
                'aout_init_struct1'      : AOUT_INIT_STRUCT1,
                'scaling_init_struct1'   : SCALING_INIT_STRUCT1}
        
    if NB_DISK_STRUCTURES > 1: dico_params.update({'rad_init_struct2': RAD_INIT_STRUCT2,
            'pa_init_struct2'        : PA_INIT_STRUCT2,
            'inc_init_struct2'       : INC_INIT_STRUCT2,
            'g1_init_struct2'        : G1_INIT_STRUCT2,
            'g2_init_struct2'        : G2_INIT_STRUCT2,
            'argperi_init_struct2'   : ARGPERI_INIT_STRUCT2,
            'ecc_init_struct2'       : ECC_INIT_STRUCT2,
            'ksi0_init_struct2'      : KSI0_INIT_STRUCT2,
            'gamma_ini_struct2'      : GAMMA_INIT_STRUCT2,
            'beta_init_struct2'      : BETA_INIT_STRUCT2,
            'alpha_init_struct2'     : ALPHA_INIT_STRUCT2,
            'ain_init_struct2'       : AIN_INIT_STRUCT2,
            'aout_init_struct2'      : AOUT_INIT_STRUCT2,
            'scaling_init_struct2'   : SCALING_INIT_STRUCT2})
    
    if NB_DISK_STRUCTURES > 2: dico_params.update({'rad_init_struct3': RAD_INIT_STRUCT3,
            'pa_init_struct3'        : PA_INIT_STRUCT3,
            'inc_init_struct3'       : INC_INIT_STRUCT3,
            'g1_init_struct3'        : G1_INIT_STRUCT3,
            'g2_init_struct3'        : G2_INIT_STRUCT3,
            'argperi_init_struct3'   : ARGPERI_INIT_STRUCT3,
            'ecc_init_struct3'       : ECC_INIT_STRUCT3,
            'ksi0_init_struct3'      : KSI0_INIT_STRUCT3,
            'gamma_ini_struct3'      : GAMMA_INIT_STRUCT3,
            'beta_init_struct3'      : BETA_INIT_STRUCT3,
            'alpha_init_struct3'     : ALPHA_INIT_STRUCT3,
            'ain_init_struct3'       : AIN_INIT_STRUCT3,
            'aout_init_struct3'      : AOUT_INIT_STRUCT3,
            'scaling_init_struct1'   : SCALING_INIT_STRUCT3})
    
    if NB_DISK_STRUCTURES > 3: dico_params.update({'rad_init_struct4': RAD_INIT_STRUCT4,
            'pa_init_struct4'        : PA_INIT_STRUCT4,
            'inc_init_struct4'       : INC_INIT_STRUCT4,
            'g1_init_struct4'        : G1_INIT_STRUCT4,
            'g2_init_struct4'        : G2_INIT_STRUCT4,
            'argperi_init_struct4'   : ARGPERI_INIT_STRUCT4,
            'ecc_init_struct4'       : ECC_INIT_STRUCT4,
            'ksi0_init_struct4'      : KSI0_INIT_STRUCT4,
            'gamma_ini_struct4'      : GAMMA_INIT_STRUCT4,
            'beta_init_struct4'      : BETA_INIT_STRUCT4,
            'alpha_init_struct4'     : ALPHA_INIT_STRUCT4,
            'ain_init_struct4'       : AIN_INIT_STRUCT4,
            'aout_init_struct4'      : AOUT_INIT_STRUCT4,
            'scaling_init_struct4'   : SCALING_INIT_STRUCT4})
    
    if NB_DISK_STRUCTURES > 4: dico_params.update({'rad_init_struct5': RAD_INIT_STRUCT5,
            'pa_init_struct5'        : PA_INIT_STRUCT5,
            'inc_init_struct5'       : INC_INIT_STRUCT5,
            'g1_init_struct5'        : G1_INIT_STRUCT5,
            'g2_init_struct5'        : G2_INIT_STRUCT5,
            'argperi_init_struct5'   : ARGPERI_INIT_STRUCT5,
            'ecc_init_struct5'       : ECC_INIT_STRUCT5,
            'ksi0_init_struct5'      : KSI0_INIT_STRUCT5,
            'gamma_ini_struct5'      : GAMMA_INIT_STRUCT5,
            'beta_init_struct5'      : BETA_INIT_STRUCT5,
            'alpha_init_struct5'     : ALPHA_INIT_STRUCT5,
            'ain_init_struct5'       : AIN_INIT_STRUCT5,
            'aout_init_struct5'      : AOUT_INIT_STRUCT5,
            'scaling_init_struct5'   : SCALING_INIT_STRUCT5})
    
    if DISPLAY_GENERAL_INFO: 
        print('\n\n(!) Check the parameters (!)')
        print('\n- System params:')
        print('The star is located at',  DISTANCE_STAR, 'pc.')
        
        print('\n- Observational params:')
        print(f'There are {NB_OBS} observation(s) of type {TYPE_OBS_ALL}.')
        if 'total_intensity' in TYPE_OBS_ALL:
            print(f"Spectral channels considered are {CHANNELS_ALL}.")
            print(f'PCA: {NB_MODES_ALL} modes are used.')
        print(f'\nThe plate scale is {PLATE_SCALE_ALL} arcsecond/pix.')
        print(f'The cropping parameters are in both direction bottom/top, left/right: {CROP_SCIENCE_ALL} (science image/cube), {CROP_NOISE_ALL} (noise map), {CROP_MASK_ALL} (mask map), and {CROP_PSF_ALL} (psf image).')  

        print('\n- Processing params:')
        print(f'The IWA used is: {IWA_ALL}')      
        print(f'The normalization factor for the science and noise cube/image is {NORM_FACTOR_SCIENCE_ALL}.')
        print(f'The exploration algo used is: {EXPLORATION_ALGO}')
        print(f'The normalization factor for the noise is {NOISE_MULTIPLICATION_FACTOR_ALL}.')

    ## Name of the backend file, where the results from the MCMC exploration are stored
    file_prefix = '' #str(np.where(NOBS ==1, EPOCH_ALL, str(NOBS) + 'epochs') )  
    name_h5 = file_prefix + 'backend_file_mcmc'

    MCMC_NB_WALKERS      = config_file['MCMC_NB_WALKERS']       # number of walkers
    MCMC_CHAINS_BURNIN   = config_file['MCMC_CHAINS_BURNIN']    # burnin value
    MCMC_CHAINS_BINNING  = config_file['MCMC_CHAINS_BINNING']   # binning the chains value
    MCMC_CHAINS_APPLY_SELECTION_CRITERION   = config_file['MCMC_CHAINS_APPLY_SELECTION_CRITERION']
    MCMC_CHAINS_LOG_PROB_CRITERION_FRACTION = config_file['MCMC_CHAINS_LOG_PROB_CRITERION_FRACTION']
    MCMC_FIG_CORNER_PLOT_SIGMA          = config_file['MCMC_FIG_CORNER_PLOT_SIGMA']
    MCMC_FIG_CORNER_PLOT_ADD_TEXT_ANNOT = config_file['MCMC_FIG_CORNER_PLOT_ADD_TEXT_ANNOT']
    MCMC_FIG_SIZE_FACTOR = config_file['MCMC_FIG_SIZE_FACTOR']


    path_mcmc_file = os.path.join(RESULTDIR_MCMC, name_h5 + '.h5')
    if not os.path.isfile(path_mcmc_file):
        raise ValueError(f"The mcmc h5 file does not exist at the location: \n {path_mcmc_file}")
    
    fig_use_tex = config_file['FIG_USE_TEX'] # True = 1, False = 0
    if fig_use_tex:
        plt.rcParams.update({
        "text.usetex": True,
        "font.serif": 'cm',
        })
    
    ## Plot MCMC results
    # Define the associated dictionary
    dico_plot_mcmc_results = {
        'display'             : DISPLAY_GENERAL_INFO,
        'params_names_flat'   : PARAMS_NAMES_FLAT,
        'params_labels_flat'  : PARAMS_LABELS_FLAT,
        'convention_unit'     : CONVENTION_UNIT,
        'nb_obs'     : NB_OBS,
        'epochs_all' : EPOCHS_ALL,
        'instru_all': INSTRU_ALL,
        'mcmc_chains_burnin': MCMC_CHAINS_BURNIN,
        'mcmc_nb_walkers'     : MCMC_NB_WALKERS,
        'mcmc_nb_dimension'   : NB_FREE_PARAMS_TOT,
        'mcmc_chains_apply_selection_criterion': MCMC_CHAINS_APPLY_SELECTION_CRITERION,
        'mcmc_chains_log_prob_criterion_fraction': MCMC_CHAINS_LOG_PROB_CRITERION_FRACTION,
        'mcmc_chains_binning'           : MCMC_CHAINS_BINNING,
        'mcmc_resultdir'      : RESULTDIR_MCMC,
        'mcmc_fig_size_factor': MCMC_FIG_SIZE_FACTOR,
        'mcmc_fig_cornerplot_sigma': MCMC_FIG_CORNER_PLOT_SIGMA,
        'mcmc_fig_cornerplot_add_text_annot': MCMC_FIG_CORNER_PLOT_ADD_TEXT_ANNOT,
        'nb_params_per_obs'   : NB_FREE_PARAMS_PER_OBS,
        # Additional parameters for the header
        'distance_star'       : DISTANCE_STAR,
        'plate_scale_all'     : PLATE_SCALE_ALL,
        'type_obs_all'        : TYPE_OBS_ALL,
    }


    # Plot the chain values
    make_chains_plot(dico_plot_mcmc_results)

    # Plot the PDFs
    make_corner_plot(dico_plot_mcmc_results)

    if DISPLAY_GENERAL_INFO: print("\n=== Compute the chisquare ===")

    ####################################################################
    ## Load the data necessary to then show the results of the models ##
    ####################################################################
    if config_file['load_data']: 

        if DISPLAY_GENERAL_INFO: print('\n- Load the PSF, (PA_ARRAY), SCIENCE, NOISE, and MASK data ')
        
        dico_psf_pa_science_ref = {
            'config_file': config_file,
            # Following data required to load PSF, PA_ARRAY, SCIENCE_DATA
            'nb_obs'            : NB_OBS,
            'type_obs_all'      : TYPE_OBS_ALL,
            'datadir_all'       : DATADIR_ALL,
            'fn_psf_all'        : FN_PSF_ALL,
            'fn_psf_1_all'      : FN_PSF_1_ALL,
            'fn_psf_2_all'      : FN_PSF_2_ALL,
            'fn_pa_all'         : FN_PA_ALL,
            'fn_science_all'    : FN_SCIENCE_ALL,
            'fn_ref_all'        : FN_REF_ALL,
            'two_psf_files_all' : TWO_PSF_FILES_ALL,
            'spatial_shift_psf_data_all': SPATIAL_SHIFT_PSF_DATA_ALL,
            'crop_psf_all'      : CROP_PSF_ALL,
            'crop_science_all'  : CROP_SCIENCE_ALL,
            'crop_ref_all'      : CROP_REF_ALL,
            'spectral_axis_all' : SPECTRAL_AXIS_ALL,
            'channels_all'      : CHANNELS_ALL,
            'inputs_resultdir'  : RESULTDIR_INPUT_FILES,
            'run_postprocessing_technique_all': RUN_POSTPROCESSING_TECHNIQUE_ALL,
            'display': DISPLAY_GENERAL_INFO,
            }

        # PSF
        PSF_ALL = load_PSF_one_or_several_datasets(dico=dico_psf_pa_science_ref, WRITETO=True)
        
        # PA (load parallactic angle files, only for pupil-stabilized mode observations which will be ADI processed, otherwise set to [..., None, ...])
        PA_ARRAY_ALL = load_PA_ARRAY_one_or_several_datasets(dico=dico_psf_pa_science_ref, WRITETO=True)

        # Science
        SCIENCE_DATA_ALL = load_SCIENCE_one_or_several_datasets(dico=dico_psf_pa_science_ref, WRITETO=False)

        # Reference
        REF_CUBE_ALL     = load_REF_CUBE_one_or_several_datasets(dico=dico_psf_pa_science_ref, WRITETO=False)

        # Noise
        dico_noise_mask = dico_psf_pa_science_ref.copy()

        dico_noise_mask.update(
            {'fn_noise_all'    : FN_NOISE_ALL,
            'path_mask_all'    : PATH_MASK_ALL,
            'crop_noise_all'   : CROP_NOISE_ALL,
            'crop_mask_all'    : CROP_MASK_ALL,
            'compute_noise_map_all'          : COMPUTE_NOISE_MAP_ALL,
            'noise_multiplication_factor_all': NOISE_MULTIPLICATION_FACTOR_ALL,
            'spatial_shift_noise_data_all'   : SPATIAL_SHIFT_SCIENCE_DATA_ALL,
            'iwa_all'          : IWA_ALL,
            'nb_modes_all'     : NB_MODES_ALL,
            'science_data_all' : SCIENCE_DATA_ALL, 
            'pa_array_all'     : PA_ARRAY_ALL,
            'run_postprocessing_technique_all': RUN_POSTPROCESSING_TECHNIQUE_ALL,
            })
        
        NOISE_MAP_ALL = load_NOISE_one_or_several_datasets(dico=dico_noise_mask, WRITETO=True)
        
        # Define/Load the mask within the residuals will be minimized in the science - model data
                # Define/Load the mask within the residuals will be minimized in the science - model data
        MASK2MINIMIZE_ALL = load_MASK2MINIMIZE_one_or_several_datasets(dico=dico_noise_mask)

        DIMENSION_ALL = [np.shape(im)[-1] for im in MASK2MINIMIZE_ALL]
        NB_RESEL_ALL = np.nansum(MASK2MINIMIZE_ALL,axis=(1,2))
        if NB_OBS == 1: NB_RESEL_ALL, DIMENSION_ALL = NB_RESEL_ALL, DIMENSION_ALL
        if DISPLAY_GENERAL_INFO: print(f'The number of pixels considered in the region to be minimized (=MASK2MINIMIZE_ALL) is {NB_RESEL_ALL}')
        
        # [important for polarized intensity observations] Computing an image with for each pixel the corresponding position angle
        IM_PA_ALL = []
        
        for i_obs in range(NB_OBS):
            if TYPE_OBS_ALL[i_obs] in ['polarised intensity', 'polarized intensity', 'polar', 'polarised_intensity', 'polarized_intensity', 'polar'] and DO_ROBUST_CONVOLUTION_ALL[i_obs]: 
                if DISPLAY_GENERAL_INFO: print('Computing an image with for each pixel the corresponding position angle.')
                dim = DIMENSION_ALL[i_obs]
                xc, yc = dim /2 + SPATIAL_SHIFT_SCIENCE_DATA_ALL[i_obs], dim /2 + SPATIAL_SHIFT_SCIENCE_DATA_ALL[i_obs]
                IM_PA_ALL.append( compute_im_pa_grid(np.zeros((dim,dim)), center='coordinates_are_given', xc=xc, yc=yc, counterclockwise=True, return_unit='rad') ) #  necessary only for polar data
            else: IM_PA_ALL.append([None])

        ######################################################
        ######## Processing the total intensity data #########
        ######################################################
        ## Stored the dataset(s) and parameters in a dico

        # Process SCIENCE data (only if we do PCA forward modelling)
        dico_process_data = dico_noise_mask
        dico_process_data['ref_cube_all']      = REF_CUBE_ALL
        dico_process_data['mask2minimize_all'] = MASK2MINIMIZE_ALL
    
        RED_DATA = postprocess_SCIENCE_DATA_PCA_several_datasets(dico=dico_process_data, WRITETO=True)

        ######################################################
        #### Compute the chisquare in some initial cases #####
        ######################################################
        if DISPLAY_GENERAL_INFO: print('\n\n=== Make test before running the MCMC or Nelder-Mead exploration algorithm ===')

        dico_chisquare = {
            'nb_obs'  : NB_OBS,
            'type_obs_all'          : TYPE_OBS_ALL,
            'display'               : DISPLAY_GENERAL_INFO,
            'do_robust_convolution_all': DO_ROBUST_CONVOLUTION_ALL,
            'im_pa_all'             : IM_PA_ALL,
            'psf_all'               : PSF_ALL,
            'mask2minimize_all'     : MASK2MINIMIZE_ALL,
            'noise_map_all'         : NOISE_MAP_ALL,
            'science_data_all'      : SCIENCE_DATA_ALL,
            'ref_cube_all'          : REF_CUBE_ALL,
            'dimension_all'         : DIMENSION_ALL,
            'nb_resel_all'          : NB_RESEL_ALL,
            'nb_free_params_tot'        : NB_FREE_PARAMS_TOT,
            'nb_free_params_per_obs_all': NB_FREE_PARAMS_PER_OBS,
            'plate_scale_all'       : PLATE_SCALE_ALL,
            'distance_star'         : DISTANCE_STAR, 
            'disk_model_polar_all'  : DISK_MODEL_POLAR_ALL,   
            'pa_array_all'          : PA_ARRAY_ALL,
            'iwa_all'               : IWA_ALL,
            'nb_modes_all'            : NB_MODES_ALL,
            'exploration_algo'      : EXPLORATION_ALGO,
            'convention_unit'       : CONVENTION_UNIT,
            'nb_disk_structures'    : NB_DISK_STRUCTURES,
            'run_postprocessing_technique_all': RUN_POSTPROCESSING_TECHNIQUE_ALL,
            # paths
            'save_some_results'     : False,
            'save_full_results'     : False,
            'inputs_resultdir'      : RESULTDIR_INPUT_FILES,
            'simu_resultdir'        : RESULTDIR_FIRSTGUESS_FILES,
            'mcmc_resultdir'        : RESULTDIR_MCMC,
            'params_names_unflat'   : PARAMS_NAMES_UNFLAT,
            'params_example_unflat' : PARAMS_INIT_UNFLAT,
            # parameters free/fixed added just afterwards
            # additional parameters for betadisk
            }
        
        dico_chisquare.update(dico_params)


        # What is the chisquare if the disk model is empty?
        chisquare_init_nodisk = do_test_disk_empty(PARAMS_INIT_UNFLAT,  dico=dico_chisquare)
            
        # What is the chisquare / likelihood of the initial model?
        ## Stored them in a dico
        dico_chisquare.update({ 
                'save_full_results': True,
                'simu_resultdir': RESULTDIR_FIRSTGUESS_FILES,
                })
        chisquare_init = do_test_disk_first_guess(PARAMS_INIT_UNFLAT,  dico=dico_chisquare)
        if DISPLAY_GENERAL_INFO: print("chisquare (no-disk) - chisquare (first-guess): {:.0f}".format((chisquare_init_nodisk-chisquare_init)))

    
        ## Processing
        # Process SCIENCE data (only if we do PCA forward modelling)
        

    ##############################################################################################
    ############# Show the results for the best model and at different percentiles ###############
    ##############################################################################################

    if config_file['show_best_model']: 
        if config_file['load_data'] == 0: raise ValueError(f"The best model can be computed only if the datasets are loaded (config_file['LOAD_DATA'] == 1 or True).")

        dF_chisquare_min = pd.read_csv( os.path.join(RESULTDIR_MCMC, 'params_chi2_min.csv'), sep=';' , index_col=0)
        params_chisquare_min = dF_chisquare_min.iloc[0,:-1].array
        chisquare_min = float(dF_chisquare_min.iloc[0,-1])

        # Change scaling and inclination values to physical
        params_chisquare_min = correct_unit(PARAMS_NAMES_FLAT, params_chisquare_min, convention_unit='user-friendly', output_unit=CONVENTION_UNIT, display=DISPLAY_GENERAL_INFO)

        ## Model for the minimal chisquare value
        if DISPLAY_GENERAL_INFO: 
            print('\n\n=== Parameters for the minimal chisquare value: ==\n')
            print_params_names_and_values_and_bounds(PARAMS_NAMES_FLAT, params_chisquare_min, PARAMS_BOUNDS_FLAT, params_flat = True)

        # Define dico
        dico_best_models = dico_chisquare
        dico_best_models['simu_resultdir'] = RESULTDIR_FINAL_FILES_CHI2MIN
        chisquare2min = chisquare_params_flat_or_unflat(params_chisquare_min, dico=dico_best_models)

        ## Model for disk parameters satifying the 50% percentile
        if DISPLAY_GENERAL_INFO: print('\n\n=== Parameters for the median disk parameters: ===')
        dico_best_models['simu_resultdir'] = RESULTDIR_FINAL_FILES_MEDIAN

        # Define dico. Basically we need parameters to compute the chisquare and to load the results from the MCMC simulations
        dico_best_models = dico_chisquare
        dico_best_models.update(dico_plot_mcmc_results)
    
        params_median = MCMC_parameters_given_percentile(dico=dico_best_models, percentile=[50], WRITETO=False)[0] # list of list [[...]]

        # Correct the unit, should be the nominal one given in the .yaml file
        params_chisquare_min = correct_unit(PARAMS_NAMES_FLAT, params_median, convention_unit='user-friendly', output_unit=CONVENTION_UNIT, display=DISPLAY_GENERAL_INFO)

        #print(theta_median, type(theta_median))
        chisquare2median = chisquare_params_flat_or_unflat(params_median, dico=dico_best_models)

        params_perc   = MCMC_parameters_given_percentile(dico=dico_best_models, percentile=[15.9, 50, 84.1, 23, 97.7, 1, 99.9])
  

