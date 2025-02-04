#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Script to load necessary files corresponding to the datasets in which the disk structure(s) will be characterized> There is also two functions to postprocess the datasets.

Functions:
    . load_PSF_one_dataset()
    . load_PSF_one_or_several_datasets() - Calls the function load_PSF_one_dataset() as many times as the number of observations
    
    . load_SCIENCE_one_dataset()
    . load_SCIENCE_one_or_several_datasets() - Calls the function load_PSF_one_dataset() as many times as the number of observations

    . load_REF_CUBE_one_dataset()
    . load_REF_CUBE_one_or_several_datasets() - Calls the function load_PSF_one_dataset() as many times as the number of observations
    
    . load_NOISE_one_dataset()
    . load_NOISE_one_or_several_datasets() - Calls the function load_PSF_one_dataset() as many times as the number of observations
    
    . load_PA_ARRAY_one_dataset()
    . load_PA_ARRAY_one_or_several_datasets() - Calls the function load_PSF_one_dataset() as many times as the number of observations

    . load_MASK2MINIMIZE_one_dataset()
    . load_MASK2MINIMIZE_one_or_several_datasets() - Calls the function load_PSF_one_dataset() as many times as the number of observations

    . postprocess_SCIENCE_DATA_PCA_one_dataset() 
    . postprocess_SCIENCE_DATA_PCA_one_or_several_datasets  - Calls the function postprocess_SCIENCE_DATA_PCA_one_dataset() as many times as the number of observations

This script would be imported by the scripts:
    ..run_modisc.py
    ..run_mcmc_results.py
'''

__author__ = 'Celia Desgrange'

# Import classic packages
from ..packages import * 

# Import MoDiSc functions
from ..functions.noise_map_derive import *
from ..functions.parameters_rearrange_print import *

##############################################
## Load stuff corresponding to the datasets ##
##############################################

def load_PSF_one_dataset(dico, WRITETO=True, ID=''):
    '''
    Returns the normalized PSF of the dataset. 
    The PSF is normalized by the sum of the pixel values, so np.nansum(PSF) = 1.
    The PSF will be used to convolve the disk model in the MCMC/Nelder-Mead simulations.

    The center of the PSF should be on a given pixel located at x = y = n // 2 with the counting starting at 0 (so NOT between four pixels).

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'datadir', 'fn_psf', 'two_psf_files', 'fn_psf_1', 'fn_psf_2', 'crop_psf', 'spatial_shift_psf_data', 'spectral_axis', 'channels', 'inputs_resultdir'. Note: Depending on the observation, some of these parameters may not be used.
        
        (optional)

        .'WRITETO' (boolean): if True (or 1), save the PSF used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the PSF file.
    
    Outputs:
        .'psf' (2D-array): the normalized PSF.
    '''

    print('\n(Load the PSF DATA)')

    # Access extra variables from dico
    display = dico.get('display', 0)
    DATADIR = dico.get('datadir', 'PROBLEM')
    FN_PSF  = dico.get('fn_psf', 'PROBLEM')
    TWO_PSF_FILES = dico.get('two_psf_files', 'PROBLEM') # indicate whether there are two different files to consider for the PSF (1 = yes, 0 = no)
    FN_PSF_1 = dico.get('fn_psf_1', 'PROBLEM')
    FN_PSF_2 = dico.get('fn_psf_2', 'PROBLEM')
    CROP_PSF = dico.get('crop_psf', 'PROBLEM')
    SPATIAL_SHIFT_PSF_DATA = dico.get('spatial_shift_psf_data','PROBLEM')
    SPECTRAL_AXIS = dico.get('spectral_axis','PROBLEM')
    CHANNELS = dico.get('channels', 'PROBLEM')
    CHANNELS_COMBINED = bool(np.where(len(CHANNELS) > 1, True, False))
    inputs_resultdir  = dico.get('inputs_resultdir', 'PROBLEM')

    if WRITETO: 
        if inputs_resultdir == 'PROBLEM': raise ValueError('There is a problem with the value of inputs_resultdir.')

    if TWO_PSF_FILES: 
        path_psf_1 = os.path.join(DATADIR, FN_PSF_1); path_psf_2 = os.path.join(DATADIR, FN_PSF_2)
        if display: print(f'The path for the PSF 1:\n {path_psf_1}\n'); print(f'The path for the PSF 2:\n {path_psf_2}\n')
        psf_1 = fits.getdata(path_psf_1); psf_2 = fits.getdata(path_psf_2)
        if display: print(f' [INFO] np.shape(psf1) = {np.shape(psf_1)} \n [INFO] np.shape(psf2) = {np.shape(psf_2)} ')
        psf = np.concatenate((psf_1, psf_2))
        if display: print(f' [INFO] np.shape(psf) = {np.shape(psf)} # after concatenation of PSF 1 and PSF 2')
        psf = np.nansum(psf, axis=0)
        if display: print(f' [INFO] np.shape(psf) = {np.shape(psf)} # after summing the PSF')
    
    else:
        path_psf = os.path.join(DATADIR + FN_PSF)
        if display: print(f'The path for the PSF data:\n {path_psf}\n')
        psf = fits.getdata(path_psf)
    
    if SPECTRAL_AXIS: 
        if display:
            print('Consider the spectral channels of interest:')
            if len(np.shape(psf)) == 4: print(f' [INFO] the axis of psf are wavelength, temporal, y, x.')
            elif len(np.shape(psf)) == 3: print(f' [INFO] the axis of psf are wavelength, y, x.')
            if display: print(f' [INFO] np.shape(psf) = {np.shape(psf)} # before selecting the spectral channels')
        psf = psf[CHANNELS] # keep the relevant spectral channels, written in CHANNELS (which is a int or a list of int)
        if display: print(f' [INFO] np.shape(psf) = {np.shape(psf)} # after selecting the spectral channels')
        if CHANNELS_COMBINED:
            psf = np.nansum(psf, axis=0) # sum channels (if CHANNELS_COMBINED set to True), otherwise consider one channel
            if display: print(f' [INFO] np.shape(psf) = {np.shape(psf)} # after summing the spectral channels')
        if len(np.shape(psf)) == 3: 
            psf = np.nansum(psf, axis=0) # sum the PSF on the temporal axis.
            if display: print(f' [INFO] np.shape(psf) = {np.shape(psf)} # after summing along the temporal axis')

    if SPATIAL_SHIFT_PSF_DATA: psf = vip.preproc.recentering.frame_shift(psf,SPATIAL_SHIFT_PSF_DATA,SPATIAL_SHIFT_PSF_DATA) 

    if CROP_PSF !=0: # crop the PSF
        if display: print(f'Crop the PSF: \n [INFO] The size of the PSF data *before* cropping is {np.shape(psf)}.')
        psf = psf[CROP_PSF:-CROP_PSF+1, CROP_PSF:-CROP_PSF+1]
        if display: print(f' [INFO] The size of the PSF data *after* cropping is {np.shape(psf)}.')
    else:
        if display: print(f'The size of the PSF data is: {np.shape(psf)}')
    
    # Normalize the PSF. The sum of the PSF should be equal to 1
    total_flux_psf = np.nansum(psf)
    if display: print(f'The total flux of the PSF image is {total_flux_psf}.')
    psf = psf/total_flux_psf

    if WRITETO: fits.writeto(os.path.join(inputs_resultdir,f'PSF{ID}.fits'), psf, overwrite=True)
    
    return psf


def load_PSF_one_or_several_datasets(dico = None, WRITETO=True):
    '''
    Returns the normalized PSF of all the observations to be matched in the MCMC/Nelder-Mead simulations.
    The PSF is normalized by the sum of the pixel values, so np.nansum(PSF) = 1.
    The PSF will be used to convolve the disk model in the MCMC/Nelder-Mead simulations.

    The center of the PSF should be on a given pixel located at x = y = n // 2 with the counting starting at 0 (so NOT between four pixels).

    Calls the function load_PSF_one_dataset().

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'nb_obs', 'datadir_all', 'fn_psf_all', 'two_psf_files_all', 'fn_psf_1_all', 'fn_psf_2_all', 'crop_psf_all', 'spatial_shift_psf_data_all', 'spectral_axis_all', 'channels_all'. 
        
        (optional)

        .'WRITETO' (boolean): if True (or 1), save the PSF used in the simulations at the path 'inputs_resultdir' when calling the function load_PSF_one_dataset()
    
    Outputs:
        .'psfs' (list of 2D-array): list of the normalized PSFs. There is one PSF per observation.

    '''
    # Access extra variables from dico
    NB_OBS       = dico.get('nb_obs', 'PROBLEM')
    DATADIR_ALL  = dico.get('datadir_all', 'PROBLEM')
    FN_PSF_ALL   = dico.get('fn_psf_all', 'PROBLEM')
    TWO_PSF_FILES_ALL = dico.get('two_psf_files_all', 'PROBLEM') # indicate whether there are two different files to consider for the PSF (1 = yes, 0 = no)
    FN_PSF_1_ALL = dico.get('fn_psf_1_all', 'PROBLEM')
    FN_PSF_2_ALL = dico.get('fn_psf_2_all', 'PROBLEM')
    CROP_PSF_ALL = dico.get('crop_psf_all', 'PROBLEM')
    SPATIAL_SHIFT_PSF_DATA_ALL = dico.get('spatial_shift_psf_data_all','PROBLEM')
    SPECTRAL_AXIS_ALL = dico.get('spectral_axis_all', 'PROBLEM')
    CHANNELS_ALL      = dico.get('channels_all', 'PROBLEM')

    if NB_OBS == 1: # only one dataset
        dico['datadir']  = DATADIR_ALL[0]
        dico['fn_psf']   = FN_PSF_ALL[0]
        dico['fn_psf_1'] = FN_PSF_1_ALL[0]
        dico['fn_psf_2'] = FN_PSF_2_ALL[0]
        dico['crop_psf'] = CROP_PSF_ALL[0]
        dico['spectral_axis']     = SPECTRAL_AXIS_ALL[0]
        dico['channels']          = CHANNELS_ALL[0]
        dico['two_psf_files']     = TWO_PSF_FILES_ALL[0]
        dico['spatial_shift_psf_data'] = SPATIAL_SHIFT_PSF_DATA_ALL[0]
        psfs = [load_PSF_one_dataset(dico=dico, ID='',WRITETO=WRITETO)]
        
    else:
        psfs = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['datadir']  = DATADIR_ALL[i]
            dico['fn_psf']   = FN_PSF_ALL[i]
            dico['fn_psf_1'] = FN_PSF_1_ALL[i]
            dico['fn_psf_2'] = FN_PSF_2_ALL[i]
            dico['crop_psf'] = CROP_PSF_ALL[i]
            dico['spectral_axis']     = SPECTRAL_AXIS_ALL[i]
            dico['channels']          = CHANNELS_ALL[i]
            dico['two_psf_files']     = TWO_PSF_FILES_ALL[i]
            dico['spatial_shift_psf_data'] = SPATIAL_SHIFT_PSF_DATA_ALL[i]
    
            psfs.append( load_PSF_one_dataset(dico=dico, ID=f'_{i}', WRITETO=WRITETO) )

    return np.array(psfs)


def load_SCIENCE_one_dataset(dico = None, WRITETO=True, ID=''):
    '''
    Returns the science data for one observation.

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'datadir', 'fn_science', 'crop_science', 'spectral_axis', 'channels', 'inputs_resultdir'. 

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the SCIENCE DATA used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the SCIENCE_DATA file.

    Outputs:
        .'science_data' (2D or 3D array): 
            .. 2D array if the SCIENCE DATA are already post-processed (e.g., polarized intensity data; or total intensity which have already been processed using the RDI technique) 
            .. 3D array if the SCIENCE DATA correspond to pre-processed total intensity data and will be post-processed using the ADI technique while running the MCMC/Nelder-Mead simulations and after having been subtracted of the disk model (to avoid biases in the disk optimization due to self-subtraction effects, Milli et al. 2012).
    '''

    print('\n(Load the SCIENCE_DATA)')

    # Access extra variables from dico
    display       = dico.get('display', 0)
    DATADIR       = dico.get('datadir', 'PROBLEM')
    FN_SCIENCE    = dico.get('fn_science', 'PROBLEM')
    CROP_SCIENCE  = dico.get('crop_science', 'PROBLEM')
    SPECTRAL_AXIS = dico.get('spectral_axis', None)
    CHANNELS      = dico.get('channels', None)
    CHANNELS_COMBINED = bool(np.where(len(CHANNELS) > 1, True, False))
    inputs_resultdir  = dico.get('inputs_resultdir', 'PROBLEM')

    if WRITETO :
        if inputs_resultdir == 'PROBLEM': raise ValueError('There is a problem with the value of inputs_resultdir.')
    
    path_science = os.path.join(DATADIR, FN_SCIENCE)
    if display: print(f'The path for the SCIENCE data:\n {path_science}\n')
    science_data = fits.getdata(path_science)

    if SPECTRAL_AXIS: # there is a spectral axis
        if display: 
            print('Consider the spectral channels of interest:')
            if len(np.shape(science_data)) == 4: print(f' [INFO] The axis of SCIENCE_DATA are wavelength, temporal, y, x.')
            elif len(np.shape(science_data)) == 3: print(f' [INFO] The axis of SCIENCE_DATA are wavelength, y, x.')
            print(f' [INFO] np.shape(science_data) = {np.shape(science_data)} # before selecting the spectral channels')
        science_data = science_data[CHANNELS] # consider the relevant spectral channels
        if display:print(f' [INFO] np.shape(science_data) = {np.shape(science_data)} # after selecting the spectral channels')
        if CHANNELS_COMBINED: 
            science_data = np.nansum(science_data, axis=0)
            if display: print(f' [INFO] np.shape(science_data) = {np.shape(science_data)} # after combining the spectral channels')

    if WRITETO: fits.writeto(os.path.join(inputs_resultdir, f'science_cube{ID}.fits'), science_data, overwrite=True)
   
    # Crop the spatial dimensions
    if  CROP_SCIENCE != 0:
        if display: print(f'Crop the SCIENCE_DATA: \n [INFO] The size of the SCIENCE_DATA *before* cropping is {np.shape(science_data)}.')
        if len(np.shape(science_data)) == 2: science_data = science_data[CROP_SCIENCE:-CROP_SCIENCE, CROP_SCIENCE:-CROP_SCIENCE]
        elif len(np.shape(science_data)) == 3: science_data = science_data[:,CROP_SCIENCE:-CROP_SCIENCE, CROP_SCIENCE:-CROP_SCIENCE]
        else: raise ValueError(f'[ERROR] The dimension of the SCIENCE_DATA is weird, it should be 2 or 3 but is {len(np.shape(science_data))}, check why!')
        if display: print(f' [INFO] The size of the SCIENCE_DATA *after* cropping is {np.shape(science_data)}.')

    elif CROP_SCIENCE == 0 and display: print(f' [INFO] The size of the SCIENCE_DATA is {np.shape(science_data)}.')

    if WRITETO and len(np.shape(science_data)) == 2: fits.writeto(os.path.join(inputs_resultdir,f'science_data{ID}.fits'), science_data, overwrite=True) # save the science data only if it has a dimension equal to 2,  to avoid saving heavy files. This file is not necessarily, but can be useful to check if the script raises some errors

    return np.array(science_data)


def load_SCIENCE_one_or_several_datasets(dico = None, WRITETO=True):
    '''
    Returns the science data of all the observations to be matched in the MCMC/Nelder-Mead simulations.

    Calls the function load_SCIENCE_one_dataset().

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'nb_obs_all', 'datadir_all', 'fn_science_all', 'crop_science_all', 'spectral_axis_all', 'channels_all'.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the SCIENCE DATA used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the SCIENCE_DATA file.

    Outputs:
        .'science_data' (list of 2D and/or 3D array): list of the SCIENCE_DATA. There is one SCIENCE_DATA per observation.
            .. 2D array if the SCIENCE DATA are already post-processed (e.g., polarized intensity data; or total intensity which have already been processed using the RDI technique) 
            .. 3D array if the SCIENCE DATA correspond to pre-processed total intensity data and will be post-processed using the ADI technique while running the MCMC/Nelder-Mead simulations and after having been subtracted of the disk model (to avoid biases in the disk optimization due to self-subtraction effects, Milli et al. 2012).
    '''

    # Access extra variables from dico
    NB_OBS    = dico.get('nb_obs', 'PROBLEM')
    DATADIR_ALL       = dico.get('datadir_all', 'PROBLEM')
    FN_SCIENCE_ALL    = dico.get('fn_science_all', 'PROBLEM')
    CROP_SCIENCE_ALL  = dico.get('crop_science_all', 'PROBLEM')
    SPECTRAL_AXIS_ALL = dico.get('spectral_axis_all', None)
    CHANNELS_ALL      = dico.get('channels_all', None)


    if NB_OBS == 1: # only one dataset
        dico['datadir'] = DATADIR_ALL[0]
        dico['fn_science'] = FN_SCIENCE_ALL[0]
        dico['crop_science'] = CROP_SCIENCE_ALL[0]
        dico['spectral_axis'] = SPECTRAL_AXIS_ALL[0]
        dico['channels'] = CHANNELS_ALL[0]

        science_data = [load_SCIENCE_one_dataset(dico=dico, ID='', WRITETO=WRITETO)]
        
    else:
        science_data = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['datadir'] = DATADIR_ALL[i]
            dico['fn_science'] = FN_SCIENCE_ALL[i]
            dico['crop_science'] = CROP_SCIENCE_ALL[i]
            dico['spectral_axis'] = SPECTRAL_AXIS_ALL[i]
            dico['channels'] = CHANNELS_ALL[i]

            science_data.append( load_SCIENCE_one_dataset(dico=dico, ID=f'_{i}', WRITETO=WRITETO) )

    return science_data


def load_REF_CUBE_one_dataset(dico = None, WRITETO=True, ID=''):
    '''
    Returns the reference cubes associated to the SCIENCE DATA for one observation. (If no reference cube, returns None.)

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'datadir', 'fn_ref_cube', 'crop_ref_cube', 'spectral_axis', 'channels'.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the REF CUBE used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the REF CUBE file.

    Outputs:
        .'ref_cube' (3D array or None):
    '''

    print('\n(Load the REF CUBE)')

    # Access extra variables from dico
    display       = dico.get('display', 0)
    DATADIR       = dico.get('datadir', 'PROBLEM')
    FN_REF        = dico.get('fn_ref', 'PROBLEM')
    CROP_REF      = dico.get('crop_ref', 'PROBLEM')
    SPECTRAL_AXIS = dico.get('spectral_axis', None)
    CHANNELS      = dico.get('channels', None)
    CHANNELS_COMBINED = bool(np.where(len(CHANNELS) > 1, True, False))
    inputs_resultdir  = dico.get('inputs_resultdir', 'PROBLEM')

    if WRITETO :
        if inputs_resultdir == 'PROBLEM': raise ValueError('There is a problem with the value of inputs_resultdir.')
    
    # Case: there is a reference cube associated to the observation
    if FN_REF not in [None, 'None']:
        path_ref = os.path.join(DATADIR, FN_REF)
        if display: print(f'The path for the REF CUBE data:\n {path_ref}\n')
        ref_cube = fits.getdata(path_ref)

        if SPECTRAL_AXIS: # there is a spectral axis
            if display: 
                print('Consider the spectral channels of interest:')
                if len(np.shape(ref_cube)) == 4: print(f' [INFO] The axis of the REF CUBE are wavelength, temporal, y, x.')
                elif len(np.shape(ref_cube)) == 3: print(f' [INFO] The axis of REF CUBE are wavelength, y, x.')
                print(f' [INFO] np.shape(ref_cube) = {np.shape(ref_cube)} # before selecting the spectral channels')
            ref_cube = ref_cube[CHANNELS] # consider the relevant spectral channels
            if display:print(f' [INFO] np.shape(ref_cube) = {np.shape(ref_cube)} # after selecting the spectral channels')
            if CHANNELS_COMBINED: 
                ref_cube = np.nansum(ref_cube, axis=0)
                if display: print(f' [INFO] np.shape(ref_cube) = {np.shape(ref_cube)} # after combining the spectral channels')

        if WRITETO: fits.writeto(os.path.join(inputs_resultdir, f'ref_cube{ID}.fits'), ref_cube, overwrite=True)
    
        # Crop the spatial dimensions
        if CROP_REF != 0:
            if display: print(f'Crop the the REF CUBE: \n [INFO] The size of the the REF CUBE *before* cropping is {np.shape(ref_cube)}.')
            if len(np.shape(ref_cube)) == 2: ref_cube = ref_cube[CROP_REF:-CROP_REF, CROP_REF:-CROP_REF]
            elif len(np.shape(ref_cube)) == 3: ref_cube = ref_cube[:,CROP_REF:-CROP_REF, CROP_REF:-CROP_REF]
            else: raise ValueError(f'[ERROR] The dimension of the the REF CUBE is weird, it should be 2 or 3 but is {len(np.shape(ref_cube))}, check why!')
            if display: print(f' [INFO] The size of the REF CUBE *after* cropping is {np.shape(ref_cube)}.')

        elif CROP_REF == 0 and display: print(f' [INFO] The size of the the REF CUBE is {np.shape(ref_cube)}.')

        if WRITETO: fits.writeto(os.path.join(inputs_resultdir,f'ref_cube{ID}.fits'), ref_cube, overwrite=True) # potentially heavy file

        return np.array(ref_cube)

    # Case: there is *no* reference cube associated to the observation
    else: return None


def load_REF_CUBE_one_or_several_datasets(dico = None, WRITETO=True):
    '''
    Returns the reference cubes associated to the SCIENCE DATA of all the observations. (If no reference cube, returns [..., None, ...])

    Calls the function load_REF_CUBE_one_dataset().

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'nb_obs_all', 'datadir_all', 'fn_ref_cube_all', 'crop_ref_cube_all', 'spectral_axis_all', 'channels_all'.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the REF CUBE used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the REF CUBE file.

    Outputs:
        .'ref_cube' (list of 3D array and/or None values): list of the REF CUBES. There is one 3D array or None value per observation.
    '''

    # Access extra variables from dico
    NB_OBS             = dico.get('nb_obs', 'PROBLEM')
    DATADIR_ALL        = dico.get('datadir_all', 'PROBLEM')
    FN_REF_ALL         = dico.get('fn_ref_all', 'PROBLEM')
    CROP_REF_ALL       = dico.get('crop_ref_all', 'PROBLEM')
    SPECTRAL_AXIS_ALL  = dico.get('spectral_axis_all', None)
    CHANNELS_ALL       = dico.get('channels_all', None)


    if NB_OBS == 1: # only one dataset
        dico['fn_ref']        = FN_REF_ALL[0]
        dico['crop_ref']      = CROP_REF_ALL[0]
        dico['spectral_axis'] = SPECTRAL_AXIS_ALL[0]
        dico['channels']      = CHANNELS_ALL[0]

        ref_cubes = [load_REF_CUBE_one_dataset(dico=dico, ID='', WRITETO=WRITETO)]
        
    else:
        ref_cubes = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['datadir']       = DATADIR_ALL[i]
            dico['fn_ref']        = FN_REF_ALL[i]
            dico['crop_ref']      = CROP_REF_ALL[i]
            dico['spectral_axis'] = SPECTRAL_AXIS_ALL[i]
            dico['channels']      = CHANNELS_ALL[i]

            ref_cubes.append( load_REF_CUBE_one_dataset(dico=dico, ID=f'_{i}', WRITETO=WRITETO) )

    return ref_cubes


def load_NOISE_one_dataset(dico = None, WRITETO=True, ID=''):
    '''
    Return the noise map for one observation. In practice, the noise map is either simply loaded similarly to the science data or computed.

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'datadir', 'fn_noise', 'crop_noise', 'noise_multiplication_factor', 'compute_noise_map', 'spatial_shift_noise_data', 'spectral_axis', 'channels', 'iwa', 'nb_modes', 'science_data', 'pa_array',  'inputs_resultdir'. Note: Depending on the observation, some of these parameters may not be used.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the NOISE MAP file.

    Outputs:
        .'noise_map' (2D array)
    '''
    print('\n(Load the NOISE_MAP)')
    # Access extra variables from dico
    display       = dico.get('display', 0)
    DATADIR       = dico.get('datadir', 'PROBLEM')
    FN_NOISE      = dico.get('fn_noise', 'PROBLEM')
    CROP_NOISE    = dico.get('crop_noise', 'PROBLEM')
    NOISE_MULTIPLICATION_FACTOR = dico.get('noise_multiplication_factor', 'PROBLEM')
    COMPUTE_NOISE_MAP           = dico.get('compute_noise_map', 'PROBLEM')
    SPATIAL_SHIFT_NOISE_DATA    = dico.get('spatial_shift_noise_data','PROBLEM')
    SPECTRAL_AXIS = dico.get('spectral_axis', None)
    CHANNELS      = dico.get('channels', None)
    CHANNELS_COMBINED = bool(np.where(len(CHANNELS) > 1, True, False)) 
    IWA           = dico.get('iwa', 7)
    NB_MODES      = dico.get('nb_modes', 10)
    SCIENCE_DATA  = dico.get('science_data','PROBLEM')
    PA_ARRAY      = dico.get('pa_array','PROBLEM')
    inputs_resultdir  = dico.get('inputs_resultdir', 'PROBLEM')

    # Load the noise map or...
    if not COMPUTE_NOISE_MAP: 
        path_noise = os.path.join(DATADIR, FN_NOISE)
        if display: print(f'The path for the NOISE data:\n {path_noise}\n')
        noise_map = fits.getdata(path_noise)

        if SPECTRAL_AXIS: # there is a spectral axis. Assumption: the noise map (so maybe not a 2D map yet) have the same dimension that the SCIENCE data indicated in the configuration file (and so before going through the functions load_SCIENCE_several_datasets, load_SCIENCE_one_dataset )
            if display: 
                print('Consider the spectral channels of interest:')
                if len(np.shape(noise_map)) == 4: print(f' [INFO] The axis of NOISE_MAP are wavelength, temporal, y, x.')
                elif len(np.shape(noise_map)) == 3: print(f' [INFO] The axis of NOISE_MAP are wavelength, y, x.')
                print(f' [INFO] np.shape(noise_map) = {np.shape(noise_map)} # before selecting the spectral channels')
            noise_map = noise_map[CHANNELS] # consider the relevant spectral channels
            if display: print(f' [INFO] np.shape(noise_map) = {np.shape(noise_map)} # after selecting the spectral channels')
            if CHANNELS_COMBINED: 
                noise_map = np.nansum(noise_map, axis=0)
                if display: print(f' [INFO] np.shape(noise_map) = {np.shape(noise_map)} # after combining the spectral channels')

    # ... compute it 
    else:
        print(f'The NOISE_MAP is computed by processing the SCIENCE_DATA with PCA using the minus parallactic angles and {NB_MODES} nmodes.')
        noise_almost = vip.psfsub.pca_fullfr.pca(SCIENCE_DATA, -PA_ARRAY, ncomp=NB_MODES, mask_center_px=IWA, imlib='opencv', full_output=False, verbose=False)
        
        if WRITETO: fits.writeto(os.path.join(inputs_resultdir,f'reduced_image_opp_pa{ID}.fits'), noise_almost, overwrite=True)
        noise_map = compute_limdet_map_ann(noise_almost, dr=2, alpha=1, center='coordinates_are_given', xc=np.shape(noise_almost)[1]//2 + SPATIAL_SHIFT_NOISE_DATA, yc=np.shape(noise_almost)[0]//2+SPATIAL_SHIFT_NOISE_DATA)

    noise_map[noise_map==0] = np.nan # later in the script, the noise_map will be used to divide something, so prevent any potential division per 0

    if CROP_NOISE != 0:
        if display: print(f'Crop the noise map: \n [INFO] The size of the NOISE_MAP *before* cropping is {np.shape(noise_map)}.')
        noise_map = noise_map[CROP_NOISE:-CROP_NOISE, CROP_NOISE:-CROP_NOISE]
        if display: print(f' [INFO] The size of the NOISE_MAP *after* cropping is {np.shape(noise_map)}.')
    elif CROP_NOISE == 0 and display: print(f'The size of the NOISE_MAP is {np.shape(noise_map)}.')

    if NOISE_MULTIPLICATION_FACTOR != 1. : 
        if display: print(f'Multiply the NOISE_MAP by {NOISE_MULTIPLICATION_FACTOR}.') 
    if WRITETO: 
        print('Save the NOISE_MAP.')
        fits.writeto(os.path.join(inputs_resultdir,f'NOISE_MAP{ID}.fits'), noise_map, overwrite=True)
        if WRITETO and len(np.shape(SCIENCE_DATA)) == 2: fits.writeto(os.path.join(inputs_resultdir,f'science_data_snr{ID}.fits'), SCIENCE_DATA/noise_map, overwrite=True) # save the S/N map of the SCIENCE_DATA only if it has a dimension equal to 2, to avoid saving heavy files. This file is not necessarily, but can be useful to check if the script raises some errors
    return  np.array(noise_map)


def load_NOISE_one_or_several_datasets(dico = None, WRITETO=True):
    '''
    Returns the noise maps of all the observations to be matched in the MCMC/Nelder-Mead simulations.

    Calls the function load_NOISE_one_dataset().
    
    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'datadir_all', 'fn_noise_all', 'crop_noise_all', 'noise_multiplication_factor_all', 'compute_noise_map_all', 'spatial_shift_noise_data_all', 'spectral_axis_all', 'channels_all', 'iwa_all', 'nb_modes_all', 'science_data_all', 'pa_array_all'.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.

    Outputs:
        .'noise_maps' (list of 2D array): There is one NOISE MAP per observation.
    '''
    
    # Access extra variables from dico
    NB_OBS    = dico.get('nb_obs', 'PROBLEM')
    DATADIR_ALL       = dico.get('datadir_all', 'PROBLEM')
    SCIENCE_DATA_ALL  = dico.get('science_data_all', 'PROBLEM')
    PA_ARRAY_ALL      = dico.get('pa_array_all', 'PROBLEM')
    FN_NOISE_ALL      = dico.get('fn_noise_all', 'PROBLEM')
    CROP_NOISE_ALL    = dico.get('crop_noise_all', 'PROBLEM')
    SPECTRAL_AXIS_ALL = dico.get('spectral_axis_all', None)
    CHANNELS_ALL      = dico.get('channels_all', 'PROBLEM')
    SPATIAL_SHIFT_NOISE_DATA_ALL    = dico.get('spatial_shift_noise_data_all', 'PROBLEM')
    COMPUTE_NOISE_MAP_ALL           = dico.get('compute_noise_map_all', 'PROBLEM')
    NOISE_MULTIPLICATION_FACTOR_ALL = dico.get('noise_multiplication_factor_all', 'PROBLEM')

    IWA_ALL =  dico.get('iwa_all', 'PROBLEM')

    if NB_OBS == 1: # only one dataset
        dico['datadir']       = DATADIR_ALL[0]
        dico['fn_noise']      = FN_NOISE_ALL[0]
        dico['crop_noise']    = CROP_NOISE_ALL[0]
        dico['science_data']  = SCIENCE_DATA_ALL[0]
        dico['pa_array']      = PA_ARRAY_ALL[0]
        dico['iwa']           = IWA_ALL[0]
        dico['spectral_axis'] = SPECTRAL_AXIS_ALL[0]
        dico['channels']      = CHANNELS_ALL[0]
        dico['compute_noise_map'] = COMPUTE_NOISE_MAP_ALL[0]
        dico['noise_multiplication_factor'] = NOISE_MULTIPLICATION_FACTOR_ALL[0]
        dico['spatial_shift_noise_data']    = SPATIAL_SHIFT_NOISE_DATA_ALL[0]

        noise_maps = [load_NOISE_one_dataset(dico, ID='', WRITETO=WRITETO)]
        
    else:
        noise_maps = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['datadir']       = DATADIR_ALL[i]
            dico['fn_noise']      = FN_NOISE_ALL[i]
            dico['crop_noise']    = CROP_NOISE_ALL[i]
            dico['science_data']  = SCIENCE_DATA_ALL[i]
            dico['pa_array']      = PA_ARRAY_ALL[i]
            dico['iwa']           = IWA_ALL[i]
            dico['spectral_axis'] = SPECTRAL_AXIS_ALL[i]
            dico['channels']      = CHANNELS_ALL[i]
            dico['compute_noise_map'] = COMPUTE_NOISE_MAP_ALL[i]
            dico['noise_multiplication_factor'] = NOISE_MULTIPLICATION_FACTOR_ALL[i]
            dico['spatial_shift_noise_data']    = SPATIAL_SHIFT_NOISE_DATA_ALL[i]

            noise_maps.append( load_NOISE_one_dataset(dico, ID=f'_{i}', WRITETO=WRITETO ) )

    return noise_maps
       

def load_PA_ARRAY_one_dataset(dico = None, ID='', WRITETO=False):
    '''
    Returns the parallactic angles array for pupil-stabilized observations.
    Otherwise, return None.

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'datadir', 'fn_pa', 'run_postprocessing_technique', 'inputs_resultdir'. Note: Depending on the observation, some of these parameters may not be used.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the NOISE MAP file.

    Outputs:
        .'pa_array' (1D array) or None

    '''
    print('\n(Load the PA -PARALLACTIC ANGLES- array)')
    # Access extra variables from dico
    display  = dico.get('display', 0)
    RUN_POSTPROCESSING_TECHNIQUE = dico.get('run_postprocessing_technique', 'PROBLEM')
    DATADIR  = dico.get('datadir', 'PROBLEM')
    FN_PA    = dico.get('fn_pa', 'PROBLEM')
    inputs_resultdir  = dico.get('inputs_resultdir', 'PROBLEM')
           
    if 'ADI' in RUN_POSTPROCESSING_TECHNIQUE:
        print('yes')
        path_pa  = os.path.join(DATADIR, FN_PA)
        if display: print(f'The path for the parallactic angles array:\n {path_pa}\n')
        pa_array = -fits.getdata(path_pa)
        if display: print(f'The shape of the parallactic angles array is {np.shape(pa_array)}.')
        if WRITETO: fits.writeto(os.path.join(inputs_resultdir,f'parallactic_angles{ID}.fits'), pa_array, overwrite=True)
        return np.array(pa_array)
    else:
        return None


def load_PA_ARRAY_one_or_several_datasets(dico = None, WRITETO=False):
    '''
    Returns the parallactic angles array (or None) of all the observations to be matched in the MCMC/Nelder-Mead simulations.

    Calls the function load_PA_ARRAY_one_dataset()

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'nb_obs', 'datadir_all', 'fn_pa_all', 'run_postprocessing_technique_all'.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the NOISE MAP file.

    Outputs:
        .'pa_arrays' (list of 1D arrays and/or None)

    '''
    # Access extra variables from dico
    NB_OBS    = dico.get('nb_obs', 'PROBLEM')
    RUN_POSTPROCESSING_TECHNIQUE_ALL = dico.get('run_postprocessing_technique_all', 'PROBLEM')
    DATADIR_ALL  = dico.get('datadir_all', 'PROBLEM')
    FN_PA_ALL    = dico.get('fn_pa_all', 'PROBLEM')

    if NB_OBS == 1: # only one dataset
        dico['run_postprocessing_technique'] = RUN_POSTPROCESSING_TECHNIQUE_ALL[0]
        dico['datadir']  = DATADIR_ALL[0]
        dico['fn_pa']    = FN_PA_ALL[0]
        pa_arrays = [load_PA_ARRAY_one_dataset(dico=dico, WRITETO=WRITETO, ID=f'')]
        
    else:
        pa_arrays = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['run_postprocessing_technique'] = RUN_POSTPROCESSING_TECHNIQUE_ALL[i]
            dico['datadir']  = DATADIR_ALL[i]
            dico['fn_pa']    = FN_PA_ALL[i]
            pa_arrays.append( load_PA_ARRAY_one_dataset(dico=dico, WRITETO=WRITETO, ID=f'_{i}') )
        
    return pa_arrays


def load_MASK2MINIMIZE_one_dataset(dico = None, WRITETO = True, ID=''):
    '''
    Returns the mask2minimize region where the disk model should match the observation.

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'path_mask', 'crop_mask', 'inputs_resultdir'. 

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the NOISE MAP file.

    Outputs:
        .'mask2minimize' (2D array)
    '''

    print('\n(Load the MASK2MINIMIZE map) \n Remark: The MASK2MINIMIZE map define the region where the disk model will be optimized to match the SCIENCE_DATA, after SCIENCE_DATA have been processed.')
    # Access extra variables from dico
    display   = dico.get('display', 0)
    PATH_MASK = dico.get('path_mask', 'PROBLEM') # the full path is already given because the mask used could be located somewhere than the science/psf(/noise/pa_angles) data
    CROP_MASK = dico.get('crop_mask', 'PROBLEM')
    inputs_resultdir = dico.get('inputs_resultdir', 'PROBLEM')
    
    if display: print(f'The path for the MASK is\n {PATH_MASK}\n')

    if PATH_MASK[-4:] == 'fits': mask2minimize = fits.getdata(PATH_MASK)

    else: raise ValueError(f' [ERROR] The extension of the MASK2MINIMIZE file should be .fits and here the path to the file is \n{PATH_MASK}\n')

    if CROP_MASK != 0: # Crop the mask
        if display: print(f'Crop spatial dimensions of MASK2MINIMIZE: \n [INFO] The shape of the MASK2MINIMIZE data *before* cropping is {np.shape(mask2minimize)}.')
        mask2minimize = mask2minimize[CROP_MASK:-CROP_MASK, CROP_MASK:-CROP_MASK]
        if display: print(f' [INFO] The shape of the MASK2MINIMIZE data *after* cropping is {np.shape(mask2minimize)}.')

    elif CROP_MASK == 0 and display: print(f'The shape of the MASK2MINIMIZE is {np.shape(mask2minimize)}.')
    
    if WRITETO: fits.writeto(os.path.join(inputs_resultdir,f'MASK2MINIMIZE{ID}.fits'), mask2minimize, overwrite=True)
    return mask2minimize


def load_MASK2MINIMIZE_one_or_several_datasets(dico = None, WRITETO = True):
    '''
    Returns the mask2minimize region where the disk model should match the observation, and this for all the observations to be matched in the MCMC/Nelder-Mead simulations.

    Calls the function load_MASK2MINIMIZE_one_dataset()

    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'nb_obs', 'path_mask_all', 'crop_mask_all'. 

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.

    Outputs:
        .'masks2minimize' (list of 2D array)
    '''

    print('\n(Load MASK data)')
    # Access extra variables from dico
    NB_OBS = dico.get('nb_obs', 'PROBLEM')
    PATH_MASK_ALL = dico.get('path_mask_all', 'PROBLEM')
    CROP_MASK_ALL = dico.get('crop_mask_all', 'PROBLEM')

    if NB_OBS == 1: # only one dataset
        dico['path_mask'] = PATH_MASK_ALL[0]
        dico['crop_mask'] = CROP_MASK_ALL[0]

        masks2minimize = [load_MASK2MINIMIZE_one_dataset(dico, ID='', WRITETO=WRITETO)]
        
    else:
        masks2minimize = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['path_mask'] = PATH_MASK_ALL[i]
            dico['crop_mask'] = CROP_MASK_ALL[i]
            masks2minimize.append( load_MASK2MINIMIZE_one_dataset(dico, ID=f'_{i}', WRITETO=WRITETO ) )

    return masks2minimize

##############################################
############## POSTPROCESS DATA ##############
##############################################

def postprocess_SCIENCE_DATA_PCA_one_dataset(dico = None, WRITETO=False, ID=''):
    '''
    Returns the post-processed data using the PCA-ADI or PCA-ADI+RDI algorithm, or if the SCIENCE DATA is already post-processed, returns None.
   
    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'run_postprocessing_technique', 'iwa', 'nb_modes', 'science_data', 'cube_ref', 'mask2minimize', 'pa_array', 'inputs_resultdir'. Note: Depending on the observation, some of these parameters may not be used.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the NOISE MAP file.

    Outputs:
        .'red_data' (2D array) or None
    '''
    # Access extra variables from dico
    display  = dico.get('display', 0)
    IWA      = dico.get('iwa', 7)
    NB_MODES = dico.get('nb_modes', 10)
    SCIENCE_DATA  = dico.get('science_data','PROBLEM')
    REF_CUBE      = dico.get('ref_cube','PROBLEM')
    MASK2MINIMIZE = dico.get('mask2minimize', 'PROBLEM')
    PA_ARRAY      = dico.get('pa_array', 'PROBLEM')
    RUN_POSTPROCESSING_TECHNIQUE = dico.get('run_postprocessing_technique', 'PROBLEM')
    inputs_resultdir = dico.get('inputs_resultdir', 'PROBLEM')

    if REF_CUBE == None: REF_CUBE = None # sometimes REF_CUBE = np.array(None) instead of None and this messes up the function vip.psfsub.pca_fullfr.pca(). So just to be sure, sets REF_CUBE = None. Note: None == np.array(None) returns True.

    if 'ADI' in RUN_POSTPROCESSING_TECHNIQUE:
        print('\n(Post-process the SCIENCE data)')
        print(f'We post-process the data by applying {RUN_POSTPROCESSING_TECHNIQUE}.')
        red_data = vip.psfsub.pca_fullfr.pca(np.copy(SCIENCE_DATA), PA_ARRAY, ncomp=NB_MODES, cube_ref=REF_CUBE, mask_center_px=IWA, imlib='opencv', full_output=False, verbose=False)

        if WRITETO:
            fits.writeto(os.path.join(inputs_resultdir,f'reduced_image{ID}.fits'), red_data, overwrite=True)
            fits.writeto(os.path.join(inputs_resultdir,f'reduced_image_mask{ID}.fits'), red_data*MASK2MINIMIZE, overwrite=True)
        return red_data
    return None


def postprocess_SCIENCE_DATA_PCA_several_datasets(dico = None, WRITETO=False):
    '''
    Returns the post-processed data using the PCA-ADI or PCA-ADI+RDI algorithm, or if the SCIENCE DATA is already post-processed, returns None.

    Calls the function postprocess_SCIENCE_DATA_PCA_one_dataset()
   
    Inputs:
        .'dico' (dictionary): dictionary storing relevant parameters to run this function and that should be defined!. 'dico' should store the following parameters: 'nb_obs', 'run_postprocessing_technique_all', 'iwa_all', 'nb_modes_all', 'science_data_all', 'cube_ref_all', 'mask2minimize_all', 'pa_array_all'. Note: Depending on the observation, some of these parameters may not be used.

        (optional)

        .'WRITETO' (boolean): if True (or 1), save the NOISE MAP used in the simulations at the path 'inputs_resultdir'.
        .'ID' (string): suffix to be added in the filename when saving the NOISE MAP file.

    Outputs:
        .'red_data' (list of 2D array and/or None)
    '''
    # Access extra variables from dico
    NB_OBS    = dico.get('nb_obs', 'PROBLEM')
    IWA_ALL      = dico.get('iwa_all', 7)
    NB_MODES_ALL = dico.get('nb_modes_all', 10)
    SCIENCE_DATA_ALL  = dico.get('science_data_all','PROBLEM')
    REF_CUBE_ALL      = dico.get('ref_cube_all', 'PROBLEM')
    MASK2MINIMIZE_ALL = dico.get('mask2minimize_all', 'PROBLEM')
    PA_ARRAY_ALL      = dico.get('pa_array_all', 'PROBLEM')
    RUN_POSTPROCESSING_TECHNIQUE_ALL = dico.get('run_postprocessing_technique_all', 'PROBLEM')

    if NB_OBS == 1: # only one dataset
        dico['science_data']  = SCIENCE_DATA_ALL[0]
        dico['ref_cube']      = REF_CUBE_ALL[0]
        dico['mask2minimize'] = MASK2MINIMIZE_ALL[0]
        dico['pa_array'] = PA_ARRAY_ALL[0]
        dico['nb_modes'] = NB_MODES_ALL[0]
        dico['iwa']      = IWA_ALL[0]
        dico['run_postprocessing_technique'] = RUN_POSTPROCESSING_TECHNIQUE_ALL[0]
        red_data = [postprocess_SCIENCE_DATA_PCA_one_dataset(dico=dico, ID='', WRITETO=WRITETO)]
        
    else:
        red_data = []
        for i in range(NB_OBS): # loop other the different datasets
            dico['science_data']  = np.array(SCIENCE_DATA_ALL[i])
            dico['ref_cube']      = np.array(REF_CUBE_ALL[i])
            dico['mask2minimize'] = np.array(MASK2MINIMIZE_ALL[i])
            dico['pa_array'] = PA_ARRAY_ALL[i]
            dico['nb_modes'] = NB_MODES_ALL[i]
            dico['iwa']      = IWA_ALL[i]
            dico['run_postprocessing_technique'] = RUN_POSTPROCESSING_TECHNIQUE_ALL[i]
            red_data.append( postprocess_SCIENCE_DATA_PCA_one_dataset(dico=dico, ID=f'_{i}', WRITETO=WRITETO) ) 

    return red_data