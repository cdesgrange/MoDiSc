#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Script with functions to fit via MCMC debris disks.

Functions:
    initialize_walkers_backend()
    from_param_to_theta_init() 
    from_param_to_theta_init_5params_1obs()
    from_param_to_theta_init_6params_1obs()
    from_param_to_theta_init_6params_2obs()
    from_param_to_theta_init_8params_2obs()
    logp()
    logl()
    logl_5params_1obs()
    logl_6params_1obs()
    logl_6params_2obs()
    logl_8params_2obs()
    lnpb_mcmc() # used when running the MCMC
    lnpb()
    generate_disk_model_5params()
    chisquare_params()
    chisquare_5params_1obs()
    chisquare_6params_1obs()
    chisquare_6params_2obs_polar()
    chisquare_8params_2obs_polar()


This script would be imported by the scripts:
    diskfit_*params_*obs.py
    plotfrombackend_*.py
    test_disks.py
'''

#from import_functions_generic import *
from ..packages import *
from ..functions.noise_map_derive import *


def print_params_names_and_values(params_names, params_values, params_flat=True, skip_first_line=False):
    '''
    Print for each parameter its name and its value.
    '''

    if not skip_first_line: print(f'\nNames and values of the disk parameters:')
    
    if params_flat:
        for j_name, name in enumerate( params_names ):
            print(f' {name} (initial value = {np.array(params_values[j_name])})')

    else:
        nb_disk_structures = len(params_names)
        for i_struct in range(nb_disk_structures):
            print(f'For disk structure #{i_struct+1}:')

            #print(params_values[i_struct])
            for j_name, name in enumerate( params_names[i_struct] ):
                #print(i_struct, j_name)
                print(f' {name} (initial value = {np.array(params_values[i_struct][j_name])})')
    print('')
    return


def print_params_names_and_values_and_bounds(params_names, params_values, params_bounds, params_flat=True, skip_first_line=False):
    '''
    Print for each parameter its name and its value.
    '''
    if not skip_first_line: print(f'\nNames, values and bounds of the disk parameters:')
    #print(params_names)
    #print(params_values)
    #print(params_bounds)

    if params_flat :
        for j_name, name in enumerate( params_names ):
            print(f' {name} (initial value = {np.array(params_values[j_name])}, bounds = {np.array(params_bounds[j_name])})')
        print('')

    else:
        nb_disk_structures = len(params_names)
        for i_struct in range(nb_disk_structures):
            print(f'\n- Disk structure #{i_struct+1}:')
            for j_name, name in enumerate( params_names[i_struct] ):
                print(f' {name} (initial value = {np.array(params_values[i_struct][j_name])}, bounds = {np.array(params_bounds[i_struct][j_name])})')
    print('')
    return



########################################################
############# INITIALIZE DISK PARAMETERS ###############
########################################################

def from_param_list_of_dict_to_param_list_of_list(PARAMS_NAMES_LIST_OF_DICT):
    """ create a initial set of parameters from the initial parameters
        store in the init yaml file
    Args:
        config_file: dic, all the parameters of the exploration algo and PCA read from yaml file

    Returns:
        initial set of parameters
    """
    nb_disk_structures = np.shape(PARAMS_NAMES_LIST_OF_DICT)[0]
   
    params_list_of_list = []
    for i_struct in range( nb_disk_structures ):
        params_names_for_struct_i = PARAMS_NAMES_LIST_OF_DICT[i_struct][f'disk_structure{i_struct+1}']
        params_list_of_list.append(params_names_for_struct_i)
    
    return params_list_of_list


def from_params_names_to_param_bounds(config_file, shape_output='unflat'):
    """ create a initial set of parameters from the initial parameters
        store in the init yaml file
    Args:
        config_file: dic, all the parameters of the exploration algo and PCA read from yaml file

    Returns:
        initial set of parameters
    """
    PARAMS_NAMES = config_file['PARAMS_NAMES']
    nb_disk_structures = len(PARAMS_NAMES)
    
    params_bounds = []
    # Dimension 
    if shape_output == 'unflat':
        for i_struct in range( nb_disk_structures ):  # loop on each disk structure
            params_bounds_for_struct_i = []

            for param in PARAMS_NAMES[i_struct]:  # for a given disk structure, loop on each parameter
                params_bounds_for_struct_i.append(config_file[f'{param}_BOUNDS_STRUCT{i_struct+1}'])

                for k_obs in range(len(params_bounds_for_struct_i[-1])):  # for a given disk structure and parameter, note that the parameter can have different values depending on the observation. In principle this should only affect the parameter named SCALING. Thus params_bounds_for_struct_i[-1] is equal to a list of one two-element list, except for the parameter named SCALING, for which params_bounds_for_struct_i[-1] is equal to a list of several two-element lists, several being equal to the number of observations.
                    
                    # Check that minimal and maximal bounds are consistent
                    min, max = params_bounds_for_struct_i[-1][k_obs][0], params_bounds_for_struct_i[-1][k_obs][1]
                    if min > max: raise ValueError(f'The upper boundary for the parameter {param} is smaller than the lower boundary... {max} < {min} !') 
            
            params_bounds.append( params_bounds_for_struct_i) 

    elif shape_output == 'flat':
        for i_struct in range( nb_disk_structures ):  # loop on each disk structure
            
            for param in PARAMS_NAMES[i_struct]:  # for a given disk structure, loop on each parameter
                param_bounds_for_struct_i = config_file[f'{param}_BOUNDS_STRUCT{i_struct+1}']
               
                #print(param_bounds_for_struct_i)

                for k_obs in range(len(param_bounds_for_struct_i)):  # for a given disk structure and parameter, note that the parameter can have different values depending on the observation. In principle this should only affect the parameter named SCALING. Thus params_bounds_for_struct_i[-1] is equal to a list of one two-element list, except for the parameter named SCALING, for which params_bounds_for_struct_i[-1] is equal to a list of several two-element lists, several being equal to the number of observations.
                    params_bounds.append( param_bounds_for_struct_i[k_obs]) 
                    #print(param, params_bounds[-1])

                    # Check that minimal and maximal bounds are consistent
                    min, max = params_bounds[-1][0], params_bounds[-1][1]
                    if min > max: raise ValueError(f'The upper boundary for the parameter {param} is smaller than the lower boundary... {max} < {min} !') 

    return params_bounds

def from_params_names_to_param_init(config_file, shape_output='unflat'):
    """ create a initial set of parameters from the initial parameters
        store in the init yaml file
    Args:
        config_file: dic, all the parameters of the exploration algo and PCA
                            read from yaml file

    Returns:
        initial set of parameters
    """
    PARAMS_NAMES = config_file['PARAMS_NAMES']
    nb_disk_structures = len(PARAMS_NAMES)

    params_init = []
    if shape_output == 'unflat':
        for i_struct in range( nb_disk_structures ):
            params_init_for_struct_i = []
        
            for param in PARAMS_NAMES[i_struct]:
                params_init_for_struct_i.append(config_file[f'{param}_INIT_STRUCT{i_struct+1}'])
            params_init.append( params_init_for_struct_i) 

    else: # flat output
        for i_struct in range( nb_disk_structures ):
            for param in PARAMS_NAMES[i_struct]:
                param_values = config_file[f'{param}_INIT_STRUCT{i_struct+1}']
                for k_obs in range(len(param_values)):
                    params_init.append(config_file[f'{param}_INIT_STRUCT{i_struct+1}'][k_obs])
    return params_init


def from_params_names_unflat_to_params_names_flat(config_file, config_file_param='PARAMS_NAMES'):
    """ create a initial set of parameters from the initial parameters
        store in the init yaml file
    Args:
        config_file: dic, all the parameters of the exploration algo and PCA
                            read from yaml file

    Returns:
        initial set of parameters
    """
    PARAMS_NAMES = config_file['PARAMS_NAMES']
    PARAMS_xx    = config_file[config_file_param]
    nb_disk_structures = len(PARAMS_NAMES)
    NB_OBS       = config_file['NB_OBS']

    params_xx_flat = []
    for i_struct in range( nb_disk_structures ):
            for j_param, param in enumerate(PARAMS_NAMES[i_struct]):
                param_values = config_file[f'{param}_INIT_STRUCT{i_struct+1}']
                for k_obs in range(len(param_values)):
                    
                    if config_file_param == 'PARAMS_NAMES':
                        params_xx_flat.append(param)

                    else:
                        params_xx_flat.append(PARAMS_xx[i_struct][j_param][param])
    return params_xx_flat


def reshape_list_of_params_from_1D_into_3D(params_to_be_reshaped, dico):
    nb_disk_structures = dico.get("nb_disk_structures")
    params_example = dico.get("params_example_unflat")

    params_reshaped = copy.deepcopy(params_example)
    
    count=0
    for i in range(len(params_example)):
        for j in range(len(params_example[i])):
            for k in range(len(params_example[i][j])):
                params_reshaped[i][j][k] = params_to_be_reshaped[count]
                count +=1

    return params_reshaped



def correct_unit(params_names_flat, params_values, convention_unit='MCMC', output_unit='user-friendly', axis_to_be_considered=-1, display=1):
    # Change scaling and inclination values to physical
    if convention_unit == 'MCMC':
        if output_unit == 'MCMC': pass
        
        elif output_unit == 'user-friendly':
            for i, name in enumerate(params_names_flat):
                if name.lower() == 'scaling':
                    if display: print(f' [{name}] Change log10 values to normal values')
                    fct = lambda x : 10**x
                    
                    if len(np.shape(params_values)) == 1:   params_values[i] = fct(params_values[i])
                    elif len(np.shape(params_values)) == 2 and axis_to_be_considered == -1 : params_values[:,i] = fct(params_values[:,i])
                    elif len(np.shape(params_values)) == 3 and axis_to_be_considered == -1 : params_values[:,:,i] = fct(params_values[:,:,i])
                    else: raise ValueError(f'Check the input parameters. The function correct_unit() is up to now coded for len(np.shape(params_values)) = 1, 2, 3 and axis_to_be_considered = -1. Here, however, len(np.shape(params_values)) = {len(np.shape(params_values))} and axis_to_be_considered = {axis_to_be_considered}.')

                elif name.lower() in ['inc', 'inclination']:
                    if display: print(f' [{name}] Change arccos values to degrees values')
                    fct = lambda x : np.arccos(x)*180/np.pi 
                    
                    if len(np.shape(params_values)) == 1:   params_values[i] = fct(params_values[i])
                    elif len(np.shape(params_values)) == 2 and axis_to_be_considered == -1 : params_values[:,i] = fct(params_values[:,i])
                    elif len(np.shape(params_values)) == 3 and axis_to_be_considered == -1 : params_values[:,:,i] = fct(params_values[:,:,i])
                    else: raise ValueError(f'Check the input parameters. The function correct_unit() is up to now coded for len(np.shape(params_values)) = 1, 2, 3 and axis_to_be_considered = -1. Here, however, len(np.shape(params_values)) = {len(np.shape(params_values))} and axis_to_be_considered = {axis_to_be_considered}.')

        else: raise ValueError(f"Check the input parameter 'output_unit'. The function correct_unit() is up to now coded only for  output_unit='MCMC' or  output_unit='user-friendly'. Here, however, output_unit = {output_unit}.")

    elif convention_unit == 'user-friendly':
        if output_unit == 'MCMC':
            for i, name in enumerate(params_names_flat):
                
                if name.lower() == 'scaling':
                    if display: print(f' [{name}] Apply log10().')
                    fct = lambda x : np.log10(x)
                    
                    if len(np.shape(params_values)) == 1:   params_values[i] = fct(params_values[i])
                    elif len(np.shape(params_values)) == 2 and axis_to_be_considered == -1 : params_values[:,i] = fct(params_values[:,i])
                    elif len(np.shape(params_values)) == 3 and axis_to_be_considered == -1 : params_values[:,:,i] = fct(params_values[:,:,i])
                    else: raise ValueError(f'Check the input parameters. The function correct_unit() is up to now coded for len(np.shape(params_values)) = 1, 2, 3 and axis_to_be_considered = -1. Here, however, len(np.shape(params_values)) = {len(np.shape(params_values))} and axis_to_be_considered = {axis_to_be_considered}.')

                elif name.lower() in ['inc', 'inclination']:
                    if display: print(f' [{name}] Apply cos().')
                    fct = lambda x : np.cos(x*np.pi/180)
                    
                    if len(np.shape(params_values)) == 1:   params_values[i] = fct(params_values[i])
                    elif len(np.shape(params_values)) == 2 and axis_to_be_considered == -1 : params_values[:,i] = fct(params_values[:,i])
                    elif len(np.shape(params_values)) == 3 and axis_to_be_considered == -1 : params_values[:,:,i] = fct(params_values[:,:,i])
                    else: raise ValueError(f'Check the input parameters. The function correct_unit() is up to now coded for len(np.shape(params_values)) = 1, 2, 3 and axis_to_be_considered = -1. Here, however, len(np.shape(params_values)) = {len(np.shape(params_values))} and axis_to_be_considered = {axis_to_be_considered}.')

        elif output_unit == 'user-friendly': pass

        else: raise ValueError(f"Check the input parameter 'output_unit'. The function correct_unit() is up to now coded only for  output_unit='MCMC' or  output_unit='user-friendly'. Here, however, output_unit = {output_unit}.")

    else: print(f" [INFO] The input parameter 'convention_unit' is different from 'MCMC' and 'user-friendly'. Therefore, no unit is going to be corrected.")

    return np.array(params_values)

