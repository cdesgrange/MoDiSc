#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Script containing functions to run MCMC or Nelder-Mead simulations to fit disks imaged in scattered-light.

Functions:
    . initialize_walkers_backend() -  Initialize the walker

    . logp() -  Log of the priors
    . logl() - Log of the Likelihood
    . lnpb_mcmc() - Sum the logs (used when running MCMC simulations)
    . lnpb() - Sum the logs

    . chisquare_params_NelderMead() - Derive the chisquare value (used when running Nelder-Mead simulations). Calls the function chisquare_params_flat_or_unflat()
    . chisquare_params_flat_or_unflat() - Derive the chisquare value (adapted to different shape -flat or unflat- of the list of the parameters). Calls the function chisquare()
    . chisquare() - Derive the chisquare value 

    . generate_disk_model() - Generate the disk model using the GRaTeR code implemented in the VIP-HCI package


This script would be imported by the scripts:
    ..run_modisc.py
    ..plot_mcmc_results.py
    ..functions.chisquare_test.py
'''

__author__ = 'Celia Desgrange'

# Import classic packages
from ..packages import * 
# Import MoDiSc functions 
from ..functions.parameters_rearrange_print import *


########################################################
################ Initialize the walkers ################
########################################################

def initialize_walkers_backend(dico, config_file, SAVINGDIR):
    # if new_backend = 0, reset the backend, if not restart the chains.
    # Be careful if you change the parameters or walkers, 
    # you have to put new_backend = 1
    
    new_backend   = dico['mcmc_new_backend']
    nb_walkers    = dico['mcmc_nb_walkers']
    nb_free_params_tot = dico['nb_free_params_tot'] # = dimension of the MCMC
    fraction_ball = dico['mcmc_fraction_ball']
    display       = dico.get('display', 1)    

    if new_backend:
        mcmcresultdir = os.path.join(SAVINGDIR, 'results_MCMC') 
        os.makedirs(mcmcresultdir, exist_ok=True)
    else:
        resultdir = config_file['RESULTDIR']
        mcmcresultdir = os.path.join(resultdir, 'results_MCMC')

    params_names = from_params_names_unflat_to_params_names_flat(config_file)
    params_init  = from_params_names_to_param_init(config_file, shape_output='flat')

    #if display: print_params_names_and_values(params_names, params_init, params_flat=True)

    # Set up the backend
    # Don't forget to clear it in case the file already exists
    filename_backend = os.path.join(mcmcresultdir,  "backend_file_mcmc.h5")
    backend_init = backends.HDFBackend(filename_backend)

    #############################################################
    # Initialize the walkers. The best technique seems to be
    # to start in a small ball around the a priori preferred position.
    # Don't worry, the walkers quickly branch out and explore the
    # rest of the space.
    if new_backend == 1:
        init_balls = []
        print(f'There are {len(params_init):.0f} parameters in the variable params_init. They should correspond to the parameters to be fitted.')
        for i in range(len(params_init)):
            param = params_names[i]

            if param in ['log_f', 'log_f0', 'log_f1', 'log_f2', 'log_f3']:
                init_balli = np.random.uniform(-10, 10, size=(nb_walkers))
            else:
                init_balli = np.random.uniform(params_init[i] * (1-fraction_ball), params_init[i] * (1+fraction_ball), size=(nb_walkers))
            init_balls.append(init_balli)

        # be careful to the shape of the array, previously done via the function np.dstack
        # p0 = np.dstack((init_ball0, init_ball1, init_ball2, init_ball3, init_ball4, init_ball5))
        # now
        p0 = np.array([np.array(init_balls).transpose()])
        
        backend_init.reset(nb_walkers, nb_free_params_tot)
        return p0[0], backend_init

    return None, backend_init

########################################################
################## Log of the priors ###################
########################################################

def logp(params, dico=None, params_flat=True):
    params_bounds = dico.get('params_bounds', 'PROBLEM')

    #if len(params_bounds)[0] != len(params): raise ValueError(f'The number of parameters listed in the list "bounds" {len(params_bounds)} should be the same than the number of parameters in "params" ({len(params)}). \n bounds = {params_bounds} \n\n params = {params}')
    nb_disk_structures = len(params_bounds)
    
    #print('params_bounds', params_bounds)
    #print('params', params)

    # Check if a parameter is beyond
    if params_flat : # the list params is flat and the list params_bounds is almost flat (dimension of list params * 2)
        n = len(params)
        params_names  = dico.get('params_names_flat', 'PROBLEM')
        for i in range(n):
            if params_bounds[i][0] != None :
                if params[i] < params_bounds[i][0] : 
                    print(f' The tested parameter {params_names[i]} ({params[i]}) is out of bounds. It should be bigger than {params_bounds[i][0]}.')
                    return -np.inf
            if params_bounds[i][1] != None :
                if params[i] > params_bounds[i][1] :
                    print(f' The tested parameter {params_names[i]} ({params[i]} is out of bounds. It should be smaller than {params_bounds[i][1]}.')
                    return -np.inf
    else :
        params_names  = dico.get('params_names_unflat', 'PROBLEM')
        for i_struct in range( nb_disk_structures ): # loop on each disk structure
            #print('params_bounds[i_struct]', params_bounds[i_struct])

            #print('params[i_struct]', params[i_struct])

            #print('')
            
            for j_param in range(len(params[i_struct])): # for a given disk structure, loop on each parameter
                #print(i_struct, j_param)
                bounds_struct_i_param_j = params_bounds[i_struct][j_param]  # for a given disk structure and parameter, parameters can have different values depending on the observation. In principle this should only affect the parameter named SCALING. Thus bounds_istruct_jparam is equal to a list of one two-element list, except for the parameter named SCALING, for which bounds_istruct_jparam is equal to a list of several two-element lists, several being equal to the number of observations.
                #print('bounds_struct_i_param_j', bounds_struct_i_param_j, params[i_struct][j_param])
                
                ##print('bounds_istruct_jparam', bounds_struct_i_param_j)
                for k_obs in range(len(bounds_struct_i_param_j)): 
                    #print('\n',k_obs,'/',len(bounds_struct_i_param_j))
                    bounds = bounds_struct_i_param_j[k_obs]
                    #print('bounds', bounds)
                    #print('params[i_struct][j_param][k_obs]', params[i_struct][j_param][k_obs])
                    if bounds[0] != None :
                        #print('bounds', bounds)
                        #print(i_struct, j_param, k_obs)
                        #print(params)
                        #time.sleep(1)
                        #print('params[i_struct][j_param]', params[i_struct][j_param])
                        #print('params[i_struct][j_param][k_obs]',params[i_struct][j_param][k_obs])
                        
                        if params[i_struct][j_param][k_obs] < bounds[0] : return -np.inf

                    if bounds[1] != None :
                        if params[i_struct][j_param][k_obs] > bounds[1] : return -np.inf

    return 0


########################################################
################ Log of the Likelihood #################
########################################################

def logl(params, dico=None, params_flat=True):
    ll = -0.5 * chisquare_params_flat_or_unflat(params, dico=dico, params_flat=params_flat)
    return ll

########################################################
#################### Sum the logs ######################
########################################################

def lnpb_mcmc(params, **kwargs):
    kwargs['params_1D'] = params
   
    lp = logp(params, dico=kwargs, params_flat=True)
    
    if not np.isfinite(lp): 
        return -np.inf

    # Log of the likelihood
    ll = logl(params, dico=kwargs, params_flat=True)
    return lp + ll


def lnpb(params, dico=None, params_flat=True):

    # Log of the priors
    lp = logp(params, dico=dico, params_flat=params_flat)
    if not np.isfinite(lp): return -np.inf

    # Log of the likelihood
    ll = logl(params, dico=dico, params_flat=params_flat)
    return lp + ll

########################################################
############## Derive the chisquare value ##############
########################################################

def chisquare_params_NelderMead(params, kwargs):
    return chisquare_params_flat_or_unflat(params, dico=kwargs, params_flat='flat')


def chisquare_params_flat_or_unflat(params, dico = None, params_flat=True):
    if not params_flat : # the list of parameters is not flat
        params3D = params
        params1D = list(deepflatten(list(params3D)))

    else: # the list of parameters is flat
        params1D = params
        params3D = reshape_list_of_params_from_1D_into_3D(params1D, dico)

    dico['params_1D'] = params1D

    return chisquare(params3D, dico = dico)


def chisquare(params_3D, dico = None):
    # Access extra keyword arguments from dictionary
    display          = dico.get('display', 0)
    nb_obs           = dico.get('nb_obs', 1)
    type_obs_all     = dico.get('type_obs_all', 'PROBLEM')
    disk_model_polar_all = dico.get('disk_model_polar_all', 'PROBLEM')
    do_robust_convolution_all = dico.get('do_robust_convolution_all', 1)
    im_pa_all        = dico.get('im_pa_all', 'PROBLEM')
    psf_all          = dico.get('psf_all', 'PROBLEM')
    mask2minimize_all= dico.get('mask2minimize_all', 'PROBLEM')
    noise_map_all    = dico.get('noise_map_all', 'PROBLEM')
    science_data_all = dico.get('science_data_all', 'PROBLEM')
    ref_cube_all     = dico.get('ref_cube_all', 'PROBLEM')
    dimension_all    = dico.get('dimension_all', 'PROBLEM')
    nb_resel_all     = dico.get('nb_resel_all', 'PROBLEM')
    nb_free_params_tot        = dico.get('nb_free_params_tot', 'PROBLEM')
    nb_free_params_per_obs_all= dico.get('nb_free_params_per_obs_all', 'PROBLEM')
    save_some_results  = dico.get('save_some_results', False)
    save_full_results  = dico.get('save_full_results', False)
    simu_resultdir     = dico.get('simu_resultdir','None')
    pa_array_all       = dico.get('pa_array_all', 'PROBLEM')
    run_postprocessing_technique_all = dico.get('run_postprocessing_technique_all', 'PROBLEM')
    iwa_all            = dico.get('iwa_all', 'MASK_RAD')
    nb_modes_all       = dico.get('nb_modes_all', 'PROBLEM')
    plate_scale_all    = dico.get('plate_scale_all', 'PROBLEM')
    params_names       = dico.get('params_names_unflat', 'PROBLEM')
    nb_disk_structures = dico.get('nb_disk_structures', 'PROBLEM')  
    params_1D          = dico.get('params_1D','PROBLEM')  
    
    if display: print('------------------------------------------------------------------------------------')
    for i_obs in range(nb_obs): # loop on each observation
        if display: print(f'Consider the observation #{i_obs+1}')
        dico['idx_obs']   = int(i_obs)

        type_obs      = type_obs_all[i_obs]
        science_data  = science_data_all[i_obs]
        ref_cube      = ref_cube_all[i_obs]
        noise_map     = noise_map_all[i_obs]
        pa_array      = pa_array_all[i_obs]
        psf           = psf_all[i_obs]
        nb_modes      = nb_modes_all[i_obs]
        iwa           = iwa_all[i_obs]
        mask2minimize = mask2minimize_all[i_obs]
        nb_resel      = nb_resel_all[i_obs]
        nb_free_params_per_obs = nb_free_params_per_obs_all[i_obs]
        do_robust_convolution = do_robust_convolution_all[i_obs]
        im_pa = im_pa_all[i_obs] # relevant only for polarized intensity data and if do_robust_convolution is set to 1
        run_postprocessing_technique = run_postprocessing_technique_all[i_obs]

        dico['plate_scale']  = plate_scale_all[i_obs]
        dico['dimension']    = dimension_all[i_obs]
        dico['nb_resel']     = nb_resel_all[i_obs]
        dico['disk_model_polar'] = disk_model_polar_all[i_obs]

        # Case 1: the data are total intensity data
        if type_obs == 'total_intensity':    
            cube_diff = np.copy(science_data)  

            models = []

            for j_struct in range(nb_disk_structures): # loop on each disk structure
                if display: print(f'Consider the disk structure #{j_struct+1}')

                dico['disk_struc'] = int(j_struct+1)

                params_values = params_3D[j_struct]
                dico['params_names_1model'] = params_names[j_struct]
                
                # Generate the disk model
                models.append(generate_disk_model(params_values, dico = dico))

                models = np.array(models)

            full_model = np.nansum(models, axis=0)    
            #print(f'\n\n CHECK WHEN SEVERAL DISK STRUCTURES THAT THE SHAPE OF full_model IS A SHAPE OF AN IMAGE { np.shape(full_model)} \n\n')
            # Rotate the disk model for different angles and convolve it by the PSF

            # Case 1a: postprocess the data (note: science data must be a 3D array)
            if 'ADI' in run_postprocessing_technique:
                if len(np.shape(science_data)) != 3: raise ValueError(f'The SCIENCE DATA must be pre-processed with one temporal axis and two spatial axis, and so a 3D array. However, their shape is {np.shape(science_data)}.')

                full_modelconvolved = vip.fm.cube_inject_fakedisk(full_model, -np.array(pa_array), psf=psf)

                # Compute the residuals and the chisquare
                # Remove the disk from the observation
                cube_diff -= full_modelconvolved

                # Reduce the observation
                if ref_cube == None: ref_cube = None # sometimes ref_cube = np.array(None) instead of None and this messes up the function vip.psfsub.pca_fullfr.pca(). So just to be sure, sets ref_cube = None. Note: None == np.array(None) returns True.
                
                try:
                    res_full_map = vip.psfsub.pca_fullfr.pca(cube_diff, pa_array, ncomp=nb_modes, mask_center_px=iwa, cube_ref=ref_cube, imlib='opencv', full_output=False)
                except: 
                    os.makedirs('PROBLEM', exist_ok=True)
                    suffix_nb_obs_nb_struct = f"_obs{i_obs+1:.0f}_struc{j_struct+1:.0f}".replace('.','p')
                    fits.writeto('PROBLEM'+f'/disk_model_{suffix_nb_obs_nb_struct}_PROBLEM.fits',full_model, overwrite=True)
                    fits.writeto('PROBLEM'+f'/disk_model_convolved_{suffix_nb_obs_nb_struct}_PROBLEM.fits',full_modelconvolved, overwrite=True)
                    fits.writeto('PROBLEM'+f'/disk_CUBE_DIFF_{suffix_nb_obs_nb_struct}_PROBLEM.fits',cube_diff, overwrite=True)

                    raise ValueError(f'The PCA reduction did not converge for the observation {i_obs+1:.0f} when after subtracting the model of the disk structure #{j_struct:.0f}. The disk parameters used are the following\n{params_3D}')

            # Case 1b: don't postprocess the data (note: science data must be a 2D array)
            else:
                if len(np.shape(science_data)) != 2: raise ValueError(f'The SCIENCE DATA must be post-processed and so a 2D array. However, their shape is {np.shape(science_data)}.')
                full_modelconvolved = convolve_fft(full_model, psf, boundary='wrap')
                res_full_map = (science_data - full_modelconvolved)
            
            res_full_map_snr = res_full_map / noise_map
            Chisquare = np.nansum(res_full_map_snr * res_full_map_snr * mask2minimize)

        # Case 2: the data are polarized intensity data
        elif type_obs in ['polarized_intensity', 'polarised_intensity', 'polarized intensity', 'polarised intensity', 'polar']:
            models = []
            for j_struct in range(nb_disk_structures): # loop on each disk structure
                if display: print(f'Consider the disk structure #{j_struct+1}')
                dico['disk_struc'] = int(j_struct+1)
                params_values = params_3D[j_struct]
                dico['params_names_1model'] = params_names[j_struct]
                
                # Generate the disk model         
                models.append(generate_disk_model(params_values, dico = dico))
                full_model = np.nansum(models, axis=0)    

                #print(f'\n\n CHECK WHEN SEVERAL DISK STRUCTURES THAT THE SHAPE OF full_model IS A SHAPE OF AN IMAGE { np.shape(full_model)} \n\n')

            if do_robust_convolution:
                # In particular important for close-in or asymmetric disks. See Heikamp & Keller 2019
                Q_neg = full_model*np.cos(2*im_pa)
                U_neg = full_model*np.sin(2*im_pa)
                Q_neg_convolved = convolve_fft(Q_neg, psf)#, mode='same')
                U_neg_convolved = convolve_fft(U_neg, psf)#, mode='same')
                full_modelconvolved = (Q_neg_convolved*np.cos(2*im_pa)+U_neg_convolved*np.sin(2*im_pa)) # corresponds to Qphi 
        
            else: full_modelconvolved = convolve_fft(full_model, psf, boundary='wrap')

            res_full_map = (science_data - full_modelconvolved)
            res_full_map_snr = (science_data - full_modelconvolved)  / noise_map
            Chisquare = np.nansum(res_full_map_snr * res_full_map_snr * mask2minimize)
            

        # Save the results
        if save_some_results or save_full_results:
            # if yes, either save only some results or all of them
            # some: include the tested model, two residuals maps (the nominal one and the one normalized by the noise) and the parameters
            suffix_chisquare = "test" # f"_{Chisquare:.0f}".replace('.','p')
            suffix_chisquare =  f"_{Chisquare:.0f}".replace('.','p')
            suffix_nb_obs = f"_obs{i_obs+1:.0f}".replace('.','p')

            if type_obs == 'total_intensity': full_modelconvolved = convolve_fft(full_model, psf, boundary='wrap') # just save the model convolved at one parallactic angle position...   

            # Save some results
            if type_obs == 'total_intensity': 
                full_modelconvolved = convolve_fft(full_model, psf, boundary='wrap') # just save the model convolved at one parallactic angle position...
            fits.writeto(simu_resultdir+f'/disk_model_convolved{suffix_chisquare}{suffix_nb_obs}.fits',full_modelconvolved, overwrite=True)  
            fits.writeto(simu_resultdir+f'/residuals{suffix_chisquare}{suffix_nb_obs}.fits', res_full_map, overwrite=True)   
            fits.writeto(simu_resultdir+f'/residuals_snr{suffix_chisquare}{suffix_nb_obs}.fits', res_full_map, overwrite=True)
            fits.writeto(simu_resultdir+f'/params{suffix_chisquare}.fits', np.array(params_1D), overwrite=True)  
            #fits.writeto(os.path.join(simu_resultdir,f'params.fits'), np.array(params_1D), overwrite=True)  
            

            if save_full_results: # save additional files
                fits.writeto(simu_resultdir+f'/disk_model{suffix_chisquare}{suffix_nb_obs}.fits',full_model, overwrite=True)  

                if nb_disk_structures!=1: # for each disk structure, save its tested model
                    for j_struct in range(nb_disk_structures): 
                        suffix_nb_obs_nb_struct = f"_obs{i_obs+1:.0f}_struc{j_struct+1:.0f}"
                        fits.writeto(simu_resultdir+f'/disk_model{suffix_chisquare}{suffix_nb_obs_nb_struct}.fits', models[j_struct], overwrite=True)

                        if type_obs == 'total_intensity': full_modelconvolved = convolve_fft(models[j_struct], psf, boundary='wrap') # just save the model convolved at one parallactic angle position...
                        fits.writeto(simu_resultdir+f'disk_model_convolved{suffix_chisquare}{suffix_nb_obs_nb_struct}.fits',full_modelconvolved, overwrite=True)  
                    
                # regarding polarized intensity, save files about the robust convolution
                if type_obs in ['polarized_intensity', 'polarised_intensity', 'polar'] and do_robust_convolution:
                    fits.writeto(simu_resultdir+f'/disk_model_Q_convolved{suffix_chisquare}{suffix_nb_obs}.fits',-Q_neg_convolved, overwrite=True)      
                    fits.writeto(simu_resultdir+f'/disk_model_U_convolved{suffix_chisquare}{suffix_nb_obs}.fits',-U_neg_convolved, overwrite=True)  

                if type_obs == 'total_intensity': fits.writeto(simu_resultdir+f'/cube_obs-model_before_processing{suffix_chisquare}{suffix_nb_obs}.fits', cube_diff, overwrite=True)   
        
        Chisquare_red = Chisquare/(nb_resel-nb_free_params_per_obs-1)    

        if display: 
            print(f'-> Summary: Considering the observation #{i_obs+1} using {nb_disk_structures} disk structure(s), the simulation gives')
            #print_params_names_and_values_for_unflatten_list(params_names, [params_values], nb_obs, nb_disk_structures)
            print(f'  chisquare = {Chisquare:.4e} i.e. {Chisquare:.0f}, reduced chisquare = {Chisquare_red:.2f}')

    if display: print('------------------------------------------------------------------------------------')

    return Chisquare

#######################################################
############### Generate the disk model ###############
#######################################################

def generate_disk_model(params_values, dico = None):

    # Access keyword arguments from kwadrgs
    params_names  = dico.get('params_names_1model','PROBLEM')
    #print('l405', params_names)
    dimension     = dico.get('dimension', 'PROBLEM')
    plate_scale  = dico.get('plate_scale', 'PROBLEM')
    distance_star = dico.get('distance_star', 'PROBLEM')
    disk_model_polar = dico.get('disk_model_polar', 'PROBLEM')
    convention_unit  = dico.get('convention_unit', 'PROBLEM')
    disk_struc = dico.get('disk_struc', 'PROBLEM') # 1 = if disk structure #1, 2 = if disk structure #2, etc
    idx_obs    = dico.get('idx_obs', 'PROBLEM') # 1 = if obs #1, 2 = obs #2, etc: useful for the value of the parameter SCALING
    display    = dico.get('display', 0)

    if type(params_names) == type(''): raise ValueError('There is a problem with the value of names.')
    
    if len(params_names) != len(params_values): 
        raise ValueError(f'The number of parameters listed in the list "NAMES" should be the same than the number of parameters in "params". \nNAMES = {params_names} \nparams = {params_values}')
    
    # Check if some extra keyword arguments were effectively given
    if disk_struc == 'PROBLEM': raise ValueError('There is a problem with the value of disk_struc.')
    if type(dimension) == type(''): raise ValueError('There is a problem with the value of dimension.')
    if type(plate_scale) == type(''): raise ValueError('There is a problem with the value of plate_scale.')
    if type(distance_star) == type(''): raise ValueError('There is a problem with the value of distance_star.')
    if type(disk_model_polar) == type(''): raise ValueError('There is a problem with the value of disk_model_polar.')
    if convention_unit == 'PROBLEM': raise ValueError('There is a problem with the value of convention_unit.')

    if display: 
        print_params_names_and_values(params_names, params_values, params_flat=True)

    # Assigned the values of the fixed parameters
    if 'RAD' not in params_names: r1  = dico.get(f'rad_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'PA' not in params_names:  pa  = dico.get(f'pa_init_struct{disk_struc}',  ['PROBLEM'])[0]
    if 'INC' not in params_names: inc = dico.get(f'inc_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'G1' not in params_names:  g1  = dico.get(f'g1_init_struct{disk_struc}',  ['PROBLEM'])[0]
    if 'G2' not in params_names:  g2  = dico.get(f'g2_init_struct{disk_struc}',  ['PROBLEM'])[0]
    if 'ARGPERI' not in params_names: argperi = dico.get(f'argperi_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'ECC' not in params_names:  ecc    = dico.get(f'ecc_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'KSI0' not in params_names: ksi0   = dico.get(f'ksi0_init_struct{disk_struc}',['PROBLEM'])[0]
    if 'GAMMA' not in params_names: gamma = dico.get(f'gamma_init_struct{disk_struc}',['PROBLEM'])[0]
    if 'BETA' not in params_names:  beta  = dico.get(f'beta_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'ALPHA' not in params_names: alpha = dico.get(f'alpha_init_struct{disk_struc}',['PROBLEM'])[0]
    if 'AIN' not in params_names:   ain   = dico.get(f'ain_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'AOUT' not in params_names:  aout  = dico.get(f'aout_init_struct{disk_struc}', ['PROBLEM'])[0]
    if 'SCALING' not in params_names:  scaling  = dico.get(f'scaling_init_struct{disk_struc}', ['PROBLEM'])[0]

    # Assigned the values of the free parameters
    for i, name in enumerate(params_names):
        if name == 'SCALING': scaling = params_values[i][idx_obs] ; #print('scaling', scaling)
        if name == 'RAD':   r1   = params_values[i][0]
        if name == 'PA':    pa   = params_values[i][0]
        if name == 'INC':   inc  = params_values[i][0]
        if name == 'G1':    g1   = params_values[i][0]
        if name == 'G2':    g2   = params_values[i][0]
        if name == 'AIN':   ain  = params_values[i][0]
        if name == 'AOUT':  aout = params_values[i][0]
        if name == 'ARGPERI': argperi = params_values[i][0]
        if name == 'ECC':   ecc   = params_values[i][0]
        if name == 'KSI0':  ksi0  = params_values[i][0]
        if name == 'GAMMA': gamma = params_values[i][0]
        if name == 'BETA':  beta  = params_values[i][0]
        if name == 'ALPHA': alpha = params_values[i][0]

        
    if 'PROBLEM' in [r1, pa, inc, g1, argperi, ecc, ksi0, gamma, beta, ain, aout, scaling] :
        raise ValueError(f'One of the disk parameter was forgotten in the parameters configuration file. \n radius = {r1} au \n PA = {pa} deg \n inc = {inc} \n g1 = {g1} \n argperi = {argperi} \n ecc = {ecc} \n ksi0 = {ksi0} \n gamma = {gamma} \n beta = {beta} \n alpha = {alpha} \n ain = {ain} \n aout = {aout}')

    if convention_unit == 'MCMC':
        scaling = 10**scaling #mt.exp(params[0])
        inc = np.arccos(inc) * 180/np.pi

    # Generate the model
    model = vip.fm.scattered_light_disk.ScatteredLightDisk(nx=dimension,
                               ny=dimension,
                               distance=distance_star,
                               itilt=inc,
                               omega=argperi,
                               pxInArcsec=plate_scale,
                               pa= pa,
                               density_dico={
                                   'name': '2PowerLaws',
                                   'ain': ain,
                                   'aout': aout,
                                   'a': r1,
                                   'e': ecc,
                                   'ksi0': ksi0,
                                   'gamma': gamma,
                                   'beta': beta
                               },
                               spf_dico={
                                   'name': 'HG',
                                   'g': g1,
                                   'polar': disk_model_polar
                               })
    
    return model.compute_scattered_light() * scaling
