#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script containing functions to make two first tests before running MCMC or Nelder-Mead simulations to fit disks imaged in scattered-light.

Functions:
    . do_test_disk_empty()
    . do_test_disk_first_guess()
  
This script would be imported by the scripts:
    ..run_modisc.py
    ..plot_mcmc_results.py
'''

__author__ = 'Celia Desgrange'

# Import classic packages
from ..packages import * 

# Import MoDiSc functions 
from ..functions.simulations import *
from ..functions.parameters_rearrange_print import *


def do_test_disk_empty(PARAMS_INIT, dico = None):

    display = dico.get('display', 0)

    if display: 
        print('\n=== Compute the likelihood and chisquare within the region to optimize (defined by MASK2MINIMIZE) *without* subtraction of any disk MODEL in the SCIENCE_DATA === \nRemark: Following disk simulations have to result in better likelihood and chisquare value than the ones computed here.')
        startTime =  datetime.now() 
    
    # Access extra keyword arguments from the dictionary
    PARAMS_NAMES = dico.get('params_names_unflat', 'PROBLEM')

    params_init_no_disk = copy.deepcopy(PARAMS_INIT)

    # Set the scaling parameter to 0. Loop on all the parameters to find it
    for i, params_names_struct_i in enumerate(PARAMS_NAMES):
        for j, param_name in enumerate(params_names_struct_i):
            if param_name[:len('SCALING')] == 'SCALING': 
                params_init_no_disk[i][j] = [0] * len(params_init_no_disk[i][j])

    if display: 
        #print(f'Check disk model parameters without flux (scaling_flux should be set to 0):\n')
        print_params_names_and_values(PARAMS_NAMES, params_init_no_disk, params_flat=False)
        print('Remark: The scaling parameter defining the flux of the disk should be equal to 0 here.')
    
    chisquare_init_nodisk = chisquare_params_flat_or_unflat(params_init_no_disk, dico = dico,  params_flat=False)
    
    if display:  print(f"\nTime for a single model: {datetime.now() - startTime}")

    return chisquare_init_nodisk


def do_test_disk_first_guess(PARAMS_INIT, dico = None):

    display = dico.get('display', 0)
    
    if display: 
        print('\n === Compute the likelihood and chisquare for the first guess disk model ===\n')
        startTime =  datetime.now() 
    
    if display: 
        # Access extra keyword arguments from the dictionary    
        PARAMS_NAMES = dico.get('params_names_unflat', 'PROBLEM')

        print(f'Check the disk model parameters for the initial guess:')
        print_params_names_and_values(PARAMS_NAMES, PARAMS_INIT, params_flat=False)
    
    chisquare_init = chisquare_params_flat_or_unflat(PARAMS_INIT, dico = dico,  params_flat=False)

    if display: print(f"\nTime for a single model: {datetime.now() - startTime}")

    return chisquare_init
    