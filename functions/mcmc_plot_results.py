#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script containing functions to plot the results of the MCMC simulations fitting disks imaged in scattered-light.

Functions:
    . open_chains() - Open the chain by loading the backend file
    . make_chains_plot() - Plot the Figure showing the chains for each walker and parameter
    . make_corner_plot() - Plot the Figure showing the cornerplot with the posteriors of the simulations
    . MCMC_parameters_given_percentile() - Save disk parameters at given percentiles
  
This script would be imported by the scripts:
    ..run_modisc.py
    ..plot_mcmc_results.py
'''

__author__ = 'Celia Desgrange'

# Import classic packages
from MoDiSc.packages import * 

# Import MoDiSc functions 
from MoDiSc.functions.parameters_rearrange_print import *

rcParams.update({
    'xtick.direction'  : "in", 
    'ytick.direction'  : "in",
    'xtick.major.width': 1,
    'xtick.major.size' : 8,
    'xtick.minor.width': 0.8,
    'xtick.minor.size' : 5,
    'ytick.major.width': 1.,
    'ytick.major.size' : 8, 
    'ytick.minor.width': 0.8,
    'ytick.minor.size' : 5,
    'ytick.right': True,
    'xtick.top'  : True,
        })


########################################################
###### Open the chain by loading the backend file ######
########################################################

def open_chains(dico, return_suffix=False, output_chains_shape='unflat'):
    '''
    Open the chain by loading the backend file
    '''

    display     = dico.get('display')
    if display: print('\n(Open the chains of the MCMC simulation)')
    mcmc_chains_binning = dico.get('mcmc_chains_binning')
    mcmc_chains_burnin  = dico.get('mcmc_chains_burnin')
    mcmc_chains_log_prob_criterion_fraction = dico.get('mcmc_chains_log_prob_criterion_fraction')
    mcmc_chains_apply_selection_criterion = dico['mcmc_chains_apply_selection_criterion']

    mcmc_resultdir = dico['mcmc_resultdir']
    name_h5 = 'backend_file_mcmc'

    # Load the MCMC chain
    reader = backends.HDFBackend(os.path.join(mcmc_resultdir, name_h5 + '.h5')) # open the file

    # Check the value of the burnin
    tau = reader.get_autocorr_time(tol=0)
    if mcmc_chains_burnin > reader.iteration - 1: raise ValueError(f"The burnin (= {mcmc_chains_burnin}) cannot be larger than the number of iterations (= {reader.iteration}).")

    # Load all the chains
    chains_full      = reader.get_chain(discard=0, thin=mcmc_chains_binning) # get the chain (without applying the burnin, this is to keep the information when we plot the chains)
    chains           = reader.get_chain(discard=mcmc_chains_burnin, thin=mcmc_chains_binning) # get the chain (by applying the burnin)
    log_prob_samples = reader.get_log_prob(discard=mcmc_chains_burnin, thin=mcmc_chains_binning) # get the log probability of the chain
    suffix = '_all_walkers'

    if display:
        chisquare_min, chisquare_max = 2 * np.nanmin( np.abs(log_prob_samples) ), 2 * np.nanmax( np.abs(log_prob_samples) )
        print(f"The first {mcmc_chains_burnin} iterations are removed for each walker.")
        print(f"The shape of the MCMC chains is {chains.shape}.")
        print(f"After the burn-in cut is applied, \n - the minimum likelihood is {np.nanmin(log_prob_samples):.0f}, corresponding to a chisquare of {chisquare_min:.0f}")
            
        print(f" - the maximum likelihood is {np.nanmax(log_prob_samples):.0f}, corresponding to a chisquare of {chisquare_max:.0f}. ")

        if mcmc_chains_binning != 1: print(f"Chains are binned by a factor {mcmc_chains_binning}.")

    # Do we want to remove some of the chains for which walkers did not converge?
    if mcmc_chains_apply_selection_criterion: # remove some of the chains for which the walkers did not converge
        log_prob_samples = reader.get_log_prob(discard=mcmc_chains_burnin, thin=mcmc_chains_binning)

        nb_walkers_tot = np.shape(log_prob_samples)[1] # total number of walkers

        # Keep all the chains for which walkers ended up at a log probability of 'mcmc_log_prob_criterion_fraction' * the minimal log probability
        log_prob_threshold  = mcmc_chains_log_prob_criterion_fraction * np.nanmin( np.abs(log_prob_samples[-1]) )
        chisquare_threshold =  2 * np.abs(log_prob_threshold)

        if display: print(f'Threshold: values below -{log_prob_threshold:.0f} (chisquare < {chisquare_threshold:.0f}) are kept.')

        cond_good_log_prob = np.abs(log_prob_samples[-1])  < log_prob_threshold  # define the condition of the values below the threshold
        chains             = chains[:,cond_good_log_prob]                        # select the according chains
        log_prob_samples   = log_prob_samples[cond_good_log_prob]
        
        nb_walkers_good = np.shape(chains)[1] # number of walkers with a log likelihood below the threshold

        if display: 
            print(f'We therefore keep {nb_walkers_good:.0f}/{nb_walkers_tot:.0f} simulations, that is to say, {nb_walkers_good/nb_walkers_tot*100:.1f}% of simulations.')
            print(f'\nRemark: the final log probability values for the {np.shape(log_prob_samples)[-1]} walkers are\n', log_prob_samples[-1] )

        suffix = f'_{nb_walkers_good:.0f}walkers_{mcmc_chains_log_prob_criterion_fraction}fraction'

        print(f"The shape of the chains after selection is {chains.shape} and without applying the burn-in.") # burnin is *not* taking into account
        print(f"The shape of log_prob_samples after selection is {np.shape(log_prob_samples)} and with applying the burn-in.") # burnin is taking into account


    if output_chains_shape == 'flat':
        chains_full = chains_full.reshape(chains_full.shape[0] * chains_full.shape[1], chains_full.shape[2])
        chains = chains.reshape(chains.shape[0] * chains.shape[1], chains.shape[2])
        log_prob_samples = log_prob_samples.flatten()

        if display: print(f'Reshape the arrays chains and log_prob_samples from 3D arrays to 2D arrays. Loose individual walker information. Their shapes are now {np.shape(chains)} for the arrays chains and {np.shape(log_prob_samples)} for the array log_prob_samples.\n')

    if return_suffix: return chains_full, chains, log_prob_samples,  suffix
    else: return chains_full, chains, log_prob_samples

########################################################
########## Plot Figures (chains, cornerplot) ###########
########################################################

def make_chains_plot(dico):
    '''
    Plot the Figure showing the chains for each walker and parameter
    '''
   
    display     = dico.get('display', 1)
    if display: print('\n(Plot the chains of the MCMC simulation)')
    mcmc_chains_burnin = dico.get('mcmc_chains_burnin')
    convention_unit = dico.get('convention_unit')

    # Figure - parameters
    mcmc_fig_size_factor = dico['mcmc_fig_size_factor']
    params_labels_flat = dico['params_labels_flat']
    params_names_flat  = dico['params_names_flat']

    # Load the MCMC chain
    mcmc_resultdir = dico['mcmc_resultdir']
    name_h5 = 'backend_file_mcmc'

    reader = backends.HDFBackend(os.path.join(mcmc_resultdir, name_h5 + '.h5')) # open the file

    chains_full, chains, log_prob, suffix = open_chains(dico, return_suffix=True, output_chains_shape='unflat') # get the log probability of the chain

    if display: print('\n(back to Plot the chains)')

    # Dimension, number of walkers
    nb_iter, nb_walkers, nb_dim_mcmc = chains_full.shape[0], chains_full.shape[1], chains_full.shape[2]

    
    # Make the unit 'user-friendly'. Currently implemented for the scaling and inclination parameters.
    chains = correct_unit(params_names_flat, chains_full, convention_unit=convention_unit, output_unit='user-friendly', axis_to_be_considered=-1, display=display)

    # Figure
    _, axarr = plt.subplots(nb_dim_mcmc, sharex=True, figsize=(6, 1 + nb_dim_mcmc))

    # Case: there is more than one parameter fitted
    if nb_dim_mcmc > 1:
        for i in range(nb_dim_mcmc): # loop on each parameter
            axarr[i].set_ylabel(params_labels_flat[i].replace('!', ' '), fontsize = 8 * mcmc_fig_size_factor) #* mcmc_chains_binning) # update labels size for the y axis
            axarr[i].tick_params(axis='y', labelsize = 8 * mcmc_fig_size_factor)#4* mcmc_chains_binning) # update params size for the y axis

            for j in range(nb_walkers): # loop on each walker
                axarr[i].plot(chains[:, j, i], linewidth= 0.5 * mcmc_fig_size_factor) #mcmc_chains_binning/2)

            axarr[i].axvline(x=mcmc_chains_burnin, color='black', linewidth=1.5 * mcmc_fig_size_factor) # add the burnin cut # * mcmc_chains_binning)

            # Customize ticks layout
            axarr[i].tick_params(axis='x',which='major',direction="in", width=1,length=5 * mcmc_fig_size_factor)
            axarr[i].tick_params(axis='x',which='minor',direction="in", width=0.8,length=4 * mcmc_fig_size_factor)
            axarr[i].tick_params(axis='y',which='major',direction="in", width=1,length=5 * mcmc_fig_size_factor)
            axarr[i].tick_params(axis='y',which='minor',direction="in", width=0.8,length=4 * mcmc_fig_size_factor)
            axarr[i].xaxis.set_minor_locator(plt.MultipleLocator(reader.iteration//20)) # add minor axis

        axarr[nb_dim_mcmc - 1].tick_params(axis='x', labelsize = 8 * mcmc_fig_size_factor)   # update params size for the x axis
        axarr[nb_dim_mcmc - 1].set_xlabel('Iterations', fontsize = 8 * mcmc_fig_size_factor) # update labels size for the x axis
        axarr[nb_dim_mcmc - 1].set_xlim([0, nb_iter])

     # Case: there is only one parameter fitted
    else: 
        axarr.set_ylabel(params_labels_flat[0].replace('!', ' '), fontsize = 8 * mcmc_fig_size_factor)
        axarr.tick_params(axis='y', labelsize = 8 * mcmc_fig_size_factor)

        for j in range(nb_walkers): # loop on each walker
            axarr.plot(chains[:, j, 0], linewidth = 1/2 * mcmc_fig_size_factor) # plot its chain 
            if params_names_flat[0].lower() == 'scaling': axarr.set_yscale('log')
        
        axarr.axvline(x=mcmc_chains_burnin, color='black', linewidth = 1.5 * mcmc_fig_size_factor) # add the burnin cut 

        axarr.set_xlabel('Iterations', fontsize=8 * mcmc_fig_size_factor)  # update labels size for the x axis
        axarr.set_xlim([0, nb_iter])

        # Customize ticks layout
        axarr.tick_params(axis='x',which='major',direction="in", width=1,length=5 * mcmc_fig_size_factor)
        axarr.tick_params(axis='x',which='minor',direction="in", width=0.8,length=4 * mcmc_fig_size_factor)
        axarr.tick_params(axis='y',which='major',direction="in", width=1,length=5 * mcmc_fig_size_factor)
        axarr.tick_params(axis='y',which='minor',direction="in", width=0.8,length=4 * mcmc_fig_size_factor)
        axarr.xaxis.set_minor_locator(plt.MultipleLocator(reader.iteration//20)) # add minor axis

    namesave = f'chains{suffix}_{mcmc_chains_burnin}burnin.pdf'
    if display: print(f'The figure is saved under the name {namesave}.')

    plt.tight_layout()

    plt.savefig(os.path.join(mcmc_resultdir, namesave)) # save the figure
    plt.close()


def make_corner_plot(dico):
    '''
    Plot the Figure showing the cornerplot with the posteriors of the simulations
    '''

    display = dico.get('display')
    if display: print('\n(Plot the cornerplot with the posteriors of the MCMC simulation)')
    mcmc_chains_burnin  = dico.get('mcmc_chains_burnin')
    mcmc_nb_walkers = dico.get('mcmc_nb_walkers')
    nb_dim_mcmc     = dico.get('mcmc_nb_dimension')
    convention_unit = dico.get('convention_unit')

    params_labels_flat = dico.get('params_labels_flat')
    params_names_flat  = dico.get('params_names_flat')
    nb_obs   = dico.get('nb_obs')
    epochs   = dico.get('epochs_all')
    instrus  = dico.get('instru_all')

    # Figure
    fig_size_factor  = dico.get('mcmc_fig_size_factor')
    cornerplot_sigma = dico.get('mcmc_fig_cornerplot_sigma')
    cornerplot_add_text_annot = dico.get('mcmc_fig_cornerplot_add_text_annot')
    
    
    #if nobs == 1: file_prefix = str(epoch)
    #else: file_prefix = f'{nobs}epochs'
    #print('file_prefix', file_prefix)
    mcmc_resultdir = dico.get('mcmc_resultdir')
    name_h5        = 'backend_file_mcmc'

    # Load the MCMC chain
    reader = backends.HDFBackend(os.path.join(mcmc_resultdir, name_h5 + '.h5'))
    chains_full_flat, chains_flat, log_prob_flat, suffix = open_chains(dico=dico, return_suffix=True, output_chains_shape='flat')
    
    # Make the unit 'user-friendly'. Currently implemented for the scaling and inclination parameters.
    chains_full_flat = correct_unit(params_names_flat, chains_full_flat, convention_unit=convention_unit, output_unit='user-friendly', axis_to_be_considered=-1, display=display)

    # Check parameters for the best chi2
    chisquare_flat = -2*np.copy(log_prob_flat)
    idx_min_chi2 = np.argmax(log_prob_flat)
    params_bestfit = chains_flat[idx_min_chi2]

    # Save min/max chi2red
    dF = pd.DataFrame([params_bestfit], columns=params_names_flat)
    dF.insert(len(dF.iloc[0]), 'chi2_min', [chisquare_flat[idx_min_chi2]])
    dF.to_csv(os.path.join(mcmc_resultdir, 'params_chi2_min' + '.csv'), sep=';')


    print('\nParameters for the best chi2:')
    if display: print_params_names_and_values(params_names=params_names_flat, params_values=params_bestfit, params_flat=True)

    # Figure
    # Title for each diagonal plot: defined via the quantiles
    # Value at 50% is the center of the Normal law ; value at 50% - value at 15.9% is -1 sigma ;value at 84.1%% - value at 50% is 1 sigma
    if cornerplot_sigma == 1: quants = (0.159, 0.5, 0.841)
    if cornerplot_sigma == 2: quants = (0.023, 0.5, 0.977)
    if cornerplot_sigma == 3: quants = (0.001, 0.5, 0.999)

    shouldweplotalldatapoints = False
    print('params_labels_flat', params_labels_flat)
    labels_hash    = [params_labels_flat[i] for i in range(nb_dim_mcmc)]
    labels_xy_axis = [params_labels_flat[i].replace('!',' ') for i in range(nb_dim_mcmc)]
    print('labels_hash', labels_hash)

    TITs, RANGE = [], []
    print('')
    for i in range(nb_dim_mcmc):
        chain_1param = chains_flat[:,i]
        values = np.percentile(chain_1param, list(np.array(quants)*100))
        q1, median_value, q2 = values[0], values[1], values[2]
        err_plus, err_moins = q2-median_value, median_value-q1
        lab = labels_hash[i].split('!')

        if params_names_flat[i].lower() == 'scaling': 
            median_value_str = r'$'   + f'{median_value:.2e}' 
            err_plus_str     = r'^{+' + f'{err_plus:.0e}'     + '}'
            err_moins_str    = r'_{-' + f'{err_moins:.0e}'    + '}$'
            TITs.append(f'{lab[0]} = \n{median_value_str}{err_plus_str}{err_moins_str} {lab[1]}')

        elif params_names_flat[i].lower() in ['rad',  'pa', 'inc', 'ain', 'beta', 'gamma', 'argperi']: 
            median_value_str = r'$'   + f'{median_value:.1f}' 
            err_plus_str     = r'^{+' + f'{err_plus:.1f}'     + '}'
            err_moins_str    = r'_{-' + f'{err_moins:.1f}'    + '}$'
            TITs.append(f'{lab[0]} = {median_value_str}{err_plus_str}{err_moins_str} {lab[1]}')

        else: 
            median_value_str = r'$'   + f'{median_value:.2f}' 
            err_plus_str     = r'^{+' + f'{err_plus:.2f}'     + '}'
            err_moins_str    = r'_{-' + f'{err_moins:.2f}'    + '}$'
            TITs.append(f'{lab[0]} = {median_value_str}{err_plus_str}{err_moins_str} {lab[1]}')

        print('For this parameter, the title will be:\n', TITs[-1])

        plusmoins, fac = (err_plus+err_moins)/2, 5
        RANGE.append([median_value-plusmoins*fac, median_value+plusmoins*fac])

    # Update the size of the labels/ticks
    plt.rcParams.update({
    "axes.titlesize"  : 16 * fig_size_factor,
    "axes.labelsize"  : 16 * fig_size_factor,
    "xtick.labelsize" : 16 * fig_size_factor,
    "ytick.labelsize" : 16 * fig_size_factor,
    })
    # Initialize the corner plot    
    fig = corner.corner(chains_flat,
                        labels=labels_xy_axis, # corresponds to the xy labels 
                        titles=TITs, # corresponds to the titles on the plots located on the diagonal
                        quantiles=quants,
                        show_titles=True,
                        title_fmt = None,
                        plot_datapoints=shouldweplotalldatapoints,
                        verbose=False,
                        range = RANGE
                        )
    

    # Extract the axes
    axes = np.array(fig.axes).reshape((nb_dim_mcmc, nb_dim_mcmc))
    
    # Add vertical lines in the diagonal plots (distribution plots) corresponding to the best theta values
    for i in range (nb_dim_mcmc):
        axes[i,i].plot([params_bestfit[i], params_bestfit[i]], axes[i,i].get_ylim(), ls='-', lw=1.5, color='teal', zorder=0)
        
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # Main title of the figure
    if nb_dim_mcmc > 1 and cornerplot_add_text_annot: # defined only if more than one parameter was fitted
        instru_epoch_str = ''
        for i in range(nb_obs):
            instru_epoch_str += f'{instrus[i]} ({epochs[i]})'
            if i != len(epochs)-1: 
                if i == len(epochs)-2: instru_epoch_str += ' and '
                else: instru_epoch_str += ', '
                
        text_annot = f"Based on {instru_epoch_str} data:\n{reader.iteration:.0f} iterations of {mcmc_nb_walkers:.0f} walkers. Burn-in = {mcmc_chains_burnin:.0f} iterations. \nIn total, {reader.iteration * mcmc_nb_walkers:.0f} simulations."
        fig.gca().annotate(text_annot, xy=(0.55, 0.99), xycoords="figure fraction", xytext=(-20, -10), textcoords="offset points", ha="center", va="top", fontsize = 18 * fig_size_factor)
    

    plt.savefig(os.path.join(mcmc_resultdir, f'cornerplot{suffix}_{mcmc_chains_burnin}burnin_{cornerplot_sigma}sigma.pdf'))
    plt.close()

########################################################
######## Save parameters at given percentiles ##########
########################################################
def MCMC_parameters_given_percentile(dico, percentile, WRITETO=True):
    '''
    Save disk parameters at given percentiles
    '''
    display = dico.get('display')
    if display: print('\n(Look for the disk parameters at specific percentile(s))')


    if len(np.shape(percentile)) != 1:
        raise ValueError(f'The dimension of the input argument percentile should be one, not {len(np.shape(percentile))} because np.shape(percentile) = {np.shape(percentile)} .')
    
    mcmc_resultdir = dico['mcmc_resultdir']
    nb_dim_mcmc    = dico.get('mcmc_nb_dimension')

    params_names_flat = dico.get('params_names_flat')
    convention_unit   = dico.get('convention_unit')
    
    # Load the MCMC chain
    chains_full_flat, chains_flat, log_prob_flat, suffix = open_chains(dico=dico, return_suffix=True, output_chains_shape='flat')
    
    # Change scaling and inclination values to physical
    chains_flat = correct_unit(params_names_flat=params_names_flat, params_values=chains_flat, convention_unit=convention_unit, output_unit='user-friendly', axis_to_be_considered=-1, display=display)
    
    # Look for the each parameter its value for several percentiles
    parameters_percentile = []

    for j in range(len(percentile)): # loop on the percentiles of interest
        L = []
        for i in range(nb_dim_mcmc): # loop on each parameter
            L.append( np.percentile(chains_flat[:,i], percentile[j]) ) # chains_flat[:,i] = all the values found out for the parameter i
        parameters_percentile.append(L)

    parameters_percentile = np.array(parameters_percentile )

    if WRITETO:
        dF = pd.DataFrame(parameters_percentile, columns=params_names_flat)
        dF.insert(len(dF.iloc[0]), 'percentile', percentile)
        path = os.path.join(mcmc_resultdir, f'params_percentile{suffix}' + '.csv')
        dF.to_csv(path, sep=';')

        if display: 
            print('\nThe following dataframe storing the values of the parameters for different percentiles...\n\n', dF.head())
            print(f' \n... have been saved at\n {path} \n')

    return np.array(parameters_percentile)


