#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script containing functions to plot the results of the MCMC simulations fitting disks imaged in scattered-light.

Functions: (version of 2025-01-15)
    . compute_im_pa_grid()   - Returns a 2D-array whose values are position angles from the top to the right (so counterclockwise) by default.
   
This script would be imported by the scripts:
    ..run_modisc.py
    ..plot_mcmc_results.py
'''

__author__ = 'Celia Desgrange'

import numpy as np

def compute_im_pa_grid(im, center='n//2', xc=None, yc=None, return_unit='deg', counterclockwise=True, display=0):
    '''
    Returns a 2D-array whose values are position angles from the top to the right (so counterclockwise) by default.

    Inputs:
        .'im' (2D array): only the shape of 'im' will be used, not its values.

        (optional)

        .'center' (string): indicates where is the center of the image. Three options: 'n//2-0.5', 'n//2-1' or 'n//2', 'coordinates_are_given'.
            Default value: center='n//2-0.5'. 
                Examples: 
                    For an even image, this means the center is between the four centered-pixels.
                    For an odd image, this means the center is the pixel at the center.
            Other possible values: 'n//2-1' or 'n//2' (but only if the image is even) or 'coordinates_are_given'.
            If 'coordinates_are_given', specify the parameters 'xc' and 'yc'.
            (!) Important: Assumption: count starts at 0.

        .'xc' (float): The coordinate in the abscisse axis of the center of the image. Default value: None.

        .'yc' (float): The coordinate in the abscisse axis of the center of the image. Default value: None.

        .'counterclockwise' (boolean or binary): if set to True (default value), position angles increase from the top to the right. If set to False, position angles increase from the top to the left.

        .'return_unit' (string): The unit of the pixels.  Default value: 'deg' (in degree).

        .'display' (boolean or binary value): prints information if sets to True or 1 (which is the default value)

    Output:
        .'thetas' (2D array): each pixel has for value its position angle from the North (up) of the image towards the East (left; astro convention) of the image. The unit is in degrees if 'return_unit' is set to 'deg', otherwise in radians. The shape of 'thetas' is the same than the input array 'im'.
    '''
    if display: print("\n(Compute 2D position angle grid for an even or odd image) \n (!) Don't forget to select the good option for the center of the image.")

    nx, ny = np.shape(im)[1], np.shape(im)[0]

    if nx % 2 == 0  and ny % 2 == 0: even_or_odd = 'even'
    elif nx % 2 == 1  and ny % 2 == 1: even_or_odd = 'odd'
    else: raise ValueError(f"The image should be in both axis either even or odd. However, its shape is {np.shape(im)}.")

    if center == 'coordinates_are_given': 
        xend, yend = nx-xc, ny-yc
        xbeg, ybeg = -xc, -yc
        
        x,y = np.arange(xbeg,xend,1), np.arange(ybeg, yend,1)
        xs,ys = np.meshgrid(x, y, sparse=True)
        zs = np.sqrt(xs**2 + ys**2)

    elif even_or_odd == 'even':
        if center == 'n//2-0.5': x,y = np.linspace(-nx//2,nx//2,nx), np.linspace(-ny//2,ny//2,ny)
        elif center == 'n//2': x,y = np.linspace(-nx//2,nx//2-1,nx), np.linspace(-ny//2,ny//2-1,ny)
        elif center == 'n//2-1': x,y = np.linspace(-nx//2+1,nx//2,nx), np.linspace(-ny//2+1,ny//2,ny)
        else: raise ValueError(f"If the image has an even shape, the value of 'center' should be 'n//2-0.5', 'n//2' or 'n//2-1' and not {center}.")

    elif even_or_odd == 'odd':
        if center == 'n//2-0.5': x,y = np.linspace(-nx//2+1,nx//2,nx), np.linspace(-ny//2+1,ny//2,ny)
        else: raise ValueError(f"If the image has an odd shape, the value of 'center' should be 'center' and not {center}.")

    xs,ys  = np.meshgrid(x, y, sparse=True)
    zs     = np.sqrt(xs**2 + ys**2)
    thetas = 2*np.arctan( xs/(ys+zs)) + np.pi 
    thetas[np.isnan(thetas)] = 0

   
    thetas = thetas[::-1,:] 
    if not counterclockwise: thetas = thetas[:,::-1] 

    if return_unit == 'deg': return thetas*180/np.pi
    else: return thetas