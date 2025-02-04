#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script containing functions to plot the results of the MCMC simulations fitting disks imaged in scattered-light.

Functions: (version of 2025-01-15)
    . compute_im_rad_grid()  -  Returns a 2D radius grid with the same shape than the input image 'im'
    . compute_limdet_map_ann -  Returns the detection limit 2D map. This map is computed by annular noise estimation with slippery annulus. Calls the function compute im_rad_grid() to derive the distance of each pixel to the center of the image.
   
This script would be imported by the scripts:
    ..functions.dataset_load_and_process.py
'''

__author__ = 'Celia Desgrange'

import numpy as np
import time

def compute_im_rad_grid(im, center='n//2-0.5', xc=None, yc=None, display=1):
    '''
    Returns a 2D radius grid with the same shape than the input image 'im'. The unit is in pixels.

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

        .'xc' (float): The coordinate in the abscisse axis of the center of the image.

        .'yc' (float): The coordinate in the abscisse axis of the center of the image.

        .'display' (boolean or binary value): prints information if sets to True or 1 (which is the default value)

    Output:
        .'zs' (2D array): each pixel has for value its distance to the center of the image. The unit is in pixels. The shape of 'zs' is the same than the input array 'im'.

    ---------
    Examples:
        .compute_im_rad_grid(np.zeros((4,4)), center='center') returns

            array( [[2.82842712, 2.10818511, 2.10818511, 2.82842712],
                    [2.10818511, 0.94280904, 0.94280904, 2.10818511],
                    [2.10818511, 0.94280904, 0.94280904, 2.10818511],
                    [2.82842712, 2.10818511, 2.10818511, 2.82842712]])

        .compute_im_rad_grid(np.zeros((4,4)), center='n//2') returns 

            array( [[2.82842712, 2.23606798, 2.        , 2.23606798],
                    [2.23606798, 1.41421356, 1.        , 1.41421356],
                    [2.        , 1.        , 0.        , 1.        ],
                    [2.23606798, 1.41421356, 1.        , 1.41421356]])

        .compute_im_rad_grid(np.zeros((4,4)), center='n//2-1')

            array( [[1.41421356, 1.        , 1.41421356, 2.23606798],
                    [1.        , 0.        , 1.        , 2.        ],
                    [1.41421356, 1.        , 1.41421356, 2.23606798],
                    [2.23606798, 2.        , 2.23606798, 2.82842712]])

        .compute_im_rad_grid(np.zeros((n,n)), 'center') returns

            array( [[1.41421356, 1.        , 1.41421356],
                    [1.        , 0.        , 1.        ],
                    [1.41421356, 1.        , 1.41421356]])
    '''
    if display: print("\n(Compute 2D radius grid for an even or odd image) \n (!) Don't forget to select the good option for the center of the image.")

    nx = np.shape(im)[1]
    ny = np.shape(im)[0]

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

    xs,ys = np.meshgrid(x, y, sparse=True)
    zs = np.sqrt(xs**2 + ys**2)
    return zs


def compute_limdet_map_ann(res_map, dr, alpha=2, center='n//2-0.5', xc=None, yc=None, display=1):
    '''
    Returns the detection limit 2D map. This map is computed by annular noise estimation with slippery annulus. Calls the function compute im_rad_grid() to derive the distance (in pixel) of each pixel to the center of the image.

    Inputs:
        .'res_map' (2D-array): residual map
        
        .'dr' (float): width of the annulus for which the detection limit is computed

        (optional)

        .'alpha' (float): factor to consider bigger annulus to derive the noise. The annulus is located in an interval [-dr-alpha; dr+alpha] at a given radius. The goal is to have a smoother transition between annuli.
        
        .'center' (string): indicates where is the center of the image. Three options: 'n//2-0.5', 'n//2-1' or 'n//2', 'coordinates_are_given'.
            Default value: center='n//2-0.5'. 
                Examples: 
                    For an even image, this means the center is between the four centered-pixels.
                    For an odd image, this means the center is the pixel at the center.
            Other possible values: 'n//2-1' or 'n//2' (but only if the image is even) or 'coordinates_are_given'.
            If 'coordinates_are_given', specify the parameters 'xc' and 'yc'.
            (!) Important: Assumption: count starts at 0.
            Additional information is available in the docstring of the function compute_im_rad_grid()

        .'xc' (float): The coordinate in the abscisse axis of the center of the image.

        .'yc' (float): The coordinate in the abscisse axis of the center of the image.

        .'display' (boolean or binary value): prints information if sets to True or 1 (which is the default value)

    Output:
        .'im_noise' (2D-array): detection limit map. 
    '''

    t0 = time.time()
    if display: print("\n(Computing the detection limit map by using the 1D-annulus method)")
    h, w = np.shape(res_map)

    im_noise, im_nan = np.zeros((h,w)), np.empty((h,w)), 
    im_radius_grid   = compute_im_rad_grid(res_map, center=center, xc=xc, yc=yc, display=display)
    im_nan[:] = np.nan
    rad       = alpha

    # Loop until the radius of the annulus considered is bigger than the size of the image divided by two
    while rad < w//2 * 1.4: #w//2 - alpha - dr :
        # One annulus is considered
        rad += dr
        cond_annulus_large = np.logical_and(im_radius_grid >= rad-alpha, rad + dr + alpha >= im_radius_grid)
        cond_annulus_thin  = np.logical_and(im_radius_grid >= rad, rad + dr >= im_radius_grid)
        im_annulus = np.where(cond_annulus_large, res_map, im_nan)
        # the noise over the annulus is computed
        noise_annulus = np.nanstd(im_annulus)
        # and the image is set at this noise for this given annulus
        im_noise[cond_annulus_thin] = noise_annulus
        rad += dr
    if display: print("Took {} seconds".format(time.time()-t0))

    return im_noise
