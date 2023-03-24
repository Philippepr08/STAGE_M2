# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 18:00:11 2023

@author: phili
"""
import numpy as np

def chi2(vis2_obs,vis2_star,vis2_obs_error,order=1): ## See Absil et al 2021
    #For linear relation of type Vis2 = k*Vis2_star
    nb_obs = vis2_obs.shape[0]
    print('Number of observations:',nb_obs)

    vis2_obs = np.array(vis2_obs)
    vis2_star = np.array(vis2_star)
    
    sum_ratio = np.sum((vis2_obs*vis2_star)/(vis2_obs_error**2))
    sum_star = np.sum(vis2_star**2/(vis2_obs_error**2))
    if order == 1:
        k = sum_ratio/sum_star
        f = (1-k)/2
        chi_2_temp = ((vis2_obs-k*vis2_star)**2)/(vis2_obs_error**2) 

    elif order == 2:
        k = np.sqrt(sum_ratio/sum_star)
        f = (1-k)
        chi_2_temp = ((vis2_obs-k**2*vis2_star)**2)/(vis2_obs_error**2) 

    
    #Calculate the chi_square for this value of k
    chi_2 = np.sum(chi_2_temp)
    print(f'The fitted flux ratio is (order={order}):',f*100, '%')
    chi2_reduced = chi_2/(nb_obs-1)
    print(f'The reduced chi squared value is:',chi2_reduced)
    return f, chi2_reduced


