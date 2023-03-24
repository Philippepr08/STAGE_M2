# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:32:30 2023

@author: Philippe Priolet
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

import get_observations as obs

import Plotting_vis as pltv
import Visibilities_new as vis


import reduced_chi_squared_fit as chi2
##############################|| SETTING PARAMETERS ||##############################

####| Observation parameters |####
data_year = '2012'
get_obs_plot = False
plot_units = 'f' # The units for the plots, either Baselines 'm' or frequencies 'f'


####| Simulation parameters |####
use_obs_baselines = True
Nimg = 3000
pix_scale_mas = 0.01 # The field of view of a single pixel in mas.px^-1

if not use_obs_baselines:
    # We create a linear distribution (with N_bas_pts points) of baselines going from B_min to B_max oriented at an angle angle_u_v_plane with respect to the u axis 
    angle_u_v_plane = 45  
    B_min = 10
    B_max = 40
    N_bas_pts = 200
    
####| Model parameters |####
angsize_star_mas = 0.707 #The angular diameter of the star in mas
ulamb = 0.27 #The limb darkening coefficient

##############################|| GETTING OBSERVATIONS ||##############################


data_paths = {'2011':'Data/Ertel2016/Beta_pic/OiXP_BETA_PIC_PIONIER_Pnat_1_5549339_1_8129009__2_D0-G1-H0-I1_2011-11-03.fits',
              '2012':'Data/Ertel2016/Beta_pic/OiXP_HD39060_PIONIER_Pnat_1_5900000_1_7679999__10_A1-B2-C1-D0_2012-10-17.fits',              
              '2014':'Data/Ertel2016/Beta_pic/OiXP_HD39060_PIONIER_Pnat_1_5900000_1_7679999__22_A1-B2-C1-D0_2014-10-12.fits',
              '2022':'Data/Latest/Beta_pic/OiXP_beta_pictoris_PIONIER_Pnat_1_5206438_1_7616950__A0-B2-C1-D0_2022-10-14.fits'
    }
if data_year in data_paths: #Checking if input date is in the data base
    obs_obj = obs.Observations(data_paths[str(data_year)])
else:
    print('No data for that year!')
    sys.exit()

## Getting the nights of observations (Some data have multiple nights)    
nights = obs_obj.get_dates()


## Getting the visibilities by night (first index is the first night and so on)    
vis2_obs = [] #The visibility from the observations
vis2_err_obs = [] #The associated errors

dic_by_night = []

for i, night in enumerate(nights):
    obs_dic_by_night = obs_obj.get_data_by_night(night)
    dic_by_night.append(obs_dic_by_night)

    
    vis2_obs.append(obs_dic_by_night['Vis2'])
    vis2_err_obs.append(obs_dic_by_night['Vis2_err'])

vis2_obs = np.array(vis2_obs)
vis2_err_obs = np.array(vis2_err_obs)


    
if get_obs_plot: #If user wants the plotted observations of the visibility
    try:
        os.makedirs("Observation_plots/")
    except FileExistsError:
        pass
    pltv.Plotting_obs(obs_obj).plot_all(units=plot_units,save=True, save_path = "Observation_plots/")
    
##############################|| GETTING BASELINES FOR THE VISIBILITY ||##############################
B_u = []
B_v = []

if use_obs_baselines:
    for i,night in enumerate(nights):
        B_u.append(dic_by_night[i]['u'])
        B_v.append(dic_by_night[i]['v'])
else:
    B = np.linspace(B_min,B_max,N_bas_pts)
    for i,night in enumerate(nights):
        B_u.append(np.cos(angle_u_v_plane)*B)
        B_v.append(np.sin(angle_u_v_plane)*B)

B_u = np.array(B_u)
B_v = np.array(B_v)

B = np.sqrt(B_u**2+B_v**2)    







##############################|| GETTING MODEL ||##############################

model_star = []

### Getting the visibilities ###
print('Getting the stellar visibilities')

for i, night in enumerate(nights):
    waves = dic_by_night[i]['wave_vis']

    vis_obj = vis.Visibilities(B_u[i], B_v[i], waves)

    vis2_stellar = vis_obj.analyt_vis_limb_dark(ulamb,angsize_star_mas)

    vis_dic_model = {'Vis2':[],'Waves':[]}
    print('Night number',i+1)
    vis_arr = []
    waves_arr = []
    # for idx, wave in enumerate(waves):
    #     print('Getting visibilities for wavelength:', round(wave*1e6,4),' µm','B:',B_u[i][idx],B_v[i][idx])
    #     vis = vis_model_obj.evaluate_vis(wave,B_u[i][idx],B_v[i][idx],plot = False) 
    #     vis_arr.append(vis**2)
    #     waves_arr.append(wave)
        
    vis_dic_model['Waves'] = waves
    vis_dic_model['Vis2'] = vis2_stellar
    # vis_dic_model['Vis2'] = np.concatenate(vis_dic_model['Vis2'])
    # vis_dic_model['Waves'] = np.concatenate(vis_dic_model['Waves'])

    # vis_dic_model['Night'] = night
    # vis_dic_model['Name'] = f'Stellar model, R_star:{R_star_mas}'
    # vis_dic_model['B_u'] = B_u[i]
    # vis_dic_model['B_v'] = B_v[i]
    # vis_dic_model['B'] = np.sqrt(B_v[i]**2+B_u[i]**2)


    # vis_dic_model['Vis2'] = np.concatenate(vis_dic_model['Vis2'])
    # vis_dic_model['Waves'] = np.concatenate(vis_dic_model['Waves'])

    model_star.append(vis_dic_model)
#%%    
##############################|| FITTING ||##############################
###|| Linear fit ||###
flux_ratios_linear = []
chi_square_linear = []
k_arr_linear = []
for i, night in enumerate(nights): #We do a fit for each night
    
    vis2_obs = dic_by_night[i]['Vis2']
    vis2_star = model_star[i]['Vis2']
    vis2_obs_error = dic_by_night[i]['Vis2_err']
    
    f, chi = chi2.chi2(vis2_obs,vis2_star,vis2_obs_error,order = 1)
    k = 1-2*f
    k_arr_linear.append(k) 
    flux_ratios_linear.append(f)
    chi_square_linear.append(chi)
    
###|| Quadratic fit ||###
flux_ratios_quadratic = []
chi_square_quadratic = []
k_arr_quadratic = []
for i, night in enumerate(nights): #We do a fit for each night

    vis2_obs = dic_by_night[i]['Vis2']
    vis2_star = model_star[i]['Vis2']
    vis2_obs_error = dic_by_night[i]['Vis2_err']

    
    f, chi = chi2.chi2(vis2_obs,vis2_star,vis2_obs_error,order = 2)
    k = 1-f
    k_arr_quadratic.append(k) 
    flux_ratios_quadratic.append(f)
    chi_square_quadratic.append(chi)
    
    
#%%   

index = 0
vis2_obs = dic_by_night[index]['Vis2']
vis2_star = model_star[index]['Vis2']
vis2_obs_error = dic_by_night[index]['Vis2_err']

plt.errorbar(vis2_star,vis2_obs,yerr=vis2_obs_error,fmt='o')
x = np.linspace(np.min(vis2_star),np.max(vis2_star),100)
plt.plot(x,k_arr_quadratic[index]**2*x,c='red',label=f'Flux_ratio:{np.round(flux_ratios_quadratic[index]*100,4)}%')
plt.legend()

#%%
plt.scatter(B[index],vis2_star,c=dic_by_night[index]['wave_vis'],cmap = 'Spectral_r')
    