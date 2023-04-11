# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 15:01:19 2023

@author: phili
"""

import plotting_observations as plt_vis
import get_observations as obs


data_paths = {'2011': 'Data/Ertel2016/Beta_pic/OiXP_BETA_PIC_PIONIER_Pnat_1_5549339_1_8129009__2_D0-G1-H0-I1_2011-11-03.fits',
              '2012': 'Data/Ertel2016/Beta_pic/OiXP_HD39060_PIONIER_Pnat_1_5900000_1_7679999__10_A1-B2-C1-D0_2012-10-17.fits',
              '2014': 'Data/Ertel2016/Beta_pic/OiXP_HD39060_PIONIER_Pnat_1_5900000_1_7679999__22_A1-B2-C1-D0_2014-10-12.fits',
              '2022': 'Data/Latest/Beta_pic/OiXP_beta_pictoris_PIONIER_Pnat_1_5206438_1_7616950__A0-B2-C1-D0_2022-10-14.fits',
              }

for year in data_paths: 
    path =   data_paths[year]         
    obs_obj = obs.Observations(path) 
    
    plot_obj = plt_vis.Plotting_obs(obs_obj)   
    plot_obj.plot_all(units='f',save=True, save_path = '')