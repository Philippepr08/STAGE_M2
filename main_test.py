# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:57:06 2023

@author: phili
"""
import numpy as np
import matplotlib.pyplot as plt
import Visibilities_new as vis
import oimodeler as oim

B_resolution = 0.1
angsize_star_mas = 0.8
ulamb = 0.27

R_int_mas = 0
R_ext_mas = 800

i_tilt_deg = 90

flux_ratio = 5e-2

B_u = np.linspace(1,100,200)
B_v = np.zeros(B_u.shape)
B = np.sqrt(B_u**2+B_v**2)

waves = 1e-6


vis_obj = vis.Visibilities(B_u, B_v, waves)

vis2_stellar = vis_obj.analyt_vis_limb_dark(ulamb,angsize_star_mas)
vis2 = vis_obj.disk_plus_star_from_img(B_resolution,angsize_star_mas,ulamb,R_int_mas,R_ext_mas,i_tilt_deg,flux_ratio,redo_FT=True)

plt.plot(B,vis2)
plt.plot(B,vis2_stellar)


f = 1-np.sqrt(vis2/vis2_stellar)
