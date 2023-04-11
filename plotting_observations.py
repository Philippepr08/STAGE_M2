# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:50:50 2023

@author: phili
"""
import numpy as np
import matplotlib.pyplot as plt

class Plotting_obs:
    def __init__(self,obs_obj):
        #INPUT:
            #obs_obj is created with the class Observations in the file get_observations.py
        self.obs_obj = obs_obj
        self.obs_dic = obs_obj.data_dic
        self.nights = obs_obj.get_dates()
        
    def plot_vis(self,axis,units,night = None):
        formatted_night = f'{night.day}/{night.month}/{night.year}'
        night_dic = self.obs_obj.get_data_by_night(night)
        
        baselines = night_dic['baselines']
        vis2 = night_dic['Vis2']
        vis2_err = night_dic['Vis2_err']
        
        if units == 'm':
            im1 = axis.scatter(baselines,vis2, c= night_dic['wave_vis'],cmap = 'Spectral_r')
            axis.set_ylim(0.85,1.15)
            axis.set_title(f'Vis2 for the night of the {formatted_night}')
            axis.errorbar(baselines,vis2, yerr=vis2_err,ecolor='darkgrey', linestyle='', marker=None, mew=0)#,**error_kwargs )
            axis.set_xlabel(r'B (m)')
            axis.set_ylabel(r'Vis2')
        else:
            freqs = baselines/np.array(night_dic['wave_vis'])
            im1 = axis.scatter(freqs,vis2, c= night_dic['wave_vis'],cmap = 'Spectral_r')
            axis.set_ylim(0.85,1.15)
            axis.set_title(f'Vis2 for the night of the {formatted_night}')
            axis.errorbar(freqs,vis2, yerr=vis2_err,ecolor='darkgrey', linestyle='', marker=None, mew=0)#,**error_kwargs )
            axis.set_xlabel(r'f ($\lambda . rad^{-1})$')
            axis.set_ylabel(r'Vis2')


        return im1
        
    def plot_uv_plane(self,axis,units,night = None):
        formatted_night = f'{night.day}/{night.month}/{night.year}'
        night_dic = self.obs_obj.get_data_by_night(night)
        
        B_u, B_v = night_dic['u'],night_dic['v']
        f_u, f_v = B_u/np.array(night_dic['wave_vis']), B_v/np.array(night_dic['wave_vis'])
        
        if units == 'm':
            im1 = axis.scatter(B_u,B_v, c= night_dic['wave_vis'],cmap = 'Spectral_r')
            im2 = axis.scatter(-B_u,-B_v, c= night_dic['wave_vis'],cmap = 'Spectral_r')
            axis.set_title(f'uv plane for the night of {formatted_night}')
            axis.set_xlabel(r'u (m)')
            axis.set_ylabel(r'v (m)')


            # axis.set_ylim(0,1.5)
        else:
            im1 = axis.scatter(f_u,f_v, c= night_dic['wave_vis'],cmap = 'Spectral_r')
            im2 = axis.scatter(-f_u,-f_v, c= night_dic['wave_vis'],cmap = 'Spectral_r')
            axis.set_title(f'uv plane for the night of {formatted_night}')
            axis.set_xlabel(r'u ($\lambda . rad^{-1})$')
            axis.set_ylabel(r'v ($\lambda . rad^{-1})$')
        
        return im1, im2
    def plot_all(self,units='m',save=True, save_path = ''):
        ## Units: can be meters (m) or spatial frequencies (f)
        
        fig, axis = plt.subplots(len(self.nights),2,gridspec_kw={'width_ratios': [1, 2]})
        fig.set_size_inches(22,6*len(self.nights))
        fig.tight_layout(pad=6)
        im = []
        if len(self.nights)==1:
            for i, night in enumerate(self.nights):
                im1 = self.plot_vis(axis[1],units,night)
                im2,im3 = self.plot_uv_plane(axis[0],units,night)
                im.append(im1)
                im.append(im2)

        else:
            for i, night in enumerate(self.nights):
                im1 = self.plot_vis(axis[i,1],units,night)
                im2,im3 = self.plot_uv_plane(axis[i,0],units,night)
                im.append(im1)
                im.append(im2)
    
            for ax in axis[:,1]:
                fig.colorbar(im[i], ax=ax)

        last_obs_date = self.nights[-1]
        fig.suptitle(f"Observations of {self.obs_dic['Target_name']}", fontsize=25)
        if save:
            plt.savefig(save_path + f"{self.obs_dic['Target_name']}_observations_{last_obs_date.month}_{last_obs_date.year}.png")
