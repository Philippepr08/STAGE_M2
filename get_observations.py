# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 13:24:48 2023

@author: phili
"""
import numpy as np
import pyhdust.oifits as oifits

class Observations:
    def __init__(self, path):
        self.data_dic = {}
        self.path = path
        self.data = oifits.open(self.path)
        self.vis2_obj = self.data.vis2
        self.create_data_dic()
        
    def create_data_dic(self):


        self.data_dic['Target_name'] = self.data.target[0].target

        self.data_dic['date_vis2'] = []
        self.data_dic['date_CP'] = []
     
        self.data_dic['wave_vis'] = []
        self.data_dic['wave_CP'] = []

        self.data_dic['Vis2'] = []
        self.data_dic['Vis2_err'] = []
        self.data_dic['u'] = []
        self.data_dic['v'] = []
        self.data_dic['freqs'] = []
        self.data_dic['baselines'] = []
        self.data_dic['CP'] = []
        self.data_dic['CP_err'] = []
        
        for obs in self.data.vis2:
            for i in range(len(obs.wavelength.eff_wave)):
                self.data_dic['Vis2'].append(obs._vis2data[i])
                self.data_dic['Vis2_err'].append(obs._vis2err[i])
                self.data_dic['wave_vis'].append(obs.wavelength.eff_wave[i])
                self.data_dic['u'].append(obs.ucoord)
                self.data_dic['v'].append(obs.vcoord)
                self.data_dic['baselines'].append(np.sqrt(obs.ucoord**2+obs.vcoord**2))
                self.data_dic['freqs'].append(np.sqrt(obs.ucoord**2+obs.vcoord**2)/obs.wavelength.eff_wave[i])
                self.data_dic['date_vis2'].append(obs.timeobs)

        for obs_cp in self.data.t3:
            for i in range(len(obs_cp.wavelength.eff_wave)):
                self.data_dic['CP'].append(obs_cp._t3phi[i])
                self.data_dic['CP_err'].append(obs_cp._t3phierr[i])
                self.data_dic['wave_CP'].append(obs_cp.wavelength.eff_wave[i])
                self.data_dic['date_CP'].append(obs_cp.timeobs)
                



        
        return self.data_dic
        
    def get_data_dic_for_wave(self,wavelength):
        idx_vis = np.where(self.data_dic['wave_vis'] == wavelength)[0]
        idx_CP = np.where(self.data_dic['wave_CP'] == wavelength)[0]


        data_dic_wave = {}
        data_dic_wave['Target_name'] = self.data_dic['Target_name']

        data_dic_wave['wave_vis'] = np.take(self.data_dic['wave_vis'], idx_vis).T
        data_dic_wave['wave_CP'] = np.take(self.data_dic['wave_CP'], idx_CP).T 
        
        data_dic_wave['date_vis2'] = np.take(self.data_dic['date_vis2'], idx_vis).T
        data_dic_wave['date_CP'] = np.take(self.data_dic['date_CP'], idx_CP).T 

        data_dic_wave['Vis2'] =  np.take(self.data_dic['Vis2'], idx_vis).T  
        data_dic_wave['Vis2_err'] = np.take(self.data_dic['Vis2_err'], idx_vis).T
        data_dic_wave['u'] = np.take(self.data_dic['u'], idx_vis).T
        data_dic_wave['v'] = np.take(self.data_dic['v'], idx_vis).T
        data_dic_wave['freqs'] = np.take(self.data_dic['freqs'], idx_vis).T 
        data_dic_wave['baselines'] = np.take(self.data_dic['baselines'], idx_vis).T
        data_dic_wave['CP'] = np.take(self.data_dic['CP'], idx_CP).T
        data_dic_wave['CP_err'] = np.take(self.data_dic['CP_err'], idx_CP).T 
        
        return data_dic_wave
    def get_possible_waves(self):
        return np.unique(self.data_dic['wave_vis'])
    
    def get_dates(self):
        return np.unique(self.data_dic['date_vis2'])
        
    def get_data_by_night(self,date_time):
     
        
        idx_vis = []
        idx_CP = []
        for i in range(len(self.data_dic['date_vis2'])):
            if self.data_dic['date_vis2'][i] == date_time:
                idx_vis.append(i)
                
        for j in range(len(self.data_dic['date_CP'])):
            if self.data_dic['date_CP'][j] == date_time:
                idx_CP.append(j)
        
        

        data_dic_night = {}
        data_dic_night['Target_name'] = self.data_dic['Target_name']
        data_dic_night['wave_vis'] = np.take(self.data_dic['wave_vis'], idx_vis).T
        data_dic_night['wave_CP'] = np.take(self.data_dic['wave_CP'], idx_CP).T 
        data_dic_night['date_vis2'] = np.take(self.data_dic['date_vis2'], idx_vis).T
        data_dic_night['date_CP'] = np.take(self.data_dic['date_CP'], idx_CP).T 
        
        data_dic_night['Vis2'] =  np.take(self.data_dic['Vis2'], idx_vis).T  
        data_dic_night['Vis2_err'] = np.take(self.data_dic['Vis2_err'], idx_vis).T
        data_dic_night['u'] = np.take(self.data_dic['u'], idx_vis).T
        data_dic_night['v'] = np.take(self.data_dic['v'], idx_vis).T
        data_dic_night['freqs'] = np.take(self.data_dic['freqs'], idx_vis).T 
        data_dic_night['baselines'] = np.take(self.data_dic['baselines'], idx_vis).T
        data_dic_night['CP'] = np.take(self.data_dic['CP'], idx_CP).T
        data_dic_night['CP_err'] = np.take(self.data_dic['CP_err'], idx_CP).T 
        
        return data_dic_night
    
if __name__ == '__main__':
    path_obs = 'Data/Latest/Beta_pic/OiXP_beta_pictoris_PIONIER_Pnat_1_5206438_1_7616950__A0-B2-C1-D0_2022-10-14.fits' 
    obs = Observations(path_obs)
    data_dic = obs.data_dic
    night_idx = 1
    nights = obs.get_dates()
    night14 = nights[0]
    dic_night = obs.get_data_by_night(nights[night_idx])
    
