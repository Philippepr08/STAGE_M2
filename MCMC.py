# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 13:40:37 2023

@author: Philippe Priolet
"""
import numpy as np
import scipy
import models_cosi as mdl
import emcee

from scipy.optimize import minimize


class MCMC:
    def __init__(self,observations,obs_err,B_u,B_v,lamb,ulamb,angsize_star_mas):
        
        self.obs = observations
        self.obs_err = obs_err
        self.model_obj = mdl.Models(B_u, B_v, lamb)
        self.B_u = B_u
        self.B_v = B_v
        self.waves = lamb
        
        self.ulamb = ulamb
        self.angsize_star_mas = angsize_star_mas
        self.stellar_model = self.model_obj.analyt_vis_limb_dark(ulamb,angsize_star_mas)
        
        self.cosi_tilt_max = 0.4
        self.cosi_tilt_min = 0.1
        
        self.ain_min = 50
        self.ain_max = 200
        
        self.da_min = 20 
        self.da_max = 55
        
        self.flux_ratio_min = 0
        self.flux_ratio_max = 0.1
        
        self.pa_max = 45
        self.pa_min = 20
        print('The limits are:','pa',(self.pa_min, self.pa_max),'cosi_tilt:',(self.cosi_tilt_min,self.cosi_tilt_max),',','a_in:',(self.ain_min,self.ain_max),',','da:',(self.da_min,self.da_max),',','f_ratio:',(self.flux_ratio_min,self.flux_ratio_max))
    
    def log_prior(self,theta):
        pa,cosi_tilt, a_in,da,flux_ratio, f = theta
        if (cosi_tilt> self.cosi_tilt_min and cosi_tilt< self.cosi_tilt_max) and (a_in>self.ain_min and a_in<self.ain_max)  and (da>self.da_min and da<self.da_max) and (flux_ratio>self.flux_ratio_min and flux_ratio<self.flux_ratio_max) and (f>0 and f<1) and (self.pa_min<pa<self.pa_max):
            log_p = 0
        else:
            log_p = -np.inf
        
        return log_p
    
    def model(self,theta):
        pa,cosi_tilt, a_in,da,flux_ratio, f = theta

        vis_disk_model = self.model_obj.inclined_disk_ccf(cosi_tilt,a_in,a_in+da,pa=pa,plot_obj = False)
        
        vis2 = (flux_ratio*vis_disk_model+(1-flux_ratio)*self.stellar_model)**2
        self.model_vis2 = vis2
        return vis2
    
    def log_likelihood(self,theta):

        pa,cosi_tilt, a_in,da,flux_ratio, f = theta
        model_vis2 = self.model(theta)
        sn2 = (self.obs_err**2)+(f**2)*(self.model_vis2**2)
        LL = -(1/2)*np.sum(((self.obs-self.model_vis2)**2/sn2+np.log(2*np.pi*sn2)))
        return LL
    
    def get_posterior(self,theta):
        if not np.isfinite(self.log_prior(theta)):
            return -np.inf

        ln_post = self.log_likelihood(theta)+self.log_prior(theta)
        return ln_post
    
    def MCMC_run(self):
        nwalkers = 30
        ndim = 6
        pos_pa = np.random.uniform(self.pa_min,self.pa_max ,nwalkers)
        pos_cosi = np.random.uniform(self.cosi_tilt_min,self.cosi_tilt_max ,nwalkers)
        pos_ain = np.random.uniform(self.ain_min,self.ain_max,nwalkers)
        pos_da = np.random.uniform(self.da_min,self.da_max,nwalkers)
        pos_fluxrat = np.random.uniform(self.flux_ratio_min,self.flux_ratio_max,nwalkers)
        pos_f = np.random.uniform(0,1,nwalkers)
        
        pos = np.array([pos_pa,pos_cosi,pos_ain,pos_da,pos_fluxrat,pos_f]).T
        print(f'Starting MCMC algorithm with {nwalkers} walkers and {ndim} parameters')        
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.get_posterior)
        sampler.run_mcmc(pos, 10000, progress=True)
        
        return sampler