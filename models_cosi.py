# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 11:31:45 2023

@author: phili
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
import oimodeler as oim

class Models:
    def __init__(self,B_u,B_v,waves):
        self.B_u = B_u
        self.B_v = B_v
        
        self.B = np.sqrt(B_u**2+B_v**2)
        self.waves = waves
        
        self.u = self.B_u/waves
        self.v = self.B_v/waves
        
    def analyt_vis_limb_dark(self,ulamb,angsize_star_mas): #angsize_star = 2*R_star_mas
        angsize_star_rad = np.radians(angsize_star_mas/(60*60*1000))
        theta_ld = np.sqrt((1-ulamb/3)*(1-7*ulamb/15))*angsize_star_rad
        
        x = np.pi*theta_ld*self.B/self.waves
        
        bessel_one = scipy.special.jv(1,x)
        bessel_3_2 = scipy.special.jv(3/2,x)

        alpha = 1 -ulamb
        
        vis2 = ((alpha*bessel_one/x+ulamb*np.sqrt(np.pi/2)*bessel_3_2/x**(3/2))**2)/(alpha/2+ulamb/3)**2
        vis = np.sqrt(vis2)
        
        return vis
    
    def inclined_disk_ccf(self,cosi_tilt,a_in,a_out,pa,plot_obj = False):
        wave = self.waves
        x = self.B_u/wave
        y = self.B_v/wave
        elong = cosi_tilt

        din = a_in*cosi_tilt
        dout = a_out*cosi_tilt

        ring = oim.oimERing(f=1,din=din,dout=dout,pa=pa,elong=elong)


        mring = oim.oimModel(ring)
        # gauss  = oim.oimGauss(fwhm=np.abs(a_out-a_in),f=1)
        # convol = oim.oimConvolutor(ring,gauss)
        # mring_pls_gauss = oim.oimModel(convol)

        if plot_obj:
            figImg=mring.showModel(512,0.1,normPow=0.2)


        ccf = mring.getComplexCoherentFlux(x,y)
        ccf = ccf/np.max(np.abs(ccf))
        
        vis = np.abs(ccf)
        
        
        return vis
    
     
    def disk_plus_star(self,cosi_tilt,a_in,a_out,flux_ratio,ulamb,angsize_star_mas):
        vis_disk = self.inclined_disk_ccf(cosi_tilt,a_in,a_out,plot_obj = False)
        vis_star = self.analyt_vis_limb_dark(ulamb,angsize_star_mas)
        vis2 = (flux_ratio*vis_disk+(1-flux_ratio)*vis_star)**2
        
        return vis2

if __name__=='__main__':
    
    Nbas = 700    
    B_u = np.linspace(-30,30,Nbas)
    B_v = np.linspace(-40,40,Nbas)
    
    x,y = np.meshgrid(B_u,B_v)
    waves = 1e-6

    mod = Models(B_u,B_v,waves)
    i_tilt_deg = 89.5
    # elong = np.sin(i_tilt_rad)
    # print(elong)
    a_in = 20
    a_out = 30
    pa = 0
    cos_i = np.cos(np.radians(i_tilt_deg))
    disk_vis = mod.inclined_disk_ccf(cos_i,a_in,a_out,pa,plot_obj = True)

    # ring = oim.oimERing(f=1,din=din,dout=dout,pa=0,elong=elong)
    # mring = oim.oimModel(ring)
    # params = mring.getParameters()
    # figImg=mring.showModel(1500,0.16,normPow=0.2)
    
    # # im=mring.getImage(512,0.1)
    
    # ccf = mring.getComplexCoherentFlux(x,y)
    
    
    # plt.imshow(im)
