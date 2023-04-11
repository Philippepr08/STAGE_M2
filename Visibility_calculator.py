# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 22:44:48 2023

@author: Philippe
"""

import numpy as np
import scipy

import vis_extractor_from_image as visext
import tools

class Visibilities:
    def __init__(self,B_u,B_v,waves):
        self.B_u = B_u
        self.B_v = B_v

        self.B = np.sqrt(B_u**2+B_v**2)
        self.waves = waves
    def analyt_vis_limb_dark(self,ulamb,angsize_star_mas): #angsize_star = 2*R_star_mas
        angsize_star_rad = np.radians(angsize_star_mas/(60*60*1000))
        theta_ld = np.sqrt((1-ulamb/3)*(1-7*ulamb/15))*angsize_star_rad
        
        x = np.pi*theta_ld*self.B/self.waves
        
        bessel_one = scipy.special.jv(1,x)
        bessel_3_2 = scipy.special.jv(3/2,x)

        alpha = 1 -ulamb
        
        vis2 = ((alpha*bessel_one/x+ulamb*np.sqrt(np.pi/2)*bessel_3_2/x**(3/2))**2)/(alpha/2+ulamb/3)**2
        return vis2
    
    def get_analyt_vis_disk_plus_star(self,vis2_star,flux_ratio,order=1):
        if order == 1:
            v2_image = (1-2*flux_ratio)*vis2_star
        if order == 2:
            v2_image = vis2_star*(1-flux_ratio)**2
         
        return v2_image
    
    def disk_not_unresolved(self,vis_star,vis_disk,flux_ratio): ##see equation 3 DiFolco et al 2007
        #Be careful, inputs are vis not vis2
        vis2 = (flux_ratio*vis_disk+(1-flux_ratio)*vis_star)**2
        return vis2

    def get_gaussian_FOV_mask(self,FOV_mas,pix_scale_mas):
        FWHM_px = FOV_mas/pix_scale_mas
        sigma = FWHM_px/(2*np.sqrt(2*np.log(2)))
        
        coord = np.arange(-self.Nimg//2,+self.Nimg//2,1)
        x, y = np.meshgrid(coord,coord)
        
        u= x**2/(2*sigma**2) + y**2/(2*sigma**2)
        
        gaussian = np.exp(-u)
        
        return gaussian

    def disk_model(self,pix_scale_mas,R_int_mas,R_ext_mas,i_tilt_deg, plot = False,redo_FT = True):
        i_tilt_rad = np.radians(i_tilt_deg)
        R_int_px = R_int_mas/pix_scale_mas
        R_ext_px = R_ext_mas/pix_scale_mas
        
        i = np.arange(-self.Nimg//2,self.Nimg//2,1)
        j = np.arange(-self.Nimg//2,self.Nimg//2,1)
        x, y = np.meshgrid(i,j)
        
        z_1 = np.sqrt(x**2+((y/np.sin(i_tilt_rad))**2))<R_ext_px 
        z_2 = np.sqrt(x**2+((y/np.sin(i_tilt_rad))**2))>R_int_px
        image = z_1*z_2
        # image = image*self.get_gaussian_FOV_mask(self.FOV_mas,pix_scale_mas)
        self.disk_image= image
        print('Vis pixscale',pix_scale_mas)
        vis_obj = visext.Visibilities(image, pix_scale_mas)
        if redo_FT:
            print('Getting Fourier Transform')
            tf_image = vis_obj.get_complex_tf(plot_bool = False)
            print('Done with Fourier Transform')
            
        vis = vis_obj.evaluate_vis(self.Nimg,self.waves,self.B_u,self.B_v,plot = False)
        return vis
    
    
    def disk_plus_star_from_img(self,B_resolution,angsize_star_mas,ulamb,R_int_mas,R_ext_mas,i_tilt_deg,flux_ratio,redo_FT):
        self.B_max = np.max(self.B)
        N_img, pix_scale_mas, FOV_mas = tools.get_pix_scale_mas(self.B_max,B_resolution, np.max(self.waves))
        self.Nimg = N_img
        self.FOV_mas = FOV_mas
        
        vis2_star = self.analyt_vis_limb_dark(ulamb,angsize_star_mas)
        vis_star = np.sqrt(vis2_star)
        vis_disk = self.disk_model(pix_scale_mas,R_int_mas,R_ext_mas,i_tilt_deg, plot = False,redo_FT=redo_FT)
        
        vis2 = self.disk_not_unresolved(vis_star,vis_disk,flux_ratio)
        return vis2
    




