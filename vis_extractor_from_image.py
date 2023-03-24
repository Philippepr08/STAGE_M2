# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:33:03 2023

@author: phili
"""
import numpy as np
import matplotlib.pyplot as plt
   
class Visibilities:
    def __init__(self,img,pix_scale_mas):
        ## INPUT:
            #lamb: The observation wavelength in meters
            #img: The intensity distribution map of the object (A 2D array)
            #FOV_mas: The field of view in mas 
            #pix_scale: The amount of sky coverage in mas that corresponds to a single pixel, in mas/px
        self.img = img
        self.N = img.shape[0]
        self.pix_scale_mas = pix_scale_mas
        self.pix_scale = np.radians(pix_scale_mas/(1000*60*60))
    def get_complex_tf(self,plot_bool = False):
        tf_img = np.fft.fftshift(np.fft.fft2(self.img))
        tf_img = tf_img/np.max(np.abs(tf_img)) #Ensures the visibility is 1 for a non-resolved object
        self.tf_img = tf_img
        self.N_tf = tf_img.shape[0]
        
        self.FOV_mas = self.N*self.pix_scale_mas
        self.FOV = self.N*self.pix_scale
        if plot_bool:
            vis = np.abs(self.tf_img)
            vis = vis/np.max(vis)
            plt.imshow(vis**0.1)
            plt.scatter([self.N//2],[self.N//2])
        
        return tf_img
    
    
    def evaluate_vis_circ_sym(self,wave,B):
        f_u = B/(wave)
        f_v = B*0

        # idx_u = np.int64(f_u*self.pix_scale*self.N)  + self.N_tf//2
        # idx_v = np.int64(f_v*self.pix_scale*self.N)  + self.N_tf//2
        
        idx_u = np.int64(f_u*self.pix_scale*self.N_tf)  + self.N_tf//2
        idx_v = np.int64(f_v*self.pix_scale*self.N_tf)  + self.N_tf//2
        
        #Evaluate the fourier transform at these points
        vis = np.abs(self.tf_img[idx_v,idx_u])
        return vis
    
    
    def evaluate_vis(self,Nimg,wave,B_u,B_v,plot = False):
        f_u = B_u/(wave)  ###### WAIT A SECOND, THE PIX_SCALE HAS TO BE IN RADIANS
        f_v = B_v/(wave)
        
        idx_u = np.int64(f_u*self.pix_scale*Nimg)  + Nimg//2
        idx_v = np.int64(f_v*self.pix_scale*Nimg)  + Nimg//2
        print('INDEX',idx_u)
        print('TF SIZE',self.tf_img.shape,'pixscale',self.pix_scale_mas,Nimg)
        #Evaluate the fourier transform at these points
        vis = np.abs(self.tf_img[idx_v.astype(int),idx_u.astype(int)])
        print('The synthetic FOV is:',self.FOV_mas,'mas','|','The pix_scale is:',self.pix_scale_mas,'mas.px^-1')
        if np.any(vis)== 'nan':
            print('NAN IN VIS')
        return vis
