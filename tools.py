# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 14:17:26 2023

@author: phili
"""
import numpy as  np
import matplotlib.pyplot as plt
def get_tf_padding(B_resolution,lamb,N_img,pix_scale_mas):
    FOV_tf = lamb/B_resolution
    pix_scale_rad = np.radians((pix_scale_mas)/(60*60*1000))
    N_tf = FOV_tf/pix_scale_rad
    padding  = (N_tf - N_img)//2 # might be wise to round up instead of //
    if padding<0:
        padding = 0
        
    return int(padding+1) #Should make it upper value, to have at least B_resolution

def get_image_size(B_resolution,lamb,pix_scale_mas):
    FOV_tf = lamb/B_resolution
    pix_scale_rad = np.radians((pix_scale_mas)/(60*60*1000))
    N_img = FOV_tf/pix_scale_rad

def get_pix_scale_mas(B_max,B_resolution,lamb):
    N_img = 2*np.ceil(B_max/B_resolution) + 1
    
    FOV_tf = lamb/B_resolution
    FOV_tf_mas = np.degrees(FOV_tf)*60*60*1000
    print('The field of view of the simulated image is:',round(FOV_tf_mas,2),'mas')
    
    pix_scale_rad = FOV_tf/N_img
    pix_scale_mas = np.degrees(pix_scale_rad)*60*60*1000
    print('The single pixel FOV of the simulated image is:',round(pix_scale_mas,2),'mas')
    print('The image size is:',N_img)
    return N_img, pix_scale_mas,FOV_tf_mas


def get_gaussian_FOV_mask(Nimg,FOV_mas,pix_scale_mas):
    FWHM_px = FOV_mas/pix_scale_mas
    sigma = FWHM_px/(2*np.sqrt(2*np.log(2)))
    
    coord = np.arange(-Nimg//2,+Nimg//2,1)
    x, y = np.meshgrid(coord,coord)
    
    u= x**2/(2*sigma**2) + y**2/(2*sigma**2)
    
    gaussian = np.exp(-u)
    
    return gaussian
#%%
# Nimg = 3000
# pix_scale_mas,FOV_tf = get_pix_scale_mas(0.1, 1e-6, Nimg)
# FOV_mas = FOV_tf
# image = np.ones(shape= (Nimg,Nimg))*get_gaussian_FOV_mask(FOV_mas,pix_scale_mas)
# tf_image = np.fft.fftshift(np.fft.fft2(image))
# vis = np.abs(tf_image)/np.max(np.abs(tf_image))
# plt.imshow(image)
# plt.colorbar()
# #%%
# zoom = 30
# zoom = zoom//2
# plt.plot(vis[Nimg//2,Nimg//2-zoom:Nimg//2+zoom])