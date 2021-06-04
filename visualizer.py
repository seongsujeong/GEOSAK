#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import grwt
import sys

#To plot the real numbered array
def plot_real(raster_in,range=None,cmap='viridis'):
    if type(raster_in)==grwt.raster:
        arr_to_plot=raster_in.z
        if raster_in.boundary==None:
            vec_extent=None
        else:
            vec_boundary=raster_in.boundary #xmin, ymin, xmax, ymax
            vec_extent=[vec_boundary[0], vec_boundary[2], vec_boundary[3], vec_boundary[1]]
    elif type(raster_in)==np.ndarray:
        arr_to_plot=raster_in
        vec_extent=None

    #determine the range
    if range==None:
        vec_clim=[np.nanpercentile(arr_to_plot,3),np.nanpercentile(arr_to_plot,97)]
    else:
        vec_clim=range

    if vec_extent==None:
        fig=plt.imshow(arr_to_plot,clim=vec_clim,cmap=cmap)

    else: #provide vec_extent
        fig=plt.imshow(arr_to_plot,clim=vec_clim,cmap=cmap,extent=vec_extent)
        
    plt.show()



        


#To plot the complex numbered data e.g. interferogram
def plot_complex(raster_in,range=None,cmap='hsv'):
    if type(raster_in)==grwt.raster:
        arr_to_plot=raster_in.z
        if raster_in.boundary==None:
            vec_extent=None
        else:
            vec_boundary=raster_in.boundary #xmin, ymin, xmax, ymax
            vec_extent=[vec_boundary[0], vec_boundary[2], vec_boundary[3], vec_boundary[1]]
    elif type(raster_in)==np.ndarray:
        arr_to_plot=raster_in
        vec_extent=None
    
    nx_arr=arr_to_plot.shape[1]
    ny_arr=arr_to_plot.shape[0]
    arr_hsv=np.zeros((ny_arr,nx_arr,3))

    if range==None:
        magmin=np.nanpercentile(arr_to_plot,3)
        magmax=np.nanpercentile(arr_to_plot,97)
    else:
        magmin=range[0]
        magmax=range[1]
    
    magn=(arr_to_plot-magmin)/(magmax-magmin)
    magn[magn>1.0]=1.0
    magn[magn<0.0]=0.0

    arr_hsv[:,:,0]=(np.arctan2(arr_to_plot.real,arr_to_plot.imag)+np.pi)/(2*np.pi)
    arr_hsv[:,:,1]=magn
    arr_hsv[:,:,2]=magn
    
    arr_rgb=matplotlib.colors.hsv_to_rgb(arr_hsv)
    
    if vec_extent==None:
        fig=plt.imshow(arr_rgb)
    else:
        fig=plt.imshow(arr_rgb,extent=vec_extent)

    plt.show()





if __name__=='__main__':
    str_usage='''

    '''
