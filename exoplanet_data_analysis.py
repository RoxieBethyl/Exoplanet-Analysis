# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:45:26 2022

@author: blybelle
"""

# import sys
import os
import numpy as np
# import astropy as aspy
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.stats import mad_std
import ccdproc as ccdp
from ccdproc import subtract_overscan
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.scale as scale
import plotly.express as px
from matplotlib.colors import Normalize 
# import pandas as pd
# from astropy import units as u
# from astropy import coordinates as coord
# from astropy.utils.data import get_pkg_data_filename
# from photutils.aperture import EllipticalAperture
#import convenience_functions as c_f
from pathlib import Path as path
import shutil
#%matplotlib inline

#def bias_master():
def fits_plot(image_path='', image_data='', cmap='afmhot', title=None, norm=None, vmin=None, vmax=None, xlim=(None, None), ylim=(None, None), details=False): #norm=colors.NoNorm(), vmin=2E3, vmax=3E3
    """
    Plots the image of the FITS file and some useful information.
    
    Parameters
    ----------
    image_path : String
        File path of the image in .FITS format to plotted
    
    colour_map : String
        Desired colour map entry for use in plot
        To see more color maps:
        https://matplotlib.org/stable/tutorials/colors/colormaps.html

    Returns
    -------
    Image_plot(plt.imshow) : 
        Image of the FITS file
        
    """
    if title is not None:
        plt.title(title)
    elif len(image_path) != 0:
        image_data = fits.getdata(image_path)
        plt.title(os.path.basename(image_path))
        
    # print(type(image_data)) # Show the Python type for image_data
    #plt.subplots(nrows=1, ncols=1)
    plt.imshow(image_data, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax) #, norm=norm, vmin=vmin, vmax=vmax, 
    plt.gca().invert_yaxis() #inverts y-axis into the right orientation
    plt.xlabel("x-pixels")
    plt.ylabel("y-pixels")
    plt.colorbar()
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()
    
    if details is True:
        print(image_data.shape) # Show the number of pixels per side in the 2-D image
        print('Min:', np.min(image_data))
        print('Max:', np.max(image_data))
        print('Mean:', np.mean(image_data))
        print('Stdev:', np.std(image_data), '\n')

def fits_hist(image_path, xlim=(None,None)):
    """

    Parameters
    ----------
    image_path : String
        File path of the image in .FITS format to plotted

    Returns
    -------
    Histogram_plot(plt.hist) : 
        Histogram of the FITS file data with bins set on auto
    
    """
    
    image_data = fits.getdata(image_path)
    print(image_data.flatten().shape)
    plt.hist(image_data.flatten(), bins='auto')
    plt.xlim(xlim)
    plt.show()

    
def get_header_data(data_dir, keywords, imagetyp='*', filter=''):
    """

    Parameters
    ----------
    data_dir : TYPE
        DESCRIPTION.
    keywords : TYPE
        DESCRIPTION.
    imagetyp : TYPE, optional
        DESCRIPTION. The default is '*'.
    filter : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    """
    files = ccdp.ImageFileCollection(data_dir)
    if filter!='' or imagetyp!='*':
        images = files.files_filtered(imagetyp=imagetyp, filter=filter, include_path=True)
    else:
        images = files.files_filtered(include_path=True)
    try: 
        print(images.summary[keywords])
    except:
        print(images.summary['file', 'imagetyp', 'filter', 'exptime', 'naxis1', 'naxis2'])
    #files.summary['imagetyp', 'biassec', 'ccdsec', 'datasec'][0]

def image_count_plot(image_path, image_label=None):
    """
    
    Parameters
    ----------
    image_paths : Matrix List
        Matrix list of image locations and name of the image

    Returns
    -------
    Plot of the counts per pixel of the image 

    """
    plt.figure(figsize=(20,10))
    image = CCDData.read(image_path, unit='count')
    #plt.plot(science_g_lfc.data.mean(axis=0), label='Science image')
    plt.plot(image.data.mean(axis=0), label=image_label)
    #plt.plot(flat_g_lfc.data.mean(axis=0), label='Flat image')
    
    plt.grid()
    #plt.axvline(x=2048, color='black', linewidth=3, linestyle='dashed', label='start of overscan')
    
    plt.legend()
    #plt.ylim(2040, 2090)
    #plt.xlim(0, 4096)
    plt.xlabel('pixel number')
    plt.ylabel('Counts')
    plt.title('Overscan region, averaged over all rows')

def copy_files(images_dir, base_dir, folder, imagetyp='*', filter=''):
    """
    Copies files into a seperate directory by filtering image type
    
    Parameters
    ----------
    image_dir : string
        Directory of images to filter
        
    imagetyp : string

        DESCRIPTION.
        
    base_dir : string
        DESCRIPTION.
        
    folder : string
        DESCRIPTION.

    Returns
    -------
    None

    """
        
    files = ccdp.ImageFileCollection(images_dir)
    calibrated_path = path(base_dir, folder)
    calibrated_path.mkdir(exist_ok=True)
    
    images = files.files_filtered(imagetyp=imagetyp, filter=filter, include_path=True)
    for image in images:
        shutil.copy(image, calibrated_path, overwrite=True)

def combine_images(base_dir, folder, filename, unit='adu', imagetyp='*', filter=''):
    """

    Parameters
    ----------
    base_dir : TYPE
        DESCRIPTION.
    folder : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.
    unit : TYPE, optional
        DESCRIPTION. The default is 'adu'.

    Returns
    -------
    None.

    """
    calibrated_path = path(base_dir, folder)
    reduced_images = ccdp.ImageFileCollection(calibrated_path, keywords=['imagetyp','filter'])
    calibrated_biases = reduced_images.files_filtered(imagetyp=imagetyp, filter=filter, include_path=True)

    combined_bias = ccdp.combine(calibrated_biases, method='average', sigma_clip=True, sigma_clip_low_thresh=5,
                                 sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, 
                                 mem_limit=350e6, unit=unit)

    combined_bias.meta['combined'] = True

    try: 
        combined_bias.write(calibrated_path / 'combined_bias.fits')
    except:
        combined_bias.write(calibrated_path / filename)
        
    print("Files combined and saved")


def main():
    # Setting the file locations for the data   
    base_dir = path("C:/Users/blysh/Documents/Telescope Group Project/")
    reduced_images_dir = path("C:/Users/blysh/Documents/Telescope Group Project/Reduced Images")
    raw_files = ccdp.ImageFileCollection(path("C:/Users/blysh/Documents/Telescope Group Project/2022-10-09 Calibration Data"))
    data_files = ccdp.ImageFileCollection(path("C:/Users/blysh/Documents/Telescope Group Project/2022-10-09 Explanets 1"))
    reduced_images = ccdp.ImageFileCollection("C:/Users/blysh/Documents/Telescope Group Project/Reduced Images")
    
    add_fil_bin = {'xbinning': 1}
    
    bias_files = raw_files.files_filtered(include_path=True, imagetyp='BIAS', **add_fil_bin)

    flats_files_V = raw_files.files_filtered(include_path=True, imagetyp='FLAT', filter='V', **add_fil_bin)
    flats_files_B = raw_files.files_filtered(include_path=True, imagetyp='FLAT', filter='B', **add_fil_bin)
    flats_files_I = raw_files.files_filtered(include_path=True, imagetyp='FLAT', filter='I', **add_fil_bin)

    exoplanet_files = data_files.files_filtered(include_path=True, imagetyp='Object', filter='V', object='TrES2b', **add_fil_bin)
    #exoplanet_files_header = data_files.hdus(imagetyp='Object', filter='V', object='TrES2b', **add_fil_bin)
    print(len(exoplanet_files))

    SA38_files_V = data_files.files_filtered(include_path=True, imagetyp='Object', filter='V', object='SA38', **add_fil_bin)
    
    exoplanet_files_names = data_files.files_filtered(imagetyp='Object', filter='V', object='TrES2b', **add_fil_bin)

    #print(bias_files, '\n')
    '''
    for bias in bias_files:
        read_image = CCDData.read(bias, unit='adu')
        trimmed_image = ccdp.subtract_overscan(read_image, overscan=read_image[:,:], median=True)
        trimmed_image = ccdp.trim_image(trimmed_image[:4076, :4061])
        trimmed_image.write('{0}/{1}_trimmed.fits'.format(os.path.join(base_dir, "Reduced Images"), os.path.basename(bias).replace('.fits', '')), overwrite=True)

    calibrated_biases = reduced_images.files_filtered(include_path=True, imagetyp='BIAS')
    print(calibrated_biases)
    combined_bias = ccdp.combine(calibrated_biases, method='average', sigma_clip=True, sigma_clip_low_thresh=5,
                                 sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,
                                 mem_limit=350e6)
    combined_bias.meta['combined'] = True
    combined_bias.write('{}/combined_bias.fits'.format(os.path.join(base_dir, "Reduced Images")), overwite=True)
    
    

    
    for flats in flats_files_V:
        read_image = CCDData.read(flats, unit='adu')
        trimmed_image = ccdp.subtract_overscan(read_image, overscan=read_image[:,:], median=True)
        trimmed_image = ccdp.trim_image(trimmed_image[:4076, :4061])
        trimmed_image.write('{0}/{1}_trimmed.fits'.format(os.path.join(base_dir, "Reduced Images"), 
                            os.path.basename(flats).replace('.fits', '')), overwrite=True)
    
    for f in reduced_images.files_filtered(include_path=True, imagetyp='FLAT'):
        filedata = fits.open(f)
        print(filedata[0].data/np.mean(filedata[0].data))
        print(filedata[0].data.shape)
    
    
    flats_trimmed = reduced_images.files_filtered(include_path=True, imagetyp='FLAT')
    for flats in flats_trimmed:
        image_data = fits.getdata(flats)
        fits_plot(image_data=image_data, cmap='afmhot',
                norm=Normalize(vmin=np.percentile(image_data, 1), vmax=np.percentile(image_data, 98)))



    EXTRA
    for ccd, file_name in raw_files.ccds(imagetyp='BIAS', ccd_kwargs={'unit': 'adu'}, return_fname=True):
        ccd = ccdp.subtract_overscan(ccd, overscan=ccd[:,:], median=True)
        ccd = ccdp.trim_image(ccd[:4076, :4061])
        ccd.write('C:/Users/blysh/Documents/Telescope Group Project/Reduced Images/{}'.format())

    for filedata in ra.hdus(imagetyp='Object', filter='V', object='TrES2b'):
        print(filedata.header['object'])

    for a_flat in bias_files.hdus(imagetyp='FLAT'):
        print(a_flat.header['EXPOSURE'])
    
    '''
    for image, image_name in zip(exoplanet_files, exoplanet_files_names):
        image_data = fits.getdata(image)
        fits_plot(image_data=image_data, title=str(image_name), cmap='afmhot', \
                  norm=Normalize(vmin=np.percentile(image_data, 35), vmax=np.percentile(image_data, 95)))
    

if __name__ == "__main__":
    main()

#for image in bias_files:
#    fits_plot(image, cmap='flag') #xlim=(-0.5, 4049.5), ylim=(-0.5,4049.5)













##combine_images(base_dir, 'Reduced Images', 'combined_bias.fits', 'bias')

#image_count_plot(images[1][0])
#fits_plot(images[0][0], cmap='flag')
#fits_plot(images[1][0], cmap='flag')

#fits_hist("E:/Telescope Group Project/Reduced Images/combined_bias.fits", (970, 990))
#fits_plot(images[0][0], cmap='flag', xlim=(-0.5, 2036.5), ylim=(-0.5, 2029.5))


#files = ccdp.ImageFileCollection(bias_dir, keywords=['imagetyp', 'filter', 'exptime', 'naxis1', 'naxis2'])


#raw_biases = files.files_filtered(include_path=True, imagetyp='BIAS')
















## Extras:
#hdu_list = fits.open(image_path)
#hdu_list.info()
# Taking info from the Primary Block [0]
#image_data = hdu_list[0].data
#print(type(image_data))
#print(image_data.shape)
#hdu_list.close()
#print(scale.get_scale_names()) # Returns the names the different scaling functions

#In case you want to filter with keyword names that cannot be used
#        as keyword argument name, you have to unpack them using a dictionary.
#        For example if a keyword name contains a space or a ``-``::

#            >>> add_filters = {'exp-time': 20, 'ESO TPL ID': 1050}
#            >>> collection.files_filtered(imagetyp='LIGHT', **add_filters)

#ccdp.read(image_path / 'ccd.001.0.fits', unit='count')
