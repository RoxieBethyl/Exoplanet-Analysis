"""
Author: Blythe Fernandes
"""

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import ccdproc as ccdp
from matplotlib.colors import Normalize
from pathlib import Path as path

def fits_plot(image_path='', image_data='', cmap='afmhot', title='' , norm=None, vmin=None, vmax=None, xlim=(None, None), ylim=(None, None), details=False): #norm=colors.NoNorm(), vmin=2E3, vmax=3E3
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
    if title != '':
        plt.title(title)
    elif len(image_path) != 0:
        image_data = fits.getdata(image_path)
        plt.title(os.path.basename(image_path))
        
    # print(type(image_data)) # Show the Python type for image_data
    plt.imshow(image_data, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax) #, norm=norm, vmin=vmin, vmax=vmax, 
    plt.gca().invert_yaxis() #inverts y-axis into the right orientation
    plt.xlabel("x-pixels")
    plt.ylabel("y-pixels")
    plt.colorbar()
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()
    
    if details == True:
        print(image_data.shape) # Show the number of pixels per side in the 2-D image
        print('Min:', np.min(image_data))
        print('Max:', np.max(image_data))
        print('Mean:', np.mean(image_data))
        print('Stdev:', np.std(image_data), '\n')

def main():
    
    add_fil_bin = {'xbinning': 1}

    data_files = ccdp.ImageFileCollection(path("C:/Users/blysh/Documents/Telescope Group Project/2022-10-09 Explanets 1"))
    exoplanet_files = data_files.files_filtered(include_path=True, imagetyp='Object', filter='V', object='TrES2b', **add_fil_bin)
    exoplanet_files_names = data_files.files_filtered(imagetyp='Object', filter='V', object='TrES2b', **add_fil_bin)

    for image, image_name in zip(exoplanet_files, exoplanet_files_names):
            image_data = fits.getdata(image)
            fits_plot(image_data=image_data, title=str(image_name), cmap='afmhot', \
                    norm=Normalize(vmin=np.percentile(image_data, 35), vmax=np.percentile(image_data, 95)))


if __name__ == "__main__":
    main()                  