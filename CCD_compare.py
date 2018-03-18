
# coding: utf-8

# # CCD Image Comparison

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import MinMaxInterval,SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib.patches import Circle
get_ipython().magic(u'matplotlib inline')

# Define the files to be used
Compare my processed files to Amber's to check for consistency!
# In[2]:

# Need to update later
dir1='images_processed_20170712_XHer'
dir2='images_processed_20170712_XHer_JRK'
filenames=['Star72_001V']


# In[3]:

# Open the files and compare
for f in filenames:
    inhdulist=fits.open(dir1+f)
    im1=inhdulist[0].data
    inhdulist.close()
    
    inhdulist=fits.open(dir2+f)
    im2=inhdulist[0].data
    inhdulist.close()
    
    plt.hist([im1,im2,im2-im1])


# In[ ]:

# Compare to MaxIm DL processed images. Must be saved UNCOMPRESSED!!!
inhdulist = fits.open('images_processed_20170712_XHer_JRK'+'/'+rawfiles[0]+'.fit',ignore_missing_end=True)
maximdl_light1=inhdulist[0].data
inhdulist.close()
# Something wrong with the files. Just get the ones we want.
print 'MaxIm DL star pix value',maximdl_light1[argmax[0],argmax[1]]
print 'MaxIm DL blank pix value',maximdl_light1[argblank[0],argblank[1]]
print 'My star pix value',image_average[0,argmax[0],argmax[1]]
print 'My blank pix value',image_average[0,argblank[0],argblank[1]]

print maximdl_light1[argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]
print image_average[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]

print (maximdl_light1[argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]-image_average[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5])/maximdl_light1[argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5].astype(np.float)


# # Perform Aperture Photometry

# In[4]:

# http://photutils.readthedocs.io/en/stable/photutils/aperture.html#local-background-subtraction
# pip install --no-deps photutils

from photutils import CircularAnnulus,CircularAperture,aperture_photometry

apertures = CircularAperture(argmax, r=5)
annulus_apertures = CircularAnnulus(argmax, r_in=7., r_out=9.)
apers = [apertures, annulus_apertures]
phot_table = aperture_photometry(image_stacked, apers)
print(phot_table)
bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
bkg_sum = bkg_mean * apertures.area()
final_sum = phot_table['aperture_sum_0'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum
print(phot_table['residual_aperture_sum'])  
print phot_table


# In[ ]:




# In[ ]:



