
# coding: utf-8

# # CCD Image Processing in Python
# Based on http://prancer.physics.louisville.edu/astrowiki/index.php/Image_processing_with_Python_and_SciPy

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import MinMaxInterval,SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib.patches import Circle
get_ipython().magic(u'matplotlib inline')


# # Define the files to be used
# This script is meant to show each step of calibration and stacking for ONE object and ONE filter at a time. In this example, we'll do the images of XHer from July 12th, 2017 in the V band. The exposure time was 3 s.

# In[2]:

dir='images_raw_20170712_XHer'
rawfiles=['XHer-001V','XHer-002V','XHer-003V','XHer-004V','XHer-005V',
          'XHer-006V','XHer-007V','XHer-008V','XHer-009V','XHer-010V']
darkfiles=['XHer-001Dark3','XHer-002Dark3','XHer-003Dark3','XHer-004Dark3','XHer-005Dark3',
           'XHer-006Dark3','XHer-007Dark3','XHer-008Dark3','XHer-009Dark3','XHer-010Dark3']
biasfiles=['XHer-001Bias','XHer-002Bias','XHer-003Bias','XHer-004Bias','XHer-005Bias',
           'XHer-006Bias','XHer-007Bias','XHer-008Bias','XHer-009Bias','XHer-010Bias']
flatfiles=['Flats-001Flat-V','Flats-002Flat-V','Flats-003Flat-V','Flats-004Flat-V','Flats-005Flat-V',
           'Flats-006Flat-V','Flats-007Flat-V','Flats-008Flat-V','Flats-009Flat-V','Flats-010Flat-V']
flatdarkfiles=['Flats-001DarkFlat10','Flats-002DarkFlat10','Flats-003DarkFlat10','Flats-004DarkFlat10','Flats-005DarkFlat10',
               'Flats-006DarkFlat10','Flats-007DarkFlat10','Flats-008DarkFlat10','Flats-009DarkFlat10','Flats-010DarkFlat10']


# Or, do the R-band instead. Exposure time was 0.5 seconds for XHer. Flats were 10 seconds? Double check.

# In[3]:

#dir='images_raw_20170712_XHer'
#rawfiles=['XHer-001R','XHer-002R','XHer-003R','XHer-004R','XHer-005R',
#          'XHer-006R','XHer-007R','XHer-008R','XHer-009R','XHer-010R']
#darkfiles=['XHer-001Dark0P5','XHer-002Dark0P5','XHer-003Dark0P5','XHer-004Dark0P5','XHer-005Dark0P5',
#           'XHer-006Dark0P5','XHer-007Dark0P5','XHer-008Dark0P5','XHer-009Dark0P5','XHer-010Dark0P5']
#biasfiles=['XHer-001Bias','XHer-002Bias','XHer-003Bias','XHer-004Bias','XHer-005Bias',
#           'XHer-006Bias','XHer-007Bias','XHer-008Bias','XHer-009Bias','XHer-010Bias']
#flatfiles=['Flats-001Flat-R','Flats-002Flat-R','Flats-003Flat-R','Flats-004Flat-R','Flats-005Flat-R',
#           'Flats-006Flat-R','Flats-007Flat-R','Flats-008Flat-R','Flats-009Flat-R','Flats-010Flat-R']
#flatdarkfiles=['Flats-001DarkFlat10','Flats-002DarkFlat10','Flats-003DarkFlat10','Flats-004DarkFlat10','Flats-005DarkFlat10',
#               'Flats-006DarkFlat10','Flats-007DarkFlat10','Flats-008DarkFlat10','Flats-009DarkFlat10','Flats-010DarkFlat10']


# 

# # Take a preliminary look at one of our "science" images, the raw file

# In[4]:

# Some useful commands for working with fits files are here: http://docs.astropy.org/en/stable/io/fits/index.html
inhdulist=fits.open(dir+'/'+rawfiles[0]+'.fit')
# Show information
inhdulist.info()


# In[5]:

# Take the image data info and keep it as an array
image=inhdulist[0].data
# Image is an array of the following dimensions
dims=image.shape
print dims
# Close the fits file; we still have our image array.
inhdulist.close()


# In[6]:

# Let's figure out some statistics about our image.
print 'Maximum pixel value: ',np.max(image)
print 'Minimum pixel value: ',np.min(image)
print 'Median pixel value: ',np.median(image)
print 'Average pixel value: ',np.average(image)
print 'Standard deviation: ',np.std(image)
argmax=np.unravel_index(np.argmax(image),dims)
print 'Location of maximum pixel value: ',argmax


# In[7]:

# Test that max pixel location:
print image[argmax]


# In[8]:

# What does the array look like?
print image


# In[9]:

# The array is very large. Let's "zoom in" on a big of blank sky.
argblank=(745,1175)
imblank=image[argblank[0]-4:argblank[0]+5,argblank[1]-4:argblank[1]+5]
print imblank


# In[10]:

# Now let's "zoom in" on just the part around our star, which is where the maximum pixel is.
immax=image[argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]
print immax


# In[11]:

# Let's view the image. NOTICE that in the following images, the color scale is DIFFERENT! Look at the colorbar.
plt.figure()
plt.imshow(image, cmap='gray')
plt.colorbar()


# In[12]:

# Now let's view zoomed in on a "blank" space
plt.figure()
plt.imshow(imblank, cmap='gray',interpolation='none')
plt.colorbar()


# In[13]:

# Finally, view zoomed in onto our star.
plt.figure()
plt.imshow(immax, cmap='gray',interpolation='none')
plt.colorbar()


# # Create "Master" Bias Frame

# In[14]:

# Create an array for all the bias files
bias_all=np.ndarray([len(biasfiles),dims[0],dims[1]])
bias_all.shape
# Open each bias file, and save it in the array.
for i, file in enumerate(biasfiles):
    inhdulist = fits.open(dir+'/'+file+'.fit')
    bias_all[i,:,:]=inhdulist[0].data
    inhdulist.close()


# In[15]:

# Examine our "star" area of the sky for the first bias image.
print bias_all[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]
plt.figure()
plt.imshow(bias_all[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5], cmap='gray',interpolation='none')
plt.colorbar()


# In[16]:

# Create the master frame
bias_median=np.median(bias_all, axis=0)
bias_average=np.average(bias_all, axis=0)


# In[17]:

# Let's see the differences between the two!
print 'Max of median combined bias frame: ',np.max(bias_median),' at pixel ',np.unravel_index(np.argmax(bias_median),dims)
print 'Min of median combined bias frame: ',np.min(bias_median)
print 'Max of average combined bias frame: ',np.max(bias_average),' at pixel ',np.unravel_index(np.argmax(bias_average),dims)
print 'Min of average combined bias frame: ',np.min(bias_average)
print 'Average combined bias frame at star pixel: ',bias_average[argmax]
print 'Average combined bias frame at blank pixel: ',bias_average[argblank]
plt.figure()
plt.imshow(bias_median-bias_average, cmap='gray',interpolation='none')
plt.colorbar()


# In[18]:

# Compare to master bias frame from MaxImDL
inhdulist = fits.open(dir+'/'+'Master_Bias1_1530x1020_Bin1x1_Temp-5C_ExpTime0ms'+'.fit')
maximdl_bias=inhdulist[0].data
inhdulist.close()
print 'Max of MaxIm DL bias frame: ',np.max(maximdl_bias),' at pixel ',np.unravel_index(np.argmax(maximdl_bias),dims)
print 'Min of MaxIm bias frame: ',np.min(maximdl_bias)
print 'MaxIm DL bias at star pixel',maximdl_bias[argmax]
print 'MaxIm DL bias at other pixel',maximdl_bias[argblank]


# I find a difference of 100 pixels! This is because of the "pedestal" setting in MaxIm DL. It adds 100 pixels to all the values.

# # Create "Master" dark

# we stack all of the original individual dark images to make a 3-d stack of 2-d arrays. Using numpy arrays we would have
# dark_stack = np.array([dark_1, dark_2, dark_3])
# where dark_1, dark_2, and dark_3 are the original dark images. We need at least 3, or any odd number in the dark stack. If the images are m rows of n columns, and if we have k images in the stack, the stack will have a shape (k,n,m): k images, each of n rows, each of m columns. A median on the first axis of this stack returns the median value for each pixel in the stack --
# dark_median = np.median(dark_stack, axis=0)
# and has a shape that is (n,m) with the elements that we wanted.

# In[19]:

# Create an array for all the dark files
dark_all=np.ndarray([len(darkfiles),dims[0],dims[1]])
dark_all.shape


# In[20]:

# Open each dark file, and save it in the array.
for i, file in enumerate(darkfiles):
    inhdulist = fits.open(dir+'/'+file+'.fit')
    dark_all[i,:,:]=inhdulist[0].data
    inhdulist.close()


# In[21]:

# Examine our "star" area of the sky for the first dark image.
print dark_all[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]
plt.figure()
plt.imshow(dark_all[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5], cmap='gray',interpolation='none')
plt.colorbar()


# In[22]:

# Usually, you create your "master" dark using a median combine. MaxIMDL seems to default to average. Let's compare!
# Be sure to bias-subtract the images! Use aveage for comparison to maximdl.
dark_median=np.median(dark_all, axis=0)-bias_average
dark_average=np.average(dark_all, axis=0)-bias_average


# In[23]:

# Let's see the differences between the two!
print 'Max of median combined dark frame: ',np.max(dark_median),' at pixel ',np.unravel_index(np.argmax(dark_median),dims)
print 'Min of median combined dark frame: ',np.min(dark_median)
print 'Max of average combined dark frame: ',np.max(dark_average),' at pixel ',np.unravel_index(np.argmax(dark_average),dims)
print 'Min of average combined dark frame: ',np.min(dark_average)
print 'Average combined bias frame at star pixel: ',dark_average[argmax]
print 'Average combined bias frame at blank pixel: ',dark_average[argblank]
plt.figure()
plt.imshow(dark_median-dark_average, cmap='gray',interpolation='none')
plt.colorbar()


# In[24]:

# Compare to master dark frame from MaxImDL
inhdulist = fits.open(dir+'/'+'Master_Dark 6_1530x1020_Bin1x1_Temp-5C_ExpTime3s'+'.fit')
maximdl_dark=inhdulist[0].data
inhdulist.close()
print 'Max of MaxIm DL dark frame: ',np.max(maximdl_dark),' at pixel ',np.unravel_index(np.argmax(maximdl_dark),dims)
print 'Min of MaxIm dark frame: ',np.min(maximdl_dark)
print 'MaxIm Dark Frame at star pixel: ',maximdl_dark[argmax]
print 'MaxIm Dark Frame at blank pixel: ',maximdl_dark[argblank]


# Again, the difference of 100 is due to the "pedestal" value.

# # Create "Master" Flat Frame

# In[25]:

# Create an array for all the flat files
flat_all=np.ndarray([len(flatfiles),dims[0],dims[1]])
# And an array for the flat-darks; those are the darks that are 10 seconds long.
flatdark_all=np.ndarray([len(flatdarkfiles),dims[0],dims[1]])
# Open each flat file, and save it in the array.
for i, file in enumerate(flatfiles):
    inhdulist = fits.open(dir+'/'+file+'.fit')
    flat_all[i,:,:]=inhdulist[0].data
    inhdulist.close()
# Open each flatdark file, and save it in the array.
for i, file in enumerate(flatdarkfiles):
    inhdulist = fits.open(dir+'/'+file+'.fit')
    flatdark_all[i,:,:]=inhdulist[0].data
    inhdulist.close()


# In[26]:

# Before we proceed, create a master flatdark, and subtract it from ALL flat files.
# Use average, for comparison to MaxIm DL. Subtract bias.
flatdark_average=np.average(flatdark_all,axis=0)-bias_average
# Can compare this to MaxImDL master dark for 10 s.


# In[27]:

# Examine our "star" area of the sky for the first flat image.
print flat_all[0,argmax[0]-3:argmax[0]+4,argmax[1]-3:argmax[1]+4] # Smaller region; floats print larger
plt.figure()
plt.imshow(flat_all[0,argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5], cmap='gray',interpolation='none')
plt.colorbar()
# Notice that the values are all very close together, but not perfectly uniform.


# In[28]:

# Examine the entire first flat frame. Make sure we are scaled from min to max.
# We want to use a "square root stretch" for our color scaling to better see the features.
norm = ImageNormalize(flat_all[0,:,:],stretch=SqrtStretch())
plt.figure()
plt.imshow(flat_all[0,:,:], cmap='gray',norm=norm)


# In[29]:

# First, median combine all the flat fields.
flat_median=np.median(flat_all, axis=0)-flatdark_average-bias_average
flat_average=np.average(flat_all,axis=0)-flatdark_average-bias_average


# In[30]:

# Let's see the differences between the two!
print 'Max of average combined flat frame: ',np.max(flat_average),' at pixel ',np.unravel_index(np.argmax(flat_average),dims)
print 'Min of average combined flat frame: ',np.min(flat_average)
print 'Average combined flat frame at star pixel: ',flat_average[argmax]
print 'Average combined flat frame at blank pixel: ',flat_average[argblank]
print 'Average flat',np.average(flat_average)


# In[31]:

# Compare to master flat frame from MaxImDL
inhdulist = fits.open(dir+'/'+'Master_Flat V 1_V_1530x1020_Bin1x1_Temp-5C_ExpTime10s.fit')
maximdl_flat=inhdulist[0].data
inhdulist.close()
print 'Max of MaxIm DL flat frame: ',np.max(maximdl_flat),' at pixel ',np.unravel_index(np.argmax(maximdl_flat),dims)
print 'Min of MaxIm flat frame: ',np.min(maximdl_flat)
print 'MaxIm Flat Frame at star pixel: ',maximdl_flat[argmax]
print 'MaxIm Flat Frame at blank pixel: ',maximdl_flat[argblank]
print 'MaxIm Flat Average: ',np.average(maximdl_flat)


# In[32]:

# We are not concerned with the actual values of the pixels. We want a NORMALIZED flat frame, where the 
# mean value corresponds to 1, and all the pixels have a value relative to 1.
# Now, divide by the average value of our combined flat field.
# The values should all be around 1, some a little higher, some a little lower
flat_normalized=flat_average / np.average(flat_average)
# Check out the result, zoomed in on our star area:
print flat_normalized[argmax[0]-4:argmax[0]+5,argmax[1]-4:argmax[1]+5]


# # Let's do a test. Can we reproduce the final image by HAND? Just try it for TWO pixels!

# In[35]:

# Commented out b/c it takes a while to run and we only need it once.
# Create a spreadsheet of all the files we need for the two pixels
#allfiles=[rawfiles,darkfiles,biasfiles,flatfiles,flatdarkfiles]
#allfiles = [item for sublist in allfiles for item in sublist]
#for i, file in enumerate(allfiles):
#    inhdulist = fits.open(dir+'/'+file+'.fit')
#    image=inhdulist[0].data
##    inhdulist.close()
#    print file,',',image[argmax],',',image[argblank]


# In[36]:

print np.mean(flat_average)
print argmax,argblank


# # Calibrate each raw image

# In[42]:

raw_all=np.ndarray([len(rawfiles),dims[0],dims[1]])
for i, file in enumerate(rawfiles):
    inhdulist = fits.open(dir+'/'+file+'.fit')
    raw_all[i,:,:]=inhdulist[0].data
    inhdulist.close()


# In[53]:

print raw_all[0,527,685]


# In[69]:

image_average=(raw_all-bias_average-dark_average)/(flat_normalized)
image_average=image_average.astype(np.int64)


# In[73]:

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


# # Stack the images
# "We would do this to create a final image that is effectively one long exposure, the sum of all the contributing image exposure times. Because of guiding errors, cosmic rays, and weather, one very long exposure is often not possible, but 10's or 100's of shorter exposures can be "co-added" after selecting the best ones and aligning them so that each pixel in the contributing image corresponds to the same place in the sky."

# In[ ]:

# Compare peak / stddev of each image.
for i in range(image_median.shape[0]):
    print np.max(image_median[i,:,:]),np.std(image_median[i,:,:]),np.max(image_median[i,:,:])/np.std(image_median[i,:,:])
image_stacked=np.sum(image_median,axis=0)
image_stacked.shape
np.max(image_stacked),np.std(image_stacked),np.max(image_stacked)/np.std(image_stacked)


# In[ ]:

argmax=np.unravel_index(np.argmax(image_stacked),dims)
subimage=image_stacked[argmax[0]-30:argmax[0]+31,argmax[1]-30:argmax[1]+31]
fig,ax=plt.subplots(1)
ax.imshow(subimage, cmap='gray',interpolation='none')
print subimage.shape
print np.unravel_index(np.argmax(subimage),subimage.shape)
circ1 = Circle((30,30),5,ec='green',fc='none')
circ2 = Circle((30,30),7,ec='blue',fc='none')
circ3 = Circle((30,30),9,ec='blue',fc='none')
ax.add_patch(circ1)
ax.add_patch(circ2)
ax.add_patch(circ3)


# # Perform Aperture Photometry

# In[ ]:

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



