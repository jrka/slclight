# Explore the results given by astrometry.net
# Assumes that all the files produced by astrometry.net are located:
#     slclight/data/lightdome_timpanogos/astrometry/ib030/

from startup_slclight import *
import glob

#####################################################################
# Get all the files
dir='data/lightdome_timpanogos/'
prefix='ib030'

origf=dir+'lights/ib030.fit'
files=glob.glob(dir+'astrometry/'+prefix+'/*.fit*')

#####################################################################
# Examine what's in each one
for f in files:
    hdu=fits.open(f)
    print f
    for h in hdu[0].header: print h
    print hdu.info() 
    hdu.close()

#####################################################################
# What is different between the original image and new-image?
# new_image has the fk5 values (RA, Dec). 
# Pixel values are preserved (the data array is the same)
orig=fits.open(origf)
print orig[0].data
new=fits.open(dir+'astrometry/'+prefix+'/new-image.fits')
print new[0].data
print new[0].data-orig[0].data

# What is different about the headers?
# fits.printdiff doesn't seem to work; try "conda update astropy"
# Going from 1.2.1 to 2.0.7
fits.printdiff(origf,dir+'astrometry/'+prefix+'/new-image.fits')
# The new file has extra keywords in the header (175 vs. 47)
# These are the ones that give the WCS values. 
# Use this for aperture photometry.
orig.close()
new.close()

#####################################################################
# What is in WCS file?
wcs=fits.open(dir+'astrometry/'+prefix+'/wcs.fits')
print wcs[0].header
# This is mostly information about the solving parameters and procedure.
# We don't necessarily need this one. 
wcs.close()

#####################################################################
# What is in "Reference star nearby" table, rdls.fits?
rdls=fits.open(dir+'astrometry/'+prefix+'/rdls.fits')
rdls[0].header
rdls[1].header
rdls[1].data
# This is just a list of the RA and Dec of the objects used as reference
# stars to find the solution.
rdls.close()

#####################################################################
# Stars detected in your image, axy.fits
axy=fits.open(dir+'astrometry/'+prefix+'/axy.fits')
axy[0].header
axy[1].header
axy[1].data
# This is a 4-column table of the stars found in the image itself.
# The columsn are x position (in pixels), y position (in pixels), flux, and background
# How were the last two calculated? Let's find out.
# The first one says it's located at 595.84, 59.88. That's a bright one 
# in the bottom, middle-rightish (with DS9 orientation; tree on top)
# Flux and background listed as 63035.109, 2499.8906
# The max pixel value there is 65535, this appears to be the sum
# of what's listed as "flux" and "background". So these columns appear
# to be pixel specific (not aperture photometry).
axy.close()

#####################################################################
# Correspondences between image and reference stars (table), corr.fits
corr=fits.open(dir+'astrometry/'+prefix+'/corr.fits')
corr[0].header
corr[1].header
corr[1].data
# This is a 13 column table. The columns are given as:
# field_x, field_y, field_ra, field_dec, index_x, index_y, index_ra, index_dec,
# index_id, field_id, match_weight, flux, background.
# First row corresponds to our identified star above in axy.fits.
# RA and Dec seem to be in degrees. Flux, background same as above.
corr.close()

