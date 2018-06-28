# Process raw light dome light images using calibration files.
# Dark subtract and flatframe
# Save images in data/lightdome_NAME/final/
#
# Example Usage:
# python lightdome_process.py lightdome_timpanogos
#
# Modification History
# 2018-06-28 JRK: Create file of lights summary table, for 
#     easy access to RA/Dec info for astrometry.net submission.
#     Also added a '/' after "dir" directory for easier string addition.
# 2018-06-27 JRK: Created file. Works on Timpanogos Cave files 
#     from NPS, but we suspect they may have already been reduced!
#     May have to change the fits header keys used for our data
#     if they are different (e.g. 'exposure' or 'imagetyp')

from startup_slclight import *
from ccdproc import CCDData
from ccdproc import ImageFileCollection
from astropy.io import ascii

if len(sys.argv)!=2:
   print('Usage:\npython lightdome_process.py [lightdome_NAME folder] \n')
   exit()
dir='./data/'+sys.argv[2]+'/'
#dir='./data/lightdome_timpanogos/'

#--------------------------------

###################### READ IN FILE INFO

print 'Processing images from '+dir

# Based off of 
# https://github.com/crawfordsm/wht_reduction_scripts/blob/master/wht_basic_reductions.py

# Find all science images in the "lights" folder.
files_l=ImageFileCollection(dir+'/lights/')
# Write summary table to file for easy reference later
ascii.write(files_l.summary,dir+'summary_lights.txt',format='fixed_width')

# Check that all light exposures have the same exposure time.
if len(np.unique(files_l.summary['exposure']))==1:
    exptime=np.unique(files_l.summary['exposure']).data.data[0]
    print 'Lights exposure time: ',exptime
else:
    print 'All light images must have the same exposure time to use this procedure.'
    exit()

# Find all calibration images in the "calib" folder.
files_c=ImageFileCollection(dir+'/calib/')

# These "ImageFileCollections" include a list of the files found
# in the folder. You can view summary information as an AstroPy table:
# print files_c.summary
# The table is very big and won't show all the columns. 
# You can also view the full names of all the columns:
# print files_c.summary.colnames

###################### CREATE MASTER FRAMES

# Create the master bias frames
print "Creating master bias frame..."
bias_list = []
if len(files_c.files_filtered(imagetyp='Bias Frame'))==0:
    print 'Zero bias frames were found. SKIPPING THIS STEP!'
else:
    for filename in files_c.files_filtered(imagetyp='Bias Frame'):
        print files_c.location + filename
    	ccd = CCDData.read(files_c.location + filename, unit = u.adu)
    	#this has to be fixed as the bias section does not include the whole section that will be trimmed
    	bias_list.append(ccd)
    master_bias = ccdproc.combine(bias_list, method='median')
    master_bias.write(dir+'/master_bias.fits', clobber=True)

# Create the flat fields (bias subtracted)
print "Creating master flat frame..."
flat_list = []
if len(files_c.files_filtered(imagetyp='Flat Frame'))==0:
    print 'Zero flat frames were found. SKIPPING THIS STEP!'
    doflat=False
else:
    doflat=True
    for filename in files_c.files_filtered(imagetyp='Flat Frame'):
        print files_c.location + filename
        ccd = CCDData.read(files_c.location + filename, unit = u.adu)
        ccd = ccdproc.subtract_bias(ccd, master_bias)
        flat_list.append(ccd)
    master_flat= ccdproc.combine(flat_list, method='median')
    master_flat.write(dir+'/master_flat.fits', clobber=True)

# Master darks will be specific to the exposure time.
# Check that we only have one!
print "Create master dark frame..."
dark_list= []
if len(files_c.files_filtered(imagetyp='Dark Frame',exptime=exptime))==0:
    print 'Zero dark frames with correct exposure time were found. SKIPPING THIS STEP!'
    dodark=False
else:
    dodark=True
    for filename in files_c.files_filtered(imagetyp='Dark Frame',exptime=exptime):
        print files_c.location + filename
        ccd = CCDData.read(files_c.location + filename, unit = u.adu)
        # Don't bias subtract; keep the bias as part of the total subtraction form image
        dark_list.append(ccd)
    master_dark = ccdproc.combine(dark_list, method='median')
    master_dark.write(dir+'/master_dark.fits', clobber=True)
    
###################### PROCESS THE SCIENCE IMAGES

print "Dark subtracting and flat-fielding the images..."

# Though all in this folder should be "Light Frame", might as well check!
# Note, ccdproc routines will update the header to indicate
# that these operations have been performed.
for filename in files_l.files_filtered(imagetyp='Light Frame'):
    hdu = fits.open(files_l.location + filename)
    print files_l.location + filename
    ccd = CCDData(hdu[0].data, header=hdu[0].header, unit = u.adu)
    if dodark: 
    	ccd = ccdproc.subtract_dark(ccd, master_dark, 
    			data_exposure=exptime*u.s,dark_exposure=exptime*u.s) # Remember, dark includes bias.
    if doflat:
        ccd = ccdproc.flat_correct(ccd, master_flat)
    ccd.write(dir+'/final/'+filename, clobber=True)

print "Image reduction completed!"

