# Process raw light dome light images using calibration files.
# Dark subtract and flatframe
# Save images in data/lightdome_NAME/final/
#
# Example Usage:
# python lightdome_process.py lightdome_timpanogos
#
# Modification History
# 2018-07-13 JRK: In the summary table of fits headers,
#     print RA and Dec in decimal degrees form for ease of
#     copy-paste into nova.astrometry.net. Does NOT check
#     if the RA is already in degrees, assumed to be in hours;
#     otherwise set ratodeg ("ra to degrees") to false.
#     Only works for sexegesimal separation of spaces.
# 2018-07-12 JRK: Modifications for more flexibility.
#     Don't require all light images to have same exposure time.
#     If they DO have same exposure time, AND you have dark frames
#     of that exposure time, do dark subtraction. Otherwise,
#     just do bias subtraction. Dark current is negligible.
#     Make directory "final" if it doesn't already exist.
# 2018-07-02 JRK: Image type for flats may be "Flat Field" 
#     or "Flat Frame". Change "clobber" to "overwrite"; deprecated.
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
dir='./data/'+sys.argv[1]+'/'
#dir='./data/lightdome_timpanogos/'

# Change to FALSE if you RA is already given in degrees and should not be
# converted from hours to degrees.
ratodeg=True

#--------------------------------

###################### READ IN FILE INFO

print 'Processing images from '+dir

# Based off of 
# https://github.com/crawfordsm/wht_reduction_scripts/blob/master/wht_basic_reductions.py

# Find all science images in the "lights" folder.
files_l=ImageFileCollection(dir+'/lights/')

# Write summary table to file for easy reference later
# See if the 'objctra' and 'objctdec' are given in sexigesimal or decimal.
# Only for separated by spaces, which is MaxImDL default.
# Assumes RA is in hours, else switch "ratodeg" above.
# Fits header may be 'ra','dec' (NPS) or 'objctra' and 'objctdec'
if np.sum([x in files_l.summary.colnames for x in ['objctra','objctdec']]):
    coordkey=['objctra','objctdec']
elif np.sum([x in files_l.summary.colnames for x in ['ra','dec']]):
    coordkey=['ra','dec']
else:
    coordkey=[]
    print 'Not sure what RA and Dec headers are; not doing a decimal conversion.'

summary=files_l.summary.copy()
if coordkey:
	colnames=summary.colnames
	for i in [0,1]:
		colnames.insert(1+i,coordkey[i]+'_decimaldeg')
		summary[coordkey[i]+'_decimaldeg']=0.0
		summary=summary[colnames]
	for i in [0,1]:
		if isinstance(summary[coordkey[i]][0],float):
		    summary[coordkey[i]+'_decimaldeg']=summary[coordkey[i]]
		elif isinstance(summary[coordkey[i]][0],basestring):
			if unicode(summary[coordkey[i]][0]).isnumeric():   # A string that is numeric
				summary[coordkey[i]+'_decimaldeg']=[np.float(x) for x in summary[coordkey[i]]]
			else:  # Assume a non-numeric string is a sexegesimal with spaces.
			    mult=-1 if np.float(summary[coordkey[i]][0].split(' ')[0])<0 else 1
			    summary[coordkey[i]+'_decimaldeg']=[np.float(str.split(x,' ')[0])+mult*np.float(str.split(x,' ')[1])/60.0+mult*np.float(str.split(x,' ')[2])/3600.0 for x in summary[coordkey[i]]]
	if ratodeg: summary[coordkey[0]+'_decimaldeg']*=360.0/24.0
ascii.write(summary,dir+'/summary_lights.txt',format='fixed_width',overwrite=True)

# Check that all light exposures have the same exposure time.
if len(np.unique(files_l.summary['exposure']))==1:
    exptime=np.unique(files_l.summary['exposure']).data.data[0]
    print 'Lights exposure time: ',exptime
    dodark=True
else:
    print 'All images do not have the same exposure time. Will '
    print 'only perform bias subtraction, not dark subtraction.'
    dodark=False

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
    master_bias.write(dir+'/master_bias.fits', overwrite=True)

# Create the flat fields (bias subtracted)
# NPS data used 'Flat Frame' whereas South Physics used 'Flat Field'
print "Creating master flat frame..."
flat_list = []
if len(files_c.files_filtered(imagetyp='Flat Frame'))==0 and len(files_c.files_filtered(imagetyp='Flat Field'))==0:
    print 'Zero flat frames were found. SKIPPING THIS STEP!'
    doflat=False
else:
    doflat=True
    tmp1=files_c.files_filtered(imagetyp='Flat Frame')
    tmp2=files_c.files_filtered(imagetyp='Flat Field')
    flatfiles = tmp1 if len(tmp1)>0 else tmp2
    for filename in flatfiles:
        print files_c.location + filename
        ccd = CCDData.read(files_c.location + filename, unit = u.adu)
        ccd = ccdproc.subtract_bias(ccd, master_bias)
        flat_list.append(ccd)
    master_flat= ccdproc.combine(flat_list, method='median')
    master_flat.write(dir+'/master_flat.fits', overwrite=True)

# Master darks will be specific to the exposure time.
# Check that we only have one!
print "Create master dark frame..."
dark_list= []
if dodark:
    if len(files_c.files_filtered(imagetyp='Dark Frame',exptime=exptime))==0:
        print 'Zero dark frames with correct exposure time were found. SKIPPING THIS STEP!'
        dodark=False
    else:
        for filename in files_c.files_filtered(imagetyp='Dark Frame',exptime=exptime):
            print files_c.location + filename
            ccd = CCDData.read(files_c.location + filename, unit = u.adu)
            # Don't bias subtract; keep the bias as part of the total subtraction form image
            dark_list.append(ccd)
        master_dark = ccdproc.combine(dark_list, method='median')
        master_dark.write(dir+'/master_dark.fits', overwrite=True)
    
###################### PROCESS THE SCIENCE IMAGES

print "Dark (or bias) subtracting and flat-fielding the images..."

# Create the "final" subdirectory if it doesn't already exist.
if not os.path.exists(dir+'/final/'): os.mkdir(dir+'/final/')

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
    else:
        ccd = ccdproc.subtract_bias(ccd, master_bias)
    if doflat:
        ccd = ccdproc.flat_correct(ccd, master_flat)
    ccd.write(dir+'/final/'+filename, overwrite=True)

print "Image reduction completed!"

