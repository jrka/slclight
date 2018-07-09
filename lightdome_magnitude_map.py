# After the instrumental zero-point has been calculated, 
# find the sky brightness in magnitudes per square arcsecond
# for each pixel (for now) in the images.
# Project this information onto an altitude-azimuth plot.
#
# Example Usage:
# python lightdome_magnitude_map.py lightdome_timpanogos
#
# MODIFICATION HISTORY
# 2018-07-09 JRK: Commented template file made.
#

# Import necessary packages
from startup_slclight import *
from ccdproc import ImageFileCollection

# Set directory. In the future, we'll use the next 4 lines of code
# so that the user can specify the directory from the command line.
# if len(sys.argv)!=2:
#    print('Usage:\npython lightdome_magnitude_map.py [lightdome_NAME folder] \n')
#    exit()
# dir='./data/'+sys.argv[1]+'/'
dir='./data/lightdome_timpanogos/' # Just hard-code for now.

# 1) Read in the zero-point magnitude (the resultant text file
#    from lightdome_calc_zeropoint, save it as a variable.

# 2) Read in a list of all the light images.
#    This is similar to what was done in lightdome_photometry.py
#    See that file, and lightdome_process.py, for examples of how you 
#    read in each file and preform operations on the data.
files=ImageFileCollection(dir+'/astrometry/',glob_exclude='*axy*')

# 3) For each file, 
for f in files.files:
#    3a) Read in the fit file so that the image itself is an array "data"

#    3b) Perform the math necessary to create a new array, in which the values
#        correspond to magnitudes per square arcsecond, instead of counts.

#    3c) Find the altitude and azimuth associate with each pixel.
#        See near the end of lightdome_photometry for how this is done
#        using cAltAz = c.transform_to(coords.AltAz(obstime = time, location = loc))

#    3d) In the end, create a table (astropy table?) that has altitude, azimuth, 
#        and sky brightness. 

#    3e) For now, I suggest writing this table as a text/ASCII file, perhaps, so we 
#        keep a copy saved for future imaging?

#    3f) Compile the results of this file's table with all the other files
#        tables. See how I did this in the start (step 2) of lightdome_calc_zeropoint.py
#        In that case, I had to read in each table from a file first.
#        Here we don't have to do that part; the table is already in memory.

# 4) Create a polar plot of sky magnitude. (The Duriscoe plots that are circle-shaped.)
#    In the case of our polar plot, the angle is the azimuth, and the radius 
#    is the zenith angle (90-altitude). See an example polar plot in 
#    plan_targets.py (but ignore the polygon/shape drawing stuff). 
#    Choose an appropriate colorbar, and be sure to include a labeled
#    colorbar on the image.
#    Search around for matplotlib help on this! 

# 5) Optional, also create a plot displaying the same information
#    in the same projection as Duriscoe. (Hammer-Aitoff Projection?)