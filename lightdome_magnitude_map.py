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
import re
from scipy.interpolate import griddata
from astropy.wcs import WCS
from astropy import coordinates as coords
from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.table import Table, Column
import glob
from astropy.table import vstack
from astropy.io import ascii
from scipy import interpolate



# Set directory. In the future, we'll use the next 4 lines of code
# so that the user can specify the directory from the command line.
# if len(sys.argv)!=2:
#    print('Usage:\npython lightdome_magnitude_map.py [lightdome_NAME folder] \n')
#    exit()
# dir='./data/'+sys.argv[1]+'/'
dir='./data/lightdome_timpanogos' # Just hard-code for now.

if not os.path.exists(dir+'/skybrightness/'): os.mkdir(dir+'/skybrightness/')

# 1) Read in the zero-point magnitude (the resultant text file
#    from lightdome_calc_zeropoint, save it as a variable.

with open(dir+'/zeropoint_information.txt') as t:
	content= t.readlines()
C = re.findall("\d+\.\d+", content[1])[0]
C = float(C)
#Bring the error in too.

# 2) Read in a list of all the light images.
#    This is similar to what was done in lightdome_photometry.py
#    See that file, and lightdome_process.py, for examples of how you 
#    read in each file and preform operations on the data.
files=ImageFileCollection(dir+'/astrometry/',glob_exclude='*axy*')


# 3) For each file, 
for f in files.files:
#    3a) Read in the fit file so that the image itself is an array "data"
	hdu=fits.open(dir+'/astrometry/'+f)
	data=hdu[0].data
	header=hdu[0].header
	 
#    3b) Perform the math necessary to create a new array, in which the values
#        correspond to magnitudes per square arcsecond, instead of counts.

	resolution = platescale(header)
	resolution_squared = resolution**2

	instrmagsky = data/(hdu[0].header['EXPTIME']*resolution_squared)
	log_instrmagsky = np.log10(instrmagsky)
	magnitude_sky = C - 2.5*log_instrmagsky
	

#    3c) Find the altitude and azimuth associated with each pixel.
#        See near the end of lightdome_photometry for how this is done
#        using cAltAz = c.transform_to(coords.AltAz(obstime = time, location = loc))
	wcs = WCS(hdu[0].header)
	dirname='lightdome_timpanogos' # Just hard-code for now.
	metad=read_metadata(dirname)
	loc = coords.EarthLocation(lat = metad['lat']*u.deg, lon = metad['lon']*u.deg, height = metad['elev']*u.m)
	time = Time(hdu[0].header['DATE-OBS'],scale='utc')
	
	xlist = np.arange(0,len(magnitude_sky[:,0]))
	ylist = np.arange(0,len(magnitude_sky[0,:]))
	xarr, yarr = np.meshgrid(xlist,ylist,indexing='ij')
	
	## Regridding
	
	# 1. Data points right above
	# 2. Values: Sky Magnitude
	# 3. xi = New smaller arrays (555,15 15 15).. points at which we want the values
	
	j = (3*3600/resolution)
	k = (2*3600/resolution)
	off_j = j/2
	off_k = k/2
	grid_x, grid_y =  np.mgrid[np.min(xlist):np.max(xlist):j, np.min(ylist):np.max(ylist):k]
	grid_x += off_j
	grid_y += off_k
	
	r = Table()
	s = Table()
	
	for i in np.arange(0,len(ylist)):
		s['x'] = xlist
		s['y'] = ylist[i]
		s['Sky_Brightness'] = magnitude_sky[:,i]
		tmp = s
		if not r:
			r = s
		else:	
			r = vstack([r, tmp])
	
	# If the magnitude is nan then it is not within the picture, delete line on the array:
	grid = griddata(np.transpose([r['x'],r['y']]), r['Sky_Brightness'], (grid_x, grid_y), method='linear')
# 	grid_nonan= (~np.isnan(grid))
# 	grid= grid[grid_nonan]
	
	
	positions_wcs = wcs.all_pix2world(grid_x,grid_y,1) #Investigate about distortion corrections.
	
	c = coords.SkyCoord(positions_wcs[0],positions_wcs[1], unit='deg')
	
	cAltAz = c.transform_to(coords.AltAz(obstime = time, location = loc))


#    3d) In the end, create a table (astropy table?) that has altitude, azimuth, 
#        and sky brightness. 

	t = Table()
	s = Table()
	
	for i in np.arange(0,len(cAltAz.alt[:,0])):
		s['Altitude'] = cAltAz.alt[i,:]
		s['Azimuth'] = cAltAz.az[i,:]
		s['Sky_Brightness'] = grid[i,:]
		tmp = s
		if not t:
			t = s
		else:	
			t = vstack([t, tmp])

# Attempt at reducing the number of data points in t:

#    3e) For now, I suggest writing this table as a text/ASCII file, perhaps, so we 
#        keep a copy saved for future imaging?

	t.write(dir+'/skybrightness/'+f.split('.')[0]+'.txt', format='ascii.fixed_width',overwrite=True)
	print 'Table created' #Put table number

#    3f) Compile the results of this file's table with all the other files
#        tables. See how I did this in the start (step 2) of lightdome_calc_zeropoint.py
#        In that case, I had to read in each table from a file first.
#        Here we don't have to do that part; the table is already in memory.

files=glob.glob(dir+'/skybrightness/*.txt')
data=[]
for f in files:
	tmp = ascii.read(f,format='fixed_width')
	tmp['file'] = f.split('/')[-1]   # Add in from which file these rows came
	if not data:
		data = tmp
	else:
		data = vstack([data,tmp])

print data


# 4) Create a polar plot of sky magnitude. (The Duriscoe plots that are circle-shaped.)
#    In the case of our polar plot, the angle is the azimuth, and the radius 
#    is the zenith angle (90-altitude). See an example polar plot in 
#    plan_targets.py (but ignore the polygon/shape drawing stuff). 
#    Choose an appropriate colorbar, and be sure to include a labeled
#    colorbar on the image.
#    Search around for matplotlib help on this!


azimuths = np.array(data['Azimuth'])
altitudes = np.array(data['Altitude'])
values= np.array(data['Sky_Brightness'])



	


 

# 5) Optional, also create a plot displaying the same information
#    in the same projection as Duriscoe. (Hammer-Aitoff Projection?)






























