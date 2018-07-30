# After the instrumental zero-point has been calculated, 
# find the sky brightness in magnitudes per square arcsecond
# for each pixel (for now) in the images.
# Project this information onto an altitude-azimuth plot.
#
# Example Usage:
# python lightdome_magnitude_map.py lightdome_timpanogos
#
# MODIFICATION HISTORY
# 2018-07-25 JRK: A method that works! 
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
import matplotlib as mpl
from scipy import ndimage

# Set directory. In the future, we'll use the next 4 lines of code
# so that the user can specify the directory from the command line.
 if len(sys.argv)!=2:
    print('Usage:\npython lightdome_magnitude_map.py [lightdome_NAME folder] \n')
    exit()
dir='./data/'+sys.argv[1]+'/'
#dir='./data/lightdome_westminster' # Just hard-code for now.
#dir='./data/lightdome_timpanogos'

if not os.path.exists(dir+'/skybrightness/'): os.mkdir(dir+'/skybrightness/')

# 1) Read in the zero-point magnitude (the resultant text file
#    from lightdome_calc_zeropoint, save it as a variable.

with open(dir+'/zeropoint_information.txt') as t:
    content= t.readlines()
C = re.findall("\d+\.\d+", content[1])[0]
C = float(C)
#C = 14.761  # 14.761 timpanogos, 14.220 Westmisnter
#Bring the error in too.

# 2) Read in a list of all the light images.
#    This is similar to what was done in lightdome_photometry.py
#    See that file, and lightdome_process.py, for examples of how you 
#    read in each file and preform operations on the data.
files=ImageFileCollection(dir+'/astrometry/',glob_exclude='*axy*')


# Setup some plotting information
cmap=mpl.cm.CMRmap_r 

skipread=False
if not skipread:
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
    
        # Plot the image, with colorbar.
        norm = mpl.colors.Normalize(vmin=np.min(magnitude_sky),vmax=np.max(magnitude_sky))
        plt.clf()
        fig,ax=plt.subplots(1,3,num=0,figsize=(10,3))
        im0=ax[0].imshow(magnitude_sky,cmap=cmap,norm=norm,aspect='equal')
        fig.colorbar(im0)
    
        ## Regridding
        #### vvvvvv Old method, interpolation.
        # 1. Data points right above
        # 2. Values: Sky Magnitude
        # 3. xi = New smaller arrays (555,15 15 15).. points at which we want the values
    #     j = (3*3600/resolution)
    #     k = (2*3600/resolution)
    #     off_j = j/2
    #     off_k = k/2
    #     grid_x, grid_y =  np.mgrid[np.min(xlist):np.max(xlist):j, np.min(ylist):np.max(ylist):k]
    #     grid_x += off_j
    #     grid_y += off_k
    #     
    #     r = Table()
    #     s = Table()
    #     
    #     for i in np.arange(0,len(ylist)):
    #         s['x'] = xlist
    #         s['y'] = ylist[i]
    #         s['Sky_Brightness'] = magnitude_sky[:,i]
    #         tmp = s
    #         if not r:
    #             r = s
    #         else:    
    #             r = vstack([r, tmp])
    #     
    #     grid = griddata(np.transpose([r['x'],r['y']]), r['Sky_Brightness'], (grid_x, grid_y), method='linear')
    #     
        #### ^^^^^^^ Old method, interpolation
    
        #### New method, downsample by taking mean
        # The great Adam Ginsburg, ladies and gentlemen...
        # https://github.com/keflavich/image_registration/blob/master/image_registration/fft_tools/downsample.py#L11
        # To determine factor by which we want to downsample... Every 1 degree?
        fov=[magnitude_sky.shape[0]*resolution/3600.0,magnitude_sky.shape[1]*resolution/3600.0] # in degrees
        factor=np.int(np.floor(np.min(fov)))
        print 'Downsampling ',magnitude_sky.shape,' by factor of ',factor
        #factor=10 # downsample by a factor of 10. Could figure out what this should be based on resolution.
        ys,xs=magnitude_sky.shape
        crarr=magnitude_sky[:ys-(ys % int(factor)),:xs-(xs % int(factor))] # Crops a bit
        xcrarr=xarr[:ys-(ys % int(factor)),:xs-(xs % int(factor))]         # Crops a bit
        ycrarr=yarr[:ys-(ys % int(factor)),:xs-(xs % int(factor))]         # Crops a bit
        grid = np.nanmean(np.concatenate([[crarr[i::factor,j::factor] for i in range(factor)] for j in range(factor)]),axis=0)
        grid_x = np.nanmean(np.concatenate([[xcrarr[i::factor,j::factor] for i in range(factor)] for j in range(factor)]),axis=0)
        grid_y = np.nanmean(np.concatenate([[ycrarr[i::factor,j::factor] for i in range(factor)] for j in range(factor)]),axis=0)
    
        # New plot, same colorscale
        im1=ax[1].imshow(grid,cmap=cmap,norm=norm,aspect='equal')
    
        # Flip?
        positions_wcs = wcs.all_pix2world(grid_y,grid_x,1) #Investigate about distortion corrections.
        c = coords.SkyCoord(positions_wcs[0],positions_wcs[1], unit='deg')
        cAltAz = c.transform_to(coords.AltAz(obstime = time, location = loc))
    
        # If the max and min are very close to 360 and 0, respectively, then we need to manually 
        #   set the axes.
        azarr=np.array(cAltAz.az)
        if np.min(cAltAz.az)<5.0*u.deg and np.max(cAltAz.az)>355.0*u.deg:
            azarr[np.where(azarr>180.0)]=azarr[np.where(azarr>180.0)]-360.0
        # We'll encounter the same thing for images that look at zenith. 
        
        # If values close to 90 are included, then we need to manually set the axes.
        #yscale=[np.min(cAltAz.alt),np.max(cAltAz.alt)]
    
        im2=ax[2].pcolor(azarr,np.array(cAltAz.alt),grid,cmap=cmap,norm=norm)
        plt.savefig(dir+'/skybrightness/image_'+f.split('.')[0]+'.png')
    
    #    3d) In the end, create a table (astropy table?) that has altitude, azimuth, 
    #        and sky brightness. 

        t = Table()
        s = Table()
    
        for i in np.arange(0,len(cAltAz.alt[:,0])):
            s['Altitude'] = cAltAz.alt[i,:]
            s['Azimuth'] = cAltAz.az[i,:]
            s['Sky_Brightness'] = grid[i,:] # Need to transpose?
            tmp = s
            if not t:
                t = s
            else:    
                t = vstack([t, tmp])

    #    3e) For now, I suggest writing this table as a text/ASCII file, perhaps, so we 
    #        keep a copy saved for future imaging?

        t.write(dir+'/skybrightness/'+f.split('.')[0]+'.txt', format='ascii.fixed_width',overwrite=True)
        print 'Table created for',f #Put table number

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

# Take column with nan sky brightness out of the data table.
numpy_skybrightness = np.array(data['Sky_Brightness'])
skybrightness_nonan = np.where(~np.isnan(numpy_skybrightness))
data = data[skybrightness_nonan]

print data

# 4) Create a polar plot of sky magnitude. (The Duriscoe plots that are circle-shaped.)
#    In the case of our polar plot, the angle is the azimuth, and the radius 
#    is the zenith angle (90-altitude). See an example polar plot in 
#    plan_targets.py (but ignore the polygon/shape drawing stuff). 
#    Choose an appropriate colorbar, and be sure to include a labeled
#    colorbar on the image.
#    Search around for matplotlib help on this!

#azimuths = np.array(data['Azimuth'])
#altitudes = np.array(data['Altitude'])
#values= np.array(data['Sky_Brightness'])

# Plot locations of altitudes, azimuths.
fig = plt.figure(num=1)
plt.clf()
ax = fig.add_subplot(111,polar='True')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
nfiles=len(np.unique(data['file']))
#scatcolors = cm.rainbow(np.linspace(0, 1, nfiles))
scatcolors=['#d7191c','#fdae61','#ffffbf','#abd9e9','#abd939','#2c7bb6']
for i,f in enumerate(np.unique(data['file'])):
    inds=np.where(data['file']==f)
    plt.scatter(data['Azimuth'][inds]*np.pi/180.0,90.0-data['Altitude'][inds],s=1,alpha=0.5,color=scatcolors[i % len(scatcolors)])
plt.ylim(0,90)
plt.title(dir+' '+str(len(data))+' downsampled points')
plt.savefig(dir+'/magnitude_map_coverage.png')

# Need to do another regridding to get this on a regular grid?
# If do the following, get a 990 x 990 array. Won't be completely filled; just diagonals.
# theta,rad=np.meshgrid(azimuths,90.0-altitudes,indexing='ij')

# Find where the azimuth jumps by at least one degree. Nope, can't do this, too continous.
#azimuths.sort()
#result=np.split(azimuths,np.where(np.diff(azimuths)>1)[0]+1)
#theta_list=[np.mean(x) for x in result]
#altitudes.sort()
#result=np.split(altitudes,np.where(np.diff(altitudes)>1)[0]+1)
#rad_list=[90.0-np.mean(x) for x in result] # Zenith angle is our radius for polar plot.

# Just create a list for new grid.
theta_list=np.arange(0,360,1) # Note we won't have equal area projection.
rad_list=np.arange(90.0-np.max(data['Altitude']),90.0-np.min(data['Altitude']),1)
theta,rad=np.meshgrid(theta_list,rad_list,indexing='ij')

values_2d = griddata(np.transpose([data['Azimuth'],90.0-data['Altitude']]), 
    data['Sky_Brightness'], (theta, rad), method='linear') # None, nans on the edges are okay.

# UNITS OF THETA ARE IN RADIANS, GENIUS!
theta=theta*np.pi/180.0

#x, y = r*np.cos(theta), r*np.sin(theta)
# Define colormap
cmap=mpl.cm.CMRmap_r # This is similar to Duriscoe... check it out? Not quite the same.
# They go pink to yellow to dark, ours is the opposite. Can look for more colormaps.
norm = mpl.colors.Normalize(vmin=np.min(values_2d[~np.isnan(values_2d)]),vmax=np.max(values_2d[~np.isnan(values_2d)]))
# Use the same normalization each time? 
#norm = mpl.colors.Normalize(vmin=17.5,vmax=22.25) Ours not nearly as dark.

fig = plt.figure(num=2)
plt.clf()
ax = fig.add_subplot(111,polar='True')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_ylim(0,90)

#plt.contour(theta,rad,values_2d)
#plt.imshow(theta,rad,values_2d)
#plt.tricontour(data['Azimuth']*np.pi/180.0,90.0-data['Altitude'],data['Sky_Brightness'],
#    norm=norm,cmap=cmap,fill=True)

plt.pcolormesh(theta,rad,values_2d,norm=norm,cmap=cmap)
plt.colorbar()
plt.title(dir)
plt.savefig(dir+'/magnitude_map.png')
#plt.show()


# 5) Optional, also create a plot displaying the same information
#    in the same projection as Duriscoe. (Hammer-Aitoff Projection?)

theta_2 = theta.copy()
theta_2[np.where(theta>np.pi)] -= 2*np.pi

# Sort all by theta values. Start negative, up to zero, and to positive.
inds=np.argsort(theta_2[:,0])
theta_2=theta_2[inds,:]
rad=rad[inds,:]
values_2d=values_2d[inds,:]

plt.figure()
plt.subplot(111,projection='mollweide')
plt.grid(True)
plt.pcolormesh(theta_2,(90.0-rad)*np.pi/180.0,values_2d,norm=norm,cmap=cmap)
plt.colorbar()
plt.savefig(dir+'/magnitude_map_mollweide.png')
