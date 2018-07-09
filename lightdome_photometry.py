# Perform aperture photometry on standard stars in our lightdome images.
# Requires the images to include WCS information (RA, Dec information),
# e.g. from astrometry.net.
#
# Example Usage:
# python lightdome_photometry.py lightdome_timpanogos
# python lightdome_photometry.py lightdome_timpanogos 2
# Default number of pixels to allow matches with catalogs is 1; set a number 
# after the directory to override this default.
#
# The CCD images must have the following headers:
#  COMMENT with scale in arcsec/pixel (from astrometry.net)
# 'DATE-OBS' in UTC, 'LONGITUD' and 'LATITUDE' in degrees,
# 'ELEVATIO' in meters, 'EXPTIME' in seconds 
#  If 'LONGITUD','LATITUDE',and 'ELEVATIO' are not present,
#  will look for those in the metadata.txt file.
#
# OUTPUTS:
# For each .fits file in the /astrometry directory, a fixed-width .txt
# file will be produced and in /photometry. The columns are as follows:
#  RA: RA of standard star, in degrees (J2000)
#  DEC: Declination of standard star in degrees (J2000)
#  Vmag: V-band magnitude of standard star, in magnitudes.
#  Source: Code for the source of the previous 3 columns in Vizier
#  idx: ID of the extracted source from astrometry.net that corresponds
#       to this row's standard star. Not useful beyond this. 
#  source_x: x-coordinate (pixel) of extracted source from astrometry.net
#  source_y: y-coordinate (pixel) of extracted source from astrometry.net
#  source_RA: conversion of x,y columns to RA using astrometry.net WCS
#  source_DEC: conversion of x,y columns to Dec using astrometry.net WCS
#  residual_aperture_sum: background-subtracted total sum (in DN) of photometry
#  instrmag: Instrumental magnitude, 2.5*log10(residual_aperture_sum / exposure time in s)
#  alt: altitude of star, in degrees (0 = horizon, 90 = zenith)
#
# Modification History
# 2018-07-02 JRK: Perform photometry, get alt/az from coordinates and other
#                 fits header information (location, elevation, date/time in UT).
#                 Write each result table to a file in /photometry.
#                 NOTE: Still need to figure out best aperture/annulus sizes, 
#                 add errors in to calculation, make plot nicer, 
#                 deal with MergeConflict warnings from Vizier results.
# 2018-06-29 JRK: Can now query Vizier and compare catalog locations
#                 to sources from astrometry.net. Not finished.
# 2018-06-28 JRK: Created file.


# Clean this up later.
from startup_slclight import *
from ccdproc import CCDData
from ccdproc import ImageFileCollection
from photutils import CircularAperture,CircularAnnulus,aperture_photometry
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.wcs import WCS
from astropy import coordinates as coords
from astropy.table import vstack, Table, join
from astropy.time import Time
from astropy.modeling import models, fitting

if len(sys.argv)<2:
    print('Usage:\npython lightdome_photometry.py lightdome_NAME_folder [pixels] \n')
    exit()
npix=np.float(sys.argv[2]) if len(sys.argv)==3 else 1.0
    
dirname=sys.argv[1]
# dirname='lightdome_timpanogos'
dir='./data/'+dirname+'/'

###################### READ IN FILE INFO

# Find all science images in the "astrometry" folder.
# Exclude the _axy files.
files=ImageFileCollection(dir+'/astrometry/',glob_exclude='*axy*')

###################### READ IN FILE INFO
for f in files.files:
    print f
    # Identify the sources. Use the ones from _axy.fits from astrometry.net
    hdu=fits.open(dir+'/astrometry/'+f.split('.')[0]+'_axy.fits')
    positions=[hdu[1].data['X'],hdu[1].data['Y']]
    source_peak=hdu[1].data['flux']
    hdu.close()
    print 'astrometry.net found ',len(positions[0]),' sources.'
    # Load in the image file
    hdu=fits.open(dir+'/astrometry/'+f)
    wcs=WCS(hdu[0].header)
    data=hdu[0].data
    # Also grab the pixel scale. This little routine is in startup_slclight.py.
    pixscale=platescale(hdu[0].header)
    # Perhaps we should narrow down to our standard stars now at this point,
    # since we don't need thousands per image...
    # but subtracting these might be nice for background estimation.
    # compare final results of sky brightness with and without doing so.
    
    # Plotting: Show ALL identified sources.
    norm=ImageNormalize(data=data,stretch=LogStretch())
    ax=plt.subplot()
    ax.imshow(data,cmap='Greys',origin='lower',norm=norm)
    apertures_allsources=CircularAperture(positions,r=10) # pixels
    apertures_allsources.plot(color='blue',lw=1.5,alpha=0.5)
    
    # Try Vizier, catalogs: II/183A is Landolt 1992, 
    # J/AJ/146/131, Landolt 2013 with +50 deg declination.
    # I/239/hip_main is Hipparcos and Tycho catalogs
    # There must be a more clever way to figure out the center coordinate
    # and width of image, but I apparently can't figure it out right now.
    center=wcs.wcs_pix2world(np.shape(data)[0]/2.0,np.shape(data)[1]/2.0,1)
    centercoord=coords.SkyCoord(ra=center[0],dec=center[1],unit="deg")
    corner=wcs.wcs_pix2world(0,0,1)
    sep1=centercoord.separation(coords.SkyCoord(ra=corner[0],dec=center[1],unit="deg"))
    sep2=centercoord.separation(coords.SkyCoord(ra=center[0],dec=corner[1],unit="deg"))
    # I don't seem to find any in the region for ib030.fits, 
    # which is centered on 186. 2-1, 70.237,
    # [<Angle 23.892974958515765 deg>, <Angle 19.945039907759938 deg>]
    # The full catalog pulls up find, and i do get results if I look in a radius of 100 degrees,
    # so I think the issue is the Dec is too high.
    # Yes, the min dec is -46 degrees, the max dec is +16 degrees
    # Add in J/AJ/146/131, which is +50 deg declination.
    # result=Vizier.query_region(centercoord,width=[sep1*2,sep2*2],
    #     catalog=['II/183A','J/AJ/146/131','I/239/hip_main'])
    
    result=query_stars(centercoord,width=[sep1*2,sep2*2])
    
    if not result:
        print 'No standard stars found for ',f,', moving on to next file'
        continue # Skips the rest of this file, continues the loop to the next file.
    
    # Limit to tables II/183A/table2 or J/AJ/146/131/standards
    # J/AJ/146/131/standards colnames: __Vmag_ and RAJ2000, DEJ2000 (sexigesimal)
    # II/183A/table2  colnames: Vmag, RAJ2000, DEJ2000 (sexigesimal). Error given, but not in other Landolt.
    # Create apertures for these standard stars, so we can visualize? 
    # How to determine which of our sources correspond?
    positions_wcs=wcs.wcs_pix2world(positions[0],positions[1],1)
    co=coords.SkyCoord(positions_wcs[0],positions_wcs[1],unit="deg")
    catalog=coords.SkyCoord(ra=result['RA'],dec=result['DEC'],unit=u.deg)
    print len(catalog),' sources in catalog'
    
    # Plotting: Add all these catalog sources to the plot
    apertures_catalog=CircularAperture(catalog.to_pixel(wcs),r=5)
    apertures_catalog.plot(color='red',lw=1.5,alpha=0.5)    
    # Save the image
    plt.savefig(dir+'/photometry/photometry_'+f.split('.')[0]+'.png') 
    
    #http://docs.astropy.org/en/stable/coordinates/matchsep.html
    # Now idx are indices into catalog that are the closest objects to each of the 
    # coordinates in c, d2d are the on-sky distances between them, and d3d are the 
    # 3-dimensional distances. Because coordinate objects support indexing, 
    # idx enables easy access to the matched set of coordinates in the catalog:
    # (3D distances only if distance to objects were input; they were not.)
    idx, d2d, d3d = catalog.match_to_catalog_sky(co)
    
    # Get the ones that are within 1 pixels of center.
    indnpix=np.where(d2d.to(u.arcsec)<(npix*pixscale)*u.arcsec)[0]
    if len(indnpix)==0:
        print 'No matched sources within ',npix*pixscale,' arcseconds for ',f
        continue
    else:
        # Cut down our "result" table to only those that match
        result['idx']=idx
        result=result[indnpix]
        catalog=catalog[indnpix]
        print 'Using ',len(result),' matches for photometry.'
        
    # Limit ourselves to these matches then do the photometry!
    # And now, restrict our list of sources to only those that also 
    # now appear in the "results" table. np.unique(result['idx'])
    result['source_x']=positions[0][result['idx']]
    result['source_y']=positions[1][result['idx']]
    result['source_peak']=source_peak[result['idx']]
    result['source_RA']=positions_wcs[0][result['idx']]*u.deg
    result['source_DEC']=positions_wcs[1][result['idx']]*u.deg
    
    # Further narrow down for reasonable values of source peak.
    # Keep only ones that are bright stars (> 5 sigma from background)
    # And likely not saturated. How to know? Assume the max pixel is saturated?
    # Allow user to set his later? THIS IS NOT RIGHT, BECAUSE FLUX IS ALREADY BG SUBTRACTED.
    idbright=np.where(np.logical_and(result['source_peak']>np.median(data)+5.0*np.std(data),
        result['source_peak']<np.max(data)*0.9))[0]
    print 'With brightness and saturation cutoff, using ',len(idbright),' for photometry.'
    result=result[idbright]
    catalog=catalog[idbright]

    # Overplot with a new color.
    apertures_catalog=CircularAperture(catalog.to_pixel(wcs),r=5)
    apertures_catalog.plot(color='yellow',lw=1.5,alpha=0.5)     
    # Note, this reveals that our 1 pixel range is actually  quite restrictive.
    # I see, by eye, many possible matches.       
    
    
    # Determine the FWHM of the sources, so that we can choose an 
    # appropriate radius for the aperture. 
    # http://docs.astropy.org/en/stable/modeling/
    # Better way to do this than looping?
    fit_g=fitting.LevMarLSQFitter()
#    g=models.Gaussian2D(amplitude=result['source_peak'],x_mean=result['source_x'],
#        y_mean=result['source_y'])
    result['source_fwhm']=1.0
    for i, row in enumerate(result):
        g_init=models.Gaussian2D(row['source_peak'],row['source_x'],row['source_y'])
        subset_inds=[np.rint(row['source_y']),np.rint(row['source_x'])]
        yi,xi=np.indices(data.shape)
        xi=xi[subset_inds[0]-5:subset_inds[0]+5,subset_inds[1]-5:subset_inds[1]+5]
        yi=yi[subset_inds[0]-5:subset_inds[0]+5,subset_inds[1]-5:subset_inds[1]+5]
        data_sub=data[subset_inds[0]-5:subset_inds[0]+5,subset_inds[1]-5:subset_inds[1]+5].copy().astype('float')
        data_sub-=np.median(data)
        g_result=fit_g(g_init,xi,yi,data_sub)
            
        # For Diagnostics
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(xi,yi,data_sub,color='blue',label='Data')
        ax.plot_wireframe(xi,yi,g_result(xi,yi),color='red',label='Model')
        
        # Get FWHM by taking average of x and y standard deviations, multiply by 2.35
        row['source_fwhm']=2.35*np.mean([g_result.x_stddev.value,g_result.y_stddev.value])

    print 'Source FWHM, Min: ',np.min(result['source_fwhm']),', Max: ',np.max(result['source_fwhm'])
    print 'Using average: ',np.mean(result['source_fwhm'])
    apr=0.6731*np.mean(result['source_fwhm'])
    # Create the apertures of sources, with local background subtraction
    # NOTE FOR FUTURE REFERENCE: See Mommert 2017 for a discussion of
    # finding the optimum aperture radius using a curve of growth analysis
    # to maximize flux AND signal-to-noise ratio simultaneously. Also Howell 2000.
    # Optimum aperture radius as 0.6731 FWHM?
    # compare this to what I was doing, using 3, 6, and 8 for circle, r_in, r_out.
    apertures=CircularAperture((result['source_x'],result['source_y']),
        r=apr) # pixels
    annulus_apertures=CircularAnnulus((result['source_x'],result['source_y']),
       r_in=apr+1.0, r_out=2.0*apr+1.0)
    
    # Can make it nice like http://docs.astropy.org/en/stable/visualization/wcsaxes/,
    #   save the files for future use.
    
    # Do aperture_photometry. Local background subtraction version,
    # http://photutils.readthedocs.io/en/stable/aperture.html#local-background-subtraction
    apers=[apertures,annulus_apertures]
    phot_table=aperture_photometry(data,apers)
    bkg_mean=phot_table['aperture_sum_1']/annulus_apertures.area()
    bkg_sum=bkg_mean*apertures.area()
    final_sum=phot_table['aperture_sum_0']-bkg_sum
    phot_table['residual_aperture_sum']=final_sum
    # Put this in instrumental magnitude.
    phot_table['instrmag']=2.5*np.log10(phot_table['residual_aperture_sum']/hdu[0].header['EXPTIME'])
    
    # Add this to the results table. 
    result['residual_aperture_sum']=phot_table['residual_aperture_sum']*u.ct
    result['instrmag']=phot_table['instrmag']*u.mag
    
    # Getting altitude --> zenith angle --> airmass. 
    c = coords.SkyCoord(result['source_RA'],result['source_DEC'], frame='icrs')
    # 2018-07-02: The below works IF your headers use this convention!
    # Read in the user-provided metadata file instead.
    #loc = coords.EarthLocation(lat = hdu[0].header['LATITUDE']*u.deg, 
    #    lon = hdu[0].header['LONGITUD']*u.deg, height = hdu[0].header['ELEVATIO']*u.m)
    metad=read_metadata(dirname)
    loc = coords.EarthLocation(lat = metad['lat']*u.deg, 
        lon = metad['lon']*u.deg, height = metad['elev']*u.m)
    # Note we DO want to use the fits header for the date and time of the observation,
    # because the sky moves over time as the observations progress!
    time = Time(hdu[0].header['DATE-OBS'],scale='utc')
    cAltAz = c.transform_to(coords.AltAz(obstime = time, location = loc))
    # Add altitude result table.
    result['alt']=cAltAz.alt.degree*u.deg
    
    # Write the results file.
    result.write(dir+'/photometry/'+f.split('.')[0]+'.txt',format='ascii.fixed_width')
    
    # Save the image
    plt.savefig(dir+'/photometry/photometry_'+f.split('.')[0]+'.png') 
    
    # Close the fits file.
    hdu.close()