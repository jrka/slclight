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
#  residual_aperture_sum_err: 1 sigma error associated with above.
#  instrmag: Instrumental magnitude, 2.5*log10(residual_aperture_sum / exposure time in s)
#  instrmag_err: 1 sigma error associated with above.
#  alt: altitude of star, in degrees (0 = horizon, 90 = zenith)
#
# Modification History
# 2018-07-16 JRK: Exclude possibly satured pixels (hardcode limit).
#                 Remove sources that are too close together (overlapping).
# 2018-07-12 JRK: If 'photometry' directory doesn't exist, make it.
#                 Set s/n cutoff in variable sncutoff (hardcode).
# 2018-07-11 JRK: Assuming each pixel has error approximately equal to sqrt(N),
#                 determine error in aperture sums and instrumental magnitudes.
# 2018-07-10 JRK: Require a peak signal/background ratio of 3.0 for sources
#                 detected by astrometry.net to be used for photometry.
#                 Plotting modifications; save FWHM modeling results to PDF.
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
from matplotlib.ticker import LogLocator
from setup_plotting import *
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv)<2:
    print('Usage:\npython lightdome_photometry.py lightdome_NAME_folder [pixels] \n')
    exit()
npix=np.float(sys.argv[2]) if len(sys.argv)==3 else 1.0 # npix = required 
# npix=1.0
#  distance between catalog and astromety.net sources to consider them a match.
    
dirname=sys.argv[1]
#dirname='lightdome_westminster'
dir='./data/'+dirname+'/'

# Some hard-coded values that user may wish to change.
sncutoff=3.0     # Minimum requirement for S/N of peak/background of course to include.
saturated=65535  # CCD limit; saturated. Peaks above 90% of this value will not be included.
npix_dup=5       # Number of pixels required between catalog sources.

###################### READ IN FILE INFO

# Find all science images in the "astrometry" folder.
# Exclude the _axy files.
files=ImageFileCollection(dir+'/astrometry/',glob_exclude='*axy*')

# If 'photometry' directory doesn't already exist, make it.
if not os.path.exists(dir+'/photometry/'): os.mkdir(dir+'/photometry/')

###################### READ IN FILE INFO
for f in files.files:
    print f,'==================================================='
    # Identify the sources. Use the ones from _axy.fits from astrometry.net
    try:
        hdu=fits.open(dir+'/astrometry/'+f.split('.')[0]+'_axy.fits')
        print 'astrometry.net found ',len(hdu[1].data['X']),' sources.'
    except:
        print 'Error finding or opening '+dir+'/astrometry/'+f.split('.')[0]+'_axy.fits.'
        continue
    
    
    # Do the cut here for sources that are very high signal to noise.
    # Note, this is just the signal of the peak pixel; not the total
    # integrated flux of the star. Use a S/N cut of 3.0 here.
    indsn=np.where(hdu[1].data['flux']/hdu[1].data['background']>sncutoff)
    positions=[hdu[1].data['X'][indsn],hdu[1].data['Y'][indsn]]
    source_peak=hdu[1].data['flux'][indsn]
    print 'Using a peak flux/background cutoff of ',str(sncutoff),' , now ',len(source_peak),' sources.'
    # Now exclude saturated pixels.
    indsat=np.where(source_peak<0.9*saturated)
    positions=[positions[0][indsat],positions[1][indsat]]
    source_peak=source_peak[indsat]
    print 'Excluding peaks above ',str(0.9*saturated),' , now ',len(source_peak),' sources.'    
    
    hdu.close()
    
    if len(source_peak)==0:
        print 'No sources in ',f,', moving on to next file.'
        continue
    
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
    fig,ax=plt.subplots(num=1)
    cax=ax.imshow(data,cmap='Greys',origin='lower',norm=norm)
    cbar=fig.colorbar(cax,ticks = LogLocator(subs=range(10)))
    apertures_allsources=CircularAperture(positions,r=10) # pixels
    apertures_allsources.plot(color='blue',lw=1.5,alpha=0.5,label='Our Sources')
    
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
    # See query_stars for more information.    
    result=query_stars(centercoord,width=[sep1*2,sep2*2])
    
    if not result:
        print 'No standard stars found for ',f,', moving on to next file'
        continue # Skips the rest of this file, continues the loop to the next file.
        
    # Put the RA and Dec from "result" into x,y position in image.
    # NOTE that catalog.to_pixel(wcs) seems to give more accurate results
    # then wcs.wcs_world2pix. 
  # Commented out by JRK for testing 7/17/18
  #  catalog=coords.SkyCoord(ra=result['RA'],dec=result['DEC'],unit=u.deg)
  #  result['x'],result['y']=catalog.to_pixel(wcs)
  #  apertures_catalog_2=CircularAperture([result['x'],result['y']],r=7)
  #  apertures_catalog_2.plot(color='green',lw=1.5,alpha=0.5)
    
    # Check for sources very close together (not distinguished in our images)
    # JRK Commented out for testing 7/17/18
    result=remove_close_sources(result,npix_dup*pixscale)
    if len(result)<1:
        print 'Only standard stars in image are too close together.'
        continue # Skips the rest of this file, continues the loop to the next file.        
    
    # JRK 7/17/18 Visualize the sources and the catalog sources.

    # Put these coordinates into WCS
    catalog=coords.SkyCoord(ra=result['RA'],dec=result['DEC'],unit=u.deg)
    print len(catalog),' sources in catalog'
    
    # Create apertures for these standard stars, so we can visualize? 
    # Plotting: Add all these catalog sources to the plot
    # We will often get thousands; only do the 100 brightest ones.
    apertures_catalog=CircularAperture(catalog[np.argsort(result['Vmag'])][0:100].to_pixel(wcs),r=5)
    apertures_catalog.plot(color='red',lw=1.5,alpha=0.5)    
    # Save the image
    plt.savefig(dir+'/photometry/photometry_'+f.split('.')[0]+'.png') 
    
    # Convert our sources to world coordinates.
    positions_wcs=wcs.wcs_pix2world(positions[0],positions[1],1)
    co=coords.SkyCoord(positions_wcs[0],positions_wcs[1],unit="deg")
    
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

    # Overplot with a new color.
    apertures_catalog=CircularAperture(catalog.to_pixel(wcs),r=5)
    apertures_catalog.plot(color='yellow',lw=1.5,alpha=0.5)     
    # Note, this reveals that our 1 pixel range is actually  quite restrictive.
    # I see, by eye, many possible matches.       
    # This generally limits us to the center of the image, 
    # where distortion is limited.
    
    # Save the image, clear for next file.
    plt.savefig(dir+'/photometry/photometry_'+f.split('.')[0]+'.png') 
    plt.clf()
    
    # Determine the FWHM of the sources, so that we can choose an 
    # appropriate radius for the aperture. 
    # http://docs.astropy.org/en/stable/modeling/
    # Better way to do this than looping?
    fit_g=fitting.LevMarLSQFitter()
    result['source_fwhm']=1.0
    # Plot this as a PDF figure.
    pdffile=dir+'/photometry/photometry_'+f.split('.')[0]+'_fwhm.pdf'
    pdf=PdfPages(pdffile)
    plt.clf()
    
    nb_plots_rows = 4 
    nb_plots_columns = 3
    nb_plots_per_page = nb_plots_rows*nb_plots_columns
    grid_size = (nb_plots_rows, nb_plots_columns)
    totplots=-1
    
    newlfd=latex_subplots(latexfd['full'],nb_plots_rows,nb_plots_columns,wspace=0,hspace=0)
    fig,axes=plt.subplots(nb_plots_rows,nb_plots_columns,figsize=(newlfd['figw'],newlfd['figh']),num=0)
    fig.subplots_adjust(top=newlfd['t'],right=newlfd['r'],bottom=newlfd['b'],left=newlfd['l'],wspace=newlfd['wspace'],hspace=newlfd['hspace'])

    for i, row in enumerate(result):
        g_init=models.Gaussian2D(row['source_peak'],row['source_x'],row['source_y'])
        subset_inds=[np.int(np.rint(row['source_y'])),np.int(np.rint(row['source_x']))]
        yi,xi=np.indices(data.shape)
        xi=xi[subset_inds[0]-5:subset_inds[0]+5,subset_inds[1]-5:subset_inds[1]+5]
        yi=yi[subset_inds[0]-5:subset_inds[0]+5,subset_inds[1]-5:subset_inds[1]+5]
        data_sub=data[subset_inds[0]-5:subset_inds[0]+5,subset_inds[1]-5:subset_inds[1]+5].copy().astype('float')
        data_sub-=np.median(data)
        g_result=fit_g(g_init,xi,yi,data_sub)
        
        totplots+=1
        if totplots % nb_plots_per_page==0:
            plt.clf()
            fig,axes=plt.subplots(nb_plots_rows,nb_plots_columns,figsize=(7,9),num=0)
            #fig.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.87,bottom=0.08) # will move right.
        axind=(np.int(np.floor((totplots % nb_plots_per_page)/np.float(nb_plots_columns))),  np.int(((totplots % nb_plots_per_page) % nb_plots_columns)))
        axes[axind].remove()
        axes[axind]=fig.add_subplot(nb_plots_rows,nb_plots_columns,totplots % nb_plots_per_page+1,projection='3d')
        axes[axind].plot_wireframe(xi,yi,data_sub,color='blue',label='Data')
        axes[axind].plot_wireframe(xi,yi,g_result(xi,yi),color='red',label='Model')
        axes[axind].set_title('FWHM = '+str(2.35*np.mean([g_result.x_stddev.value,g_result.y_stddev.value])))
        
        # Get FWHM by taking average of x and y standard deviations, multiply by 2.35
        row['source_fwhm']=2.35*np.mean([g_result.x_stddev.value,g_result.y_stddev.value])
    pdf.savefig()
    pdf.close()
    # print pdffile

    print 'Source FWHM, Min: ',np.min(result['source_fwhm']),', Max: ',np.max(result['source_fwhm'])
    print 'Using median: ',np.median(result['source_fwhm'])
    apr=3.0*np.mean(result['source_fwhm'])

    # Create the apertures of sources, with local background subtraction
    # NOTE FOR FUTURE REFERENCE: See Mommert 2017 for a discussion of
    # finding the optimum aperture radius using a curve of growth analysis
    # to maximize flux AND signal-to-noise ratio simultaneously. Also Howell 2000.
    # Optimum aperture radius as 0.6731 FWHM? Use 3.0 to be sure to get all flux.
    # compare this to what I was doing, using 3, 6, and 8 for circle, r_in, r_out.
    apertures=CircularAperture((result['source_x'],result['source_y']),
        r=apr) # pixels
    annulus_apertures=CircularAnnulus((result['source_x'],result['source_y']),
       r_in=apr+1.0, r_out=2.0*apr+1.0)
    
    # Can make it nice like http://docs.astropy.org/en/stable/visualization/wcsaxes/,
    #   save the files for future use.
    
    # Do aperture_photometry. Local background subtraction version,
    # http://photutils.readthedocs.io/en/stable/aperture.html#local-background-subtraction
    # For an estimate of the error, assume all errors are Poisson and that
    # the sky/source signal dominates the error. 
    apers=[apertures,annulus_apertures]
    phot_table=aperture_photometry(data,apers,error=np.sqrt(data.copy()))
    bkg_mean=phot_table['aperture_sum_1']/annulus_apertures.area()
    bkg_sum=bkg_mean*apertures.area()
    final_sum=phot_table['aperture_sum_0']-bkg_sum
    phot_table['residual_aperture_sum']=final_sum
    # Based on error propagation calculation by Nicole.
    phot_table['residual_aperture_sum_err']=np.sqrt(phot_table['aperture_sum_err_0']**2+(apertures.area()/annulus_apertures.area()*phot_table['aperture_sum_err_1'])**2)
    # Put this in instrumental magnitude.
    phot_table['instrmag']=2.5*np.log10(phot_table['residual_aperture_sum']/hdu[0].header['EXPTIME'])
    # Based on error propagation calculation by Nicole.
    phot_table['instrmag_err']=1.085*phot_table['residual_aperture_sum_err']/phot_table['residual_aperture_sum']
    
    # Add this to the results table. Add errors as well.
    result['residual_aperture_sum']=phot_table['residual_aperture_sum']*u.ct
    result['residual_aperture_sum_err']=phot_table['residual_aperture_sum_err']*u.ct
    result['instrmag']=phot_table['instrmag']*u.mag
    result['instrmag_err']=phot_table['instrmag_err']*u.mag
    
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
    result.write(dir+'/photometry/'+f.split('.')[0]+'.txt',
        format='ascii.fixed_width',overwrite=True)
    
    # Close the fits file.
    hdu.close()