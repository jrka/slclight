# Perform aperture photometry on standard stars in our lightdome images.
# Requires the images to include WCS information (RA, Dec information),
# e.g. from astrometry.net.
#
# Example Usage:
# python lightdome_photometry.py lightdome_timpanogos
#
# Modification History
# 2018-06-28 JRK: Created file.
# 2018-06-29 JRK: Can now query Vizier and compare catalog locations
#                 to sources from astrometry.net. Not finished.

# Clean this up later.
from startup_slclight import *
from ccdproc import CCDData
from ccdproc import ImageFileCollection
from photutils import CircularAperture,CircularAnnulus
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.wcs import WCS
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.table import vstack, Table, join
from astroquery.vizier import Vizier

# if len(sys.argv)!=2:
#    print('Usage:\npython lightdome_process.py [lightdome_NAME folder] \n')
#    exit()
# dir='./data/'+sys.argv[2]+'/'
dir='./data/lightdome_timpanogos/'

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
    hdu.close()
    print 'astrometry.net found ',len(positions[0]),' sources.'
    # Load in the image file
    hdu=fits.open(dir+'/astrometry/'+f)
    wcs=WCS(hdu[0].header)
    data=hdu[0].data
    # Perhaps we should narrow down to our standard stars now at this point,
    # since we don't need thousands per image...
    # but subtracting these might be nice for background estimation.
    # compare final results of sky brightness with and without doing so.
    
    # Use astroquery to get SDSS magnitudes. Is this the right origin? 0 or 1?
   # positions_wcs=wcs.wcs_pix2world(positions[0],positions[1],1)
    # An example with ib030.fits: the first position is 594.84, 59.88. In DS9, 
    # that is, in degrees, about 222.63, 74.184. That is quite close to what is
    # produced.
    
    # This queries SDSS, but this gets galaxies, not stars? 
    # Try VizieR next?
#     match=False
#     for i in np.arange(len(positions_wcs[0])):
#         co=coords.SkyCoord(positions_wcs[0][i],positions_wcs[1][i],unit="deg")
#         result=SDSS.query_crossid(co,photoobj_fields=['psfMag_z','psfMagerr_z','type'])
#         if result:
#             result['ourid']=i
#             result['x']=positions[0][i]
#             result['y']=positions[1][i]
#             print i
#             print result
#         if result and not match: # First time we find a match.
#             match=True
#             standard_tbl=result
#         elif result: # Subsequent times we find a match.
#         	standard_tbl=vstack([standard_tbl,result]) # Cant' get this to work
        	
        	
    # Try Vizier, catalogs: II/183A is Landolt 1992, 
    # J/AJ/146/131, Landolt 2013 with +50 deg declination.
    # I/239/hip_main is Hipparcos and Tycho catalogs
    Vizier.ROW_LIMIT=-1  # Row limit is 50 by default.
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
    result=Vizier.query_region(centercoord,width=[sep1*2,sep2*2],
         catalog=['II/183A','J/AJ/146/131','I/239/hip_main'])
    
        
    result=query_stars(centercoord,width=[sep1*2,sep2*2])
    
    if not result:
        print 'No standard stars found for ',f,', moving on to next file'
        continue # Skips the rest of this file, continues the loop to the next file.
    
    print result
    
    # Limit to tables II/183A/table2 or J/AJ/146/131/standards
    # J/AJ/146/131/standards colnames: __Vmag_ and RAJ2000, DEJ2000 (sexigesimal)
    # II/183A/table2  colnames: Vmag, RAJ2000, DEJ2000 (sexigesimal). Error given, but not in other Landolt.
    # Create apertures for these standard stars, so we can visualize? 
    # How to determine which of our sources correspond?
    positions_wcs=wcs.wcs_pix2world(positions[0],positions[1],1)
    co=coords.SkyCoord(positions_wcs[0],positions_wcs[1],unit="deg")
    
    catalog=coords.SkyCoord(ra=result['RA'],dec=result['DEC'],unit=u.deg)
    
    
    #http://docs.astropy.org/en/stable/coordinates/matchsep.html
    # Now idx are indices into catalog that are the closest objects to each of the 
    # coordinates in c, d2d are the on-sky distances between them, and d3d are the 
    # 3-dimensional distances. Because coordinate objects support indexing, 
    # idx enables easy access to the matched set of coordinates in the catalog:
    # (3D distances only if distance to objects were input; they were not.)
    idx, d2d, d3d = co.match_to_catalog_sky(catalog)
    
    # This is about 3 pixels to be off by for NPS data; IN FUTURE
    # CALCULATE FROM WCS OR WHATEVER, PLATE SCALE.
    ##### PICK UP FROM HERE MONDAY!!!!!
    
    print len(catalog),' sources in catalog'
    
    if len(np.where(d2d.to(u.arcsec)<270.0*u.arcsec)[0])==0:
        print 'No matched sources within 270 arcseconds for ',f
        continue
        
    # Limit ourselves to these matches then do the photometry!
    # idx is the length of our positions data (our sources)
    
    
    
    print 'MADE IT'
    stop
    apertures_catalog=CircularAperture(catalog.to_pixel(wcs),r=5)
    
    # Create the apertures of sources, with local background subtraction
    # NOTE FOR FUTURE REFERENCE: See Mommert 2017 for a discussion of
    # finding the optimum aperture radius using a curve of growth analysis
    # to maximize flux AND signal-to-noise ratio simultaneously. Also Howell 2000.
    apertures=CircularAperture(positions,r=3) # pixels
    annulus_apertures=CircularAnnulus(positions,r_in=6., r_out=8.)
    # Plot them?
    norm=ImageNormalize(stretch=LogStretch())
    ax=plt.subplot()
    ax.imshow(data,cmap='Greys',origin='lower',norm=norm)
    apertures.plot(color='blue',lw=1.5,alpha=0.5)
    apertures_catalog.plot(color='red',lw=1.5,alpha=0.5)
    
    
    # Can make it nice like http://docs.astropy.org/en/stable/visualization/wcsaxes/,
    #   save the files for future use.
    
    # Do aperture_photometry. Local background subtraction version,
    # http://photutils.readthedocs.io/en/stable/aperture.html#local-background-subtraction
    apers=[apertures,annulus_apertures]
    phot_table=aperture_photometry(data,apers)
    print(phot_table)
    bkg_mean=phot_table['aperture_sum_1']/annulus_apertures.area()
    bkg_sum=bkg_mean*apertures.area()
    final_sum=phot_table['aperture_sum_0']-bkg_sum
    phot_table['residual_aperture_sum']=final_sum
    print(phottable['residual_aperture_sum'])
    # Put this in instrumental magnitude
    
    # Want to output a list of the stars, their catalog magnitudes, airmass, 
    # total counts, 
    
    hdu.close()