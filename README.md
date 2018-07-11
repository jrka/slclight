# slclight

Create a full-sky mosaic of sky brightness (mag/arcsec^2) given all-sky 
CCD imagery.

Assumes you have some of the standard packages installed with e.g. Anaconda (numpy, scipy,
matplotlib...) as well as astropy and 3 related astropy packages: ccdproc, photutils, astroquery.

May need to do: 

`conda install -c astropy ccdproc`

`conda install -c astropy photutils`

`conda install -c astropy astroquery`

## Setting up the data

1. For each set of observations, create a folder in /data/, 
for example, slclight/data/lightdome\_timpanogos. 
In the examples that follow, "lightdome_timpanogos" would be DIRNAME.

2. In the folder, create a metadata.txt file with important information 
about the observing site and conditions. See example file.

3. Place all calibration frames (bias, dark, flat)in a subfolder "calib"

4. Place all sky frames in a subfolder "lights"

5. Also create subfolders "final," "astrometry," and "photometry."

## Doing the processing

Each script can be called from the command line, and should be done
from the main slclight directory.

1. `python lightdome_process DIRNAME` will dark-subtract and flat-field the images.
New set of .fit files will be in directory "final".

2. Use astrometry.net to plate-solve the images. This can be downloaded and run
locally, but can be tricky to install. Instead, we used the web interface.
Upload all the .fit files in directory "final" to an album on astrometry.net
Provide some Advanced Settings using the web interface before uploading
 * Set public visibility to "no" for now
 * Very wide field for images > 10 deg wide
 * Input the RA and Dec from the fits file header (in degrees), if applicable, radius 30 deg. In general,
the horizon images can be solved by astrometry.net blind, but inputting the RA and Dec will
significantly decrease the amount of time it takes.
 * Don't downsample
 * Batch upload as .tar or .tar.gz file doesn't seem to be working for now.
 
3. Download the "new-image.fits" files, one at time, into a subfolder "astrometry", 
renaming each one to match its original filename (e.g. ib030.fits). 
These new files have "World Coordinate System" (WCS) information, which allows
 us to identify stars by RA and Dec. Also download the "axy.fits" files, an rename as 
 e.g. ib030_axy.fits
    
4. `python lightdome_photometry DIRNAME` will do the following:
 * Read in the sources in the image detected by astrometry.net
 * Limit those numerous sources to only those of high signal/noise.
 * Use astroquery and the WCS information to query the magnitudes of known 
 stellar sources (from Hipparchos and Tycho) in the RA/Dec box contained in the image.
 * Match those sources to those corresponding to within 1 pixel of their 
 expected location (by default). 
 * Model the FWHM of the image sources.
 * Perform aperture photometry on the image sources using a radius of 3 times the 
 median FWHM.
 * Save the aperture photometry results, along with catalog source information, into
 a text file in "photometry" folder. Also saves a png file with the image and sources identified,
 and a PDF with the FWHM models.
 
5. `python lightdome_calc_zeropoint.py DIRNAME` will read in all the text files created
in the previous step, plot (catalog magnitude + instrumental magnitude) vs. airmass,
and find a linear best fit to the result. The intercept is the instrumental zero point magnitude.

6. `python lightdome_magnitude_map.py DIRNAME`

## Other Scripts Available

1. `plan_targets.py`

2. `altaz_to_radec.py`

