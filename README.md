# slclight

Create a full-sky mosaic of sky brightness (mag/arcsec^2) given all-sky 
CCD imagery. Contact Julia Kamenetzky at Westminster College for inquiries.

## Background

The goal of this work is to follow the analysis procedure of Duriscoe et al. 2007 
(http://iopscience.iop.org/article/10.1086/512069/meta)
to create all-sky maps of sky brightness (in magnitudes per square arcsecond) 
given a set of raw CCD images (fits files), using only freely available
data analysis tools (python, astropy) and minimal observing equipment.

This code is made freely available to other researchers, especially those 
at small colleges. The procedures require the use of CCD techniques that
are typically part of the undergraduate astronomy curriculum (calibration of 
CCD images, aperture photometry, determination of zero-point magnitude, 
familiarity with alt-az coordinate system).

At this time, the code is presented "as-is" for other researchers to use
as possible templates for their own work. It is designed to work following 
the procedure described below, but may currently contain bugs and 
content specific to our work. It is not actively maintained during the 
academic year.

## Citing

In addition to the citation for the original analysis procedure given
above (Durscoe et al. 2007), if using this Python code in 
whole or in part, please cite this DOI: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1471633.svg)](https://doi.org/10.5281/zenodo.1471633)


## Future Goals, Ongoing Work
- Calculate and present basic statistics of resulting magnitude maps 
(surface brightness at zenith, brightest location, and darkest location, 
total integrated sky background, and total integrated brightness of 
a city’s light dome).
- Subtraction of natural sky background (ala Duriscoe 2014, 
http://iopscience.iop.org/article/10.1086/673888/meta)
- More robust, user friendly interface. Admittedly, the initial version 
of the code is accurately described by this XKCD comic: https://xkcd.com/2054/

## Python Requirements 

Assumes you have some of the standard packages installed with e.g. Anaconda 
(numpy, scipy,matplotlib...) as well as astropy and 3 related astropy 
packages: ccdproc, photutils, astroquery.

May need to do: 

`conda install -c astropy ccdproc`

`conda install -c astropy photutils`

`conda install -c astropy astroquery`

## Observing Procedure

CCD image files must contain an accurate timestamp and an identification
as dark, bias, flat or light frames.
No other information is required in the fits headers (e.g. lat, long, elev).

The routine `plan_targets.py` can create a map of the sky broken up
into sections depending on the size of your CCD field of view (which will
depend on the combination of camera and lens). Good coverage can be 
achieved manually using an inexpensive alt-az mount if degrees are marked
on the mount. If using a computerized mount system that requires RA and Dec
for pointing, use `altaz_to_radec.py'.

Somehow, you want to observe CCD images of large portions of the sky,
with as much coverage as reasonably possible, if not the whole sky.

## Setting up the data

1. For each set of observations, create a folder in /data/, 
for example, slclight/data/lightdome\_timpanogos. 
In the examples that follow, "lightdome_timpanogos" would be DIRNAME.

2. In the folder, create a metadata.txt file with important information 
about the observing site and conditions. See example file.

3. Place all calibration frames (bias, dark, flat) in a subfolder "calib"

4. Place all sky frames in a subfolder "lights"

5. Other subfolders shoudl be created automatically: 
"final," "astrometry," "photometry," and "skybrightness"

## Doing the processing

Each script can be called from the command line, and should be done
from the main slclight directory.

1. `python lightdome_process DIRNAME` will dark-subtract and flat-field the images.
New set of .fit files will be in directory "final".

2. Use astrometry.net to plate-solve the images. This can be downloaded and run
locally, but can be tricky to install. Instead, we used the web interface.
Upload all the .fit files in directory "final" to an album on astrometry.net
Provide some Advanced Settings using the web interface before uploading
 * Very wide field for images > 10 deg wide
 * Input the RA and Dec from the fits file header (in degrees), if applicable, radius 30 deg. In general,
the horizon images can be solved by astrometry.net blind, but inputting the RA and Dec will
significantly decrease the amount of time it takes.
 * Don't need to downsample
 * Batch upload as .tar or .tar.gz file doesn't seem to be working at the moment.
 
3. Download the "new-image.fits" files, one at time, into a subfolder "astrometry", 
renaming each one to match its original filename (e.g. ib030.fits). 
These new files have "World Coordinate System" (WCS) information, which allows
 us to identify stars by RA and Dec. Also download the "axy.fits" files, an rename as 
 e.g. ib030_axy.fits. The prefix of image files and \_axy files must match.
    
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
and find a linear best fit to the result. Outliers will be excludeed.
The intercept is the instrumental zero point magnitude.

6. `python lightdome_magnitude_map.py DIRNAME` will read in each image, use the 
zeropoint magnitude to calculate sky brightness in mag/arcsec^2 for each pixel, 
and then downsample the image to 1 deg x 1 deg sections. All images will then 
be combined and presented in an alt-az map. Currently, color scale is hard-coded.

## Other Scripts Available

1. `plan_targets.py`: enter width and height of CCD FOV in arcseconds, produce a map
and list of target locations to map the sky. Not a lot of overlap is currently built in.

2. `altaz_to_radec.py`: Given an intended altitude and azimuth for pointing, return 
the current RA and Dec at that location given your location and the current date/time.
Helpful if observing with a computerized mount control that wants RA and Dec input.
