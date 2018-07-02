Instructions for the "data" folder in slclight

Each site/evening's observations will be collected in one folder, 
with the name lightdome_NAME, where "Name" is some identifying name.

Within each folder, you'll setup the following:
- metadata.txt. Follow the template to record observing information.
- lights, a directory for the raw science images.
- calib, a directory for the calibration images (bias, darks, flats)

After following procedures below, you'll produce 
- final, a directory for the calibrated science images.
- astrometry, a directory for results from astrometry.net
- photometry, a directory which includes results from lightdome_photometry.py
- results.txt, a file which summarizes the results, like the "photometric indicators"
  on the NPS Google Earth kml file. 
  
.fit or .fits image files should not be included in the repo.


PROCEDURES AFTER SETTING UP THE 3 INITIAL THINGS ABOVE.
1) Run lightdome_process to dark-subtract and flat-field the images.
	Creates new set of .fit files in directory "final"
2) Upload all the .fit files in directory "final" to an album on astrometry.net
	Provide some Advanced Settings using the web interface before uploading
		Set public visibility to "no" for now
		Very wide field for images > 10 deg wide
		Input the RA and Dec from the fits file header, if applicable, radius 30 deg
		Don't downsample
	Batch upload as .tar or .tar.gz file doesn't seem to be working for now.
3) Download the "new-image.fits" files, one at time, into a subfolder "astrometry",
    renaming each one to match its original filename (e.g. ib030.fits).
    These new files have "World Coordinate System" (WCS) information, which allows
    us to identify stars by RA and Dec. 
    Also download the "axy.fits" files, an rename as e.g. ib030_axy.fits