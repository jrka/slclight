Instructions for the "data" folder in slclight

Each site/evening's observations will be collected in one folder, 
DIRNAME.

Within each folder, you'll setup the following:
- metadata.txt. Follow the template to record observing information.
- lights, a directory for the raw science images.
- calib, a directory for the calibration images (bias, darks, flats)

.fits image files are not included in the repository, but example
images are available upon request.

After following procedures below, you'll produce 
- final, a directory for the calibrated science images.
- astrometry, a directory for results from astrometry.net
- photometry, a directory which includes results from lightdome_photometry.py
- skybrightness, a directory which includes one text file and one .png file per image
- zeropoint.png and zeropoint_information.txt
- magnitude_map_coverage.png, magnitude_map.png, magnitude_map_mollweide.png
- Future: results.txt, a file which summarizes the results, like the "photometric indicators"
  on the NPS Google Earth kml file. 
  
.fit or .fits image files should not be included in the repo.

See the main readme file in previous directory for more info.
