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

See the main readme file in previous directory for more info.