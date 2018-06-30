# The purpose of this file is to load in all of the packages and 
# variables that we are always using to complete this project.
# By loading in this file at the start of each script, we reduce
# the need to copy/paste and continually update each individual
# script.
#
# Revision History
# 2018-06-29 JRK: Added query_stars to query a few catalogs on 
#                 Vizier for star photometry given a center
#                 coordinate and box width.
#
# Data directory: This is the ONE line that will be custom to 
# each individual person. Update it do reflect where your data is stored.
#datadir='/Users/jkamenetzky3/data/'
#
# Import packages
# General packages 
import numpy as np
import matplotlib.pyplot as plt
# Fits and photometry
from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits
import ccdproc


###########
# Query various catalogs for standard stars using Vizier,
# and combine the results in a way that makes sense.
def query_stars(centercoord,width):

    Vizier.ROW_LIMIT=-1

    # Try Vizier, catalogs: II/183A is Landolt 1992, 
    # J/AJ/146/131, Landolt 2013 with +50 deg declination.
    # I/239/hip_main is Hipparcos and Tycho catalog
    catalogs=['II/183A/table2',         # Landolt 1992
              'J/AJ/146/131/standards',    # Landolt 2013 (> +50 dec)
              'I/239/hip_main',#  Hipparcos and Tycho
         #     'II/336/apass9', # AAVSO Photometric All Sky Survey (APASS) DR9 - HUGE
              ] 

    result=Vizier.query_region(centercoord,width=[sep1*2,sep2*2],
         catalog=catalogs)
    
    tbl=Table(names=('RA','DEC','Vmag','Source'),
        dtype=('float64','float64','float32','str'))
         
    for i in np.arange(np.shape(result)[0]):
        if result.keys()[i]=='I/239/hip_main':
            # Probably a smarter way to do this.
            tmp=result[i]['RAICRS','DEICRS','Vmag']
            tmp.rename_column('RAICRS','RA')
            tmp.rename_column('DEICRS','DEC')
            tmp['Source']='I/239/hip_main'
            tbl=vstack([tbl,tmp])
        elif result.keys()[i]=='II/183A/table2':
            # I think these are too dim anyway... shrug.
            tmp=result[i]['RAJ2000','DEJ2000','Vmag']
            # Convert to decimal degrees
            tmp['RA']=coords.Angle(tmp['RAJ2000'],unit=u.hour).degree
            tmp['DEC']=coords.Angle(tmp['DEJ2000'],unit=u.degree).degree
            tmp['Source']='II/183A/table2'
            tbl=vstack([tbl,tmp['RA','DEC','Vmag','Source']])
        elif result.keys()[i]=='J/AJ/146/131/standards':
            # I think these are too dim anyway... shrug.
            tmp=result[i]['RAJ2000','DEJ2000','__Vmag_']
            # Convert to decimal degrees
            tmp['RA']=coords.Angle(tmp['RAJ2000'],unit=u.hour).degree
            tmp['DEC']=coords.Angle(tmp['DEJ2000'],unit=u.degree).degree
            tmp.rename_column('__Vmag_','Vmag')
            tmp['Source']='J/AJ/146/131/standards'
            tbl=vstack([tbl,tmp['RA','DEC','Vmag','Source']])            
        else: 
            print 'Table type unrecognized by query_stars.'
         
    return tbl
    
    



#######
# A function which calculates altitude, azimuth, and airmass
# given a pixel location x, y, 
# from image header information (requires WCS information).
