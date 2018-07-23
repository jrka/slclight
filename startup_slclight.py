# The purpose of this file is to load in all of the packages and 
# variables that we are always using to complete this project.
# By loading in this file at the start of each script, we reduce
# the need to copy/paste and continually update each individual
# script.
#
# Revision History
# 2018-07-16 JRK: In order to get variability and error 
#                 columns from Vizier, query all 3 catalogs 
#                 separately. Exclude any Hipparchos/Tycho sources
#                 with any variability flag. Add e_Vmag to table.
# 2018-07-10 JRK: Suppress warnings for metadata_conflicts when 
#                 stacking tables in query_stars.
# 2018-07-03 JRK: Added read_metadata to return a dictionary of 
#                 info from the metadata file; helpful to use this
#                 instead of fits headers for lat, long, and elevation,
#                 which may use different conventions or be missing. 
# 2018-07-02 JRK: Added platescale to get the platescale in arcsec/pix
#                 from the header comments added by astrometry.net
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
import sys
from matplotlib.backends.backend_pdf import PdfPages
import os


###########
# Query various catalogs for standard stars using Vizier,
# and combine the results in a way that makes sense.
def query_stars(centercoord,width):
    from astroquery.vizier import Vizier
    from astropy.table import vstack, Table, join
    from astropy import coordinates as coords
    
    # Try Vizier, catalogs: II/183A is Landolt 1992, 
    # J/AJ/146/131, Landolt 2013 with +50 deg declination.
    # I/239/hip_main is Hipparcos and Tycho catalog
    catalogs=['II/183A/table2',         # Landolt 1992
              'J/AJ/146/131/standards',    # Landolt 2013 (> +50 dec)
              'I/239/hip_main',#  Hipparcos and Tycho
         #     'II/336/apass9', # AAVSO Photometric All Sky Survey (APASS) DR9 - HUGE
              ] 
     ##### vvvvvvvvv OLD VERSION    
#      Vizier.ROW_LIMIT=-1
#
#     result=Vizier.query_region(centercoord,width=width,
#          catalog=catalogs) # When queried all at once, miss error and variability.
#     
#     tbl=Table(names=('RA','DEC','Vmag','Source'),
#         dtype=('float64','float64','float32','str'))
#          
#     for i in np.arange(np.shape(result)[0]):
#         if result.keys()[i]=='I/239/hip_main':
#             # Probably a smarter way to do this.
#             tmp=result[i]['RAICRS','DEICRS','Vmag']
#             tmp.rename_column('RAICRS','RA')
#             tmp.rename_column('DEICRS','DEC')
#             tmp['Source']='I/239/hip_main'
#             tbl=vstack([tbl,tmp],metadata_conflicts='silent')
#         elif result.keys()[i]=='II/183A/table2':
#             # I think these are too dim anyway... shrug.
#             tmp=result[i]['RAJ2000','DEJ2000','Vmag']
#             # Convert to decimal degrees
#             tmp['RA']=coords.Angle(tmp['RAJ2000'],unit=u.hour).degree
#             tmp['DEC']=coords.Angle(tmp['DEJ2000'],unit=u.degree).degree
#             tmp['Source']='II/183A/table2'
#             tbl=vstack([tbl,tmp['RA','DEC','Vmag','Source']],metadata_conflicts='silent')
#         elif result.keys()[i]=='J/AJ/146/131/standards':
#             # I think these are too dim anyway... shrug.
#             tmp=result[i]['RAJ2000','DEJ2000','__Vmag_']
#             # Convert to decimal degrees
#             tmp['RA']=coords.Angle(tmp['RAJ2000'],unit=u.hour).degree
#             tmp['DEC']=coords.Angle(tmp['DEJ2000'],unit=u.degree).degree
#             tmp.rename_column('__Vmag_','Vmag')
#             tmp['Source']='J/AJ/146/131/standards'
#             tbl=vstack([tbl,tmp['RA','DEC','Vmag','Source']],metadata_conflicts='silent')            
#         else: 
#             print 'Table type unrecognized by query_stars.'
    
    ##### ^^^^^^^^^ OLD VERSION    

    ##### vvvvvvvvv NEW VERSION    
    tbl=Table(names=('RA','DEC','Vmag','e_Vmag','Source'),
        dtype=('float64','float64','float32','float32','str'))
        
    for cat in catalogs:
        if cat=='I/239/hip_main':
            v=Vizier(columns=['RAICRS','DEICRS','Vmag','VarFlag','e_VTmag'],row_limit=-1)
            tmp=v.query_region(centercoord,width=width,catalog=cat)
            # Remove any variability (even though a flag of 1 is <0.06 mag)
            # Unvariable ones are those whose column value is masked.
            # Sometimes an error happens in the parsing by Vizier.
            if not tmp or len(tmp[0].colnames)!=5: 
                print 'Error or no sources from I/239/hip_main'
            else:
                print 'Removing ',len(np.where(tmp[0]['VarFlag'].mask==False)[0]),' sources due to variability'
                tmp=tmp[0][tmp[0]['VarFlag'].mask==True]['RAICRS','DEICRS','Vmag','e_VTmag']
                tmp.rename_column('RAICRS','RA')
                tmp.rename_column('DEICRS','DEC')
                tmp.rename_column('e_VTmag','e_Vmag')
                tmp['Source']='I/239/hip_main'
                tbl=vstack([tbl,tmp],metadata_conflicts='silent')            
        elif cat=='II/183A/table2':
            # I think these are too dim anyway... shrug.
            v=Vizier(columns=['RAJ2000','DEJ2000','Vmag','e_Vmag'],row_limit=-1)
            tmp=v.query_region(centercoord,width=width,catalog=cat)
            if tmp and len(tmp[0].colnames)==4:
                tmp=tmp[0]
                # Convert to decimal degrees
                tmp['RA']=coords.Angle(tmp['RAJ2000'],unit=u.hour).degree
                tmp['DEC']=coords.Angle(tmp['DEJ2000'],unit=u.degree).degree
                tmp['Source']='II/183A/table2'
                tbl=vstack([tbl,tmp['RA','DEC','Vmag','e_Vmag','Source']],metadata_conflicts='silent')
        elif cat=='J/AJ/146/131/standards':
            # I think these are too dim anyway... shrug.
            v=Vizier(columns=['RAJ2000','DEJ2000',r'<Vmag>',r'e_<Vmag>'],row_limit=-1)
            tmp=v.query_region(centercoord,width=width,catalog=cat)
            if tmp and len(tmp[0].colnames)==4:
                tmp=tmp[0]       
                # Convert to decimal degrees
                tmp['RA']=coords.Angle(tmp['RAJ2000'],unit=u.hour).degree
                tmp['DEC']=coords.Angle(tmp['DEJ2000'],unit=u.degree).degree
                tmp.rename_column('__Vmag_','Vmag')
                tmp.rename_column('e__Vmag_','e_Vmag')
                tmp['Source']='J/AJ/146/131/standards'
                tbl=vstack([tbl,tmp['RA','DEC','Vmag','e_Vmag','Source']],metadata_conflicts='silent')            
        else: 
            print 'Table type unrecognized by query_stars.'
            
    ##### ^^^^^^ NEW VERSION
         
    return tbl
    
###########
# Given a table with RA and DEC (in decimal degrees), remove
# rows that are within "dist" arcseconds of each other.
def remove_close_sources(result,dist):
    from astropy import coordinates as coords
    import numpy as np
    
    corows=coords.SkyCoord(ra=result['RA'],dec=result['DEC'],unit="deg")
    include=np.ones(len(result),dtype=bool)
    for i,row in enumerate(result):
        coi=coords.SkyCoord(row['RA']*u.deg,row['DEC']*u.deg)
        di=coi.separation(corows)
        if len(np.where(di.arcsec<dist)[0])>1: # Will be 1; closest to itself
            include[i]=False
    result=result[include]
    print 'Removed ',len(np.where(include==False)[0]),' sources due to being too close together.'
    return result

###########
# Return the pixel scale reported by astrometry.net in the FITS header.
def platescale(hdr):
    # One of the COMMENT fields will look like this:
    # COMMENT scale: 97.0455 arcsec/pix
    # You don't want the one that says "Field scale"
    pixscale=np.nan
    for i,line in enumerate(hdr['COMMENT']):
        if 'arcsec/pix' in line and 'Field' not in line: 
            pixscale=np.float(line.split(' ')[1])
    return pixscale




#######
# Return a dictionary of the contents of metadata.txt file
# for a given lightdome folder. See an example file
# for how to create it with the right format, units, etc.

def read_metadata(dir):
    keys=['sitename','lat','lon','elev','date','time','camera','lens',
          'observers','airtmp','sqm','comments']
    flt=[1,2,3,9,10] # The indices of above keys to return as floats.
                       # The rest will be strings.
    d={}
    i=0
    with open('./data/'+dir+'/metadata.txt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                if i in flt and line.rstrip('\n')=='Unknown':
                    d[keys[i]]=np.nan
                else:
                    d[keys[i]]=np.float(line.rstrip('\n')) if i in flt else line.rstrip('\n')
                i+=1
    return d

#######
# Visualize (zoom-in) to our sources.
# Input f (filename),result is 
def plot_sources(f,positions,result):

    # Plot this as a PDF figure.
    pdffile=dir+'/photometry/photometry_'+f.split('.')[0]+'_sources.pdf'
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
    

