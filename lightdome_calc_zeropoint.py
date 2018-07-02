# Given the results of lightdome_photometry, calculate
# the zero-point magnitude of the instrument ("photometric
# calibration") and extinction, k, in magnitudes per airmass.
# This is Falchi 2011 equations 1, 2, 3, and Figure 4.
# See lightdome_photometry.py for a full description of the 
# text files produced, which will be read in here.
#
# Example Usage:
# python lightdome_photometry.py lightdome_timpanogos
#
# MODIFICATION HISTORY
# 2018-07-02 JRK: Commented template file made.
#
#

# Import necessary packages
from astropy.io import ascii
from astropy.table import vstack
import glob

# Set directory. In the future, we'll use the next 4 lines of code
# so that the user can specify the directory from the command line.
# if len(sys.argv)!=2:
#    print('Usage:\npython lightdome_calc_zeropoint.py [lightdome_NAME folder] \n')
#    exit()
# dir='./data/'+sys.argv[2]+'/'
dir='./data/lightdome_timpanogos/' # Just hard-code for now.

# 1) Create a list of all the files in dir+'/photometry/'
#    using a tool like glob. Don't forget to import it.
files=glob.glob(dir+'/photometry/*.txt')


# 2) Create a loop that reads in each table, and then stacks them 
#    all together to create one full table.
data=[]
for f in files:
    tmp=ascii.read(f,format='fixed_width')
    tmp['file']=f.split('/')[-1]   # Add in from which file these rows came
    if not data:
        data=tmp
    else:
        data=vstack([data,tmp])


# 3) Create a new column that calculates airmass. Note 
#    that usually you see a formula for airmass as a function of zenith angle,
#    and our tables have altitude. Also, trigonometric functions are likely
#    expecting angles in radians.
#
# 4) Plot: on horizontal axis, airmass. On vertical axis, the sum 
#    of the catalog magnitude + instrumental magnitude.
# 
# 5) Use a built-in package to find the best-fit line: we need to know
#    the parameters (slope, intercept) and their errors. 
#    Look for how to do this using numpy or scipy.
# 
# 6) Add the best-fit line onto the plot.
#
# 7) Save the plot as a .png file.
#
# 8) Somehow save the fit results (slope, intercept) into a text file that can then be 
#    added to git and read in with future python scripts.