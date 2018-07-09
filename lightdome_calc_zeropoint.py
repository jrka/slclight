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
# 2018-07-05 NRC: Calculated airmass and created scatter plot.
#
#

# Import necessary packages
from startup_slclight import *
from astropy.io import ascii
from astropy.table import vstack
import glob
from math import pi
from scipy import stats

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
numpy_alt=np.array(data['alt'])
zenith_ang=(90-(numpy_alt))*(pi/180)
cos_zenith=np.cos(zenith_ang)
sec_zenith=1/(cos_zenith)
airmass= sec_zenith
data['airmass']= airmass
airmass_threshold= np.where(airmass<=4)
data=data[airmass_threshold]
airmass= data['airmass']



# 4) Plot: on horizontal axis, airmass. On vertical axis, the sum 
#    of the catalog magnitude + instrumental magnitude.
numpy_vmag=np.array(data['Vmag'])
numpy_instrmag=np.array(data['instrmag'])
magnitudes = (numpy_vmag + numpy_instrmag)
data['magnitudes']= magnitudes
plt.scatter(airmass,magnitudes, color='c', s=10)
plt.xlabel('Airmass')
plt.ylabel('Magnitudes')
plt.title('Instrumental Zero Point')


# 5) Use a built-in package to find the best-fit line: we need to know
#    the parameters (slope, intercept) and their errors. 
# 	We might still need code that takes all the nan values, if any.
slope, intercept, r_value, p_value, std_err = stats.linregress(airmass,magnitudes)
r_squared = r_value**2

 
# 6) Add the best-fit line onto the plot.
# Why when our x range is till five it graphs till 4?
x= np.array(range(0,5))
y= slope*x+intercept


plt.plot(x,y, color='k', label='y= %fx + %f\nR^2=%f' %(slope,intercept,r_squared))
plt.scatter(airmass,magnitudes, color='c', s=10)
plt.xlabel('Airmass')
plt.ylabel('Magnitudes')
plt.title('Instrumental Zero Point')
plt.legend()
plt.show()


# 7) Save the plot as a .png file.
plt.plot(x,y, color='k', label='y= %fx + %f\nR^2=%f' %(slope,intercept,r_squared))
plt.scatter(airmass,magnitudes, color='c', s=10)
plt.xlabel('Airmass')
plt.ylabel('Magnitudes')
plt.title('Instrumental Zero Point')
plt.legend()
plt.savefig(dir+'zeropoint.png')


# 8) Somehow save the fit results (slope, intercept) into a text file that can then be 
#    added to git and read in with future python scripts.










