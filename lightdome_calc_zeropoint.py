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
# 2018-07-23 JRK: Commented out some unnecessary plotting lines.
#                 Impose a signal/noise cut (of only 1.0) to remove
#                 points with extreme error bars.
#                 Do first once, remove 2 sigma outliers, then fit again.
# 2018-07-05 NRC: Calculated airmass and created scatter plot.
# 2018-07-02 JRK: Commented template file made.
#
#

# Import necessary packages
from setup_plotting import *
from startup_slclight import *
from astropy.io import ascii
from astropy.table import vstack
import glob
from math import sqrt
from math import pi
from scipy import stats

# Set directory. In the future, we'll use the next 4 lines of code
# so that the user can specify the directory from the command line.
# if len(sys.argv)!=2:
#    print('Usage:\npython lightdome_calc_zeropoint.py [lightdome_NAME folder] \n')
#    exit()
dir='./data/'+sys.argv[1]+'/'
#dir='./data/lightdome_westminster/' # Just hard-code for now.

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

print len(data)

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
airmass_threshold= np.where(abs(airmass)<=4)
data=data[airmass_threshold]

print len(data),' after airmass threshold'

# 4) Plot: on horizontal axis, airmass. On vertical axis, the sum 
#    of the catalog magnitude + instrumental magnitude.
#  	 We still need code that takes all the nan values, if any.

# Include a signal/noise threshold?
data=data[np.where(data['instrmag']/data['instrmag_err']>1.0)]

numpy_instrmag= np.array(data['instrmag'])
instrmag_nonan= np.where(~np.isnan(numpy_instrmag))
data= data[instrmag_nonan]
numpy_vmag= np.array(data['Vmag'])
numpy_instrmag= data['instrmag']
magnitudes= (numpy_vmag + numpy_instrmag)
data['magnitudes']= magnitudes
# mag_threshold= np.where(10<magnitudes)  # Removed; this was trying out a special case?
# data= data[mag_threshold]
# magnitudes= data['magnitudes']
instrmag_err = data['instrmag_err']
airmass= data['airmass']

#plt.scatter(airmass,magnitudes, color='c', s=10, label=None)
lfd=latexfd['full']
fig = plt.figure(num=0,figsize=(lfd['figw'],lfd['figh']))
plt.clf()
fig.subplots_adjust(right=lfd['r'])
plt.errorbar(airmass,magnitudes,yerr=instrmag_err, ls='None', capsize=2, label=None,
   marker='o',markersize=2,elinewidth=1.0) # I added markers to "errorbar", so no need to call "scatter"
plt.xlabel('Airmass')
plt.ylabel('Instrumental Magnitude + Literature Magnitude')
plt.title('Instrumental Zero Point')
#plt.show()


# 5) Use a built-in package to find the best-fit line: we need to know
#    the parameters (slope, intercept) and their errors. 

#warnings.simplefilter('ignore', np.RankWarning)
best_fit, cov = np.polyfit(airmass,magnitudes,1,full=False, w=(1.0/instrmag_err), cov=True)
slope= best_fit[0]
intercept= best_fit[1]
err_slope= sqrt(cov[0][0])
err_intercept= sqrt(cov[1][1])

# Add this original line to graph, dotted.
x= np.array(range(0,4))
y= slope*x+intercept
plt.plot(x,y, color='k', marker='None',linestyle=':',label='$\mathbf{y= %.3fx }\pm %.3f \mathbf{+ %.3f }\pm %.3f$' %(slope, err_slope, intercept, err_intercept))

# 5b) Remove 2 sigma outliers, redo the fit.
model=airmass*slope+intercept
diff=magnitudes-model
outlier_mask=np.logical_or(diff>np.median(diff)+2.0*np.std(diff),diff<np.median(diff)-2.0*np.std(diff))
plt.errorbar(airmass[outlier_mask],magnitudes[outlier_mask],yerr=instrmag_err[outlier_mask],
    ls='None',capsize=2,label='Excluded as Outliers',marker='o',markersize=2,color='r',elinewidth=1)

best_fit, cov = np.polyfit(airmass[~outlier_mask],magnitudes[~outlier_mask],1,full=False, w=(1.0/instrmag_err[~outlier_mask]), cov=True)
slope= best_fit[0]
intercept= best_fit[1]
err_slope= sqrt(cov[0][0])
err_intercept= sqrt(cov[1][1])

 
# 6) Add the best-fit line onto the plot.
# Why when our x range is till five it graphs till 4?
x= np.array(range(0,5))
y= slope*x+intercept

plt.plot(x,y, color='k', marker=None,label='$\mathbf{y= %.3fx }\pm %.3f \mathbf{+ %.3f }\pm %.3f$' %(slope, err_slope, intercept, err_intercept))
#plt.scatter(airmass,magnitudes, color='c', s=10, label=None) # No need to plot, relabel again.
#plt.errorbar(airmass,magnitudes,yerr=instrmag_err, ls='None', capsize=2, label=None)
#plt.xlabel('Airmass')
#plt.ylabel('Magnitudes')
#plt.title('Instrumental Zero Point with Error')
plt.legend()
#plt.show()


# 7) Save the plot as a .png file.
# plt.plot(x,y, color='k', label='$\mathbf{y= %.3fx }\pm %.3f \mathbf{+ %.3f }\pm %.3f$' %(slope, err_slope, intercept, err_intercept))
# plt.scatter(airmass,magnitudes, color='c', s=10, label=None)
# plt.errorbar(airmass,magnitudes,yerr=instrmag_err, ls='None', capsize=2, label=None)
# plt.xlabel('Airmass')
# plt.ylabel('Magnitudes')
# plt.title('Instrumental Zero Point with Error')
# plt.legend()
plt.savefig(dir+'zeropoint.png')


# 8) Somehow save the fit results (slope, intercept) into a text file that can then be 
#    added to git and read in with future python scripts.

text = open(dir+'zeropoint_information.txt','w+')
text.write("Atmospheric extinction coefficient in magnitudes per airmass (k) = %.3f +/- %.3f\n"%(slope, err_slope))
text.write("Instrument photometric calibration constant (C) = %.3f +/- %.3f "%(intercept, err_intercept))
text.close()









