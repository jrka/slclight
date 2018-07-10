# Given an image width and height in degrees, where should the telescope be 
# pointed to image the whole sky?
#
# To figure out the Field of View in degrees, you can use:
# https://astronomy.tools/calculators/field_of_view/
# (Click on Imaging Mode)
# 
# Use: edit "width" and height" below, then run the script
# to produce a text file and image with pointings in altitude and azimuth.

from startup_slclight import *
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle,Wedge,Polygon
from collections import OrderedDict

# SET IMAGE PROPERTIES: WIDTH IN DEGREES AND HEIGHT IN DEGREES
width=28.14
height=18.76

# Determine number of images required around the horizon, and their overlap.
# Round up.
nhorizon=np.ceil(360.0/width)
# Determine number of images to cover horizon to zenith.
# Round up.
nzenith=np.ceil(90.0/height)

# Determine the center pointings for our horizon images. (in azimuth, degrees)
centers_horiz=np.arange(0,360.0,360.0/nhorizon)
# Determine center pointings for horizon to zenith images (in altitude, degrees)
centers_zenith=np.arange(0.5*(90.0/nzenith),90.0,90.0/nzenith)

# http://burro.case.edu/Academics/Astr306/Coordinates.pdf
# For small separations, Delta(Dec) = (Dec2-Dec1)
#     and Delta(RA)=(RA2-RA1)*Cos(Dec)
# Determine the pointings required as one goes higher in altitude.
centers_layers=OrderedDict()
for i,alt in enumerate(centers_zenith):
    adjwidth=width/np.cos(alt*3.14/180.0)
    nlayer=np.ceil(360.0/adjwidth)
    nudge=(nlayer*adjwidth-360.0)/nlayer
    print nlayer
    centers_layers.update({str(alt):np.arange(0,360.0,adjwidth-nudge)})
    print len(centers_layers[str(alt)])
    
# Report results
print 'Calculation: '
print len(centers_horiz),' exposures around horizon.'
print len(centers_zenith), 'layers of exposures at increasing altitude.'
print np.sum([len(centers_layers[x]) for x in centers_layers]),' total images.'

# Setup a polar plot.
fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.set_theta_zero_location=('N')
ax.set_ylim(0,90)
plt.grid(True)
ax.set_yticks(range(0,90+10,10))
#yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
yLabel = ['', '', '', '', '', '', '', '', '', '']
ax.set_yticklabels(yLabel)
# The plot appears to be taking x-coordinate = AZ in RADIANS, 
# and y-coordinate=ZENITH angle in degrees
# Don't bother too much about direction/orientation, where zero is...
for i,key in enumerate(centers_layers):
    ax.plot([centers_layers[key]*3.14/180.0],[90.0-centers_zenith[i]],'r*')
    print [centers_layers[key]*3.14/180.0],[90.0-centers_zenith[i]]

    for j in np.arange(len(centers_layers[key])):
        theta1=np.linspace((centers_layers[key][j]-0.5*width/np.cos(centers_zenith[i]*3.14/180.0))*3.14/180.0,
             (centers_layers[key][j]+0.5*width/np.cos(centers_zenith[i]*3.14/180.0))*3.14/180.0,10)
        theta2=np.flipud(theta1)
        r1=np.ones(10)*(90.0-centers_zenith[i]+0.5*height)
        r2=np.ones(10)*(90.0-centers_zenith[i]-0.5*height)
        theta=np.concatenate((theta1,theta2))
        r=np.concatenate((r1,r2))
        polygon = Polygon(zip(theta,r),fill=False)
        ax.add_line(polygon)

# Save the image.
ax.set_title('Image Width = '+str(width)+' deg x '+str(height)+' deg')        
plt.savefig('plan_targets_'+str(width)+'_'+str(height)+'.png')

# Write the altitude, azimuth in a file.
file='plan_targets_'+str(width)+'_'+str(height)+'.txt'
with open(file,'w') as file:
    file.write('alt,az\n')
    for i,key in enumerate(centers_layers):
        for j in centers_layers[key]: file.write(str(key)+','+str(j)+'\n')