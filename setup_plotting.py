# Some plotting preferences taken from previous work (ftssurvey/pysurvey_plotting.py)
# Based on 
# http://blog.dmcdougall.co.uk/publication-ready-the-first-time-beautiful-reproducible-plots-with-matplotlib/
#
# How to use:
#   from setup_plotting import *
# This automatically sets some matplotlib rcParams, and
# loads a dictionary called "latexfd" with some plotting specifications, 
# including the sizes of single-column and two-column plots for aastex6.
# These dimensions can be modified for subplots using latex_subplots 
# function, described in more detail below.


from matplotlib import rcParams  
from matplotlib.ticker import MaxNLocator  
import numpy as np

# Fonts. 9 is a good size for half-width (one-column) images. 14 for full.
rcParams['axes.labelsize'] = 14
rcParams['xtick.labelsize'] = 14  
rcParams['ytick.labelsize'] = 14
rcParams['legend.fontsize'] = 14
rcParams['axes.titlesize'] = 14
rcParams['font.family'] = 'sans-serif'  
rcParams['font.sans-serif'] = ['Verdana']  
#rcParams['text.usetex'] = True  # Get rid of this for plotting that uses underscores for proposal IDs

#print rcParams.keys()

# Line thickness
rcParams['axes.linewidth'] = 2
rcParams['lines.linewidth']=3.0
rcParams['lines.markeredgewidth']=1
rcParams['lines.marker']='o'
rcParams['lines.markersize'] = 3
# Need to set elinewidth=1.0 elsewhere.

# Figure size (and tick location)
# In survey3.tex, using \documentclass[twocolumn]{aastex6}, 
#  I put \showthe\textwidth, and got back 513.11743 pt
#  So use this values for full-width figures. Use 0.45x that for single column.
latexfd={'half':{'width':242,
                't':0.95,  # subplots_adjust. Note if multiple panels, divide these.
                'l':0.15,  # EG for 2 rows, divide bottom by 2. Top would be 1-(1-t)/2
                'r':0.95,  
                'b':0.18},
         'inchesperpt':1.0/72.27,
         'goldenratio': (np.sqrt(5)-1.0)/2.0, # Not used.
         'ratio':0.8,                         # I like this better.
         #'tickloc':MaxNLocator(6),
         }
# Preserve same AMOUNT IN INCHES of space
latexfd['full']={'width':510}
latexfd['full']['l']=latexfd['half']['l']*latexfd['half']['width']/latexfd['full']['width']
latexfd['full']['b']=latexfd['half']['b']*latexfd['half']['width']/latexfd['full']['width']
latexfd['full']['t']=1.0-(1.0-latexfd['half']['t'])*latexfd['half']['width']/latexfd['full']['width']
latexfd['full']['r']=1.0-(1.0-latexfd['half']['r'])*latexfd['half']['width']/latexfd['full']['width']

for s in ['full','half']:
    latexfd[s]['figw']=latexfd[s]['width']*latexfd['inchesperpt']
    latexfd[s]['figh']=latexfd[s]['figw']*latexfd['ratio']
# fig=plt.figure(figsize=(latexfd['figw_half'],latexfd['figh_half']))
    
def latex_subplots(lfd,nrow,ncol,wspace=0,hspace=0):
    #================================
    # Calculate new figure properties, given lfd is a dictionary with:
    # figw = figure width in inches (or whatever actual unit)
    # l, b, t, r, as percentage locations for subplots_adjust
    # wspace and hspace as INCHES. 
    
    # ==== Example inputs:
    # # This dictionary defines a figure that is 10x5 inches, and has 
    # 0.1x margins all around (1 inch left/right, 0.5 inch top/bottom).
    # lfd={'figw':10,'figh':5,'l':0.1,'b':0.1,'t':0.9,'r':0.9}
    # # Now we want to divide that same figure up into 3x2 figures, 
    # # with a bit of space inbetween.
    # nrow=3
    # ncol=2
    # wspace=lfd['l']*lfd['figw']*0.34
    # hspace=lfd['b']*lfd['figh']*0.5
    # # newlfd returns a dictionary with the new numbers we wnat to use.
    # newlfd=latex_subplots(lfd,nrow,ncol,wspace=wspace,hspace=hspace)
    # 
    # # Here's what the original figure would look like.
    # fig,ax=plt.subplots(1,1,figsize=(lfd['figw'],lfd['figh']),num=0)
    # fig.subplots_adjust(top=lfd['t'],bottom=lfd['b'],left=lfd['l'],right=lfd['r'])
    # ax.plot([1,2,3],[8,10,12])
    #
    # # And here's our new figure. The outside margins are retained in actual size,
    # # and the width/height of each subfigure (actual axis size) is retained.
    # fig2,ax2=plt.subplots(nrow,ncol,figsize=(newlfd['figw'],newlfd['figh']),num=1)
    # fig2.subplots_adjust(top=newlfd['t'],bottom=newlfd['b'],left=newlfd['l'],right=newlfd['r'],
    #     wspace=newlfd['wspace'],hspace=newlfd['hspace'])
    # ax2[1,1].plot([1,2,3],[8,10,12])
    
    #================================
    
    # First, maintain all margins in actual inches.
    margins={'l':lfd['l']*lfd['figw'],
                'b':lfd['b']*lfd['figw'],
                't':(1.0-lfd['t'])*lfd['figw'],
                'r':(1.0-lfd['r'])*lfd['figw'],
                'wspace':wspace,
                'hspace':hspace}
                
    wa=lfd['figw']-margins['l']-margins['r'] # Actual axis width, original 
    ha=lfd['figh']-margins['b']-margins['t'] # Actual axis height, original   
    
    # Keep width same, now how wide is each plot?
    newwa=(lfd['figw']-margins['l']-margins['r']-(ncol-1)*margins['wspace'])/ncol
    # So now how tall is each plot, keeping original dimensions
    newha=newwa*lfd['figh']/lfd['figw']
    # Total new height will be...
    newh=(newha*nrow)+margins['b']+margins['t']+(nrow-1)*margins['hspace']
    
    newlfd={'figw':lfd['figw'],'figh':newh}
    # Now the inches back to percentages
    newlfd['l']=margins['l']/newlfd['figw']
    newlfd['b']=margins['b']/newlfd['figh']
    newlfd['r']=1.0-margins['r']/newlfd['figw']
    newlfd['t']=1.0-margins['t']/newlfd['figh']
    
    # hspace and wspace are fractions of AXES height
    newlfd['hspace']=margins['hspace']/newha
    newlfd['wspace']=margins['wspace']/newwa
    
    return newlfd
    
                