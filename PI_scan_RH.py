#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to vary ETM radii of curvature to mitigate PI risk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import constants
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
matplotlib.rcParams.update({'savefig.dpi':250})
import pickle
from PI_tools import *



# Set up some parameters for the IFO

P = 120e3
lamd = 1064.0e-9
M = 40.0 # in kilograms!
Q = 1e8

min_ROC = -10
max_ROC = 1.5
ROC_step = 0.2

ROC = arange(min_ROC,max_ROC,ROC_step)


### If you just want to plot an earlier results file, add a comment here ###

# Load the mechanical mode frequencies
data = genfromtxt('modi.txt')
mode_num = data[:,1]
f1 = data[:,0]

# Load the overlaps between mechanical modes and optical modes
# NOTE - the mode overlap results are probably wrong!
pkl_file = open('new_mode_overlaps.pkl', 'rb')
modes = pickle.load(pkl_file)
pkl_file.close()


unstable = zeros((len(ROC),len(ROC)))

# scan NE and WE ROC and count unstable modes
for ii in range(len(ROC)):

    for jj in range(len(ROC)):

        print ROC[ii], ROC[jj]

        counter=0
        # Loop over the mechanical modes
        for i in range(len(f1)):

            ww = 2*pi*f1[i]
        
            list_of_overlaps = modes[i]

            Gn_real_sum = 0
            for j in range(len(list_of_overlaps)):

                type = list_of_overlaps[j][0]
                n = list_of_overlaps[j][1]
                m = list_of_overlaps[j][2]
                Bmn = list_of_overlaps[j][3]

                if type=='HG':
                    order = n+m
                if type=='LG':
                    order = 2*n+m
        
                # sanity check
                #Gn_real = calc_optical_TF_fullIFO(ww,order,'V',[-6.4, 1.0, 1.0]) * Bmn**2
                
                Gn_real = calc_optical_TF_fullIFO(ww,order,'V',[ROC[ii],ROC[jj],1.0]) * Bmn**2
                Gn_real_sum += Gn_real

            # calculate the parametric gain for this mode
            R = 8*pi*Q*P / (M * ww**2 * constants.c * lamd) * Gn_real_sum

            # if the mode is unstable, increment this pixel by one count
            if R > 1.0:

                #print i+1, f1[i], R
                #print list_of_overlaps
                unstable[ii,jj] += 1



# save the data in case you want to plot it later
output = open('unstable_modes_ROC.pkl', 'wb')
pickle.dump(unstable, output)
output.close()


# Make a plot
fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.imshow(unstable,aspect='equal',origin='lower',cmap=pylab.cm.jet_r, interpolation='nearest')
cbar = pylab.colorbar()
cbar.ax.tick_params(labelsize=10)
cbar.set_label('Number of Unstable Modes', size=10)

labels = arange(min_ROC,max_ROC,2.0)
locs = [argmin(abs(ROC-label)) for label in labels]

pylab.yticks(locs, labels, fontsize=10)
pylab.xticks(locs, labels, fontsize=10)

pylab.xlabel('ETMX ROC Change',fontsize=10)
pylab.ylabel('ETMY ROC Change',fontsize=10)

pylab.title('ETM Ring Heater ROC Scan at ' + str(P/1e3) + ' kW',fontsize=12)

pylab.savefig('120kW_unstable.png',bbox_inches='tight')


"""

# retrieve the data in case you want to plot it later
pkl_file = open('unstable_modes_ROC.pkl', 'rb')
unstable = pickle.load(pkl_file)
pkl_file.close()

# Make a plot
fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.imshow(unstable,aspect='equal',origin='lower',cmap=pylab.cm.jet_r, interpolation='nearest')
cbar = pylab.colorbar()
cbar.ax.tick_params(labelsize=10)
cbar.set_label('Number of Unstable Modes', size=10)

labels = arange(min_ROC,max_ROC,2.0)
locs = [argmin(abs(ROC-label)) for label in labels]

pylab.yticks(locs, labels, fontsize=10)
pylab.xticks(locs, labels, fontsize=10)

pylab.xlabel('ETMX ROC Change',fontsize=10)
pylab.ylabel('ETMY ROC Change',fontsize=10)

pylab.title('ETM Ring Heater ROC Scan at ' + str(P/1e3) + ' kW',fontsize=12)

pylab.savefig('120kW_unstable.png',bbox_inches='tight')

"""
