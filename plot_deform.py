#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
matplotlib.rcParams.update({'savefig.dpi':250})
from scipy import interpolate

# get the frequencies
data = genfromtxt('Maps_NEW/Modes1-400.txt')
num = data[:,0]
f1 = data[:,1]


fignum=0

#for i in range(1,1201):
for i in range(1,400):
    print i

    #data = genfromtxt('surface_maps/shape_' + str(i) + '.txt',skip_header=1)
    #data = genfromtxt('Modes_newmodel_nolayer/shape_' + str(i) + '.txt',skip_header=1)
    data = genfromtxt('Maps_NEW/shape' + str(i) + '_v9.txt',skip_header=1)
    
    # old data (previous to 2017)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    # new data
    x = data[:,0]
    y = data[:,1]
    z = data[:,4]
    points = vstack((x.T,y.T)).T
    
    xx, yy = mgrid[min(x):max(x):0.0008, min(y):max(y):0.0008]
    
    zz = interpolate.griddata(points, z, (xx, yy), method='cubic')
    
    fignum=fignum+1
    pylab.figure(fignum)
    
    pylab.imshow(zz.T,aspect='equal',origin='lower',cmap=pylab.cm.jet_r, interpolation='nearest', extent=(min(x),max(x),min(y),max(y)))
    #pylab.plot(points[:,0], points[:,1], 'k.', markersize=0.8)

    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.title('Shape ' +str(i) + ' - ' + str(int(f1[i-1])) + ' Hz',fontsize=10)
    pylab.colorbar()

    #pylab.savefig('images/shape_' + str(i) + '.png')
    #pylab.savefig('images_newmodel_nolayer/shape_' + str(i) + '.png')
    pylab.savefig('images_Maps_NEW/shape_' + str(i) + '.png')
    pylab.close()
