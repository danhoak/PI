#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
#from scipy import constants
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
matplotlib.rcParams.update({'savefig.dpi':250})
from numpy.polynomial import hermite
from scipy.special import genlaguerre
from scipy import interpolate
#from math import factorial

# Function to return either HG or LG field of given order
# x and y are coorinates for field
# waist is waist size in meters
# aperture is radius of field, in meters - outside is NaN
# xm and ym are orders (radial and azimuthal for LG modes)
# type is a string, HG or LG, to return
def HOM_field(aa,bb,waist,aperture,xm,ym,type):

    if type=='HG':

        xcoeff = zeros(xm+1)
        xcoeff[xm]=1
        ycoeff = zeros(ym+1)
        ycoeff[ym]=1

        Hx = hermite.hermval(sqrt(2)*aa/waist,xcoeff)
        ex = exp(-aa**2/waist**2)
        Hy = hermite.hermval(sqrt(2)*bb/waist,ycoeff)
        ey = exp(-bb**2/waist**2)

        HermiteGauss = zeros(len(aa))
        for i in range(len(aa)):
        
            if (aa[i]**2 + bb[i]**2) > aperture:
                HermiteGauss[i] = 0.0
                continue

            HermiteGauss[i] = Hx[i] * Hy[i] * ex[i] * ey[i]


        return HermiteGauss

    if type=='LG':

        L = genlaguerre(xm,ym)
        LaguerreGauss = zeros(len(aa))
        for i in range(len(aa)):

            r2 = aa[i]**2 + bb[i]**2
        
            if r2 > aperture:
                LaguerreGauss[i] = 0.0
                continue

            rho = 2*r2/waist**2
            phi = arctan2(bb[i],aa[i])+2*pi

            LaguerreGauss[i] = rho**(ym/2.) * L(rho) * cos(ym*phi) * exp(-rho/2)

        return LaguerreGauss



# ligo waists: ETM=6.2cm, ITM=5.3cm
# virgo waists: ETM=5.3cm, ITM=4.6cm
#waist = 0.053
waist = 0.062
#waist = 0.046

# virgo TM radius
lim = 0.175
aperture = lim**2

xx, yy = mgrid[-lim:lim:0.001, -lim:lim:0.001]

idx = 35#9#159
type = 'LG'
n = 1
m = 1

# get mode data
data = genfromtxt('Maps_NEW_v2/shape_' + str(idx) + '.txt',skip_header=1)
y = data[:,0]
x = data[:,1]
z = data[:,2]

points = vstack((x.T,y.T)).T

data = genfromtxt('Maps_NEW_v2/freqs.txt')
num = data[:,0]
f1 = data[:,1]

# Get the pump field mode
E_LG00 = HOM_field(x,y,waist,aperture,0,0,'LG')

f0 = interpolate.griddata(points, E_LG00, (xx, yy), method='cubic', fill_value = NaN)
f0_norm = f0 / sqrt(nansum(f0**2))

# interpolate surface map and optical mode to meshgrid for plotting
U = interpolate.griddata(points, z, (xx, yy), method='cubic', fill_value=NaN)
U_norm = U * sqrt(40.0)  # mass matrix normalization

# pick the optical mode
if type=='HG':
    E_HOM = HOM_field(x,y,waist,aperture,n,m,'HG')
if type=='LG':
    E_HOM = HOM_field(x,y,waist,aperture,n,m,'LG')

fn = interpolate.griddata(points, E_HOM, (xx, yy), method='cubic', fill_value = NaN)
fn_norm = fn / sqrt(nansum(fn**2))

print nansum(U_norm**2)
print nansum(f0_norm**2)
print nansum(fn_norm**2)

overlap = abs(nansum(U_norm*f0_norm*fn_norm))

print 'Overlap integral:', overlap


fignum=0

fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(1,3,1)
#pylab.imshow(U_norm,aspect='equal',origin='lower',cmap=pylab.cm.cubehelix, interpolation='nearest')
pylab.imshow(U_norm,aspect='equal',origin='lower',cmap=pylab.cm.jet_r, interpolation='nearest')
pylab.xticks(visible=False)
pylab.yticks(visible=False)
pylab.title('Mechanical Mode (' + str(int(f1[idx-1])) + ' Hz)',fontsize=8)
#pylab.colorbar()

pylab.subplot(1,3,2)
pylab.imshow(fn_norm,aspect='equal',origin='lower',cmap=pylab.cm.jet, interpolation='nearest')
pylab.xticks(visible=False)
pylab.yticks(visible=False)
pylab.title('Optical Mode, ' + type + str(n) + str(m),fontsize=8)
#pylab.colorbar()

pylab.subplot(1,3,3)
#pylab.imshow(fn,aspect='equal',origin='lower',cmap=pylab.cm.cubehelix, interpolation='nearest')
pylab.imshow(U_norm * f0_norm * fn_norm,aspect='equal',origin='lower',cmap=pylab.cm.jet, interpolation='nearest')
pylab.xticks(visible=False)
pylab.yticks(visible=False)
pylab.title('Overlap = ' + str(round(overlap,3)),fontsize=8)
#pylab.colorbar()

pylab.savefig('mode_overlap.png',bbox_inches='tight')
