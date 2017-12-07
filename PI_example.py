#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Examples that illustrate PI code, following Evans et al.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import constants
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
matplotlib.rcParams.update({'savefig.dpi':250})
#import pickle
from PI_tools import *
import sys

# Parameters for IFO

P = 1000e3
lamd = 1064.0e-9
M = 40.0 # in kilograms!
Q = 1e7
Bmn = 0.21 # overlap for a 2nd order mode



fignum=0


# Plot single F-P cavity parametric gain for a single optical mode; similar to Fig 4 of Evans et al.
f1 = arange(20e3,50e3,3.0)

fignum=fignum+1
pylab.figure(fignum)

R_FP = zeros(len(f1))
for i in range(len(f1)):

    Gn_real = calc_optical_TF(2*pi*f1[i],2,'L')

    # There is an error in the Evans et al. Eq. 9 of a factor of two!
    R_FP[i] = 8*pi*Q*P / (M * (2*pi*f1[i])**2 * constants.c * lamd) * (Gn_real * Bmn**2)

pylab.plot(f1/1e3,R_FP,'red',marker='.',markersize=5,markeredgecolor=None,label='LIGO')

pylab.grid(True, which='both', alpha=0.8)
pylab.yticks(fontsize=12)
pylab.xticks(fontsize=12)
pylab.ylabel(r'Parametric Gain, R$_{m,HG11}$',fontsize=14)
pylab.xlabel('Mode Frequency [kHz]',fontsize=14)

pylab.savefig('FP_paramgain_aLIGO.png',bbox_inches='tight')




# Plot single F-P cavity parametric gain for a single optical mode; similar to Fig 4 of Evans et al.
# Note, the calc_optical_TF_fullIFO function uses the as-built IFO parameters, which are a little
# different than the ones used in the 2009 paper. So the peak position is a little different.
f2 = arange(47.5e3,47.8e3,1.0)

fignum=fignum+1
pylab.figure(fignum)

R_IFO = zeros(len(f2))
for i in range(len(f2)):

    Gn_real = calc_optical_TF_fullIFO(2*pi*f2[i],2,'L')

    R_IFO[i] = 8*pi*Q*P / (M * (2*pi*f2[i])**2 * constants.c * lamd) * (Gn_real * Bmn**2)

pylab.semilogy(f1/1e3,R_FP,'blue',marker='.',markersize=5,markeredgecolor=None,label='LIGO')
pylab.semilogy(f2/1e3,R_IFO,'red',marker='.',markersize=5,markeredgecolor=None,label='LIGO')

pylab.xlim(47.5,47.8)
pylab.ylim(1e-2,1e2)
pylab.grid(True, which='both', alpha=0.8)
pylab.yticks(fontsize=12)
pylab.xticks(fontsize=12)
pylab.ylabel(r'Parametric Gain, R$_{m,HG11}$',fontsize=14)
pylab.xlabel('Mode Frequency [kHz]',fontsize=14)

pylab.savefig('IFO_paramgain_aLIGO.png',bbox_inches='tight')


