#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
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


# Set up our IFO: 40kg masses, 10kW arm power


P = 10e3
lamd = 1064.0e-9
M = 40.0 # in kilograms!
Q_start = 50e6  # central Q for lognormal distribution - based on conservative estimate from logbook:37588


# get the mechanical mode frequencies
data = genfromtxt('modi.txt')
mode_num = data[:,1]
f1 = data[:,0]

# get the mode overlap calculations - this is done previously and written to a file, to save time
# NOTE - the mode overlap results are probably wrong!
pkl_file = open('new_mode_overlaps.pkl', 'rb')
modes = pickle.load(pkl_file)
pkl_file.close()

# array to store monte carlo results for each mechanical mode
# change the value of N to increase the number of Monte Carlo trials
N = 300
MC_gains = zeros((len(f1),N))

# get random numbers for jittering Gouy phases

# meters
sigma_ETM_ROC = 2.0

# fraction
sigma_PRC_Gouy = 0.1

ROC_jitter = random.randn(N,2)*sigma_ETM_ROC
PRC_jitter = random.randn(N)*sigma_PRC_Gouy + 1.0

# now start looping over different trials
# in each trial, we will calculate the parametric gain for each mechanical mode
# we will vary the ETM ROCs, the PRC Gouy phase, the mechanical mode Q, and the mode overlap
for ii in range(N):

    print ii

    jitter = array([ROC_jitter[ii,0], ROC_jitter[ii,1], PRC_jitter[ii]])

    # Loop over the mechanical modes
    for i in range(len(f1)):

        # mechanical mode angular frequency
        ww = 2*pi*f1[i]
        
        # jitter the mode Q following a lognormal distribution
        Q = Q_start * random.lognormal(-0.3,0.3)

        list_of_overlaps = modes[i]

        Gn_real_sum = 0
        for j in range(len(list_of_overlaps)):

            type = list_of_overlaps[j][0]
            n = list_of_overlaps[j][1]
            m = list_of_overlaps[j][2]
            Bmn = list_of_overlaps[j][3]

            """
            if Bmn > 0.5:
                print i+1, type, n, m, Bmn
            """

            if type=='HG':
                order = n+m
            if type=='LG':
                order = 2*n+m
        
            # jitter the mode overlap Bmn by a few 10s of percent
            Bmn_sigma = 0.2
            Bmn_jitter = Bmn * (random.randn()*Bmn_sigma+1)

            # calculate the optical gain for this mechanical mode and this optical mode
            # sum the optical gain over optical modes
            Gn_real = calc_optical_TF_fullIFO(ww,order,'V',jitter) * Bmn_jitter**2
            Gn_real_sum += Gn_real
            

        MC_gains[i,ii] = 8*pi*Q*P / (M * ww**2 * constants.c * lamd) * Gn_real_sum



# We now have a list of parametric gains for each optical mode, at a given input power
# We want to calculate when these modes will become unstable (as a function of input power)



# calculate the number of unstable modes across a vector of input powers

# factor to multiply gain
pwr = logspace(0,2,4000)

#counts = array([1.,2.,3.,5.,10.,20.,40.,60.])
counts = array([1.,2.,3.,5.,10.,20.])

pwr_counts = zeros((len(pwr),N))
pwr_prob = zeros((len(counts),len(pwr)))


for i in range(N):
    for j in range(len(pwr)):
        # for each trial, calculate number of unstable modes for a given circulating power
        pwr_counts[j,i] = sum(MC_gains[:,i] * pwr[j] >= 1.0)


for i in range(len(counts)):
    for j in range(len(pwr)):
    
        # fraction of trials with >N unstable modes at given circulating power
        pwr_prob[i,j] = sum(pwr_counts[j,:] >= counts[i])/float(N)



"""
print '10W input:', sum(5*MC_gains >= 1.0)
print '25W input:', sum(12.2*MC_gains >= 1.0)
print '35W input:', sum(17.1*MC_gains >= 1.0)
#print '1MW:', sum(100*MC_gains >= 1.0)
"""


fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.semilogx([50, 50],[0.0, 1.0],'k--',linewidth=2.0,alpha=0.3)
pylab.semilogx([125, 125],[0.0, 1.0],'k--',linewidth=2.0,alpha=0.3)

for i in range(len(counts)):
    pylab.semilogx(pwr*P/1e3,pwr_prob[i,:],'-',label=str(counts[i]),linewidth=1.2)

pylab.grid(True, which='both', alpha=0.8)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylabel(r'Probability of at least N Unstable Modes',fontsize=10)
pylab.xlabel('Circulating Power in Each Arm Cavity [kW]',fontsize=10)
pylab.legend(loc=2,fontsize=10)
#pylab.xlim(10,15)
#pylab.ylim(1e-3,1e3)

pylab.savefig('Virgo_mode_probabilities.png',bbox_inches='tight')





# And, we want to calculate the parametric gain for each mechanical mode, marginalizing over the jittered values (ROC, phase, Q, overlap, etc)
# For each mechanical mode, we plot the 90% confidence value for the parametric gain at 50kW (5x baseline)

pwr = 5.0
maxP = percentile(pwr*MC_gains,90.0,1)

idx = argsort(maxP)[::-1]


# print the top 30 modes
for i in range(30):
    list_of_overlaps = modes[idx[i]]
    #maxmode = 0
    #Bmn = 0.0
    """
    for j in range(len(list_of_overlaps)):
        #print list_of_overlaps[j], abs(list_of_overlaps[j][3])
        #if abs(list_of_overlaps[j][3]) > Bmn:
        #    maxmode = j
        #    Bmn = abs(list_of_overlaps[j][3])
    """
    #type = list_of_overlaps[maxmode][0]
    #n = list_of_overlaps[maxmode][1]
    #m = list_of_overlaps[maxmode][2]
    #Bmn = list_of_overlaps[maxmode][3]

    print mode_num[idx[i]], f1[idx[i]], maxP[idx[i]]#, type, n, m, Bmn



fignum=fignum+1
pylab.figure(fignum)

pylab.semilogy(f1/1e3,maxP,'ro',markersize=5)

pylab.grid(True, which='both', alpha=0.8)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylabel(r'Parametric Gain (90% confidence level)',fontsize=10)
pylab.xlabel('Mechanical Mode Frequency [kHz]',fontsize=10)
#pylab.legend(loc=2,fontsize=10)
#pylab.xlim(10,15)
pylab.ylim(1e-3,max(1e1,1.1*max(maxP)))

pylab.savefig('Virgo_mode_gains.png',bbox_inches='tight')
