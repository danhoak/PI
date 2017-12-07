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
#from numpy.polynomial import hermite
#from scipy.special import genlaguerre



# Implements Eqs.10-13 of Evans et al arXiv:0910.2716
# w   - angular frequency of mechanical mode
# order - order of optical mode
def calc_optical_TF(w,order,IFO,jitter=array([0,0,1.0])):

    if IFO=='L':
        # LIGO numbers
        phi = 2.7169
        L = 3994.5
        TA = 0.014
        TB = 1e-5

    if IFO=='V':
        # Virgo numbers
        #phi = 2.7456
        L = 2999.8
        TA = 0.0137
        TB = 1e-5

        # calculate the arm gouy phase independently, using the jitter parameters
        R1 = 1424.57
        R2 = 1695.5 + jitter[0]
        g1 = 1 - L/R1
        g2 = 1 - L/R2
        phi = arccos(-1*sqrt(g1*gN))


    tA = sqrt(TA)
    tB = sqrt(TB)

    rA = sqrt(1-TA)
    rB = sqrt(1-TB)

    e4 = array([0,0,0,1,0])

    def prop_op(phi0, phiG, w, L, order, sign):
        if sign=='p':
            return exp(1j*((phi0 - order*phiG) + w*L/constants.c) )
        if sign=='m':
            return exp(1j*((phi0 - order*phiG) - w*L/constants.c) )

    pL_plus = prop_op(0.0,phi,w,L,order,'p')
    pL_minus = prop_op(0.0,phi,w,L,order,'m')

    Sn_plus = array([[0,0,0,0,0],
                     [tA,0,0,0,-1*rA],
                     [0,pL_plus,0,0,0],
                     [0,0,-1*rB,0,0],
                     [0,0,0,pL_plus,0]])

    Sn_minus = array([[0,0,0,0,0],
                     [tA,0,0,0,-1*rA],
                     [0,pL_minus,0,0,0],
                     [0,0,-1*rB,0,0],
                     [0,0,0,pL_minus,0]])

    Gn_plus = dot(e4.T,dot(linalg.inv(identity(5)-Sn_plus),e4))
    Gn_minus = dot(e4.T,dot(linalg.inv(identity(5)-Sn_minus),e4))

    Gn = Gn_minus - Gn_plus.conj()

    #print w/(2*pi), Gn_plus, Gn_minus, Gn.real

    return Gn.real




# Implements Eqs.15-17 of Evans et al arXiv:0910.2716
# w   - angular frequency of mechanical mode
# order - order of optical mode
# IFO - string, either 'L' or 'V'
# jitter - 3-element array, NE and WE change to ROC (meters) and change to PRC gouy phase (fraction)
def calc_optical_TF_fullIFO(w,order,IFO,jitter=array([0,0,1.0])):

    if IFO=='V':
        
        # Virgo numbers
        Tix = 0.0137
        Tiy = 0.0136
        Tex = 1e-5
        Tey = 1e-5
        Tpr = 0.04835
        Tsr = 1.0
        Tbs = 0.5

    if IFO=='L':

        # LIGO numbers
        Tix = 0.014
        Tiy = 0.014
        Tex = 1e-5
        Tey = 1e-5
        Tpr = 0.03
        Tsr = 0.37
        Tbs = 0.5


    tix = sqrt(Tix)
    tiy = sqrt(Tiy)
    tex = sqrt(Tex)
    tey = sqrt(Tey)
    tpr = sqrt(Tpr)
    tsr = sqrt(Tsr)
    tbs = sqrt(Tbs)

    rix = sqrt(1-Tix)
    riy = sqrt(1-Tiy)
    rex = sqrt(1-Tex)
    rey = sqrt(1-Tey)
    rpr = sqrt(1-Tpr)
    rsr = sqrt(1-Tsr)
    rbs = sqrt(1-Tbs)

    e4 = zeros(12)
    e4[0] = 1.0

    M = zeros((12,12))

    M[0,1] = -rex
    M[1,0] = -rix
    M[1,5] = tix
    M[2,3] = -rey
    M[3,2] = -riy
    M[3,7] = tiy
    M[4,0] = tix
    M[4,5] = rix
    M[5,9] = tbs
    M[5,11] = rbs
    M[6,2] = tiy
    M[6,7] = riy
    M[7,9] = -rbs
    M[7,11] = tbs
    M[8,4] = tbs
    M[8,6] = -rbs
    M[9,8] = -rpr
    M[10,4] = rbs
    M[10,6] = tbs
    M[11,10] = -rsr

    if IFO=='V':

        # Virgo numbers
        L1234 = 2999.8
        L56 = 6.142
        L78 = 5.912
        L910 = 5.925
        L1112 = 5.925

        # calculate the arm gouy phase independently, using the jitter parameters
        # for now use a constant ROC for the input mirrors
        R1 = 1424.57
        RN = 1695.5 + jitter[0]
        RW = 1695.5 + jitter[1]
        
        g1 = 1 - L1234/R1
        gN = 1 - L1234/RN
        gW = 1 - L1234/RW

        phiGN = arccos(-1*sqrt(g1*gN))
        phiGW = arccos(-1*sqrt(g1*gW))

        phi018 = 0.0
        phi0910 = pi/2
        #phiG1234 = 2.7456
        phiG5678 = 0.0
        phiG910 = 0.002 * jitter[2]
        phiG1112 = 0.002

        fHOM = (1-phiGN/pi) * constants.c / (2*L1234)

    if IFO=='L':
        
        # LIGO numbers
        L1234 = 3994.5
        L56 = 4.85
        L78 = 4.9
        L910 = 52.3
        L1112 = 50.6

        # calculate the arm gouy phase independently, using the jitter parameters
        R1 = 1939.3
        RN = 2240.0 + jitter[0]
        RW = 2240.0 + jitter[1]
        
        g1 = 1 - L1234/R1
        gN = 1 - L1234/RN
        gW = 1 - L1234/RW

        phiGN = arccos(-1*sqrt(g1*gN))
        phiGW = arccos(-1*sqrt(g1*gW))
        
        phi018 = 0.0
        phi0910 = pi/2
        #phiG1234 = 2.7169
        phiG5678 = 0.0
        phiG910 = 0.44 * jitter[2]
        phiG1112 = 0.35

        #ep = 1e-3
        #phiGN = phiG1234
        #phiGW = phiG1234 * (1-ep)
        
        fHOM = (1-phiGN/pi) * constants.c / (2*L1234)


    def prop_op(phi0, phiG, w, L, order, sign):
        if sign=='p':
            return exp(1j*((phi0 - order*phiG) + w*L/constants.c) )
        if sign=='m':
            return exp(1j*((phi0 - order*phiG) - w*L/constants.c) )

    p_p = array([prop_op(phi018, phiGN, w, L1234, order, 'p'),
                    prop_op(phi018, phiGN, w, L1234, order, 'p'),
                    prop_op(phi018, phiGW, w, L1234, order, 'p'),
                    prop_op(phi018, phiGW, w, L1234, order, 'p'),
                    prop_op(phi018, phiG5678, w, L56, order, 'p'),
                    prop_op(phi018, phiG5678, w, L56, order, 'p'),
                    prop_op(phi018, phiG5678, w, L78, order, 'p'),
                    prop_op(phi018, phiG5678, w, L78, order, 'p'),
                    prop_op(phi0910, phiG910, w, L910, order, 'p'),
                    prop_op(phi0910, phiG910, w, L910, order, 'p'),
                    prop_op(phi018, phiG1112, w, L1112, order, 'p'),
                    prop_op(phi018, phiG1112, w, L1112, order, 'p')])
    
    p_plus = diag(p_p)

    p_m = array([prop_op(phi018, phiGN, w, L1234, order, 'm'),
                   prop_op(phi018, phiGN, w, L1234, order, 'm'),
                   prop_op(phi018, phiGW, w, L1234, order, 'm'),
                   prop_op(phi018, phiGW, w, L1234, order, 'm'),
                   prop_op(phi018, phiG5678, w, L56, order, 'm'),
                   prop_op(phi018, phiG5678, w, L56, order, 'm'),
                   prop_op(phi018, phiG5678, w, L78, order, 'm'),
                   prop_op(phi018, phiG5678, w, L78, order, 'm'),
                   prop_op(phi0910, phiG910, w, L910, order, 'm'),
                   prop_op(phi0910, phiG910, w, L910, order, 'm'),
                   prop_op(phi018, phiG1112, w, L1112, order, 'm'),
                   prop_op(phi018, phiG1112, w, L1112, order, 'm')])

    p_minus = diag(p_m)
    
    Sn_plus = dot(M,p_plus)
    Sn_minus = dot(M,p_minus)

    Gn_plus = dot(e4.T,dot(linalg.inv(identity(12)-Sn_plus),e4))
    Gn_minus = dot(e4.T,dot(linalg.inv(identity(12)-Sn_minus),e4))

    Gn = Gn_minus - Gn_plus.conj()

    #print w/(2*pi), Gn_plus, Gn_minus, Gn.real

    return Gn.real, fHOM




P = 120e3
lamd = 1064.0e-9
M = 40.0 # in kilograms!

# pick a mode
f1 = 44061.0
Bmn = 0.3
order = 7
Q = 6e7

f2 = 45308.0
f3 = 45251.0
f4 = 43506.0
f5 = 44159.0
f6 = 43111.0

f = arange(40e3,50e3,2.0)

#data = genfromtxt('TM_modes.txt')
#mode_num = data[:,0]
#f1 = data[:,1]

data = genfromtxt('modi.txt')
mode_num = data[:,1]
f11 = data[:,0]

# Figure out which modes near 482 overlap with 7th order...
pkl_file = open('new_mode_overlaps.pkl', 'rb')
modes = pickle.load(pkl_file)
pkl_file.close()
for i in range(440,530):

    list_of_overlaps = modes[i]

    for j in range(len(list_of_overlaps)):

        type = list_of_overlaps[j][0]
        n = list_of_overlaps[j][1]
        m = list_of_overlaps[j][2]
        Bmn = list_of_overlaps[j][3]

        if type=='HG':
            order = n+m
        if type=='LG':
            order = 2*n+m

        if order==7 and Bmn > 0.05:
            print f11[i], i+1, type, n, m, Bmn


# pick a mode
f1 = 44061.0
Bmn = 0.3
order = 7
Q = 6e7

w = 2*pi*f

R = zeros(len(f))
for i in range(len(f)):
    [Gn_real, fHOM] = calc_optical_TF_fullIFO(w[i],order,'V')
    R[i] = 8*pi*Q*P / (M * w[i]**2 * constants.c * lamd) * Gn_real * Bmn**2

R_RH = zeros(len(f))
for i in range(len(f)):
    [Gn_real, fHOM] = calc_optical_TF_fullIFO(w[i],order,'V',[-8.0,0.0,1.0])
    R_RH[i] = 8*pi*Q*P / (M * w[i]**2 * constants.c * lamd) * Gn_real * Bmn**2




fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.semilogy(array([f1,f1])/1e3,[1e-6,1e2],'k--',linewidth=1.2)
pylab.semilogy(array([f2,f2])/1e3,[1e-6,1e2],'k--',linewidth=1.2)
pylab.semilogy(array([f3,f3])/1e3,[1e-6,1e2],'k--',linewidth=1.2)
pylab.semilogy(array([f4,f4])/1e3,[1e-6,1e2],'k--',linewidth=1.2)
pylab.semilogy(array([f5,f5])/1e3,[1e-6,1e2],'k--',linewidth=1.2)
#pylab.semilogy([f5,f5],[1e-6,1e2],'k--',linewidth=1.2)
#pylab.semilogy([order*fHOM,order*fHOM],[1e-6,1e2],'r--',linewidth=1.2)
pylab.semilogy(f/1e3,R,'orange',linestyle='-',linewidth=1.8,alpha=1.0,label='Nominal Optical TF')

#pylab.semilogy(f/1e3,R,'orange',linestyle='--',linewidth=1.8,alpha=0.6,label='Nominal Optical TF')
#pylab.semilogy(f/1e3,R_RH,'orange',linestyle='-',linewidth=1.6,alpha=1.0,label='ETM ROC - 8m')

pylab.grid(True, which='both', alpha=0.8)
pylab.yticks(fontsize=12)
pylab.xticks(fontsize=12)
pylab.ylabel(r'Parametric Gain, |R$_m$|',fontsize=14)
pylab.xlabel('Mechanical Mode Frequency[kHz]',fontsize=14)
pylab.legend(loc=2,fontsize=12)
pylab.xlim(42,46)
pylab.ylim(1e-4,1e2)

pylab.savefig('Mode482.png',bbox_inches='tight')

