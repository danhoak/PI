#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A small library of functions for Parametric Instability calculations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import constants
#import matplotlib
#from matplotlib.font_manager import FontProperties
#matplotlib.use("Agg")
#import pylab
#matplotlib.rcParams.update({'savefig.dpi':250})
#import pickle
#from numpy.polynomial import hermite
#from scipy.special import genlaguerre



# calc_optical_TF -- function to return optical transfer coefficient for a single Fabry-Perot arm cavity,
# for a particular optical mode frequency. Implements Eqs.10-13 of Evans et al arXiv:0910.2716
#
# Inputs:
# w   - angular frequency of mechanical mode (float)
# order - order of optical mode (int)
# IFO - 'L' or 'V' (advanced LIGO or advanced Virgo arm cavities) (string)
# jitter - parameters for small random variation of arm (3-element array)
#
# jitter is of the form array([a,b,c])
# a is change to ETMX radius of curvature in meters
# b is change to ETMY ROC in meters
# c is fractional change to PRC Gouy phase; c=1.0 leaves PRC Gouy phase unchanged from as-built value
# only a=jitter[0] is used in this function
#
#
# Output:
# real(Gn) - real part of optical transfer coefficient, Eq. 11 from Evans et al.
#
def calc_optical_TF(w,order,IFO,jitter=array([0,0,1.0])):

    # First, define optical parameters for the arm cavities.
    # Arm length, mirror ROC, power transmissivities
    # Use these to calculate Gouy phase

    if IFO=='L':
        # LIGO numbers
        #phi = 2.7169 # this is calculated below
        L = 3994.5
        TA = 0.014
        TB = 1e-5

        # calculate the arm gouy phase independently, using the jitter parameters
        R1 = 1939.3
        R2 = 2241.54 + jitter[0]
        g1 = 1 - L/R1
        g2 = 1 - L/R2
        phi = arccos(-1*sqrt(g1*g2))

        #print phi

    if IFO=='V':
        # Virgo numbers
        #phi = 2.7456 # this is calculated below
        L = 2999.8
        TA = 0.0137
        TB = 1e-5

        # calculate the arm gouy phase independently, using the jitter parameters
        R1 = 1424.57
        R2 = 1695.5 + jitter[0]
        g1 = 1 - L/R1
        g2 = 1 - L/R2
        phi = arccos(-1*sqrt(g1*g2))

        #print phi

    tA = sqrt(TA)
    tB = sqrt(TB)

    rA = sqrt(1-TA)
    rB = sqrt(1-TB)

    e4 = array([0,0,0,1,0])

    
    # calculate the propagation operator for the optical mode in the cavity (depends on Gouy phase and mode frequency)
    def prop_op(phi0, phiG, w, L, order, sign):
        if sign=='p':
            return exp(1j*((phi0 - order*phiG) + w*L/constants.c) )
        if sign=='m':
            return exp(1j*((phi0 - order*phiG) - w*L/constants.c) )

    pL_plus = prop_op(0.0,phi,w,L,order,'p')
    pL_minus = prop_op(0.0,phi,w,L,order,'m')


    # populate the scattering matrices
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

    # calculate the optical transfer coefficients
    Gn_plus = dot(e4.T,dot(linalg.inv(identity(5)-Sn_plus),e4))
    Gn_minus = dot(e4.T,dot(linalg.inv(identity(5)-Sn_minus),e4))

    Gn = Gn_minus - Gn_plus.conj()

    #print w/(2*pi), Gn_plus, Gn_minus, Gn.real

    return Gn.real





# calc_optical_TF_fullIFO -- function to return optical transfer coefficient for a complete DRFPMI,
# for a particular optical mode frequency. Implements Section 4 of Evans et al arXiv:0910.2716 (Eqs 15-17)
#
# Inputs:
# w   - angular frequency of mechanical mode (float)
# order - order of optical mode (int)
# IFO - 'L' or 'V' (advanced LIGO or advanced Virgo arm cavities) (string)
# jitter - parameters for small random variation of arm (3-element array)
#
# jitter is of the form array([a,b,c])
# a is change to ETMX radius of curvature in meters
# b is change to ETMY ROC in meters
# c is fractional change to PRC Gouy phase; c=1.0 leaves PRC Gouy phase unchanged from as-built value
# For now we don't change the SRC Gouy phase
#
#
# Output:
# real(Gn) - real part of optical transfer coefficient, Eq. 11 from Evans et al.
#
def calc_optical_TF_fullIFO(w,order,IFO,jitter=array([0,0,1.0])):

    # define the mirror power transmissivities

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
        Tsr = 0.37  # 0.2, someday
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

    # set up the mirror matrix for the full IFO
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

    # Now, calculate the propagation operators, using the cavity Gouy phases
    if IFO=='V':

        # Virgo numbers
        L1234 = 2999.8
        L56 = 6.142
        L78 = 5.912
        L910 = 5.925
        L1112 = 5.925

        # calculate the arm gouy phase independently, using the jitter parameters
        # for now use a constant ROC for the input mirrors (should jitter this eventually)
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

    if IFO=='L':
        
        # LIGO numbers
        L1234 = 3994.5
        L56 = 4.85
        L78 = 4.9
        L910 = 52.3
        L1112 = 50.6

        # calculate the arm gouy phase independently, using the jitter parameters
        # 1934, 2245

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
        phiG910 = -0.44 * jitter[2]
        phiG1112 = -0.35

        # Note: the Gouy phase of the aLIGO PRC and SRC is negative (?)
        # This returns the same full-IFO response function as in the Evans et al paper (Fig 7)
        # Since AdVirgo has such a small Gouy phase the sign shouldn't matter (?)


    # Now calculate the propagation operators and populate the scattering matrix
    def prop_op(phi0, phiG, w, L, order, sign):
        if sign=='p':
            return exp(1j*((phi0 - order*phiG) + w*L/constants.c) )
        if sign=='m':
            return exp(1j*((phi0 - order*phiG) - w*L/constants.c) )


    # Is it possible to pre-calculate these arrays and call them from a lookup table?
    # Might speed up the calculations

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
    

    # Should do some code profiling to see how much these steps dominate the calculation time

    Sn_plus = dot(M,p_plus)
    Sn_minus = dot(M,p_minus)

    Gn_plus = dot(e4.T,dot(linalg.inv(identity(12)-Sn_plus),e4))
    Gn_minus = dot(e4.T,dot(linalg.inv(identity(12)-Sn_minus),e4))

    Gn = Gn_minus - Gn_plus.conj()

    #print w/(2*pi), Gn_plus, Gn_minus, Gn.real

    return Gn.real

