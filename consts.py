# this file defines constants, the material properties of the liquid crystal,
# how far to set the bounds of computation, the starting point, and the zero round-off point

#import eulerfuncs
from sympy import *
import numpy as np

# material parameters
#LAMBDA = 0.6
#NU = 0.8

# values of actuation:
Lambda = [0.0, 0.1, 0.2]  # [0.0, 1.0/3.33, 2.0/3.33]

# actual expansion/contraction factors depend on u,v:

# shouldn't LAMBDA2 i.e. the v expansion be the one to have target strain? whereas u curves don't expand


def LAMBDA1(Lambda, oldpoint):
    # epsilon ratio material to cavity
    # target_strain_plus_one(Lambda, oldpoint.epsilon)
    return 1 + Lambda*oldpoint.epsilon


def LAMBDA2(Lambda, oldpoint):
    return 1


E_val = 1.5
nu_val = 0.5  # at 1/2 no expansion along air channels
# d_width = 1.2 not needed for epsilon as ratio


def PHI_var(epsilon):
    # d/d+dw i.e. ratio of cavity to total spacing
    return epsilon


# looking at the measurement tables, both were measured/acceptable around this value
PSI_var = 0.66  # 0.8  # 0.66  # 0.75  # 0.66


def target_strain_plus_one(Lambda, epsilon):  # lambda -> pressure
    return 1.0 + (Lambda/E_val)*(2.0-PHI_var(epsilon))*PHI_var(epsilon)*(
        (PSI_var/(1-PSI_var)) - (nu_val*PHI_var(epsilon)*PSI_var/(1-PSI_var*PHI_var(epsilon)))*(
            1 + nu_val*((1-PHI_var(epsilon))/(PHI_var(epsilon)*(1-PSI_var)))))


# initial phi (i.e. v) and theta (i.e. u) on the starting surface map
# that is, we set the u,v coordinates to start corresponding to theta,phi given below
V0 = 0  # np.pi / 2
U0 = 0

# This is the bounds for distances on the surface in terms of theta,phi.
# For example, to cover a quadrant of the sphere we can have both vary pi/2 (from pi/2,0 to pi,pi/2, to pi,-pi/2, to 0,pi/2, to 0,-pi/2 such that there is no overlap)
# pi is just arbitrary, and for parabaloids there is no bound on theta,phi
BOUNDS = pi * 0.4  # 0.7

# length in u,v to take per step for the euler approximation
STEP = 0.02

# how many euler steps to perform in order to move BOUNDS distance
SIZE = int(BOUNDS/STEP)  # ie in uv

# the actual size of a step in u and in v yields a different distance in theta and phi
# no longer can exactly cover each surface since multiple, but it's ok
# int(BOUNDS/(STEP*LAMBDA))
# # make u and v different length?
THETASIZE = int(0.7*BOUNDS/(STEP))
PHISIZE = int(1.3*BOUNDS/(STEP))  # int(BOUNDS/(STEP*CONTRACT))

# 0 approximation
DELTA = 0.2
