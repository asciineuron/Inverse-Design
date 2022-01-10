# this file takes the 3 supporting files and computes for a desired start
# and end surface, along with u v curves on the surface from which to start (note the code only supports "straight" lines in theta/phi)

import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from mayavi import mlab
import math

from consts import *
import eulerfuncs
import mathfuncs
import plotting

u, v, alpha, beta, phi, theta, t, s = symbols(
    'u v alpha beta phi theta t s')


# --- CONTROLS ---
# define start and end surfaces, here from a wide to narrow parabaloid:
# here we follow the paper, but *starting* flat and then going to the sphere and anisotropic gaussian:
R0_symb = Matrix([theta, phi, 0])
# R1_symb = Matrix([sin((1*phi + np.pi/2))*cos(1*theta),
#                   sin(1*theta)*sin((1*phi + np.pi / 2)),
#                   cos((1*phi + np.pi / 2))])
R1_symb = Matrix([sin((1*theta + np.pi/2))*cos(1*phi),
                  sin(1*phi)*sin((1*theta + np.pi / 2)),
                  cos((1*theta + np.pi / 2))])

# at p=0.1, intermediate sphere radius 1.3269 from curvature=0.5679 (assuming single surface calculations correct??)
# R2_symb = Matrix([1.3269*sin((1*theta + np.pi/2))*cos(1*phi),
#                   1.3269*sin(1*phi)*sin((1*theta + np.pi / 2)),
#                   1.3269*cos((1*theta + np.pi / 2))])
R2_symb = Matrix([1.5*sin((1*theta + np.pi/2))*cos(1*phi),
                  1.4*sin(1*phi)*sin((1*theta + np.pi / 2)),
                  1.6*cos((1*theta + np.pi / 2))])
# R1_symb = Matrix([sin((theta + np.pi/2))*cos(phi),
#                   sin(phi)*sin((theta + np.pi / 2)),
#                   cos((theta + np.pi / 2))])
# R1_symb = Matrix([sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)])
# R2_symb = Matrix([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])
# R2_symb = Matrix([theta, phi, exp((theta*theta)/4 + phi*phi)])
# R2_symb = Matrix([2*sin((1*phi + np.pi/2))*cos(1*theta),
#                   2*sin(1*theta)*sin((1*phi + np.pi / 2)),
#                   2*cos((1*phi + np.pi / 2))])  # Matrix([phi, theta, exp(-((1*phi*phi)/4 + 2*theta*theta))])
#R2_symb = Matrix([theta, phi, theta**3 - 3*theta*phi**2])

# define initial curves in omega space:
curve1u = t
curve1v = 0*t

curve2u = 0*t
curve2v = t

# set starting point in omega space, found in consts:
# --- END CONTROLS ---

# create the math object from these data:
mathobj = mathfuncs.MathObj(
    [R0_symb, R2_symb, R1_symb], curve1u, curve1v, curve2u, curve2v)


# calculate for each quadrant, then plot

Data1 = eulerfuncs.euler_procedure(mathobj, 1)  # eulerfuncs.gen_empty_Data(1)
# # Data1  # eulerfuncs.euler_procedure(mathobj, 2)
#Data2 = eulerfuncs.euler_procedure(mathobj, 2)
# # # Data1  # eulerfuncs.euler_procedure(mathobj, 3)
#Data3 = eulerfuncs.euler_procedure(mathobj, 3)
# # # Data1  # eulerfuncs.euler_procedure(mathobj, 4)
#Data4 = eulerfuncs.euler_procedure(mathobj, 4)
# # eulerfuncs.json_out(mathobj, Data1)
# # display:
plotting.plot_3d_4_quadrants_2_surfaces_epsilon(
    mathobj, Data1, Data1, Data1, Data1)
# # # write to csv:
#plotting.diagnostic_epsilon_csv_xy(Data1, 1)
# # plotting.spaced_epsilon_csv_xy(Data2, 2)
# # plotting.spaced_epsilon_csv_xy(Data3, 3)
# # plotting.spaced_epsilon_csv_xy(Data4, 4)
plotting.balanced_epsilon_csv_xy(Data1, 1, PHISIZE)
# plotting.balanced_epsilon_csv_xy(Data2, 2, PHISIZE)
# plotting.balanced_epsilon_csv_xy(Data3, 3, PHISIZE)
# plotting.balanced_epsilon_csv_xy(Data4, 4, PHISIZE)
# plotting.v_curve_balanced_epsilon_csv_xy(Data1, 1, THETASIZE)
