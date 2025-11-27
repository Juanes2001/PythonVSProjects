import functions as fun

import sympy as sp

import sympy.vector as vec
from sympy.vector import CoordSys3D

from IPython.display import display, Math
import scipy.integrate as integ
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math



Ax = sp.Integer(1)
Ay = -sp.I

# Onda incidente y onda reflejada
E_vec_i = 1/math.sqrt(2)*fun.Evec(Ax,Ay,0,omega = 2*math.pi*fun.c/fun.lambda0,dir = 1)
E_vec_r = 1/math.sqrt(2)*fun.Evec(Ax,-Ay,0,omega = 2*math.pi*fun.c/fun.lambda0, dir = -1)

# 

print(fun.find_amplitudes(1.41,1))

