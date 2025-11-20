"""
Se definiran aqui las funciones que se necesiten para el trabajo para ser luego usadas en en la terminal principal
"""


# Librerias a usar
import numpy as np
import sympy as sp

import sympy.vector as vec
from sympy.vector import CoordSys3D

from IPython.display import display, Math
import scipy.integrate as integ
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

# Se definiran los parametros iniciales

u0r = 1
u0 = 4*math.pi* 10e-7
e0r = 1
e0 = 8.854e-12 


r = CoordSys3D('r')
k = CoordSys3D('k')
E = CoordSys3D('E')
Ax,Ay,phix,phiy,t = sp.symbols('Ax Ay phix phiy t')

# Definición en forma simbólica de la forma incial de la onda electromagnética la cual tendremos que solucionar en sus variables.

def Evec(t,omega,k_vec):

    phase = sp.exp(-sp.I*(r.dot(k_vec)+omega*t))
    E0_vec = Ax*sp.exp(sp.I*(phix))*E.i + Ay*sp.exp(sp.I*(phiy))*E.j 

    return E0_vec * phase

# Definición en forma simbólica del operador curl. curl

def curl_curl(E_vec):

    return vec.curl(vec.curl(E_vec))




