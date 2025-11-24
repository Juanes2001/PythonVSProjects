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
p = 10e-9
q = 2*math.pi/p

## Parametros de entrada
delta_ab = 0.001
epsilon_ab = 1.41**2


r = CoordSys3D('r')
E = CoordSys3D('E')
z,t = sp.symbols('z t')
Ax = sp.Function("Ax")
Ay = sp.Function("Ay")

er_tensor = sp.Matrix([[epsilon_ab + delta_ab*sp.cos(2*q*z), -delta_ab*sp.sin(2*q*z)],
                       [-delta_ab*sp.sin(2*q*z), epsilon_ab - delta_ab*sp.cos(2*q*z)]])

# Definición en forma simbólica de la forma incial de la onda electromagnética la cual tendremos que solucionar en sus variables.

def Evec(t,omega,r_vec,k_vec):

    temporal_phase = sp.exp(sp.I*(omega*t))
    E0_vec = (Ax(z)*E.i + Ay(z)*E.j)*sp.exp(-sp.I*(r_vec.dot(k_vec)))

    return E0_vec * temporal_phase

# Definición en forma simbólica del operador rotacion S(theta)

def S(phi):# al ser un input negativo obtenemos S-1,no necesitamos definirlo aparte

    matr = sp.Matrix( [[sp.cos(phi), sp.sin(phi)],
                       [-sp.sin(phi), sp.cos(phi)]] )
    return matr

def Matx_vec_mult(M,V):

    M = sp.sympify(M)
    V = sp.sympify(V)

    vec_matr =  sp.Matrix(V.to_matrix(r))

    mult = M * vec_matr
    re_vec = mult[0]*r.i + mult[1]*r.j + mult[2]*r.k 

    return re_vec

def solve(A,):




