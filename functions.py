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
c = 299792458

## Parametros de entrada
delta_ab = 0.001
epsilon_ab = 1.41**2
lambda0 = 500E-9
k0 = 2*math.pi / lambda0
alf = lambda0/p


nx = math.sqrt(epsilon_ab+delta_ab)
ny = math.sqrt(epsilon_ab-delta_ab)


r = CoordSys3D('r')
z,t = sp.symbols('z t')

# Con esta funcion definimos los autovalores dependientes del valor de lambda
def N_(lam,type,dir):

    if type == 1 and dir == 1:
        auto_vals = + ((lam/p)**2 + epsilon_ab + 
                          math.sqrt(4*(lam/p)**2 * epsilon_ab + delta_ab**2))**(0.5) + 1j*0
    elif type == 1 and dir == -1:
        auto_vals = - ((lam/p)**2 + epsilon_ab + 
                          math.sqrt(4*(lam/p)**2 * epsilon_ab + delta_ab**2))**(0.5) + 1j*0
    elif type == 2 and dir == 1:
        auto_vals = + ((lam/p)**2 + epsilon_ab - 
                          math.sqrt(4*(lam/p)**2 * epsilon_ab + delta_ab**2))**(0.5) + 1j*0
    elif type == 2 and dir == -1:
        auto_vals = - ((lam/p)**2 + epsilon_ab - 
                          math.sqrt(4*(lam/p)**2 * epsilon_ab + delta_ab**2))**(0.5) + 1j*0
    return auto_vals


# Definimos con esta funcion los autovectores representativos como base
# dependientes del valor que adopte n como entrada
def Aprime_0(n):

    A0_vec = ((1*r.i + -sp.I * (nx**2-n**2-(alf)**2)/(2*alf*n)*r.j)
                *sp.exp(-sp.I*(z*k0*n)))

    return A0_vec 

# Definición en forma simbólica de la forma incial de la onda electromagnética 
# incidente

def Evec(Ax,Ay,t,omega,dir):

    temporal_phase = sp.exp(sp.I*(omega*t))
    if dir == +1:    
        E0_vec = (Ax*r.i + Ay*r.j)*sp.exp(-sp.I*(z*k0))
    elif dir == -1:
        E0_vec = (Ax*r.i + Ay*r.j)*sp.exp(sp.I*(z*k0))
    return E0_vec * temporal_phase

# Definimos una función que halle las amplitudes de entrada con respecto a las amplitudes de los 
#Eigenmodos

def find_amplitudes(n_cnc,n_medium):

    w = 1j*((epsilon_ab+delta_ab)-n_cnc**2-alf**2)/(2*alf*n_cnc)

    u,r,v1,v2 = sp.symbols('u r v1 v2')

    f1 = sp.Eq(u+r                       ,v1 + v2)
    f2 = sp.Eq(u-r                       ,-1j*w*v1 + 1j*w*v2)
    
    sol = sp.solve((f1,f2),(u,r))
    
    return sol


# Definición en forma simbólica del operador rotacion S(theta)

def S(phi):
    # al ser el input negativo obtenemos S-1,no necesitamos definirlo aparte

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



