import functions as fun

import sympy as sp

import sympy.vector as vec
from sympy.vector import CoordSys3D

from IPython.display import display, Math
import scipy.integrate as integ
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math



A_vec_in = fun.Evec(0,omega = 2*math.pi*fun.c/fun.lambda0)