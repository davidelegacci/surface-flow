from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

# Vector field tangent to SCALAR FIELD i.e. surface that can be written as image of z = f(u,v) so that ( x, y, z ) = ( u, v, f(u,v) )
# This is a particular case of a generic parametric surface ( x, y, z ) = F(u,v)

# This works, but it is a particular case of the parametric surface.
# To achieve a graph surface in the parametric surface just set x = u, y = v, z = f(u,v)
# Hence -----> this is not updated anymore <------; the plotting routines are more advanced in the parametric surface version 

#################### PARAMETERS ###########################################

# Range in domain plane
u_inf = -2
u_sup = 2

# Risoluzione della superficie
surf_res = 30

# Numero di vettori plottati = field_res^2
field_res = 7

# queste curve integrali corrispondono a un vettore tangente del campo plottato
number_initial_points = 6

# queste curve integrali hanno un punto iniziale a caso sulla superficie
number_extra_initial_points = 15

# dominio di curve integrali e accuratezza di ODE solver
time=np.linspace(0,8,300)

###########################################################################
#################### DEFINE THE GEOMERTY  #################################

# These are the only 3 things conceptually to be changed: the surface, and the components of the vector field on it
# These are indeed the components; it is tangent because built as linear combination of local basis
# of vector fields, spanning the tangent space at evert point

###########################################################################
#################### DEFINE SURFACE #######################################

# symbolic variables, they could as well be called u,v, whatever
a,b  = symbols('a b')

# surface as graph of function z(u,v)
# must be a scalar function of the variables a, b
f_s = a#0.2* sin(a) * cos(b)    # <---------------------------------------------------------------------- SURFACE HERE

###########################################################################
#################### DEFINE COMPONENTS OF VECTOR FIELD ####################
# Scalar functions of u,v

def cu(u,v):
	return v  # <--------------------------------------------------------------------------------- VECTOR FIELD HERE

def cv(u,v):
	return -u  # <--------------------------------------------------------------------------------- VECTOR FIELD HERE


#################### (almost) NOTHING TO EDIT BELOEW THIS POINT, CAN RUN  ##########
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Things that can be edited below
# 1. graphs layout (color title ....)
# 2. plot or not basis vector fields, see arrow below

#################### sympy procedure, do not change #######################

# Symbolic partial derivatives
da_f_s = diff(f_s,a)
db_f_s = diff(f_s,b)

# Define usual Python functions from symbolic functions with lambdify

# Surface
f = lambdify((a,b), f_s)
# Partial derivatives
fu = lambdify((a,b), da_f_s)
fv = lambdify((a,b), db_f_s)

#################### CHECK ###########################################
if number_initial_points > field_res:
    raise Exception('number_initial_points should not exceed field_res')


#################### PLOT SURFACE ###########################################

# domain
u_surface = np.outer(np.linspace(u_inf, u_sup, surf_res), np.ones(surf_res))
v_surface = u_surface.copy().T # transpose
# evaluate f on domain
z_surface = f(u_surface,v_surface)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(u_surface, v_surface, z_surface,cmap='viridis', edgecolor='none', alpha = 0.8)
ax.set_title('Surface with tangent vector field and its flow')

#################### VECTOR FIELD ###########################################

# Basis vector fields spanning tangent space

# Domain of vector field: points on the surface
u_vector = np.outer(np.linspace(u_inf, u_sup, field_res), np.ones(field_res))
v_vector = u_vector.copy().T # transpose
z_vector = f(u_vector,v_vector)

# Basis vector field as derivatives of inverse of chart

def DU(u,v):
	return np.array([1, 0, fu(u,v)])

def DV(u,v):
	return np.array([0, 1, fv(u,v)])

# Unpack values of Basis vector field to plot
DU1, DU2, DU3 = DU(u_vector, v_vector)
DV1, DV2, DV3 = DV(u_vector, v_vector)

# Plot two basis fields # <--------------------------------------------------------------------------------- 2. Plot or not basis vector fields
ax.quiver(u_vector, v_vector, z_vector, DU1, DU2, DU3, length=0.2, color='k')
ax.quiver(u_vector, v_vector, z_vector, DV1, DV2, DV3, length=0.2, color='k')

##############################################################
#DISTRIBUTION
##############################################################
# Chosen vector field, linear combination of two generators
# This is THE vector fiels, tangen to the manifold, whose flow is studied
# the product has to be done like this and not simply cu * DU + cv * DV
def Z(u,v): 
	return np.array([cu(u,v), 0, fu(u,v) * cu(u,v)]) + np.array([0, cv(u,v), fv(u,v) * cv(u,v)])

# Unpack values of Z to plot
Z1, Z2, Z3 = Z(u_vector, v_vector)
ax.quiver(u_vector, v_vector, z_vector, Z1, Z2, Z3, length=0.1, color='red')

# Curva integrale ######### ODE SOLVER NUMERICO
# set of initial conditions for integral curves

# one single initial point is of the kind [ x, y, f(u,v) ] 

# Exactly same function as above, just with different structure to be used as input in the ode solver
def Z_ode(VAR, t):
	u,v,z = VAR
	return np.array([cu(u,v), 0, fu(u,v) * cu(u,v)]) + np.array([0, cv(u,v), fv(u,v) * cv(u,v)])

# Takes only INITIAL POINT as argument and plots integral curve of Z through that point
def plot_integral_curve(P, col = 'purple', lw = 2):
	data=odeint(Z_ode, P, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	ax.plot(curve1, curve2, curve3, color = col, linewidth = lw)
	ax.plot([P[0]], [P[1]], [P[2]], 'ko', ms = 5)

# curve integrali con punti iniziali dove c'è plottato un vettore
for i in np.linspace(0, field_res-1, number_initial_points):
	i = int(i)
	P = [ u_vector[i][i], v_vector[i][i], z_vector[i][i] ]
	plot_integral_curve(P)

# altre curve integrali con punti iniziali a caso sulla superficie

for i in range(number_extra_initial_points):
	r1 = rm.uniform(u_inf, u_sup)
	r2 = rm.uniform(u_inf, u_sup, )
	P = [r1, r2, f(r1, r2)]
	plot_integral_curve(P)

ax.set_zlim(-1,1)
plt.show()


