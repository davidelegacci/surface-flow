from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

'''

Plotta campo vettoriale on 2-dimensional symplex. 

'''

# x parameter interval domain
u_inf = 0
u_sup = 1

# number of plotted points = this number^2
surface_res = 3

# the y parameter range depends on the x one
def y_bound(x):
	return 1-x # 1-x for 2-simplex


plot_basis = True

# Numero di vettori plottati = field_res^2
field_res = 3

# queste curve integrali corrispondono a un vettore tangente del campo plottato
# numero punti iniziali plottati = questo parametro --> AL QUADRATO <--
# distribuiti su griglia regolare, dove c'è un vettore del flusso plottato
number_initial_points = 3

# queste curve integrali hanno un punto iniziale a caso sulla superficie
# numero di punti iniziali con relativa curva plottate = questo parametro
number_extra_initial_points = 5

# dominio di curve integrali e accuratezza di ODE solver
time=np.linspace(0,15,1000)

#########################################################
u, v = symbols('u v')
#### Symbolic parametric surface <--------------------------------------------------------------- HERE SIMPLEX EQUATION
X_s = u 
Y_s = v 
Z_s = 1-(u+v)

#########################################################

###########################################################################
#################### DEFINE COMPONENTS OF VECTOR FIELD ####################
# Scalar functions of u,v

# def cu(u,v):
# 	return 0.1*u*(1-u) # <--------------------------------------------------------------------------------- VECTOR FIELD HERE

# def cv(u,v):
# 	return 0.2*v*(1-v)  # <--------------------------------------------------------------------------------- VECTOR FIELD HERE

#NICE ONE
def cu(u,v):
	return v-u 
def cv(u,v):
	return -2*(u-1)

###########################################################################

# there is symbolic jacobian but ok adesso a mano <-------------------------------------------- 1. Improve with Jacobian here

du_X_s = diff(X_s, u)
du_Y_s = diff(Y_s, u)
du_Z_s = diff(Z_s, u)

dv_X_s = diff(X_s, v)
dv_Y_s = diff(Y_s, v)
dv_Z_s = diff(Z_s, v)

# make functions python
FX = lambdify((u,v), X_s)
FY = lambdify((u,v), Y_s)
FZ = lambdify((u,v), Z_s)

duFX = lambdify((u,v), du_X_s)
duFY = lambdify((u,v), du_Y_s)
duFZ = lambdify((u,v), du_Z_s)

dvFX = lambdify((u,v), dv_X_s)
dvFY = lambdify((u,v), dv_Y_s)
dvFZ = lambdify((u,v), dv_Z_s)



# ------------ PLOT SIMPLEX
x_surface = np.linspace(u_inf, u_sup, surface_res)

X_surface = [x*np.ones(len(x_surface)) for x in x_surface]
Y_surface = [ np.linspace(0, y_bound(x), surface_res) for x in x_surface ]
Z_surface = [ [FZ (x,y) for y in Y ] for x,Y in zip(x_surface, Y_surface) ]

X_plot = []
for el in X_surface:
	X_plot += list(el)

Y_plot = []
for el in Y_surface:
	Y_plot += list(el)

Z_plot = []
for el in Z_surface:
	Z_plot += list(el)


fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.plot_trisurf(X_plot, Y_plot, Z_plot, antialiased=True, alpha = 0.5)

# --------------------------------------------------------------------------------------------------------------------------------------------

# basis tangent vector fields
def DU(u,v):
	return np.array([duFX(u,v), duFY(u,v), duFZ(u,v)])

def DV(u,v):
	return np.array([dvFX(u,v), dvFY(u,v), dvFZ(u,v)])

# vectors domain on surface

# domain
u_vector = np.linspace(u_inf, u_sup, field_res)
v_vector = np.linspace(u_inf, u_sup, field_res)
#------------------------------------------------------------------  START DIRTY TRICK 
for i in range(len(u_vector)):
	if u_vector[i] + v_vector[i] > 1:
		u_vector = np.delete(u_vector,i)
		v_vector = np.delete(v_vector,i)
#------------------------------------------------------------------ END DIRTY TRICK 
u_vector, v_vector = np.meshgrid(u_vector, v_vector)
# evaluate
x_vector = FX(u_vector,v_vector)
y_vector = FY(u_vector,v_vector)
z_vector = FZ(u_vector,v_vector)

# evaluate basis vector fields to plot span of tangent space at every point
DU_1, DU_2, DU_3 = DU(u_vector, v_vector)
DV_1, DV_2, DV_3 = DV(u_vector, v_vector)

# plot basis vectors  # <--------------------------------------------------------------------------------- PLOT BASIS
if plot_basis == True:
	ax.quiver(x_vector, y_vector, z_vector, DU_1, DU_2, DU_3, length=0.4, color='k')
	ax.quiver(x_vector, y_vector, z_vector, DV_1, DV_2, DV_3, length=0.2, color='k')

# THE vector field whose flow is studied
# cannot be written just as cu * DU + cv * DV for techincal plotting reasons (need some matrix product)
# but conceptually it is exactly the linear combination of the basis vector fields DU, DV with coefficients cu, cv
def Z(u,v):
	return np.array([cu(u,v) * duFX(u,v), cu(u,v) * duFY(u,v), cu(u,v) * duFZ(u,v)]) + np.array([cv(u,v) * dvFX(u,v), cv(u,v) * dvFY(u,v), cv(u,v) * dvFZ(u,v)])

Z1, Z2, Z3 = Z(u_vector, v_vector)
ax.quiver(x_vector, y_vector, z_vector, Z1, Z2, Z3, length=0.2, color='red')


# input for ode solver, exactly as Z, just formatted to be used in odenit
def Z_ode(VAR, t):
	u,v,w = VAR
	return np.array([cu(u,v) * duFX(u,v), cu(u,v) * duFY(u,v), cu(u,v) * duFZ(u,v)]) + np.array([cv(u,v) * dvFX(u,v), cv(u,v) * dvFY(u,v), cv(u,v) * dvFZ(u,v)])


# <-------------------------------------------------------------------------------------------------------  numeric ODE solver

# Takes only INITIAL POINT as argument and plots integral curve of Z through that point
def plot_integral_curve(P, col = 'purple', lw = 2):
	data=odeint(Z_ode, P, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	ax.plot(curve1, curve2, curve3, color = col, linewidth = lw)
	ax.plot([P[0]], [P[1]], [P[2]], 'ko', ms = 5)

# curve integrali con punti iniziali dove c'è plottato un vettore
for i in np.linspace(0, len(x_vector)-1, number_initial_points):
	for j in np.linspace(0, len(x_vector)-1, number_initial_points):
		i = int(i)
		j = int(j)
		P = [ x_vector[i][j], y_vector[i][j], z_vector[i][j] ]
		plot_integral_curve(P, col = 'orange')

# altre curve integrali con punti iniziali a caso sulla superficie
for i in range(number_extra_initial_points):
	j = rm.randint(0, len(X_plot)-1)
	r1 = X_plot[j]
	r2 = Y_plot[j]
	r3 = Z_plot[j]
	P = [r1, r2, r3]
	plot_integral_curve(P)

# ax.set_xlim(-range_param,range_param)
# ax.set_ylim(-range_param,range_param)
# ax.set_zlim(-1,1)
plt.show()

