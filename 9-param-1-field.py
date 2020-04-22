from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

'SOLVE ODE IN PARAMETER SPACE, NOT ON MANIFOLD!!!!!!'

# parametrization range

range_param = 3

u_inf = 0
u_sup = 2*3.14
v_inf = 0
v_sup = 6.28

# Risoluzione della superficie
surface_res = 50

plot_basis = False

# Numero di vettori plottati = field_res^2
field_res = 0

# queste curve integrali corrispondono a un vettore tangente del campo plottato
# numero punti iniziali plottati = questo parametro --> AL QUADRATO <--
# distribuiti su griglia regolare, dove c'è un vettore del flusso plottato
number_initial_points = 0

# queste curve integrali hanno un punto iniziale a caso sulla superficie
# numero di punti iniziali con relativa curva plottate = questo parametro
number_extra_initial_points = 3

# dominio di curve integrali e accuratezza di ODE solver
time=np.linspace(0,15,1000)

#########################################################
u, v = symbols('u v')
#### Symbolic parametric surface <--------------------------------------------------------------- HERE PARAM SURFACE

'Waves'
# X_s = u 
# Y_s = v 
# Z_s = 0.1*sin(u)*cos(v)

'Sphere'
# X_s = sin(u) * cos(v)
# Y_s = sin(u) * sin(v)
# Z_s = cos(u)

'Torus'
X_s = ( 2 + cos(u) ) * cos(v)
Y_s = ( 2 + cos(u) ) * sin(v)
Z_s = sin(u)

# derivatives
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


# Generate surface mesh

# surface domain
u_surface = np.linspace(u_inf, u_sup, surface_res)
v_surface = np.linspace(v_inf, v_sup, surface_res)
u_surface, v_surface = np.meshgrid(u_surface, v_surface)
# evaluate
x_surface = FX(u_surface,v_surface)
y_surface = FY(u_surface,v_surface)
z_surface = FZ(u_surface,v_surface)
# Display the mesh
fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.plot_surface(x_surface, y_surface, z_surface, cmap = 'viridis', rstride = 1, cstride = 1, alpha = 0.5)



# basis tangent vector fields
def DU(u,v):
	return np.array([duFX(u,v), duFY(u,v), duFZ(u,v)])

def DV(u,v):
	return np.array([dvFX(u,v), dvFY(u,v), dvFZ(u,v)])

# domain
u_vector = np.linspace(u_inf, u_sup, field_res)
v_vector = np.linspace(u_inf, u_sup, field_res)
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
# vectors domain on surface
# <--------------------------------------------------------------------------------- VECTOR FIELD HERE
#NICE ONE
def cu(u,v):
	return v-u 
def cv(u,v):
	return -2*(u-1)

# cannot be written just as cu * DU + cv * DV for techincal plotting reasons (need some matrix product)
# but conceptually it is exactly the linear combination of the basis vector fields DU, DV with coefficients cu, cv
def Z(u,v):
	return np.array([cu(u,v) * duFX(u,v), cu(u,v) * duFY(u,v), cu(u,v) * duFZ(u,v)]) + np.array([cv(u,v) * dvFX(u,v), cv(u,v) * dvFY(u,v), cv(u,v) * dvFZ(u,v)])

Z1, Z2, Z3 = Z(u_vector, v_vector)
ax.quiver(x_vector, y_vector, z_vector, Z1, Z2, Z3, length=0.2, color='red')

# <------------------------------------------------------------------------------------------------------- 2. Symbolic ODE solver?

# input for ode solver, exactly as Z, just formatted to be used in odenit
def Z_ode(VAR, t):
	u,v = VAR
	return [cu(u,v), cv(u,v)]

# <------------------------------------------------------------------------------------------------------- 3. Check numeric ODE solver

# Takes only INITIAL POINT as argument and plots integral curve of Z through that point
def plot_integral_curve(P, col = 'purple', lw = 2):
	data=odeint(Z_ode, P, time)
	u_curve, v_curve = data[:,0],data[:,1]
	curve1, curve2, curve3 = FX(u_curve, v_curve), FY(u_curve, v_curve), FZ(u_curve, v_curve)
	ax.plot(curve1, curve2, curve3, color = col, linewidth = lw)
	ax.plot([FX(P[0],P[1])], [FY(P[0],P[1])], [FZ(P[0],P[1])], 'ko', ms = 5)

# curve integrali con punti iniziali dove c'è plottato un vettore
for u in np.linspace(u_inf, u_sup, field_res):
	for v in np.linspace(v_inf, v_sup, field_res):
		P = [u,v]
		plot_integral_curve(P, col = 'orange')

uu = np.linspace(u_inf, u_sup, surface_res)
vv = np.linspace(v_inf, v_sup, surface_res)
# altre curve integrali con punti iniziali a caso sulla superficie
for i in range(number_extra_initial_points):
	j = rm.randint(0, surface_res-1)
	k = rm.randint(0, surface_res-1)
	P = [uu[i], vv[j]]
	plot_integral_curve(P)

ax.set_xlim(-range_param,range_param)
ax.set_ylim(-range_param,range_param)
ax.set_zlim(-range_param,range_param)
plt.show()

