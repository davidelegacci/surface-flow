from sympy import *
from sympy.diffgeom import *

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d

import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

######### Plotting region 
plot_range = 2
param_range = 5

# risoluzione delle foglie
surface_res = 40

# number of plotted vectors = this number ^ 3
vector_res = 2
# number of initial points for integral curves  = this number ^ 3
curves_res = 1

time=np.linspace(0,6,500)

# Set up manifold

M = Manifold('M', 3)
patch = Patch('P', M)
phi = CoordSystem('phi', patch, ['x','y', 'z'])

# Test point

# coord function
x,y,z = phi.coord_functions()


# Basis vector fields
Dx, Dy, Dz = phi.base_vectors()

# Basis oneforms
dx, dy, dz = phi.base_oneforms()


# <------------------------------------------------------------------------------- CAMPI

X1 = y+z
X2 = -x+z
X3 = -x-y

Y1 = y+z
Y2 = -x
Y3 = -x
################

X_s = X1 * Dx + X2 * Dy + X3 * Dz
Y_s = Y1 * Dx + Y2 * Dy + Y3 * Dz

Z = Commutator(X_s,Y_s)
# pprint(Z)
# if Z  == 0:
# 	print('COMMUTANO')
# else:
# 	raise Exception('NON COMMUTANO')

C1 = [X1, Y1]
C2 = [X2, Y2]
C3 = [X3, Y3]

colorX = 'red'
colorY = 'green'

color_vf = [colorX, colorY]

fig = plt.figure()
ax = fig.gca(projection = '3d')

def plot_integral_curve(P, lw = 3):
	data=odeint(CV_ode, P, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	ax.plot(curve1, curve2, curve3, color = used_color, linewidth = lw)
	ax.plot([P[0]], [P[1]], [P[2]], 'ko', ms = 5)


for i in range(len(C1)):

	used_color = color_vf[i]

	c1 = lambdify((x,y,z), C1[i])
	c2 = lambdify((x,y,z), C2[i])
	c3 = lambdify((x,y,z), C3[i])

	a, b, c = np.meshgrid(np.linspace(-plot_range, plot_range, vector_res),
						  np.linspace(-plot_range, plot_range, vector_res),
						  np.linspace(-plot_range, plot_range, vector_res))
	
	#  < -------------------------------------------------------------------------------------------------- VECTOR FIELD PLOT
	ax.quiver(a,b,c, c1(a,b,c), c2(a,b,c), c3(a,b,c), length = 0.2, color = color_vf[i])

	def CV_ode(VAR,t):
		u, v, w = VAR
		return [c1(u,v,w), c2(u,v,w), c3(u,v,w)]

	for i in np.linspace(0, vector_res-1, curves_res):
		for j in np.linspace(0, vector_res-1, curves_res):
			for k in np.linspace(0, vector_res-1, curves_res):
				i = int(i)
				j = int(j)
				k = int(k)
				P = [ a[i][j][k], b[i][j][k], c[i][j][k] ]
				plot_integral_curve(P) #  < ------------------------------------------------------------------ INT CURVE  PLOT

# now P is the last used point

# symbolic point, FIXED! determina foglia
P_s = phi.point([P[0], P[1], P[2]])

from sympy.abc import t,s

def INT_CURVE(VF,init_point, variable = t):

	flow_system, init_cond = intcurve_diffequ(VF, variable, init_point)

	f_0 = Function('f_0')
	f_1 = Function('f_1')
	f_2 = Function('f_2')
	ic = {f_0(0):x(init_point), f_1(0):y(init_point), f_2(0):z(init_point)}
	
	# Symbolic solution
	solution = dsolve(flow_system, ics = ic)
	
	x_flow_symb = solution[0].rhs
	y_flow_symb = solution[1].rhs
	z_flow_symb = solution[2].rhs
	
	# Python functions solution to use as part of leaf parametrization
	# x_flow = lambdify(t, x_flow_symb)
	# y_flow = lambdify(t, y_flow_symb)
	# z_flow = lambdify(t, z_flow_symb)


	return x_flow_symb, y_flow_symb, z_flow_symb

# foo, bar, lol = INT_CURVE(X_s, P_s)


'''
Funziona, ma riceve T e S da subito, quindi per ogni valutazione (x,y,z) = PSI(T,S) deve ricalcolare le curve integrali! Impensabile plottare una surface
INVECE la versione sotto, ad oggetti, fa lo stesso tenendo T,S come placeholders :) quindi sputa fuori una funzione di T,S che si può valutare al volo!
'''
# def PSI(T, S):

# 	'Parametrizzazione delle foglie di Frobenius (x, y, z) = PSI (T,S)'

# 	'T, S numeri, P punto simbolico'

# 	theta_X_x, theta_X_y, theta_X_z = INT_CURVE(X_s, P_s) # these are functions of t

# 	# the new initial point is the point along the integral curve of X_s at parameter distance T from P
# 	Q_s = phi.point([theta_X_x.subs(t,T), theta_X_y.subs(t,T), theta_X_z.subs(t,T)])

# 	theta_Y_x, theta_Y_y, theta_Y_z = INT_CURVE(Y_s, Q_s) # these are functions of t

# 	final_point = phi.point([theta_Y_x.subs(t,S), theta_Y_y.subs(t,S), theta_Y_z.subs(t,S)])

# 	return x(final_point).evalf(), y(final_point).evalf(), z(final_point).evalf()

# Parametrizzazione delle foglie di Frobenius (x, y, z) = PSI (T,S)'
def PSI():

	theta_X_x, theta_X_y, theta_X_z = INT_CURVE(X_s, P_s, variable = t) # these are functions of t

	Q_s = phi.point([theta_X_x, theta_X_y, theta_X_z])

	theta_Y_x, theta_Y_y, theta_Y_z = INT_CURVE(Y_s, Q_s, variable = s) # these are functions of s

	final_point = phi.point([theta_Y_x, theta_Y_y, theta_Y_z])

	print('')
	print('---PSI() CALLED')
	print('')

	return x(final_point), y(final_point), z(final_point) # queste sono le tre funzioni cercate (x,y,z) = PSI(T,S) !!!

#RUN PSI() ONLY ONCE

frobenius_param = PSI()

FX_sym, FY_sym, FZ_sym = frobenius_param

'''
A questo punto due modi equivalenti per maneggiarle come vere funzioni Python

1. 
foo = lambdify((t,s), FX_sym)
print(foo(1,2))

2.
print(FX_sym.subs(t,1).subs(s,2).evalf())
'''
pprint(FX_sym)
pprint(FY_sym)
pprint(FZ_sym)
# PARAMETRIZZAZIONE DI FROBENIUS
############################
FX = lambdify((t,s), FX_sym)
FY = lambdify((t,s), FY_sym)
FZ = lambdify((t,s), FZ_sym)
############################

# Parametric surface

# Generate surface mesh

# Per costruzione il punto P_s da cui tutto si è costruito
# di coordinate x(P_s), y(P_s), z(P_s)
# corrisponde al parametro (t,s) = (0,0), essendo punto iniziale di
# composizione di flussi
# quindi faccio una mesh grid di t,s intorno al punto (0,0)
# domain
t_surface = np.linspace(-param_range, param_range, surface_res)
s_surface = np.linspace(-param_range, param_range, surface_res)
t_surface, s_surface = np.meshgrid(t_surface, s_surface)
# evaluate
x_surface = FX(t_surface,s_surface)
y_surface = FY(t_surface,s_surface)
z_surface = FZ(t_surface,s_surface)

# Display the leaf 	#  < -------------------------------------------------------------------------------------------------- SURFACE PLOT

ax.plot_surface(x_surface, y_surface, z_surface, cmap = 'viridis', rstride = 1, cstride = 1, alpha = 0.1)


plt.show()







