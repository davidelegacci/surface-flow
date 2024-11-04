import matplotlib.pyplot as plt
import numpy as np

import sys

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d


import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

# ROTAZIONE, per teoria vedi

# /Users/davidelegacci/coding/python/sympy/LinearAlgebra/3d-rotation.html 

# grid where vectors are plotted, chill, can grow, non-critical
vector_range = 1

# number of plotted vectors = this number ^ 3. chill, can grow, non-critical
vector_res = 0

surface_res = 100

ode_res = 2000

# not important, just time of 2 integral curves plotted
integral_curve_end_time = 1000

# Parametric surface

fig = plt.figure()
ax = fig.add_subplot(projection = '3d')

t_surface = np.linspace(0, np.pi, surface_res)
s_surface = np.linspace(0, 2*np.pi, surface_res)

t_surface, s_surface = np.meshgrid(t_surface, s_surface)

x_surface = np.sin(t_surface) * np.cos(s_surface)
y_surface = np.sin(t_surface) * np.sin(s_surface)
z_surface = np.cos(t_surface)

ax.plot_surface(x_surface, y_surface, z_surface, cmap = 'viridis', rstride = 1, cstride = 1, alpha = 0.5)


# <------------------------------------------------------------------------------- CAMPI, UNICA COSA CONCETTUALE DA EDITARE

A = 1
B = 1
C = 1


def X(x,y,z):
	# return [A*y+B*z, -A*x+C*z, -B*x-C*y]
	return [x*y*z, x*x*z, -2*x*x*y]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ 
# PLOT AXIS

# axis_time = np.linspace(-1.1,1.1,100)
# x_axis_generator = C 
# y_axis_generator = -B 
# z_axis_generator = A 

# x_axis = x_axis_generator * axis_time
# y_axis = y_axis_generator * axis_time
# z_axis = z_axis_generator * axis_time

# ax.plot(x_axis, y_axis, z_axis)

# # axis points on unit sphere, fp = fixpoint
# z_fp_1 = (C**2 / A**2 + B**2 / A**2 + 1)**(-0.5)
# x_fp_1 = C * z_fp_1 / A
# y_fp_1 = -B * z_fp_1 / A

# z_fp_2 = -z_fp_1
# x_fp_2 = C * z_fp_2 / A
# y_fp_2 = -B * z_fp_2 / A

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ 
# PLOT SCATTER DI VETTORI
x_vector, y_vector, z_vector = np.meshgrid(	np.linspace(-vector_range, vector_range, vector_res),
											np.linspace(-vector_range, vector_range, vector_res),
						  					np.linspace(-vector_range, vector_range, vector_res))

X_plot = X(x_vector, y_vector, z_vector)


X_color = 'red'


ax.quiver(x_vector, y_vector, z_vector, X_plot[0], X_plot[1], X_plot[2], length = 0.2, color = X_color)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# make point on unit 2-sphere
def mk_p(a,b):
	cq = 1 - a**2 - b**2

	if cq < 0:
		raise Exception('nope')

	c = np.sqrt(cq)
	return [a,b,c]

POINTS = [mk_p(0.5,0.8), mk_p(0.3,0.2), mk_p(-0.1,0.4), mk_p(0.2,0.3), mk_p(0.3,0.1), mk_p(-0.4,0.6)] # some free points
# POINTS = [mk_p(1,0), mk_p(0,1), mk_p(0,0)] # some free points
# POINTS.append( [x_fp_1, y_fp_1, z_fp_1] )
# POINTS.append( [x_fp_2, y_fp_2, z_fp_2] )

# add axis points


def X_ode(VAR, t):
	x,y,z = VAR
	return X(x,y,z)

def plot_integral_curve(field, P, col = 'red'):
	time=np.linspace(0,integral_curve_end_time,ode_res)
	data=odeint(field, P, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	ax.plot(curve1, curve2, curve3, color = col)
	ax.plot([P[0]], [P[1]], [P[2]], 'ko', ms = 5)

for p in POINTS:
	plot_integral_curve(X_ode, p, X_color)


print(POINTS)
plt.show()





