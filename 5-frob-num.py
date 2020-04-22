import matplotlib.pyplot as plt
import numpy as np

import sys

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d


import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

# NUMERIC Frobenius leaves, cfr. 3-frobenius.py

'''
Given two vector fields plots the surface obtained composing the two flows.
It does not check if they commute or not, this has to be done somewhere else.
It follows X first and Y after; if they commute, it does not matter. 
'''

#########

# grid where vectors are plotted, chill, can grow, non-critical
vector_range = 2

# number of plotted vectors = this number ^ 3. chill, can grow, non-critical
vector_res = 5

# Range of the double flow, noncritical.
surface_range = 2

# risoluzione delle foglie !!--- CRITICAL: number od ODE solved = this param^2 ---!!
surface_res = 100

# Resolution of every ODE solution !!--- CRITICAL: resolution of every single ODE solution ---!!
ode_res = 2000

# not important, just time of 2 integral curves plotted
integral_curve_end_time = 10


# <------------------------------------------------------------------------------- CAMPI, UNICA COSA CONCETTUALE DA EDITARE
# <------------------------------------------------------------------------------- 
# <------------------------------------------------------------------------------- 

def X(x,y,z):
	return [np.sin(x), -x+z, -x-y] #[y+z, -x+z, -x-y]

def Y(x,y,z):
	return [y+z,-x,-x]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------ 


fig = plt.figure()
ax = fig.gca(projection = '3d')


x_vector, y_vector, z_vector = np.meshgrid(	np.linspace(-vector_range, vector_range, vector_res),
											np.linspace(-vector_range, vector_range, vector_res),
						  					np.linspace(-vector_range, vector_range, vector_res))

X_plot = X(x_vector, y_vector, z_vector)
Y_plot = Y(x_vector, y_vector, z_vector)

X_color = 'red'
Y_color = 'green'

ax.quiver(x_vector, y_vector, z_vector, X_plot[0], X_plot[1], X_plot[2], length = 0.2, color = X_color)
ax.quiver(x_vector, y_vector, z_vector, Y_plot[0], Y_plot[1], Y_plot[2], length = 0.2, color = Y_color)

i, j, k = rm.randint(0, vector_res-1), rm.randint(0, vector_res-1), rm.randint(0, vector_res-1)

POINT = [x_vector[i][j][k], y_vector[i][j][k], z_vector[i][j][k]] # <---------------------------------------------- POINT

def X_ode(VAR, t):
	x,y,z = VAR
	return X(x,y,z)

def Y_ode(VAR, t):
	x,y,z = VAR
	return Y(x,y,z)

def plot_integral_curve(field, P, col = 'red'):
	time=np.linspace(0,integral_curve_end_time,ode_res)
	data=odeint(field, P, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	ax.plot(curve1, curve2, curve3, color = col)
	ax.plot([P[0]], [P[1]], [P[2]], 'ko', ms = 5)

plot_integral_curve(X_ode, POINT, X_color)
plot_integral_curve(Y_ode, POINT, Y_color)

def PSI(t,s):

	# parametrizzazione
	# parti da P, segui un flusso per t, l'altro per s, raggiungi R
	# concetto cruciale: qui t è la variabile

	time = np.linspace(0, t, ode_res)
	data=odeint(X_ode, POINT, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	Q = [curve1[-1], curve2[-1], curve3[-1]]

	time = np.linspace(0, s, ode_res)
	data=odeint(Y_ode, Q, time)
	curve1, curve2, curve3 = data[:,0],data[:,1], data[:,2]
	R = [curve1[-1], curve2[-1], curve3[-1]]

	return R


# Parametric surface

t_surface = np.linspace(-surface_range, surface_range, surface_res)
s_surface = np.linspace(-surface_range, surface_range, surface_res)

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CRUCIAL HEAVY STEP, COMPUTE A POINT FOR EVERY (t,s) PAIR

# This is faster and more elegant but does not allow to show progress
# points = [PSI(t,s) for t in t_surface for s in s_surface] # < --------------- GOOD and self contained

# SO normal cycle for progress
points = []
TOT = len(t_surface) * len(s_surface)
counter = 0
for t in t_surface:
	for s in s_surface:
		points.append(PSI(t,s))
		if counter % 100 == 0:
			print(f'counter = {counter} out of {TOT}')
		counter +=1
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------

x_surface = [el[0] for el in points]
y_surface = [el[1] for el in points]
z_surface = [el[2] for el in points]

# https://matplotlib.org/examples/mplot3d/trisurf3d_demo.html
# here the surface is plotted from points, there is no meshgrid
ax.plot_trisurf(x_surface, y_surface, z_surface, antialiased=True, alpha = 0.5)


plt.show()







