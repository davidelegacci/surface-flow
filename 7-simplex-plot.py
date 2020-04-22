import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d

'''

Plots a surface being the graph of a function on a domain wich is not just [a, b] x [c, d], instead with y depending on x (in particular a simplex)
For example, triangular 2D domain to plot simplex

'''

# x parameter interval domain
x_inf = 0
x_sup = 1

# number of plotted points = this number^2
surface_res = 10

# the y parameter range depends on the x one
def y_bound(x):
	return 1-x # 1-x for 2-simplex

# graph of z = f(x,y) making surface
def FZ(x,y):
	return 1-(x+y) # 1-(x+y) for 2-simplex


x_surface = np.linspace(x_inf, x_sup, surface_res)

X_surface = [x*np.ones(len(x_surface)) for x in x_surface]
Y_surface = [ np.linspace(0, y_bound(x), surface_res) for x in x_surface ]
Z_surface = [ [FZ (x,y) for y in Y ] for x,Y in zip(x_surface, Y_surface) ]

# Now X_surface, Y_surface, Z_surface are surface_res lists list of surface_res elements each
# Every X list contains copies of one element
# They must be concatenated
X_plot = []
for el in X_surface:
	X_plot += list(el)

Y_plot = []
for el in Y_surface:
	Y_plot += list(el)

Z_plot = []
for el in Z_surface:
	Z_plot += list(el)

# Now we have three normal lists of coordinates.

fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.plot_trisurf(X_plot, Y_plot, Z_plot, antialiased=True, alpha = 0.5)


plt.show()




