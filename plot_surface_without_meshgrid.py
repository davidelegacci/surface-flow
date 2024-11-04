import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np

X = np.arange(1, 2, 0.25)
Y = np.arange(1, 2, 0.25)
X, Y = np.meshgrid(X, Y)

def f(x,y):
	return x+y


Z = f(X,Y)

myZ = []
for i in range(len(X)):
	myZ.append([])
	for j in range(len(Y)):
		myZ[-1].append( f(X[i][j], Y[i][j]) )


myZ = np.array([ np.array(z) for z in myZ])

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Plot the surface.
surf = ax.plot_surface(X, Y, myZ, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)

# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()