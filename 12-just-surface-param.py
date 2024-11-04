import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d

# parametrization range
range_param = 3

u_inf = -10
u_sup = 10
v_inf = -10
v_sup = 10

# Risoluzione della superficie
surface_res = 100

#### <--------------------------------------------------------------- HERE PARAM SURFACE

'Waves'
# X = u 
# Y = v 
# Z = 0.1*np.sin(u)*np.cos(v)

'Sphere'
# X = np.sin(u) * np.cos(v)
# Y = np.sin(u) * np.sin(v)
# Z = np.cos(u)

'Torus'

P = [  (0,0,0)  ]

def FX(u,v, p):
    point1, point2, point3 = p
    y1 = u + point1
    y3 = v + point3
    y2 = 3*(y1-y3) + point2
    return np.exp(y1) / ( 1 + np.exp(y1) + np.exp(y2) + np.exp(y3) ) 

def FY(u,v, p):
    point1, point2, point3 = p
    y1 = u + point1
    y3 = v + point3
    y2 = 3*(y1-y3) + point2
    return np.exp(y2) / ( 1 + np.exp(y1) + np.exp(y2) + np.exp(y3) ) 

def FZ(u,v, p):
    point1, point2, point3 = p
    y1 = u + point1
    y3 = v + point3
    y2 = 3*(y1-y3) + point2
    return np.exp(y3) / ( 1 + np.exp(y1) + np.exp(y2) + np.exp(y3) )  

# Generate surface mesh

# surface domain
u_surface = np.linspace(u_inf, u_sup, surface_res)
v_surface = np.linspace(v_inf, v_sup, surface_res)
u_surface, v_surface = np.meshgrid(u_surface, v_surface)

fig = plt.figure()
ax = fig.gca(projection = '3d')
for p in P:
# evaluate
    x_surface = FX(u_surface,v_surface, p)
    y_surface = FY(u_surface,v_surface, p)
    z_surface = FZ(u_surface,v_surface, p)
    # Display the mesh
    ax.plot_surface(x_surface, y_surface, z_surface, cmap = 'viridis', rstride = 1, cstride = 1, alpha = 0.5)
    ax.scatter(1/4, 1/4, 1/4, color = 'red')


ax.set_xlim(-0,1)
ax.set_ylim(-0,1)
ax.set_zlim(-0,1)
plt.show()

