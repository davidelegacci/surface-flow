from sympy import *
from sympy.diffgeom import *

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d

import scipy.integrate as inte
from scipy.integrate import odeint
import random as rm

#########################################


# pessima notazioen tra x y u v ...

# Given two vector fields X and Y, pushes one along the flow of the other (Z = Y pushed along x)
# and plots 
# 1. the integral curve of X used to push
# 2. the vector fields Y and Z

range_param = 5
plot_basis = False
time_of_push = pi/2

u_inf = - range_param
u_sup = range_param
v_inf = -range_param
v_sup = range_param

# Risoluzione della superficie
surface_res = 50

# Numero di vettori plottati = field_res^2
field_res = 4

# dominio di curve integrali e accuratezza di ODE solver
time=np.linspace(0,10,100)

######################################### DIFFERENTIAL GEOMETRY

x_, y_, t = symbols('x_ y_ t')

M = Manifold('M', 2)
patch = Patch('P', M)
A = CoordSystem('A', patch, ['x_','y_'])

x,y  = A.coord_functions()
Dx, Dy = A.base_vectors()
dx, dy = A.base_oneforms()

def flow(X):
    
    P = A.point([x_, y_])
    
    flow_system, init_cond = intcurve_diffequ(X, t, P)
    f_0 = Function('f_0')
    f_1 = Function('f_1')
    ic = {f_0(0):x(P), f_1(0):y(P)}
    
    solution = dsolve(flow_system, ics = ic)
    
    x_flow = solution[0].rhs
    y_flow = solution[1].rhs
    
    return [x_flow, y_flow]

def push(X, func):
    
    # Push X along func
    
    # func is a list with dim(M) functions of the coordintes of M, e.g.
    # func = [x_ + y_, x_ - y_]
    # NB it's crucial that func is given ad function of the coordinates x_ and y_ WITH THE _
    # x and y without the _ are coodrinate functions
    
    u_, v_, = symbols('u_ v_')
    F = CoordSystem('F', patch, ['u_','v_'])
    
    A.connect_to(F, [x_, y_], func) # this means dummy = func(x,y)
    
    u,v = F.coord_functions()
    Du, Dv = F.base_vectors()
    
    point = F.point([u,v])
    
    X_tmp = X.subs(x,x(point)).subs(y,y(point)).subs(Dx,vectors_in_basis(Dx,F)).subs(Dy,vectors_in_basis(Dy,F))
    X_push = X_tmp.subs(u,x).subs(v,y).subs(Du,Dx).subs(Dv,Dy)
    
    return X_push

def push_flow(X,Y):
    # push Y along the flow of X
    # order = same of Lie Der
    f = flow(X)
    return f, push(Y, f)

#--------------------------------------------------------------------------- VECTORS HERE
X1 = y
X2 = -x

Y1 = y*sin(x)
Y2 = -x
#--------------------------------------------------

X_s = X1 * Dx + X2 * Dy
Y_s = Y1 * Dx + Y2 * Dy

Flow_of_X, Z_s = push_flow(X_s, Y_s)
Z_s = Z_s.subs(t, time_of_push)

Z1 = Z_s.rcall(x)
Z2 = Z_s.rcall(y)

# vectors to plot
C1 = [Y1, Z1]
C2 = [Y2, Z2]
color_vf = ['red', 'green']


##########################------------------------------------############### SURFACE
u, v = symbols('u v')
#--------------------------
x_surf_param = u #(u + v)
y_surf_param = v 
z_surf_param = 0.1 * sin(u) * cos(v)
#--------------------------
du_x_surf_param = diff(x_surf_param, u)
du_y_surf_param = diff(y_surf_param, u)
du_z_surf_param = diff(z_surf_param, u)
#--------------------------
dv_x_surf_param = diff(x_surf_param, v)
dv_y_surf_param = diff(y_surf_param, v)
dv_z_surf_param = diff(z_surf_param, v)
#--------------------------
FX = lambdify((u,v), x_surf_param)
FY = lambdify((u,v), y_surf_param)
FZ = lambdify((u,v), z_surf_param)
#--------------------------
duFX = lambdify((u,v), du_x_surf_param)
duFY = lambdify((u,v), du_y_surf_param)
duFZ = lambdify((u,v), du_z_surf_param)
#--------------------------
dvFX = lambdify((u,v), dv_x_surf_param)
dvFY = lambdify((u,v), dv_y_surf_param)
dvFZ = lambdify((u,v), dv_z_surf_param)

# Generate surface mesh

# domain
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
ax.plot_surface(x_surface, y_surface, z_surface, cmap = 'viridis', rstride = 1, cstride = 1, alpha = 0.8)

# basis tangent vector fields
def DU(u,v):
	return np.array([duFX(u,v), duFY(u,v), duFZ(u,v)])

def DV(u,v):
	return np.array([dvFX(u,v), dvFY(u,v), dvFZ(u,v)])

# vectors domain on surface

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

##########################################################################
def plot_integral_curve(P, lw = 2):

	# key: the flow of X is a curve in parameter space - 2 space

	u_func = lambdify(t,Flow_of_X[0].subs(x_, P[0]).subs(y_, P[1]))
	v_func = lambdify(t,Flow_of_X[1].subs(x_, P[0]).subs(y_, P[1]))
	u, v, = u_func(time), v_func(time) # curve in parameter space, to bring to the surface

	if (isinstance(u,int)):
		u = np.ones(len(time)) * u
	if (isinstance(v, int)):
		v = np.ones(len(time)) * v
	curve1, curve2, curve3 = FX(u,v), FY(u,v), FZ(u,v)
	ax.plot(curve1, curve2, curve3, color = 'k', linewidth = lw)
	ax.plot([FX(P[0],P[1])], [FY(P[0],P[1])], [FZ(P[0],P[1])], 'ko', ms = 5)

plot_integral_curve([1,1])

############################################################################

for i in range(len(C1)):

	used_color = color_vf[i]

	cu = lambdify((x,y), C1[i])
	cv = lambdify((x,y), C2[i])

	u_vector = np.linspace(u_inf, u_sup, field_res)
	v_vector = np.linspace(u_inf, u_sup, field_res)
	u_vector, v_vector = np.meshgrid(u_vector, v_vector)
	# evaluate
	x_vector = FX(u_vector,v_vector)
	y_vector = FY(u_vector,v_vector)
	z_vector = FZ(u_vector,v_vector)

	test_u = cu(u_vector,v_vector)
	test_v = cv(u_vector,v_vector)

	def dummy(a,b):
		return a+b

	if (isinstance(test_u,int)):
		test_u = np.ones_like(dummy(u_vector,v_vector)) * test_u

	if (isinstance(test_v,int)):
		test_v = np.ones_like(dummy(u_vector,v_vector)) * test_v


	Z1, Z2, Z3 = np.array([test_u * duFX(u_vector, v_vector), test_u * duFY(u_vector, v_vector), test_u * duFZ(u_vector, v_vector)]) + np.array([test_v * dvFX(u_vector, v_vector), test_v * dvFY(u_vector, v_vector), test_v * dvFZ(u_vector, v_vector)])
	
	ax.quiver(x_vector, y_vector, z_vector, Z1, Z2, Z3, length=0.1, color=used_color)


ax.set_zlim(-1,1)

plt.show()


##### OOOOOOK HO ESAGERATO. QUESTO FUNZIONA ma è scritto da cani e non è per nulla illuminante ed è notte fonda anzi mattina







