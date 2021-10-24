#################################
#                               #
# FEniCS Example Comet Equation #
#                               #
#################################


# Relevant libraries

from fenics import *

import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd


# Settings

mu    = 5                 # diffusion parameter
theta = 1                 # advection angle parameter

n_segments_h = 64         # number of segments on the horizontal axis
n_segments_v = 64         # number of segments on the vertical   axis

x_0 = np.array([0.5,0.5]) # center of the bump source-force

save_to_csv = True        # 


# Create mesh and define function space

mesh = UnitSquareMesh(n_segments_h, n_segments_v)
V = FunctionSpace(mesh, 'P', 1)


# Define boundary condition

u_D = Expression('0', degree=1)     # Homogeneous (null) boundary conditions

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)


# Define variational problem

u = TrialFunction(V)
v = TestFunction(V)

f = Expression('10*exp(-100*pow( pow(x[0]-x_00, 2) + pow(x[1]-x_01, 2), 0.5))', 
               degree=1, x_00 = x_0[0], x_01 = x_0[1])


b = Constant((np.cos(theta), np.sin(theta)))
a = mu * dot(grad(u), grad(v)) * dx + 10 * dot(b,grad(u)) * v * dx
L = f*v*dx


# Compute solution

u = Function(V)
solve(a == L, u, bc)


# Plot solution and mesh

fig1 = plt.figure(figsize=(8, 8), dpi=300)

plot(u)
plot(mesh)


# Get mesh vertices and corresponding solution

vertex_values    = u.compute_vertex_values()
vertex_positions = mesh.coordinates()


# Plot the solution on the vertices

fig2 = plt.figure(figsize=(8, 8), dpi=300)

plt.scatter(x = vertex_positions[:,0], 
            y = vertex_positions[:,1],
            c = vertex_values)
plt.colorbar()
plt.show()


# To debug, save solution to csv

if(save_to_csv):
    
    df = pd.DataFrame({"x" : vertex_positions[:,0], 
                       "y" : vertex_positions[:,1], 
                       "u" : vertex_values})
    df.to_csv("comet_solution.csv", index=False)