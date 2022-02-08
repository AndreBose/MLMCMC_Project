####################
#                  #
# FEniCS Example 1 #
#                  #
####################

##############################
#   PDE definition       

#  - delta(u) = f  in Omega
#           u = g  in d_Omega
# domain: Omega = [0,1]^2

##############################


# Settings

n_segments_h = 32
n_segments_v = 32


# Relevant libraries

from fenics import *

import numpy             as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 8), dpi=1200)


# Create mesh and define function space

mesh = UnitSquareMesh(n_segments_h, n_segments_v)   # create a 2D square mesh in [0,1]^2
V = FunctionSpace(mesh, 'P', 1)                     # create a Galerkin functional space with elements of degree 1


# Define Dirichlet boundary condition

u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)   #it is an user defined mathematical expression to describe the function over the BC

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)


# Define variational problem: find u in V st a(u,v) = F(v) for all v in V

u = TrialFunction(V)    # set it as the unknown function
v = TestFunction(V)     # set it as test function
f = Constant(-6.0)      # set it as the force function

a = dot(grad(u), grad(v))*dx    # it represents a(u,v) = integral( grad(u)' grad(v) dx )
F = f*v*dx                      # it represents F(v) = integral( f(x) * v(x) * dx )
# dx is a special character to say it integrates

# Compute solution

u = Function(V)
solve(a == F, u, bc)            # solve variational problem a(u,v) = f(v) with u as unknown and given BC


# Plot solution and mesh

# plot(u)
# plot(mesh)


# Compute error in L2 norm

error_L2 = errornorm(u_D, u, 'L2')


# Compute maximum error at vertices

vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u   = u.compute_vertex_values(mesh)
error_max         = np.max(np.abs(vertex_values_u_D - vertex_values_u))


# Print errors

print('error_L2 =', error_L2)
print('error_max =', error_max)


# Get mesh vertices and corresponding solution

vertex_values    = u.compute_vertex_values()
vertex_positions = mesh.coordinates()


# Plot the solution on the vertices

plt.scatter(x = vertex_positions[:,0], 
            y = vertex_positions[:,1],
            c = vertex_values)
plt.colorbar()
plt.show()

