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
import time


# Settings

x_0 = np.array([0.5,0.5]) # center of the bump source-force


n_segments_hv = [4, 8, 16, 32, 64, 128, 256, 512]  # number of segments on horizontal and vertical axis

n_attempts = 100                                   # number of iterations per refinement
n_burn     =  30                                   # number of times that are discarded for the mean


# Perform iterations

times_list = []


for n_segments in n_segments_hv:
    
    current_times = []
    
    for _ in range(n_attempts):
    
        mu    = np.random.uniform(low=0.1, high=5.0    )
        theta = np.random.uniform(low=0. , high=2*np.pi)
    
    
        # Create mesh and define function space
        
        initial_time = time.time()
        
        mesh = UnitSquareMesh(n_segments, n_segments)
        V    = FunctionSpace(mesh, 'P', 2)
    
    
        # Define boundary condition
        
        u_D = Constant(0)         # Homogeneous (null) boundary conditions
        
        def boundary(x, on_boundary):
            return on_boundary
        
        bc = DirichletBC(V, u_D, boundary)
        
        
        # Define variational problem
        
        u = TrialFunction(V)
        v = TestFunction (V)
        
        f = Expression('10*exp(-50*pow( pow(x[0]-x_00, 2) + pow(x[1]-x_01, 2), 0.5))', 
                       degree=2, x_00 = x_0[0], x_01 = x_0[1])
        
        
        b = Constant((np.cos(theta), np.sin(theta)))
        a = mu * dot(grad(u), grad(v)) * dx + 10 * dot(b,grad(u)) * v * dx
        L = f*v*dx
    
    
        # Compute solution
        
        u = Function(V)
        solve(a == L, u, bc)
        
        current_times.append( time.time() - initial_time )

    
    times_list.append(current_times)


for index in range(len(n_segments_hv)):
    
    print(n_segments_hv[index],'refinement takes on average ', np.mean(np.array(times_list[index])[(n_burn+1):]),'seconds \n')


