"""
Compute eigenvalues of Laplacian on unit square. using Lagrange multiplier for
boundary condition.

Exact eigenvalues are

eig/pi^2 = -2, -5, -5, -8, -10, -10, -13, -13, -17, -17, -18, ...

After running this program, run eigval.m in matlab.
"""

import scipy.io as sio
from dolfin import *
from param import *

parameters.linear_algebra_backend = "Eigen"

mesh = Mesh('mesh.xml')
V = FunctionSpace(mesh, 'CG', degree)

u = TrialFunction(V)
v = TestFunction(V)

M = assemble(u*v*dx)
M = as_backend_type(M).sparray()

A = assemble(-inner(grad(u),grad(v))*dx + Constant(shift)*u*v*dx)
A = as_backend_type(A).sparray()

# Boundary mass matrix
N = assemble(u*v*ds)

bd = boundary()
bc = DirichletBC(V, 0.0, bd)
binds = bc.get_boundary_values().keys()
N = as_backend_type(N).sparray()
Mb= N[binds,:][:,binds] # boundary mass matrix
N = N[:,:][:,binds]

print 'Saving matrices into linear.mat'
sio.savemat('linear.mat', mdict={'M':M, 'A':A, 'N':N, 'Mb':Mb}, oned_as='column')
