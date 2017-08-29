"""
Compute eigenvalues of Laplacian on unit square. using Lagrange multiplier for
boundary condition.

Exact eigenvalues are

eig/pi^2 = -2, -5, -5, -8, -10, -10, -13, -13, -17, -17, -18, ...

After running this program, run eigval.m in matlab.
"""

from scipy.sparse.linalg import cg
import scipy.io as sio
from dolfin import *

parameters.linear_algebra_backend = "Eigen"

class boundary(SubDomain):
   def inside(self,x,on_boundary):
      return on_boundary


degree = 1
n      = 50
mesh = UnitSquareMesh(n,n)
V = FunctionSpace(mesh, 'CG', degree)

u = TrialFunction(V)
v = TestFunction(V)

M = assemble(u*v*dx)
M = as_backend_type(M).sparray()

A = assemble(-inner(grad(u),grad(v))*dx)
A = as_backend_type(A).sparray()

# Boundary mass matrix
N = assemble(u*v*ds)

bd = boundary()
bc = DirichletBC(V, 0.0, bd)
binds = bc.get_boundary_values().keys()
N = as_backend_type(N).sparray()
N = N[:,:][:,binds]

print "Saving matrices into linear.mat"
sio.savemat('linear.mat', mdict={'M':M, 'A':A, 'N':N}, oned_as='column')
