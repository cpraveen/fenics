"""
Steady flow over cylinder in external flow

No stress BC on top, bottom and outflow boundaries
Picard iteration on convective term
"""
from dolfin import *
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as la

parameters.linear_algebra_backend = "uBLAS"

class inlet_velocity(Expression):
   def __init__(self, t=0.0):
      self.t = t
   def eval(self, value, x):
      value[0] = 1.0
      value[1] = 0.0
      return value
   def value_shape(self):
      return (2,)

class initial_condition(Expression):
   def eval(self, value, x):
      value[0] = 1.0
      value[1] = 0.0
      value[2] = 0.0
      return value
   def value_shape(self):
      return (3,)

mesh = Mesh("cylinder.xml")
boundaries = MeshFunction("size_t", mesh, "cylinder_facet_region.xml")

# Function space
udeg = 2
pdeg = udeg - 1

V = VectorFunctionSpace(mesh, 'CG', udeg)
W = FunctionSpace(mesh, 'CG', pdeg)
X = V * W

print "Velocity dofs = ", V.dim()
print "Pressure dofs = ", W.dim()
print "Total    dofs = ", X.dim()

# Solution variables
ups = Function(X)  # steady solution
File("steady.xml") >> ups.vector()
us = as_vector((ups[0],ups[1]))

# Trial functions
up  = TrialFunction(X)
u   = as_vector((up[0],up[1]))
p   = up[2]

# Test functions
vq  = TestFunction(X)
v   = as_vector((vq[0],vq[1]))
q   = vq[2]

# Boundary condition
vin  = inlet_velocity()
bcin = DirichletBC(X.sub(0), vin,   boundaries, 1)
bccyl= DirichletBC(X.sub(0), (0,0), boundaries, 2)
bcs  = [bcin, bccyl]

# Same as above, but all zero dirichlet bc.
# Used in residual computation
bcin0 = DirichletBC(X.sub(0), (0,0), boundaries, 1)
bccyl0= DirichletBC(X.sub(0), (0,0), boundaries, 2)
bcsid0= DirichletBC(X.sub(0).sub(1), 0, boundaries, 4)
bcs0  = [bcin0, bccyl0, bcsid0]

Ur = 1.0                # Reference velocity
D  = 0.1                # dia of cylinder
Re = 50.0               # Reynolds number
nu = Constant(Ur*D/Re)  # viscosity coefficient


# Form for linearized NS equation
F  = - inner(grad(u)*us, v)*dx      \
     - inner(grad(us)*u, v)*dx      \
     + p*div(v)*dx                   \
     - nu*inner(grad(u), grad(v))*dx \
     + q*div(u)*dx

Aa = assemble(F)

# Convert to sparse format
rows, cols, values = Aa.data()
Aa = sps.csc_matrix((values, cols, rows))
print "Size of Aa =",Aa.shape

m  = inner(u,v)*dx
Ma = assemble(m)

# Convert to sparse format
rows, cols, values = Ma.data()
Ma = sps.csc_matrix((values, cols, rows))
print "Size of Ma =",Ma.shape

bcinds = []
for bc in bcs0:
   bcinds.extend(bc.get_boundary_values().keys())

print len(bcinds)
N = X.dim()
freeinds = np.setdiff1d(range(N),bcinds,assume_unique=True).astype(np.int32)

A = Aa[freeinds,:][:,freeinds]
print "Size of A =",A.shape

M = Ma[freeinds,:][:,freeinds]
print "Size of M =",M.shape

# Compute eigenvalues/vectors of (A,M)
print "Computing eigenvalues/vectors ..."
sigma = 0.0
vals, vecs = la.eigs(A, k=100, M=M, sigma=sigma, which='LM', ncv=400, tol=1.0e-8)
ii = np.argsort(vals)[::-1]
fv = File("eig.pvd")
up = Function(X)
fe = open("eig.dat","w")
for i in ii:
    vr, vi = np.real(vals[i]), np.imag(vals[i])
    print vr, vi
    fe.write(str(vr)+"  "+str(vi)+"\n")

    up.vector()[freeinds] = np.array(np.real(vecs[:,i]))
    u,p = up.split()
    fv << u

    up.vector()[freeinds] = np.array(np.imag(vecs[:,i]))
    u,p = up.split()
    fv << u
