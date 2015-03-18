"""
This demo program solves the steady incompressible Navier-Stokes equations
for cylinder in channel problem using Taylor-Hood elements.
   Author: Praveen. C
   www   : http://math.tifrbng.res.in/~praveen
"""
import sys, math

if len(sys.argv) < 2:
   sys.exit("Must specify Reynolds number; restart is optional")

from dolfin import *
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as la

parameters.linear_algebra_backend = "uBLAS"

# Set parameter values
Re   = float(sys.argv[1])
D    = 0.1
Uinf = 1.0
nu   = D * Uinf / Re

# Load mesh from file
mesh = Mesh("cylinder_in_channel.xml")
sub_domains = MeshFunction("size_t", mesh, "subdomains.xml")
dss = Measure("ds")[sub_domains]

# Define function spaces (P2-P1)
udeg = 2
pdeg = udeg - 1
V = VectorFunctionSpace(mesh, "CG", udeg)
Q = FunctionSpace(mesh, "CG", pdeg)
W = V * Q

ups = Function(W)
print "Reading stationary solution from file ..."
File("steady.xml") >> ups.vector()
us = as_vector((ups[0],ups[1]))

# Define test functions
(v,q) = TestFunctions(W)

# Define trial functions
(u,p) = TrialFunctions(W)

# Define boundary conditions
cyl    = DirichletBC(W.sub(0), (0, 0), sub_domains, 0)
inlet  = DirichletBC(W.sub(0), (0, 0), sub_domains, 1)
noslip = DirichletBC(W.sub(0), (0, 0), sub_domains, 3)
bcs    = [noslip, inlet, cyl]

# Weak form
F = - inner(grad(us)*u, v)*dx        \
    - inner(grad(u)*us, v)*dx        \
    - nu*inner(grad(u), grad(v))*dx \
    + p*div(v)*dx                   \
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
for bc in bcs:
   bcinds.extend(bc.get_boundary_values().keys())

print len(bcinds)
N = W.dim()
freeinds = np.setdiff1d(range(N),bcinds,assume_unique=True).astype(np.int32)

A = Aa[freeinds,:][:,freeinds]
print "Size of A =",A.shape

M = Ma[freeinds,:][:,freeinds]
print "Size of M =",M.shape

# Compute eigenvalues/vectors of (A,M)
print "Computing eigenvalues/vectors ..."
sigma = 10.0
vals, vecs = la.eigs(A, k=200, M=M, sigma=sigma, which='LM', ncv=400, tol=1.0e-8)
ii = np.argsort(vals)[::-1]
fv = File("eig.pvd")
up = Function(W)
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
