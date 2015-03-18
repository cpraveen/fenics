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

A = PETScMatrix()
assemble(F, tensor=A)

m  = inner(u,v)*dx
M  = PETScMatrix()
assemble(m, tensor=M)

for bc in bcs:
    (bc.apply(A) and bc.zero(M))

#create eigensolver
eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters['spectrum'] = 'smallest magnitude'
#eigensolver.parameters['solver'] = 'arnoldi'
#eigensolver.parameters['problem_type'] = 'gen_non_hermitian'
eigensolver.parameters["spectral_transform"] = "shift-and-invert"
eigensolver.parameters["spectral_shift"] = 0.0
#eigensolver.parameters['verbose'] = True

print 'solving: start'
num = 100
eigensolver.solve(num)
print 'solving: end'

print "Number of converged eigenvalues = %d" % eigensolver.get_number_converged()

fe = File("eig.pvd")
ue = Function(W)
for i in range(num):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    print "Eigenvalue: %5d %20.10e %20.10e" % (i, r, c)
    ue.vector()[:] = rx
    u,p = ue.split()
    fe << u
