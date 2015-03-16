"""
This demo program solves the steady incompressible Navier-Stokes equations
for cylinder in channel problem using Taylor-Hood elements.
   Author: Praveen. C
   www   : http://math.tifrbng.res.in/~praveen
"""

import sys, math

if len(sys.argv) < 2:
   sys.exit("Must specify Reynolds number; restart is optional")

# Whether to use previous solution as initial condition
restart = "no"
if len(sys.argv) == 3:
   restart = sys.argv[2]

# Set parameter values
Re   = float(sys.argv[1])
D    = 0.1
Uinf = 1.0
nu   = D * Uinf / Re

from dolfin import *

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

# Define test functions
(v,q) = TestFunctions(W)

# Define trial functions
w     = Function(W)
(u,p) = (as_vector((w[0], w[1])), w[2])

# Define boundary conditions
uinlet = Expression(("1.5*(1.0 - (x[1]/0.2)*(x[1]/0.2))", "0"))
cyl    = DirichletBC(W.sub(0), (0, 0), sub_domains, 0)
inlet  = DirichletBC(W.sub(0), uinlet, sub_domains, 1)
noslip = DirichletBC(W.sub(0), (0, 0), sub_domains, 3)
bc     = [noslip, inlet, cyl]

# Weak form
F =   inner(grad(u)*u, v)*dx        \
    + nu*inner(grad(u), grad(v))*dx \
    - p*div(v)*dx                   \
    - q*div(u)*dx

# Derivative of weak form
dw = TrialFunction(W)
dF = derivative(F, w, dw)

problem = NonlinearVariationalProblem(F, w, bc, dF)
solver  = NonlinearVariationalSolver(problem)
# Set linear solver parameters
itsolver = solver.parameters["newton_solver"]
itsolver["absolute_tolerance"] = 1.0e-10
itsolver["relative_tolerance"] = 1.0e-10

# To see various solver options, uncomment following line
#info(solver.parameters, True); quit()

# If you want to initialize solution from previous computation
if restart == "yes":
   print "Setting initial condition from file ..."
   File("steady.xml") >> w.vector()

# Solve the problem
solver.solve()

# Save steady solution
File("steady.xml") << w.vector()

# Save vtk for visualization
(u,p) = w.split()
File("velocity.pvd") << u
File("pressure.pvd") << p

# Compute and save vorticity in vtk format
r = TrialFunction(Q)
s = TestFunction(Q)
a = r*s*dx
L = (u[0].dx(1) - u[1].dx(0))*s*dx
vort = Function(Q)
solve(a == L, vort)
File("vorticity.pvd") << vort

# Stress tensor
T = nu*(grad(u) + grad(u).T) - p*Identity(2)
# Face normals
n = FacetNormal(mesh)

# Compute force on cylinder
drag = -T[0,j]*n[j]*dss(0)
lift = -T[1,j]*n[j]*dss(0)
print "Drag =", assemble(drag)
print "Lift =", assemble(lift)
