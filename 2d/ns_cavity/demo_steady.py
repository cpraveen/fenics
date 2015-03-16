"""
This demo program solves the steady incompressible Navier-Stokes equations
for lid-driven cavity problem using Taylor-Hood elements.
Author: Praveen. C
www   : http://math.tifrbng.res.in/~praveen
"""

from dolfin import *

n = 100
Re= 100

# Load mesh from file
mesh = UnitSquare(n,n,"crossed")

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define test functions
(v,q) = TestFunctions(W)

# Define trial functions
w     = Function(W)
(u,p) = (as_vector((w[0], w[1])), w[2])

# Set parameter values
nu = 1.0/Re

# Define boundary conditions
noslip  = DirichletBC(W.sub(0), (0, 0), "x[0] < DOLFIN_EPS || x[0] > 1.0 - DOLFIN_EPS || x[1] < DOLFIN_EPS")
lid  = DirichletBC(W.sub(0), (1,0), "x[1] > 1.0 - DOLFIN_EPS")
pref = DirichletBC(W.sub(1), 0, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")

bc = [noslip, lid, pref]

# Tentative velocity step
F =   inner(grad(u)*u, v)*dx \
    + nu*inner(grad(u), grad(v))*dx \
    - div(v)*p*dx \
    - q*div(u)*dx

dw = TrialFunction(W)
dF = derivative(F, w, dw)

nsproblem = NonlinearVariationalProblem(F, w, bc, dF)
solver = NonlinearVariationalSolver(nsproblem)
solver.solve()

(u,p) = w.split()
File("velocity.pvd", "compressed") << u
File("pressure.pvd", "compressed") << p

# Compute vorticity by L2 projection
r = TrialFunction(Q)
s = TestFunction(Q)
a = r*s*dx
L = (u[0].dx(1) - u[1].dx(0))*s*dx
vort = Function(Q)
solve(a == L, vort)
File("vorticity.pvd", "compressed") << vort

# Compute stream function
# Laplace(psi) = -vort
a = inner(grad(r), grad(s))*dx
L = vort*s*dx
psi = Function(Q)
wall = DirichletBC(Q, 0, "on_boundary")
solve(a == L, psi, wall)
File("stream.pvd", "compressed") << psi
