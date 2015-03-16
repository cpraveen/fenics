"""
This demo program solves Poisson's equation

    - (u_xx + u_yy) = f(x)

on the unit square with source f given by

    f(x) = sin(pi*x[0])*cos(pi*x[1])

and boundary conditions given by

    u = 0
"""

from dolfin import *

# Read mesh
mesh = UnitSquareMesh(20,20)

# FE function space
V    = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
def Boundary(x, on_boundary):
     return on_boundary

# Define variational problem
f = Expression("sin(pi*x[0])*cos(pi*x[1])")

u  = TrialFunction(V)
v  = TestFunction(V)

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, Boundary)

a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Solve PDE
w = Function(V)
problem = LinearVariationalProblem(a, L, w, bc)
solver  = LinearVariationalSolver(problem)
solver.solve()

# Plot solution
plot(w, interactive=True)

# Save solution to file
File("u.pvd") << w
