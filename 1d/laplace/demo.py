"""
This demo program solves Poisson's equation

    - u_xx = f(x)

on the unit interval with source f given by

    f(x) = sin(pi*x[0])

and boundary conditions given by

    u(0) = u(1) = 0
"""

from dolfin import *

# Read mesh
mesh = UnitInterval(20)

# FE function space
V    = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
def Boundary(x, on_boundary):
     return on_boundary

# Define variational problem
f = Expression("sin(pi*x[0])")

u  = TrialFunction(V)
v  = TestFunction(V)

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, Boundary)

a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Solve PDE
problem = VariationalProblem(a, L, bc)
u = problem.solve()

# Save solution to file
File("u.pvd") << u
