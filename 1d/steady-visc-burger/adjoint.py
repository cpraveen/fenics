"""
This demo program solves Poisson's equation

    (u^2/2)_x - u_xx = f(x)

on the unit interval with source f given by

    f(x) = 9*pi^2*sin(3*pi*x[0])

and boundary conditions given by

    u(0) = u(1) = 0
"""

from dolfin import *

# Read mesh
mesh = Mesh("mesh.xml")

# FE function space
V    = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
def DirichletBoundary(x):
     return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

u  = Function(V)
v  = TestFunction(V)

# Read primal solution from file
File("primal.xml") >> u.vector()

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, DirichletBoundary)

# Adjoint problem
phi = TrialFunction(V)
a = (u*phi.dx(0)*v - inner(grad(phi), grad(v)))*dx
L = u*v*dx

adjoint = VariationalProblem(a, L, bc)
phi = adjoint.solve()

# Save solution to file
File("adjoint.pvd") << phi
File("adjoint.xml") << phi.vector()
