""" 
Steady state advection-diffusion equation,
discontinuous formulation using full upwinding.

Implemented in python from cpp demo by Johan Hake.
"""

from dolfin import *

# Sub domain for Dirichlet boundary condition
class Boundary(SubDomain):
   def inside(self, x, on_boundary):
      return (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS) and on_boundary

   def map(self, x, y):
      y[0] = x[0] - 1.0

# Initial condition
u0 = Expression("sin(2*pi*x[0])")

# Load mesh
mesh = UnitInterval(100)

# Defining the function spaces
V = FunctionSpace(mesh, "DG", 1)

# Periodic boundary condition
bc = PeriodicBC(V, Boundary())

# Test and trial functions
phi = TestFunction(V)
uh  = TrialFunction(V)

uold = Function(V)  # previous solution
u    = Function(V)  # current solution

# Set initial condition
u = project(u0, V)
uold.assign(u)

# Mesh-related functions
n = FacetNormal(mesh)
h = CellSize(mesh)
h_avg = (h('+') + h('-'))/2

un = dot(u, n)

# Numerical flux function - Roe flux
Hn = 0.25*(un('+')**2 + un('-')**2) - 0.5*abs(un('+')+un('-'))*(u('-') - u('+'))

# Mass matrix
M = uh*phi*dx
MMat = assemble(M)

dt = 0.005 # time step
T  = 0.35  # final time
t  = 0.0   # time counter
beta = 1.0/10.0
Ceps = 0.001

# Viscosity of Jaffre
eps = Ceps * h**(2.0-beta) * abs( (u-uold)/dt + u*u.dx(0) )

# RHS
R = uold*phi*dx + dt*0.5*u**2*phi.dx(0)*dx - dt*Hn*jump(phi)*dS \
      - dt*eps*u.dx(0)*phi.dx(0)*dx

ufile = File("u.pvd")
ufile << u

while t < T + DOLFIN_EPS:
   rhs = assemble(R)
   uold.assign(u)
   solve(MMat, u.vector(), rhs)
   t += dt
   print t
   ufile << u
