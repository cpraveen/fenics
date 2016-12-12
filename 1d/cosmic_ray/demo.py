import numpy as np
import matplotlib.pyplot as plt
from dolfin import *

#PETScOptions().set("snes_monitor")
#PETScOptions().set("snes_view")
#PETScOptions().set("snes_converged_reason")

nc = 101
T  = 0.02  # final time
h  = 1.0/(nc-1)
ep = Constant(h)
dt = h
idt= Constant(1.0/dt)

print "h  = ", h
print "dt = ", dt

# Sub domain for Dirichlet boundary condition
class PeriodicBoundary(SubDomain):
   def inside(self, x, on_boundary):
      return (near(x[0],0.0) and on_boundary)

   def map(self, x, y):
      y[0] = x[0] - 1.0

# Initial condition
u0 = Expression("1.0+0.5*sin(2*pi*x[0])",degree=1)

# Load mesh
mesh = UnitIntervalMesh(nc)

# Defining the function spaces
V = FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicBoundary())

# Test and trial functions
v    = TestFunction(V)
u    = Function(V)
du   = TrialFunction(V)
uold = Function(V)  # previous solution

ux = u.dx(0)
vx = v.dx(0)
K0 = u/sqrt(u.dx(0)**2 + ep**2)
K1 = uold/sqrt(uold.dx(0)**2 + ep**2)
R0 = idt*(u-uold)*v*dx + K0*ux*vx*dx
R1 = idt*(u-uold)*v*dx + K1*ux*vx*dx 
J0  = derivative(R0, u, du)  # Gateaux derivative in dir. of du
J1  = derivative(R1, u, du)  # Gateaux derivative in dir. of du

R, J = R0, J0

problem = NonlinearVariationalProblem(R, u, J=J)
solver  = NonlinearVariationalSolver(problem)
prm = solver.parameters
prm['nonlinear_solver'] = 'snes'
prm['snes_solver']['line_search'] = 'bt'
prm['snes_solver']['sign'] = 'nonnegative'
prm['snes_solver']['linear_solver']= 'lu'
prm['snes_solver']['maximum_iterations'] = 200

# Set initial condition
u.interpolate(u0)
uinit = u.vector().array()
x = V.tabulate_dof_coordinates()
i = np.argsort(x)
plt.plot(x[i],uinit[i])
plt.xlabel('x')
plt.ylabel('u')

it, t = 0, 0.0
while t < T:
    uold.assign(u)
    solver.solve()
    it += 1; t += dt

print "it, t =", it, t

# Plot final solution
unew = u.vector().array()
plt.plot(x[i],unew[i]);
plt.savefig("sol.pdf")
