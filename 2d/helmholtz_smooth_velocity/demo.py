"""
Solves
   -eps*Laplace(wx) + wx = wx0
   -eps*Laplace(wy) + wy = wy0
with Dirichlet bc given by rhs (wx0,wy0). The two equations are decoupled but 
we solve them together.
"""
from dolfin import *

# domain
xmin, xmax = -5.0, +5.0
ymin, ymax = -5.0, +5.0

# mesh size
nx, ny = 50, 50

# Smoothing coefficient
eps = 0.01

# Velocity of isentropic vortex
class Vortex(Expression):
    def eval(self, value, x):
        value[0] =  x[1]
        value[1] = -x[0]
        return value
    def value_shape(self):
        return(2,)

# mesh
p0 = Point(xmin, ymin)
p1 = Point(xmax, ymax)
mesh = RectangleMesh(p0, p1, nx, ny)

# FE function space
V    = VectorFunctionSpace(mesh, "CG", 1)

wx, wy  = TrialFunctions(V)
tx, ty  = TestFunctions(V)

# Define boundary condition
w0 = Vortex(degree=1)
bc = DirichletBC(V, w0, 'on_boundary')

a = eps*dot(grad(wx), grad(tx))*dx + eps*dot(grad(wy), grad(ty))*dx \
        + (wx*tx + wy*ty)*dx
L = (w0[0]*tx + w0[1]*ty)*dx

# Solve PDE
w = Function(V)
problem = LinearVariationalProblem(a, L, w, bc)
solver  = LinearVariationalSolver(problem)
solver.solve()

# Save solution to file
File("w.pvd") << w
