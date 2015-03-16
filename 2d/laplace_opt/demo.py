from dolfin import *
import numpy
import sys

if len(sys.argv) < 2:
   sys.exit("Specify gradient method: linear or adjoint")

mode = sys.argv[1]

set_log_level(100)

# Boundary
def Boundary(x, on_boundary):
   return on_boundary

# Characteristic function
class CHAR_A(Expression):
   def eval(self, value, x):
      if x[0]**2 + x[1]**2 < 1.0:
         value[0] = 1.0
      else:
         value[0] = 0.0

class CHAR_B(Expression):
   def eval(self, value, x):
      if (x[0]+3)**2 + x[1]**2 < 1.0:
         value[0] = 1.0
      else:
         value[0] = 0.0

class CHAR_C(Expression):
   def eval(self, value, x):
      if x[0]**2 + (x[1]+3)**2 < 1.0:
         value[0] = 1.0
      else:
         value[0] = 0.0

# Diffusion coefficient
class KAPPA(Expression):
   def __init__(self, A, B, C):
      self._A = A
      self._B = B
      self._C = C
   def eval(self, value, x):
      value[0] = 1.0
      if x[0]**2 + x[1]**2 < 1.0:
         value[0] += self._A
      elif (x[0]+3)**2 + x[1]**2 < 1.0:
         value[0] += self._B
      elif x[0]**2 + (x[1]+3)**2 < 1.0:
         value[0] += self._C

# Load the mesh
mesh = Mesh('circle.xml')

# Create function space and trial/test functions
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

# Define problem
kappa = KAPPA(A=2, B=3, C=4)
a = kappa*inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx
g = Expression('pow(x[0],3) - pow(x[1],3)')
bc= DirichletBC(V, g, Boundary)

# Compute target solution
ud = Function(V)
solve( a == L, ud, bc )
plot(kappa, mesh=mesh, title='Optimized kappa')
plot(ud, title='Desired solution')

File('sol_desired.pvd')   << ud
File('kappa_desired.pvd') << project(kappa, V)

# Variable to store new solution
un = Function(V)

# Objective function
def J(K):
   kappa._A = K[0]
   kappa._B = K[1]
   kappa._C = K[2]
   solve(a==L, un, bc)
   return assemble(0.5*(un-ud)**2*dx)

# BC and linear functionals for linearized problems
bc_du = DirichletBC(V, 0, Boundary)
LA = -CHAR_A()*inner(grad(un), grad(v))*dx
LB = -CHAR_B()*inner(grad(un), grad(v))*dx
LC = -CHAR_C()*inner(grad(un), grad(v))*dx

# Linear functional for adjoint
L_adj = - (un - ud)*v*dx

# Gradient of objective function
def dJ_lin(K):
   kappa._A = K[0]
   kappa._B = K[1]
   kappa._C = K[2]
   dua = Function(V)
   dub = Function(V)
   duc = Function(V)
   solve(a==LA, dua, bc_du)
   solve(a==LB, dub, bc_du)
   solve(a==LC, duc, bc_du)
   G = numpy.array([0.0, 0.0, 0.0])
   G[0] = assemble((un-ud)*dua*dx)
   G[1] = assemble((un-ud)*dub*dx)
   G[2] = assemble((un-ud)*duc*dx)
   return G

# Gradient of objective function using adjoint equation
def dJ_adj(K):
   kappa._A = K[0]
   kappa._B = K[1]
   kappa._C = K[2]
   phi = Function(V)
   solve(a==L_adj, phi, DirichletBC(V, 0, Boundary))
   G = numpy.array([0.0, 0.0, 0.0])
   G[0] = assemble( CHAR_A()*inner(grad(un), grad(phi))*dx )
   G[1] = assemble( CHAR_B()*inner(grad(un), grad(phi))*dx )
   G[2] = assemble( CHAR_C()*inner(grad(un), grad(phi))*dx )
   return G

# Compute gradient
def dJ(K):
   if mode == 'linear':
      return dJ_lin(K)
   elif mode == 'adjoint':
      return dJ_adj(K)
   else:
      print "Unknown option = ", mode

# Initial guess
K=numpy.array([0, 0, 0])

obj = J(K)
G   = dJ(K)
print "Control variable   =", K
print "Objective function =", obj
print "Gradient           =", G

# Initial step length
step = 1.0e-3

# Optimization iterations
# We use simple steepest descent and step size is slowly
# increased as iterations progress
iter = 1
while numpy.linalg.norm(G) > 1.0e-5:
   K   = K - step*G
   obj = J(K)
   G   = dJ(K)
   print "Iteration = ", iter
   print "Control variable   =", K
   print "Objective function =", obj
   print "Gradient           =", G
   print "Step size          =", step
   step *= 1.1
   iter += 1


# Plot and save optimized solution
plot(kappa, mesh=mesh, title='Optimized kappa')
plot(un, title='Optimized solution')
File('sol_opt.pvd')   << un
File('kappa_opt.pvd') << project(kappa, V)
interactive()
