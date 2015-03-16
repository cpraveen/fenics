"""
Steady flow over cylinder in external flow

No stress BC on top, bottom and outflow boundaries
Picard iteration on convective term
"""
from dolfin import *

class inlet_velocity(Expression):
   def __init__(self, t=0.0):
      self.t = t
   def eval(self, value, x):
      value[0] = 1.0
      value[1] = 0.0
      return value
   def value_shape(self):
      return (2,)

class initial_condition(Expression):
   def eval(self, value, x):
      value[0] = 1.0
      value[1] = 0.0
      value[2] = 0.0
      return value
   def value_shape(self):
      return (3,)

mesh = Mesh("cylinder.xml")
boundaries = MeshFunction("size_t", mesh, "cylinder_facet_region.xml")

# Function space
udeg = 2
pdeg = udeg - 1

V = VectorFunctionSpace(mesh, 'CG', udeg)
W = FunctionSpace(mesh, 'CG', pdeg)
X = V * W

print "Velocity dofs = ", V.dim()
print "Pressure dofs = ", W.dim()
print "Total    dofs = ", X.dim()

# Solution variables
ups = Function(X)  # steady solution
File("steady.xml") >> ups.vector()
us = as_vector((ups[0],ups[1]))

# Trial functions
up  = TrialFunction(X)
u   = as_vector((up[0],up[1]))
p   = up[2]

# Test functions
vq  = TestFunction(X)
v   = as_vector((vq[0],vq[1]))
q   = vq[2]

# Boundary condition
vin  = inlet_velocity()
bcin = DirichletBC(X.sub(0), vin,   boundaries, 1)
bccyl= DirichletBC(X.sub(0), (0,0), boundaries, 2)
bcs  = [bcin, bccyl]

# Same as above, but all zero dirichlet bc.
# Used in residual computation
bcin0 = DirichletBC(X.sub(0), (0,0), boundaries, 1)
bccyl0= DirichletBC(X.sub(0), (0,0), boundaries, 2)
bcsid0= DirichletBC(X.sub(0).sub(1), 0, boundaries, 4)
bcs0  = [bcin0, bccyl0, bcsid0]

Ur = 1.0                # Reference velocity
D  = 0.1                # dia of cylinder
Re = 100.0              # Reynolds number
nu = Constant(Ur*D/Re)  # viscosity coefficient


# Form for linearized NS equation
F  = - inner(grad(u)*us, v)*dx      \
     - inner(grad(us)*u, v)*dx      \
     + p*div(v)*dx                   \
     - nu*inner(grad(u), grad(v))*dx \
     + q*div(u)*dx

A = PETScMatrix()
assemble(F, tensor=A)

m  = inner(u,v)*dx
M  = PETScMatrix()
assemble(m, tensor=M)

for bc in bcs0:
    (bc.apply(A) and bc.zero(M))

#create eigensolver
eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters['spectrum'] = 'largest real'
#eigensolver.parameters['spectrum'] = 'smallest magnitude'
#eigensolver.parameters['solver'] = 'arnoldi'
#eigensolver.parameters['problem_type'] = 'gen_non_hermitian'
eigensolver.parameters["spectral_transform"] = "shift-and-invert"
eigensolver.parameters["spectral_shift"] = 0.0
#eigensolver.parameters['verbose'] = True

print 'solving: start'
num = 20
eigensolver.solve(num)
print 'solving: end'

print "Number of converged eigenvalues = %d" % eigensolver.get_number_converged()

fe = File("eig.pvd")
ue = Function(X)
for i in range(num):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    print "Eigenvalue: %5d %20.10e %20.10e" % (i, r, c)
    ue.vector()[:] = rx
    u,p = ue.split()
    fe << u
    ue.vector()[:] = cx
    u,p = ue.split()
    fe << u
