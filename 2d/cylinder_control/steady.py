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
up0 = Function(X)  # u^{n-1}
up1 = Function(X)  # u^{n}

# Trial functions
u,p = TrialFunctions(X)

# Test functions
v,q = TestFunctions(X)

# Boundary condition
vin  = inlet_velocity()
bcin = DirichletBC(X.sub(0), vin,   boundaries, 1)
bccyl= DirichletBC(X.sub(0), (0,0), boundaries, 2)
bcsid= DirichletBC(X.sub(0).sub(1), 0, boundaries, 4)
bcs  = [bcin, bccyl, bcsid]

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

up0.interpolate(initial_condition())
u0 = as_vector((up0[0], up0[1]))
u1 = as_vector((up1[0], up1[1]))

fu = File("steady.pvd")

# Linear form
F  = inner(grad(u)*u0, v)*dx      \
   - p*div(v)*dx                   \
   + nu*inner(grad(u), grad(v))*dx \
   - q*div(u)*dx

# Right hand side is zero
L = Constant(0)*inner(u0,v)*dx \
  + Constant(0)*q*dx

# Picard iteration
for i in range(20):
    A = assemble(F)
    b = assemble(L)
    solver = LUSolver(A)
    [bc.apply(A,b) for bc in bcs]
    solver.solve(up1.vector(), b)
    up0.assign(up1)
    # save solution to file
    u,p = up1.split()
    fu << u
    # compute residual norm
    res = assemble(action(F, up1))
    [bc.apply(res) for bc in bcs0]
    res_norm = norm(res)/sqrt(X.dim())
    print "Iter = %d,  res norm = %e" % (i, res_norm)
    if res_norm < 1.0e-14:
        break

File("steady.xml") << up1.vector()
