"""
Steady Flow over cylinder in external flow

Using selective frequency damping

THIS IS NOT WORKING YET.

No stress BC on top, bottom and outflow boundaries
Picard iteration on convective term
BDF1 in first step, BDF2 subsequently
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

uf0 = Function(V)
uf1 = Function(V)
uf  = TrialFunction(V)
vf  = TestFunction(V)

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

# These are used to estimate cfl number
DG = FunctionSpace(mesh, 'DG', 0)
vdg= TestFunction(DG)
h    = [cell.diameter() for cell in cells(mesh)]
area = [cell.volume()   for cell in cells(mesh)]

Ur = 1.0                # Reference velocity
D  = 0.1                # dia of cylinder
Re = 100.0              # Reynolds number
nu = Constant(Ur*D/Re)  # viscosity coefficient
dt = 0.002
idt= Constant(1.0/dt)
chi= 0.1
delta = 20.0
idelta = Constant(1.0/delta)

# Initial condition for velocity
up0.interpolate(initial_condition())
u0 = as_vector((up0[0], up0[1]))
u1 = as_vector((up1[0], up1[1]))

t  = 0.0
Tf = 500.0
it = 0
fu = File("u.pvd")
fw = File("vorticity.pvd")

# First time step: BDF1
F1 = idt*inner(u - u0, v)*dx       \
   + inner(grad(u0)*u0, v)*dx      \
   - p*div(v)*dx                   \
   + nu*inner(grad(u), grad(v))*dx \
   + chi*inner(u - uf0,v)*dx       \
   - q*div(u)*dx

a1 = lhs(F1)
L1 = rhs(F1)
A1 = assemble(a1)
solver1 = LUSolver(A1)
[bc.apply(A1) for bc in bcs]

# Filtering equation
F2 = idt*inner(uf - uf0, vf)*dx \
   - idelta*inner(u1 - uf, vf)*dx
a2 = lhs(F2)
L2 = rhs(F2)
A2 = assemble(a2)
solver2 = LUSolver(A2)
[bc.apply(A2) for bc in bcs]


while t < Tf:
    # estimate cfl number
    uavg = assemble(sqrt(u0[0]**2+u0[1]**2)*vdg*dx)
    uavg = uavg.array()/area
    cfl  = dt * max(uavg/h)

    # Picard iteration
    for i in range(1):
        b1 = assemble(L1)
        vin.t = t + dt
        [bc.apply(b1) for bc in bcs]
        res= A1 * up1.vector() - b1
        res_norm = norm(res)/sqrt(X.dim())
        print "%3d %12.4e" % (i, res_norm)
        solver1.solve(up1.vector(), b1)

        b2 = assemble(L2)
        [bc.apply(b2) for bc in bcs]
        solver2.solve(uf1.vector(), b2)

    up0.assign(up1)
    uf0.assign(uf1)
    t += dt
    it+= 1
    print "it = %6d,   t = %12.6e,   cfl = %e" % (it,t,cfl)
    if cfl > 10.0:
        print "cfl is too large !!!"
        break
    if it%100 == 0:
        u,p = up1.split()
        fu << u
