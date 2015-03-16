"""
Flow over cylinder

BDF1 in first step, BDF2 subsequently
Extrapolation for convection term: uext . grad(u)
"""
from dolfin import *

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
up0 = Function(X)  # u^{n-2}
up1 = Function(X)  # u^{n-1}
up2 = Function(X)  # u^{n}

# Trial functions
u,p = TrialFunctions(X)

# Test functions
v,q = TestFunctions(X)

# Boundary condition
bcin = DirichletBC(X.sub(0), (1.0,0),   boundaries, 1)
bccyl= DirichletBC(X.sub(0), (0,  0), boundaries, 2)
bcsid= DirichletBC(X.sub(0).sub(1), 0, boundaries, 4)
bcs  = [bcin, bccyl, bcsid]

# These are used to estimate cfl number
DG = FunctionSpace(mesh, 'DG', 0)
vdg= TestFunction(DG)
h    = [cell.diameter() for cell in cells(mesh)]
area = [cell.volume()   for cell in cells(mesh)]

Ur = 1.0                # Reference velocity
D  = 0.1                # dia of cylinder
Re = 100.0              # Reynolds number
nu = Constant(Ur*D/Re)  # viscosity coefficient
dt = 0.01
idt= Constant(1.0/dt)

#up0.interpolate(initial_condition())
File("steady.xml") >> up0.vector()
u0 = as_vector((up0[0], up0[1]))
u1 = as_vector((up1[0], up1[1]))
u2 = as_vector((up2[0], up2[1]))

# Compute force on cylinder
dss = Measure("ds")[boundaries]
# Stress tensor
T    = nu*(grad(u2) + grad(u2).T) - up2[2]*Identity(2)
# Face normals
n = FacetNormal(mesh)
drag = -T[0,j]*n[j]*dss(2)
lift = -T[1,j]*n[j]*dss(2)
ffile = open('force.dat', 'w')

t  = 0.0
Tf = 50.0
it = 0
fu = File("u.pvd")

# First time step: BDF1
F1 = idt*inner(u - u0, v)*dx       \
   + inner(grad(u)*u0, v)*dx      \
   - p*div(v)*dx                   \
   + nu*inner(grad(u), grad(v))*dx \
   - q*div(u)*dx

a  = lhs(F1)
L  = rhs(F1)

A  = PETScMatrix()
assemble(a, tensor=A)
solver = LUSolver(A)

b  = assemble(L)
[bc.apply(A,b) for bc in bcs]
solver.solve(up1.vector(), b)
t += dt
it+= 1

# Now switch to BDF2
uext = 2.0*u1 - u0
F2 = idt*inner(1.5*u - 2.0*u1 + 0.5*u0, v)*dx  \
   + inner(grad(u)*uext, v)*dx                  \
   - p*div(v)*dx                               \
   + nu*inner(grad(u), grad(v))*dx             \
   - q*div(u)*dx

a  = lhs(F2)
L  = rhs(F2)

A  = PETScMatrix()
solver = LUSolver(A)
solver.parameters['same_nonzero_pattern'] = True

while t < Tf:
    # estimate cfl number
    uavg = assemble(sqrt(u1[0]**2+u1[1]**2)*vdg*dx)
    uavg = uavg.array()/area
    cfl  = dt * max(uavg/h)

    if it==1:
        assemble(a, tensor=A)
    else:
        assemble(a, tensor=A, reset_sparsity=False)
    assemble(L, tensor=b)
    [bc.apply(A,b) for bc in bcs]
    solver.solve(up2.vector(), b)
    up0.assign(up1)
    up1.assign(up2)
    t += dt
    it+= 1
    print "it = %6d,   t = %12.6e,   cfl = %12.3e" % (it,t,cfl)
    # Compute lift/drag and store in arrays
    cl = assemble(lift)
    cd = assemble(drag)
    force=str(it)+" "+str(t)+" "+str(cl)+" "+str(cd)+"\n"
    ffile.write(force)
    ffile.flush()
    if cfl > 100.0:
        print "cfl is too large !!!"
        break
    if it%50 == 0:
        u,p = up2.split()
        fu << u
