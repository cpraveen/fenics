"""
Flow over cylinder in channel
Test case from Turek

GRPC algorithm from dolfin/nsbench, also see fenics book
"""
from dolfin import *

tau1, tau2 = 2.0, 2.0
maxiter = 100

Ur = 1.0                # Reference velocity
D  = 0.1                # dia of cylinder
Re = 100.0              # Reynolds number
Tf = 40.0

def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, nu):
    "Return stress tensor."
    return 2*nu*epsilon(u) - p*Identity(2)

class inlet_velocity(Expression):
   def __init__(self, t=0.0):
      self.t = t
   def eval(self, value, x):
      yc = x[1]/0.41
      ff = 1.0 - exp(-5.0*self.t)
      value[0] = 6.0*yc*(1.0-yc)
      value[1] = 0.0
      return value
   def value_shape(self):
      return (2,)

mesh = Mesh("cylinder.xml")
boundaries = MeshFunction("size_t", mesh, "cylinder_facet_region.xml")

# Function space
udeg = 2
pdeg = udeg - 1

V = VectorFunctionSpace(mesh, 'CG', udeg)
Q = FunctionSpace(mesh, 'CG', pdeg)

print "Velocity dofs = ", V.dim()
print "Pressure dofs = ", Q.dim()

u0 = interpolate(inlet_velocity(), V)
p0 = interpolate(Constant(0), Q)

u1 = interpolate(u0, V)
p01= interpolate(p0, Q)

# Trial functions
u  = TrialFunction(V)
p  = TrialFunction(Q)

# Test functions
v  = TestFunction(V)
q  = TestFunction(Q)

# Boundary condition
vin  = inlet_velocity(0.0)
bcin = DirichletBC(V, vin,   boundaries, 1)
bccyl= DirichletBC(V, (0,0), boundaries, 2)
bcwal= DirichletBC(V, (0,0), boundaries, 4)
bcu  = [bcin, bccyl, bcwal]

pbar = Constant(0.0)
bcout= DirichletBC(Q, pbar, boundaries, 3)
bcp  = [bcout]

nu = Constant(Ur*D/Re)  # viscosity coefficient
dt = 0.001
k  = Constant(dt)
n  = FacetNormal(mesh)

U = 0.5*(u0 + u1)
P = p01
Ru = inner(v, u1 - u0)*dx + k*inner(v, (grad(U)*U))*dx \
     + k*inner(epsilon(v), sigma(U, P, nu))*dx \
     - k*nu*inner(v, grad(U).T*n)*ds + k*inner(v, pbar*n)*ds
Rp = k*q*div(U)*dx

# Assemble preconditioners
ax  = inner(v, u)*dx + 0.5*k*inner(v, (grad(u)*u0))*dx \
    + 0.5*k*2*nu*inner(epsilon(v), epsilon(u))*dx \
    - 0.5*k*nu*inner(v, grad(u).T*n)*ds

ay1 = k**2*(inner(grad(q), grad(p)))*dx
ay2 = k**2*((1.0/(nu*k))*q*p)*dx

Kx  = assemble(ax)
Ky1 = assemble(ay1)
Ky2 = assemble(ay2)
[bc.apply(Kx) for bc in bcu]
[bc.apply(Ky1) for bc in bcp]

# Get solution vectors
x = u1.vector()
y = p01.vector()
delta_x = Vector(x)
delta_y = Vector(y)

fu = File("u.pvd")

t, it = 0.0, 0
while t < Tf:
    for iit in range(maxiter):
        # Velocity update
        rx = assemble(Ru)
        [bc.apply(rx, x) for bc in bcu]
        delta_x.zero()
        solve(Kx, delta_x, rx, 'gmres', 'ilu')
        x.axpy(-1.0, delta_x)

        # Pressure update
        ry = assemble(Rp)
        delta_y.zero()
        solve(Ky1, delta_y, ry, 'cg', 'petsc_amg')
        y.axpy(-tau1, delta_y)

        delta_y.zero()
        solve(Ky2, delta_y, ry, 'cg', 'jacobi')
        y.axpy(-tau2, delta_y)

        r = sqrt(norm(rx)**2 + norm(ry)**2)
        print iit, r
        if r < 1.0e-7: break
    t += dt; it += 1
    print "Time = ", t
    print "------------------------------------------------------------------"
    u0.assign(u1)
    if it%100 == 0:
        fu << u1
