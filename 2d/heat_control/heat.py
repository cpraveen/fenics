import scipy.io as sio
from dolfin import *
from param import *

with_control = False
Tf = 1.0
dt = 0.01

# Read mesh from file
mesh = Mesh('mesh.xml')

V = FunctionSpace(mesh, 'CG', degree)
u, v = TrialFunction(V), TestFunction(V)

# Load feedback
if with_control:
    print('Reading gain matrix from file')
    gain = sio.loadmat('gain.mat')

uc = Function(V) # used in boundary condition
uc.vector()[:] = 0.0
bd = boundary()
bc = DirichletBC(V, uc, bd)
binds = bc.get_boundary_values().keys()

idt = Constant(1/dt)

u0 = Function(V) # sol at n-2
u1 = Function(V) # sol at n-1
u2 = Function(V) # sol at n

# Set initial condition
eigvec = Expression('eps*sin(pi*x[0])*sin(pi*x[1])',degree=degree,eps=1e-2)
u0 = interpolate(eigvec,V)
energy0 = sqrt(assemble(u0**2*dx))
print('Initial energy = %12.6e' % energy0)

# Open file to save some info
flog = open('log.txt','w')

# Time counter
t, it = 0.0, 0
flog.write('%5d %12.6e %12.6e\n' % (it,t,energy0))

# First time step: use BDF1
F1 = idt*(u - u0)*v*dx + inner(grad(u),grad(v))*dx - Constant(shift)*u*v*dx
a, L = lhs(F1), rhs(F1)

if with_control:
    uc.vector()[binds] = -np.dot(gain['Kt'], u0.vector().array())

solve(a == L, u1, bc)
t += dt; it += 1
energy = sqrt(assemble(u1**2*dx))
print('it,t,energy = %5d %12.6e %12.6e' % (it,t,energy))
flog.write('%5d %12.6e %12.6e\n' % (it,t,energy))

# Now define BDF2 for remaining steps
F2 = idt*(1.5*u - 2.0*u1 + 0.5*u0)*v*dx \
        + inner(grad(u),grad(v))*dx - Constant(shift)*u*v*dx
a, L = lhs(F2), rhs(F2)

while t < Tf:
    if with_control:
        uc.vector()[binds] = -np.dot(gain['Kt'], u1.vector().array())
    solve(a == L, u2, bc)
    energy = sqrt(assemble(u2**2*dx))
    t += dt; it += 1
    print('it,t,energy = %5d %12.6e %12.6e' % (it,t,energy))
    flog.write('%5d %12.6e %12.6e\n' % (it,t,energy))
    u0.assign(u1)
    u1.assign(u2)

flog.close()
