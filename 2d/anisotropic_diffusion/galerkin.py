from dolfin import *

xmin, xmax =-1.0, 1.0
ymin, ymax =-1.0, 1.0
chi_par = 0.01

degree = 1
nx = 50
ny = 50
Cip = 10.0

def qflux(B, u):
    return -chi_par * B * dot(B,grad(u))

class InitialCondition(Expression):
    def eval(self, value, x):
        if abs(x[0]) < 0.25 and abs(x[1]) < 0.25:
            value[0] = 10.0
        else:
            value[0] = 0.1
        return value

p0 = Point(xmin, ymin)
p1 = Point(xmax, ymax)
mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'CG', degree)

# Magnetic field
B = Expression(("1.0/sqrt(2.0)","-1.0/sqrt(2.0)"))
# Dirichlet bc
g = Constant(0.0)

u = TrialFunction(V)
v = TestFunction(V)

u0 = interpolate(InitialCondition(), V)

# Set initial condition

dt = 0.001
idt = Constant(1.0/dt)
n = FacetNormal(mesh)
qu = qflux(B,u)
F = idt*(u - u0)*v*dx       \
    - inner(qu, grad(v))*dx

a, L = lhs(F), rhs(F)

A  = PETScMatrix(); assemble(a, tensor=A)
solver = LUSolver(A)
solver.parameters['reuse_factorization'] = True

u = Function(V)
u.assign(u0)

Tf = 10.0
t  = 0.0
it = 0
sol = File('sol.pvd')
sol << u
while t < Tf:
    rhs = assemble(L)
    solver.solve(u.vector(), rhs)
    u0.assign(u)
    t += dt; it += 1
    print "it, dt, t = ", it, dt, t
    if it%100 == 0:
        sol << u
