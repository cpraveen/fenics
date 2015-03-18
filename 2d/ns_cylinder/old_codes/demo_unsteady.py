"""
This demo program solves the unsteady incompressible Navier-Stokes equations
for cylinder in channel problem using Taylor-Hood elements.
   Author: Praveen. C
   www   : http://math.tifrbng.res.in/~praveen
"""

# Set parameter values
Re   = 200
D    = 0.1
Uinf = 1.0
nu   = D * Uinf / Re

theta= 0.5    # implicit scheme
cfl  = 1.0   # cfl number
tf   = 1e-0   # final time

# Functions to perturb the boundary velocity
def goo(t):
   if t <= -1.0:
      return 0.0
   elif t > 0.0:
      return 1.0
   else:
      return 1.0 - 3*t**2 - 2*t**3

def foo(t):
   return 0.1*(goo((t-0.2)/0.1) + goo(-(t-0.3)/0.1) - 1)

from dolfin import *

# Load mesh from file
mesh = Mesh("cylinder_in_channel.xml")
sub_domains = MeshFunction("size_t", mesh, "subdomains.xml")
dss = Measure("ds")[sub_domains]

hmin = mesh.hmin()
print "hmin =", hmin
dt = cfl * hmin / Uinf
print "dt   =", dt
tf = 2000*dt

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define test functions
(v,q) = TestFunctions(W)

# Define trial functions
w    = Function(W)
wold = Function(W)

# Initial condition is steady solution from file
File("steady.xml") >> w.vector()
wold.assign(w)

(u,p) = (as_vector((w[0], w[1])), w[2])
(uold,pold) = (as_vector((wold[0], wold[1])), wold[2])

utheta = (1.0-theta)*uold + theta*u
ptheta = (1.0-theta)*pold + theta*p

# Define boundary conditions
uinlet = Expression(("F*(1.0 - (x[1]/0.2)*(x[1]/0.2))", "0"), F=1.5)
cyl    = DirichletBC(W.sub(0), (0, 0), sub_domains, 0)
inlet  = DirichletBC(W.sub(0), uinlet, sub_domains, 1)
noslip = DirichletBC(W.sub(0), (0, 0), sub_domains, 3)
bc     = [noslip, inlet, cyl]

# Weak form
F =   (1.0/dt)*inner(u - uold, v)*dx     \
    + inner(grad(utheta)*utheta, v)*dx   \
    + nu*inner(grad(utheta), grad(v))*dx \
    - ptheta*div(v)*dx                   \
    - q*div(u)*dx

# Derivative of weak form
dw = TrialFunction(W)
dF = derivative(F, w, dw)

problem = NonlinearVariationalProblem(F, w, bc, dF)
solver  = NonlinearVariationalSolver(problem)
# Set linear solver parameters
itsolver = solver.parameters["newton_solver"]
itsolver["absolute_tolerance"] = 1.0e-10
itsolver["relative_tolerance"] = 1.0e-10

# To see various solver options, uncomment following line
#info(solver.parameters, True); quit()

# Compute force on cylinder
# Stress tensor
T    = nu*(grad(u) + grad(u).T) - p*Identity(2)
# Face normals
n = FacetNormal(mesh)
drag = -T[0,j]*n[j]*dss(0)
lift = -T[1,j]*n[j]*dss(0)

fv   = File("velocity.pvd", "compressed")
fp   = File("pressure.pvd", "compressed")
t    = 0.0
iter = 0

ffile = open('force.dat', 'w')
while t < tf:
   uinlet.F = 1.0 + foo(t+dt)
   wold.assign(w)
   solver.solve()
   if iter % 10 == 0:
      (u1,p1) = w.split()
      fv << u1
      fp << p1
   t    = t + dt
   iter = iter + 1
   # Compute lift/drag and store in arrays
   cl = assemble(lift)
   cd = assemble(drag)
   force=str(iter)+" "+str(t)+" "+str(cl)+" "+str(cd)+"\n"
   ffile.write(force)
   ffile.flush()
   print "iter =", iter, ",  t =", t, ",cl =", cl, ", cd =", cd
