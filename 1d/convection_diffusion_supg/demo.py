"""
Solution of convection-diffusion equation using SUPG and Crank-Nicholson 

PDE: u_t + c u_x = mu u_xx in (-1,+1)
IC : u(x,0) = -sin(pi*x)
BC : periodic
"""

from dolfin import *

xmin = -1.0
xmax =  1.0

c      = 1.0
mu     = 0.001

nc     = 50 # no. of elements
degree = 1  # degree of FEM
T      = 5  # final time

mesh = Interval(nc, xmin, xmax)
h    = CellSize(mesh)

hmin = mesh.hmin()
dt   = min(hmin/c, hmin**2/mu)

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return bool(x[0]-xmin <  DOLFIN_EPS and \
                    on_boundary)

    def map(self, x, y):
        y[0] = x[0] - xmax + xmin

Vh = FunctionSpace(mesh, "CG", degree)

u = Function(Vh)
uo= Function(Vh)
v = TrialFunction(Vh)
w = TestFunction(Vh)

uinit = Expression("-sin(pi*x[0])")
u = interpolate(uinit, Vh)

ut = 0.5*(uo + v)
tau= h/(2*abs(c))
s  = mu*Dx(ut,0)

R =   (1/dt)*(v-uo)*w*dx \
    + c*Dx(ut,0)*w*dx             \
    + mu*Dx(ut,0)*Dx(w,0)*dx      \
    - mu*Dx(ut,0)*w*ds(1) + mu*Dx(ut,0)*w*ds(0) \
    + c*Dx(w,0)*tau*((v-uo)/dt + c*Dx(ut,0) - Dx(s,0))*dx

a = lhs(R)
L = rhs(R)

pbc= PeriodicBoundary()
bc = PeriodicBC(Vh, pbc)

ff = File("u.pvd", "compressed")
ff << u

t = 0
it= 0
while t < T:
   # Adjust dt so that we exactly reach T
   if t+dt > T:
      dt = T - t
   uo.assign(u)
   solve(a == L, u, bc)
   t = t + dt
   it= it + 1
   print "Iter=", it, ", t=", t
   ff << u
