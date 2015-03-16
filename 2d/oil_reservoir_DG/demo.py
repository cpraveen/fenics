"""
2-D oil reservoir problem using DGFEM and DFLU flux
Author: Praveen. C
DOES NOT WORK YET
"""

from dolfin import *

# Material properties
mu_o   = 1.0
mu_w   = 1.0
pinlet = 1.0
poutlet= 0.0
sinlet = 1.0
cinlet = 0.0

class Inlet(SubDomain):
   def inside(self, x, on_boundary):
      return ((x[0] < DOLFIN_EPS and x[1]-0.1 < DOLFIN_EPS) or \
              (x[1] < DOLFIN_EPS and x[0]-0.1 < DOLFIN_EPS)) and \
             on_boundary

class Outlet(SubDomain):
   def inside(self, x, on_boundary):
      return ((x[0]-1 > -DOLFIN_EPS and x[1]-0.9 > -DOLFIN_EPS) or \
              (x[1]-1 > -DOLFIN_EPS and x[0]-0.9 > -DOLFIN_EPS)) and \
             on_boundary

def mobility_water(s):
   return s**2/mu_w

def mobility_oil(s):
   return (1-s)**2/mu_o

def dflu1(sl, sr, cl, cr, vn):
   if vn > 0.0:
      m_w = mobility_water (sl)
      m_o = mobility_oil (sl)
   else:
      m_w = mobility_water (sr)
      m_o = mobility_oil (sr)
   return vn * m_w / (m_w + m_o)

def dflu2(sl, sr, cl, cr, vn):
   flux = dflu1(sl, sr, cl, cr, vn)
   if vn > 0.0:
      return cl*flux
   else:
      return cr*flux

# Numerical parameters
np    = 50
theta = 0.5
dt    = 0.0001

mesh = UnitSquare(np, np)
n    = FacetNormal(mesh)

sub_domains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(100)
inlet = Inlet()
inlet.mark(sub_domains, 0)
outlet = Outlet()
outlet.mark(sub_domains, 1)

V = FunctionSpace(mesh, "DG", 0)
W = V * V

ts, tc = TestFunctions(W)

z = Function(W)
zold = Function(W)

# Set initial condition
ic = Expression(("0.0", "0.0"))
z  = project(ic, W)

# Saturation and concentration
s    = z[0]
c    = z[1]
sold = zold[0]
cold = zold[1]

stheta = (1-theta)*sold + theta*s
ctheta = (1-theta)*cold + theta*c

# Space for pressure
Q = FunctionSpace(mesh, "CG", 1)
r = TestFunction(Q)
q = TrialFunction(Q)
p = Function(Q)

lw_old = sold**2/mu_w
lo_old = (1-sold)**2/mu_o
lt_old = lw_old + lo_old

pa = lt_old*inner(grad(q), grad(r))*dx
pL = Constant(0)*r*dx
pbc_inlet  = DirichletBC(Q, pinlet,  sub_domains, 0)
pbc_outlet = DirichletBC(Q, poutlet, sub_domains, 1)
pbc = [pbc_inlet, pbc_outlet]

v = -lt_old*grad(p)
vn= avg(dot(v, n))

lw = stheta**2/mu_w
lo = (1-stheta)**2/mu_o
f  = lw/(lo + lw)
F  = v*f
cF = ctheta*F
H  = dflu1(stheta('+'), stheta('-'), ctheta('+'), ctheta('-'), vn)
cH = dflu2(stheta('+'), stheta('-'), ctheta('+'), ctheta('-'), vn)

vnb= dot(v,n)
Hi = dflu1(sinlet, sinlet, cinlet, cinlet, vnb)
cHi= dflu2(sinlet, sinlet, cinlet, cinlet, vnb)

dss = Measure("ds")[sub_domains]

L = (1/dt)*(s - sold)*ts*dx + (1/dt)*(s*c - sold*cold)*tc*dx \
    - inner(F, grad(ts))*dx - inner(cF, grad(tc))*dx \
    + H*jump(ts)*dS + cH*jump(tc)*dS \
    + Hi*ts*dss(0) + cHi*tc*dss(0) \
    + inner(F, n)*ts*dss(1) + inner(cF, n)*tc*dss(1)

#   + inner(F, n)*ts*dss(0) + inner(cF, n)*tc*dss(0) \

sbc = DirichletBC(W.sub(0), sinlet, sub_domains, 0, "geometric")
cbc = DirichletBC(W.sub(1), cinlet, sub_domains, 0, "geometric")
bc  = [sbc, cbc]

dz  = TrialFunction(W)
dL  = derivative(L, z, dz)

problem = NonlinearVariationalProblem(L, z, [], dL)
solver  = NonlinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "gmres"
#itsolver= solver.parameters["newton_solver"]
#info(solver.parameters, True)
#quit()

fp = File("p.pvd", "compressed")
fs = File("s.pvd", "compressed")
fc = File("c.pvd", "compressed")

t = 0
T = 100*dt
iter = 0
while t < T:
   zold.assign(z)
   solve(pa == pL, p, pbc)
   solver.solve()
   t = t + dt
   iter = iter + 1
   print "iter =", iter, ", t =", t
   if iter % 10 == 0:
      fp << p
      s1, c1 = z.split()
      fs << s1
      fc << c1
