from dolfin import *
from param import *

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
class Left(SubDomain):
   def inside(self, x, on_boundary):
      return near(x[0],0.0) and on_boundary

class Right(SubDomain):
   def inside(self, x, on_boundary):
      return near(x[0],1.0) and on_boundary

#------------------------------------------------------------------------------
# Stationary solution
#------------------------------------------------------------------------------
class StationarySolution(Expression):
   def eval(self, values, x):
      values[0] = rho_stat
      u = u_stat
      values[1] = values[0] * u

   def value_shape(self):
      return (2,)

#------------------------------------------------------------------------------
# Initial condition
#------------------------------------------------------------------------------
class InitialCondition(Expression):
   def eval(self, values, x):
      values[0] = 1.0
      u = 0.01*sin(2*pi*x[0])
      values[1] = values[0] * u

   def value_shape(self):
      return (2,)

#------------------------------------------------------------------------------
# VFRoe flux
#------------------------------------------------------------------------------
def VFRoe(Ul, Ur):
   Ua  = 0.5*(Ul + Ur)
   rho = Ua[0]
   u   = Ua[1]/rho
   c   = sqrt(gamma*rho**(gamma-1.0))
   l1  = u - c
   l2  = u + c
   dU  = Ur - Ul
   a1  = (l2*dU[0] - dU[1])/(2.0*c)
   a2  = dU[0] - a1
   l1m = Min(l1,0.0)
   fm1 = Ul[1] + l1m*a1
   fm2 = Ul[0]**gamma + Ul[1]**2/Ul[0] + l1m*a1*l1
   fm  = as_vector([fm1, fm2])
   l2p = Max(l2,0.0)
   fp1 = Ur[1] - l2p*a2
   fp2 = Ur[0]**gamma + Ur[1]**2/Ur[0] - l2p*a2*l2
   fp  = as_vector([fp1, fp2])
   return conditional( ge(u,0.0), fm, fp )

#------------------------------------------------------------------------------
# Lax-Friedrich flux
def LxF(Ul, Ur):
   Ua  = 0.5*(Ul + Ur)
   rho = Ua[0]
   u   = Ua[1]/rho
   c   = sqrt(gamma*rho**(gamma-1.0))
   lam = abs(u) + c
   f1  = 0.5*(Ul[1] + Ur[1])
   f2  = 0.5*(Ul[0]**gamma + Ul[1]**2/Ul[0] + Ur[0]**gamma + Ur[1]**2/Ur[0])
   res = as_vector([f1, f2]) - 0.5*lam*(Ur - Ul)
   return res
#------------------------------------------------------------------------------
mesh = IntervalMesh(nc, xmin, xmax)
h = (xmax-xmin)/nc

sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(100)
left = Left()
left.mark(sub_domains, 0)
right = Right()
right.mark(sub_domains, 1)
ds = Measure("ds")[sub_domains]

Vh = FunctionSpace(mesh, "DG", degree)
Qh = MixedFunctionSpace([Vh, Vh])

U = Function(Qh)
W = TestFunction(Qh)

rho = U[0]
m   = U[1]
u   = m/rho
p   = rho**gamma

# Matrix for viscous flux
D = mu *                  \
      as_matrix([         \
      [0,      0       ], \
      [-u/rho, 1.0/rho ]  \
      ])

Dl = mu *                  \
      as_matrix([         \
      [0,      0       ], \
      [-ul/rho, 1.0/rho ] \
      ])

Dr = mu *                  \
      as_matrix([         \
      [0,      0       ], \
      [-ur/rho, 1.0/rho ] \
      ])


# Shear stress
tau = mu*Dx(u,0)
Ux  = Dx(U,0)
Wx  = Dx(W,0)

# Total flux vector
Finv = as_vector([m, p + m*u])
Fvis = as_vector([0, tau])

Fi = LxF(U('+'), U('-'))
Fl = as_vector([rho*ul, p + rho*ul**2])
Fr = as_vector([rho*ur, p + rho*ur**2])
delta = Cip * avg(D) * jump(U) / h
#delta = Cip * mu/avg(rho) * jump(U) / h

Ul     = as_vector([rho, rho*ul])
Ur     = as_vector([rho, rho*ur])
deltal = Cip * Dl * (U - Ul) / h
deltar = Cip * Dr * (U - Ur) / h
#deltal = Cip * mu/rho * (U - Ul) / h
#deltar = Cip * mu/rho * (U - Ur) / h

# Weak formulation
# Volume terms
B_int_inv = - inner(Finv, Wx)*dx            \
            + inner(Fi, jump(W))*dS

B_int_vis =   inner(Fvis, Wx)*dx            \
            - inner(avg(D*Ux), jump(W))*dS   \
            - inner(avg(D.T*Wx), jump(U))*dS \
            + inner(delta, jump(W))*dS

# Boundary terms
B_bdy_inv = - inner(Fl, W)*ds(0)            \
            + inner(Fr, W)*ds(1)

B_bdy_vis =   inner(deltal, W)*ds(0)        \
            + inner(deltar, W)*ds(1)        \
            + inner(Dl*Ux, W)*ds(0)         \
            - inner(Dr*Ux, W)*ds(1)         \
            + inner(Dl.T*Wx, U-Ul)*ds(0)    \
            - inner(Dr.T*Wx, U-Ur)*ds(1)

B2 = B_int_inv + B_int_vis
B3 = B_bdy_inv + B_bdy_vis
