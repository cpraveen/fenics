from dolfin import *

nc   = 20
xmin = 0
xmax = 1

R      = 0.5
gamma  = 1.4
Pr     = 0.72
mu_ref = 0.0005
T_ref  = 300

Cp     = gamma*R/(gamma-1)

degree = 1
dt     = 0.001

mesh = Interval(nc, xmin, xmax)

Vh = FunctionSpace(mesh, "DG", degree)
Qh = MixedFunctionSpace([Vh, Vh, Vh])

V = Function(Qh)
Vo= Function(Qh)
W = TestFunction(Qh)

T   = -1/V[2]
u   = V[1]*V[2]
rho = ln(V[0])

p   = rho*R*T
E   = p/(gamma-1) + 0.5*rho*u**2

A0 = as_matrix([ \
      [rho,   rho*u,      E                   ], \
      [rho*u, p+rho*u**2, (E+p)*u             ], \
      [E,     (E+p)*u,    Cp*T*E+0.25*rho*u**4] \
      ])
A0 = replace(A0, {V:Vo})

# Coefficient of viscosity and heat conduction
mu = mu_ref*(T/T_ref)**0.8
k  = mu*Cp/Pr

# Shear stress and heat flux
tau = (4/3)*mu*Dx(u,0)
q   = -k*Dx(T,0)

# Total flux vector
F = as_vector([            \
      rho*u,              \
      p + rho*u**2 - tau, \
      (E+p)*u - tau*u + q \
      ])

def Finv(rhop, up, pp, rhom, um, pm):
   betap = 0.5*rhop/pp
   sp    = up*sqrt(betap)
   Ap    = 0.5*(1 + erf(sp))
   return as_vector([0, 0, 0])

Fi = Finv(rho('+'), u('+'), p('+'), rho('-'), u('-'), p('-'))

# Weak formulation
B =   (1/dt)*inner(A0*(V-Vo), W)*dx \
    - inner(F, Dx(W,0))*dx          \
    + inner(Fi, jump(W))*dS

problem = NonlinearVariationalProblem(B, V)
solver  = NonlinearVariationalSolver(problem)
