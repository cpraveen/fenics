"""
Based on primitive variables (p, u, T)
Ref: Lars Pesch, "Discontinuous Galerkin Finite Element Methods for the
Navier-Stokes Equations in Entropy Variable Formulatio", PhD Thesis
"""
from dolfin import *

degree = 1
parameters['form_compiler']['quadrature_degree'] = 2*degree

R     = 287.0
gamma = 1.4
mu    = 0.01
Pr    = 0.72

Cp    = gamma*R/(gamma-1)
Cv    = R/(gamma-1)
k     = mu*Cp/Pr
gamma1= gamma/(gamma-1)

omg   = 10  # angular speed in rad/sec
Tbc   = 300 # Boundary temperature

# Initial pressure
p_init = 1.0*R*Tbc

dt    = 0.01
Tf    = 1000*dt

def Boundary(x, on_boundary):
   return on_boundary

# Load mesh file
mesh = Mesh("annulus.xml")

n    = FacetNormal(mesh)
h    = CellSize(mesh)
hmin = mesh.hmin()

print "Viscous dt =", hmin**2/mu

V = FunctionSpace(mesh, "CG", degree)
X = VectorFunctionSpace(mesh, "CG", degree)
Vh= MixedFunctionSpace([V, X, V, V])

v = Function(Vh)
vo= Function(Vh)
w = TestFunction(Vh)

# Initial condition
v_init = Expression(("p", "0", "0", "omg*x[0]", "T"), \
                     p=p_init, omg=omg, T=Tbc)
v = interpolate(v_init, Vh)
vo.assign(v)

p   = v[0]
ur  = v[1]
uy  = v[2]
ut  = v[3]
T   = v[4]

rho = p/(R*T)
e   = p/(gamma-1) + 0.5*rho*(ur**2 + uy**2 + ut**2)
H   = gamma1*R*T + 0.5*(ur**2 + uy**2 + ut**2)
ke  = 0.5*(ur**2 + uy**2 + ut**2)

alpha_p = 1.0/T
beta_T  = 1.0/p

U   = as_vector([rho, rho*ur, rho*uy, rho*ut, e])
Uo  = replace(U, {v:vo})

r   = Expression("x[0]", element=V.ufl_element())

# Inviscid flux matrix
F = as_matrix([[rho*ur,      rho*uy     ], \
               [p+rho*ur**2, rho*ur*uy  ], \
               [rho*ur*uy,   p+rho*uy**2], \
               [rho*ur*ut,   rho*uy*ut  ], \
               [(e+p)*ur,    (e+p)*uy   ]])

# Inviscid boundary flux
Fb = as_vector([0, p*n[0], p*n[1], 0, 0])

# Divergence
rdivu   = ur   + r*ur.dx(0) + r*uy.dx(1)
divu    = ur/r + ur.dx(0)   + uy.dx(1)

# Stress tensor
rtau_rr = 2.0*mu*r*ur.dx(0) - (2.0/3.0)*mu*rdivu
rtau_ry = mu*r*(ur.dx(1) + uy.dx(0))
rtau_yr = rtau_ry
rtau_yy = 2.0*mu*r*uy.dx(1) - (2.0/3.0)*mu*rdivu
rtau_tt = 2.0*mu*ur - (2.0/3.0)*mu*rdivu
rtau_tr = mu*r*ut.dx(0) - mu*ut
rtau_rt = rtau_tr
rtau_ty = mu*r*ut.dx(1)
rtau_yt = rtau_ty

tau_tt = 2.0*mu*ur/r - (2.0/3.0)*mu*divu
tau_tr = mu*ut.dx(0) - mu*ut/r

# Heat flux vector
rq_r = -k*r*T.dx(0)
rq_y = -k*r*T.dx(1)

# Energy flux
er = rtau_rr*ur + rtau_ry*uy + rtau_rt*ut - rq_r
ey = rtau_yr*ur + rtau_yy*uy + rtau_yt*ut - rq_y

# Viscous flux matrix, includes radius
G = as_matrix([[0,       0,     ], \
               [rtau_rr, rtau_ry], \
               [rtau_yr, rtau_yy], \
               [rtau_tr, rtau_ty], \
               [er,      ey     ]])

# Source term
S = as_vector([0,                      \
               p + rho*ut**2 + tau_tt, \
               0,                      \
               -rho*ur*ut + tau_tr,    \
               0])

# Jacobian of conserved variable wrt v = dU/dV
Ao = as_matrix([ \
     [rho*beta_T,             0,      0,      0,      -rho*alpha_p         ], \
     [rho*beta_T*ur,          rho,    0,      0,      -rho*alpha_p*ur      ], \
     [rho*beta_T*uy,          0,      rho,    0,      -rho*alpha_p*uy      ], \
     [rho*beta_T*ut,          0,      0,      rho,    -rho*alpha_p*ut      ], \
     [beta_T*(e+p)-T*alpha_p, rho*ur, rho*uy, rho*ut, -alpha_p*(e+p)+rho*Cp]
     ])
Ao = replace(Ao, {v:vo})

# Weak form
B_GAL = (1.0/dt)*r*inner(U-Uo,w)*dx \
        - r*F[i,j]*Dx(w[i], j)*dx \
        + G[i,j]*Dx(w[i], j)*dx   \
        - S[i]*w[i]*dx
B_BCS = r*Fb[i]*w[i]*ds - w[i]*G[i,j]*n[j]*ds

# Jacobian of flux wrt v
g1 = beta_T*(e+p) - alpha_p*T + 1
g2 = -alpha_p*(e+p) + rho*Cp

Ar = as_matrix([ \
      [rho*beta_T*ur,      rho,            0,              0,         -rho*alpha_p*ur   ], \
      [rho*beta_T*ur*ur+1, 2.0*rho*ur,     0,              0,         -rho*alpha_p*ur*ur], \
      [rho*beta_T*ur*uy,   rho*uy,         rho*ur,         0,         -rho*alpha_p*ur*uy], \
      [rho*beta_T*ur*ut,   rho*ut,         0,              rho*ur,    -rho*alpha_p*ur*ut], \
      [g1*ur,              rho*ur*ur+e+p,  rho*ur*uy,      rho*ur*ut, g2*ur             ]
      ])

Ay = as_matrix([ \
      [rho*beta_T*uy,      0,              rho,            0,         -rho*alpha_p*uy   ], \
      [rho*beta_T*uy*ur,   rho*uy,         rho*ur,         0,         -rho*alpha_p*uy*ur], \
      [rho*beta_T*uy*uy+1, 0,              2.0*rho*uy,     0,         -rho*alpha_p*uy*uy], \
      [rho*beta_T*uy*ut,   0,              rho*ut,         rho*uy,    -rho*alpha_p*uy*ut], \
      [g1*uy,              rho*uy*ur,      rho*uy*uy+e+p,  rho*uy*ut, g2*uy             ]
      ])

# Stabilization matrix

# limit(Re)=Re if 0 <= Re < 1, else limit(Re)=1
def limit(Re):
   return conditional(lt(Re,1.0), Re, 1.0)

Reh  = rho*sqrt(ur**2 + uy**2 + ut**2)*h/(3.0*mu)
tauc = h*sqrt(ur**2 + uy**2 + ut**2)/2.0
taum = h*limit(Reh)/(2.0*rho*sqrt(ur**2 + uy**2 + ut**2))
taue = taum/Cv

om = 0.0

f1 = rho*(om+1.0)*(taum+(H-ke)*alpha_p*taue)
f2 = om*(taum + (H - ke)*alpha_p*taue)
f3 = (alpha_p*T - 1.0)*taue

tauY = as_matrix([\
      [tauc,         f1*ur,  f1*uy,  f1*ut,  2.0*rho*alpha_p*taue*ke], \
      [f2*ur,        taum,   0,      0,      alpha_p*taue*ur        ], \
      [f2*uy,        0,      taum,   0,      alpha_p*taue*uy        ], \
      [f2*ut,        0,      0,      taum,   alpha_p*taue*ut        ], \
      [-(H-ke)*taue, f3*ur,  f3*uy,  f3*ut,  taue                   ]  \
      ])

RES   = as_vector(r*Ao[i,j]*(v[j]-vo[j])/dt + Ar[i,j]*Dx(r*v[j],0) \
                  + Ay[i,j]*Dx(r*v[j],1) - Dx(G[i,j],j) - S[i], i)

# Ar and Ay must be transposed
LW    = as_vector(Ar[j,i]*Dx(w[j],0) + Ay[j,i]*Dx(w[j],1), i)
PSUP  = as_vector(LW[i]*tauY[i,j], j)
PSUP  = replace(PSUP, {v:vo})

# For derivative, we consider SUPG terms as constant
B_SUP = PSUP[i]*RES[i]*dx
B     = B_GAL + B_BCS + B_SUP

dv = TrialFunction(Vh)
dB = derivative(B, v, dv)

# Boundary condition
ur_bc_value = Expression("0")
uy_bc_value = Expression("0")
ut_bc_value = Expression("omg*x[0]", omg=omg)
T_bc_value  = Expression("Tbc", Tbc=Tbc)
ur_bc = DirichletBC(Vh.sub(1).sub(0), ur_bc_value, Boundary)
uy_bc = DirichletBC(Vh.sub(1).sub(1), uy_bc_value, Boundary)
ut_bc = DirichletBC(Vh.sub(2), ut_bc_value, Boundary)
T_bc  = DirichletBC(Vh.sub(3), T_bc_value,  Boundary)
bc    = [ur_bc, uy_bc, ut_bc, T_bc]

problem = NonlinearVariationalProblem(B, v, bc, dB)
solver  = NonlinearVariationalSolver(problem)

solver.parameters["linear_solver"] = "gmres"
itsolver = solver.parameters["newton_solver"]
itsolver["absolute_tolerance"] = 1.0e-8
itsolver["relative_tolerance"] = 1.0e-3

fp  = File("p.pvd",   "compressed")
fry = File("vel.pvd", "compressed")
fut = File("ut.pvd",  "compressed")
fT  = File("T.pvd",   "compressed")

iter = 0
t    = 0
while t < Tf:
   vo.assign(v)
   solver.solve()
   t    = t + dt
   iter = iter + 1
   print "Iter=", iter, ", t=", t
   if iter % 1 == 0:
      p,velry,ut,T = v.split()
      p.rename("p", "Pressure")
      fp  << p
      fry << velry
      fut << ut
      fT  << T
