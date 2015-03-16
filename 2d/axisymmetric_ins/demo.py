"""
Axisymmetric incompressible NS equations with swirl
See Batchelor appendix for the equations
"""
from dolfin import *

mu    = 0.01

omg   = 10   # angular speed in rad/sec

mesh = Mesh("annulus.xml")
sub_domains = MeshFunction("uint", mesh, "subdomains.xml")

n    = FacetNormal(mesh)
h    = CellSize(mesh)

# Taylor-Hood element
udeg = 3
pdeg = udeg - 1

Vry = VectorFunctionSpace(mesh, "CG", udeg)
Vt  = FunctionSpace(mesh, "CG", udeg)
Q   = FunctionSpace(mesh, "CG", pdeg)
X   = MixedFunctionSpace([Vry, Vt, Q])

V = Function(X)
W = TestFunction(X)

# Initial condition
V_init = Expression(("0", "0", "omg*x[0]", "0"), omg=omg)
V = interpolate(V_init, X)

ur  = V[0]
uy  = V[1]
ut  = V[2]
p   = V[3]

wr  = W[0]
wy  = W[1]
wt  = W[2]
q   = W[3]

r   = Expression("x[0]", element=Vt.ufl_element())

# Divergence
rdivu = ur + r*ur.dx(0) + r*uy.dx(1)
rdivw = wr + r*wr.dx(0) + r*wy.dx(1)

# Weak form
B =  r*(ur*Dx(ur,0) + uy*Dx(ur,1))*wr*dx \
   - ut**2*wr*dx                         \
   + mu*r*inner(grad(ur), grad(wr))*dx   \
   + (mu/r)*ur*wr*dx                     \
   + r*(ur*Dx(uy,0) + uy*Dx(uy,1))*wy*dx \
   + mu*r*inner(grad(uy), grad(wy))*dx   \
   + r*(ur*Dx(ut,0) + uy*Dx(ut,1))*wt*dx \
   + ur*ut*wt*dx                         \
   + mu*r*inner(grad(ut), grad(wt))*dx   \
   + (mu/r)*ut*wt*dx                     \
   - p*rdivw*dx                          \
   - q*rdivu*dx

dV = TrialFunction(X)
dB = derivative(B, V, dV)

# Boundary condition
ut_bc_value = Expression("omg*x[0]", omg=omg)

ury_inner_bc  = DirichletBC(X.sub(0), (0,0), sub_domains, 0)
ury_outer_bc  = DirichletBC(X.sub(0), (0,0), sub_domains, 1)
ury_bottom_bc = DirichletBC(X.sub(0), (0,0), sub_domains, 2)
ury_top_bc    = DirichletBC(X.sub(0), (0,0), sub_domains, 3)

ut_inner_bc  = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 0)
ut_outer_bc  = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 1)
ut_bottom_bc = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 2)
ut_top_bc    = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 3)

p_bc  = DirichletBC(X.sub(2), 0, "x[0]-1<DOLFIN_EPS && x[1]<DOLFIN_EPS", "pointwise")

bc    = [ury_inner_bc, ury_outer_bc, ury_bottom_bc, ury_top_bc, \
         ut_inner_bc, ut_outer_bc, ut_bottom_bc, ut_top_bc, \
         p_bc]

problem = NonlinearVariationalProblem(B, V, bc, dB)
solver  = NonlinearVariationalSolver(problem)

#solver.parameters["linear_solver"] = "gmres"
#itsolver = solver.parameters["newton_solver"]
#itsolver["absolute_tolerance"] = 1.0e-8

fury = File("ury.pvd", "compressed")
fut  = File("ut.pvd",  "compressed")
fp   = File("p.pvd",   "compressed")

solver.solve()
ury,ut,p = V.split()
fury << ury
fut  << ut
fp   << p
