"""
Scalar advection on circles using SUPG
a.grad(u) = 0 in [0,1]x[0,1],  a = (y, -x)
u = uin on x=0
"""
from dolfin import *

def Inlet(x, on_boundary):
   return x[0] < DOLFIN_EPS and on_boundary

np = 100

mesh = UnitSquareMesh(np, np)
h = CellSize(mesh)

X = FunctionSpace(mesh, "CG", 1)

v = TestFunction(X)
u = TrialFunction(X)

# Note: moda = |a|
# If you change "a", then change "moda" also.
a = Expression(("x[1]", "-x[0]"))
moda = Expression("sqrt(x[0]*x[0] + x[1]*x[1])")

L = inner(a, grad(u))*v*dx + (h/moda)*inner(a, grad(u))*inner(a, grad(v))*dx
b = Constant(0)*v*dx # Dummy zero rhs, remove in new version

uin = Expression("exp(-50*(x[1]-0.5)*(x[1]-0.5))")
bc  = DirichletBC(X, uin, Inlet)

U = Function(X)
solve(L == b, U, bc)
File("sol.pvd") << U
