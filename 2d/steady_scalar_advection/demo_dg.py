"""
Scalar advection on circles using DG
a.grad(u) = 0 in [0,1]x[0,1],  a = (y, -x)
u = uin on x=0

Since div(a), we write DG scheme for 
div(a*u) = 0
"""
from dolfin import *

def Inlet(x, on_boundary):
   return x[0] < DOLFIN_EPS and on_boundary

np = 100

mesh = UnitSquareMesh(np, np)

X = FunctionSpace(mesh, "DG", 1)

v = TestFunction(X)
u = TrialFunction(X)

class BoundaryValue(Expression):
   def eval(self, values, x):
      if x[0] < DOLFIN_EPS:
         values[0] = exp(-50*(x[1]-0.5)*(x[1]-0.5))
      else:
         values[0] = 0

uin = BoundaryValue()

a = Expression(("x[1]", "-x[0]"))

n = FacetNormal(mesh)
an= dot(a,n)
anp = 0.5*(an+abs(an))
anm = 0.5*(an-abs(an))
H = anp('+')*u('+') + anm('+')*u('-')

L = -inner(a, grad(v))*u*dx + H*jump(v)*dS + v*anp*u*ds
b = -anm*v*uin*ds

U = Function(X)
solve(L == b, U, bcs=None)
File("sol.pvd") << U
plot(U, interactive=True)
