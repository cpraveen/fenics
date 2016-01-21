"""
Induction equation using DG
B_t + div(uB-Bu) + udiv(B) = 0 
in [-1,1]x[-1,1] or [0,1]x[0,1], u = (-y, x)
B0 = Initial data, g = 0 bdy data
Backward differentiation formula
y_n+2 - 4/3 y_n+1 + 1/3 y_n = 2/3 h f(t_n+2,y_n+2) 
Authors: Tanmay Sarkar, Praveen C

To get help
   python ./bdf2.py -h
"""
from dolfin import *
import math
import numpy
import argparse

def BForm(B,v,u,g,n):
   Fl = as_tensor(B[i]*u[j]-B[j]*u[i],(i,j))
   un = dot(u,n)
   Bn = dot(B,n)
   unp = 0.5*(un+abs(un))
   unm = 0.5*(un-abs(un))
   H =  unp('+')*B('+') + unm('+')*B('-') - u('+')*inner(avg(B),n('+'))
   Hb = B*unp + g*unm - u*Bn
   F1 = -inner(Fl, grad(v))*dx + dot(u,v)*div(B)*dx \
      + inner(H,jump(v))*dS - 2*dot(avg(u),avg(v))*avg(dot(B,n))*dS + inner(Hb,v)*ds
   return F1

def solve_induction(degree,np,itsave):
   mesh = UnitSquareMesh(np, np)
   X = VectorFunctionSpace(mesh, "DG", degree)

   B = TrialFunction(X)
   v = TestFunction(X)

   B1 = Function(X)
   # Velocity field
   u = Expression(("-x[1]", "x[0]"))
   # Exact solution
   ge = (("4.0*(-x[1]+0.5*sin(t))*exp(-20*(x[0]*x[0]+x[1]*x[1]-(x[0]*cos(t)+x[1]*sin(t))+0.25))",
         "4.0*(x[0]-0.5*cos(t))*exp(-20*(x[0]*x[0]+x[1]*x[1]-(x[0]*cos(t)+x[1]*sin(t))+0.25))"))
   g = Expression(ge,t=0.0)

   # Set initial condition
   B0 = interpolate(g, X)
   B1 = Function(X)

   # Save initial condition to file
   fsol = File("sol.pvd")
   B1.assign(B0)
   fsol << B1

   T = 0.5*pi
   h = 1.0/np
   dt = 0.5 * h
   N = int(T/dt)
   dt = T/N
   n = FacetNormal(mesh)

   Bt= 0.5*(B + B0)
   F = inner(B-B0,v)*dx + dt*BForm(Bt,v,u,g,n)
   
   a = lhs(F)
   L = rhs(F)
   A  = PETScMatrix(); assemble(a, tensor=A)
   solver = LUSolver(A)
   solver.parameters['reuse_factorization'] = True

   it, t = 0, 0.0
   while t < T:
      g.t = t + 0.5*dt
      b = assemble(L)
      solver.solve(B1.vector(), b)
      B0.assign(B1)
      it += 1; t += dt
      print "it, dt, t = ", it, dt, t
      if it%itsave == 0:
         fsol << B1

   # Compute error norms
   Be = Expression(ge,t=t)
   err_l2 = errornorm(Be, B1, 'l2')
   Bd = div(B1)**2*dx
   div_l2 = sqrt(assemble(Bd))
   # Save error into file
   Berr = Function(X)
   Bex  = interpolate(Be, X)
   Berr.vector()[:] = B1.vector() - Bex.vector()
   File("Berr.pvd") << Berr
   return div_l2, err_l2

if __name__ == "__main__" :
   parser = argparse.ArgumentParser()
   parser.add_argument('-deg',type=int,help='Degree of polynomial space',required=True)
   parser.add_argument('-N',type=int,nargs='+',help='No. of cells e.g., 20 40 80',required=True)
   parser.add_argument('-s',type=int,help='Interval to save results',default=1000000)
   args = parser.parse_args()

   err_l2 = numpy.zeros(len(args.N))
   div_l2 = numpy.zeros(len(args.N))
   for m,np in numpy.ndenumerate(args.N):
      (div_l2[m],err_l2[m]) = solve_induction(args.deg, np, args.s)
      print "np, div, err = ", np, div_l2[m], err_l2[m]

   print "l2 error = ", err_l2
   print "div error= ", div_l2

   print "L2 rate,        Div rate"
   for m in range(len(args.N)-1):
      rt = err_l2[m]/err_l2[m+1]
      p = math.log(rt)/math.log(2)
      rt = div_l2[m]/div_l2[m+1]
      q = math.log(rt)/math.log(2)
      print p, q
