from dolfin import *

class boundary(SubDomain):
   def inside(self,x,on_boundary):
      return on_boundary


shift  = 3.0*pi**2
degree = 1
n      = 50
