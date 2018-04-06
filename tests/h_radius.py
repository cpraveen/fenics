'''
Construct a piecewise constant function which gives inradius and circumradius
Works in parallel also.
'''
import numpy as np
from dolfin import *

mesh = UnitSquareMesh(5,5)

# Iterate over cells to get circumradius/area/...
h1 = []; h2 = []
for cell in cells(mesh):
    h1.append(cell.circumradius())
    h2.append(cell.inradius())

V = FunctionSpace(mesh,'DG',0)
H1 = Function(V)
H1.vector().set_local(np.array(h1))
H1.vector().apply('insert')
H2 = Function(V)
H2.vector().set_local(np.array(h2))
H2.vector().apply('insert')

File('H1.pvd') << H1
File('H2.pvd') << H2
