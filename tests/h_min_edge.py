'''
Construct a piecewise constant function which is equal to min edge length on
the cell.  Works in parallel also.
'''
import numpy as np
from dolfin import *

mesh = UnitSquareMesh(5,5)
hmin = np.array([min(edge.length() for edge in edges(cell)) for cell in         cells(mesh)])

V = FunctionSpace(mesh,'DG',0)
H = Function(V)
H.vector().set_local(hmin)
H.vector().apply('insert')

File('h.pvd') << H
