"""
Compute area of unit square
Try:
    python parallel_assemble1.py
    mpirun -np 4 python parallel_assemble1.py
"""
from fenics import *

mesh = UnitSquareMesh(10,10)

area = assemble(Constant(1)*dx(mesh))

print "Area = ", area
