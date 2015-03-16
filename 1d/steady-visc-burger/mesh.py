"""
This demo program solves Poisson's equation

    (u^2/2)_x - u_xx = f(x)

on the unit interval with source f given by

    f(x) = 9*pi^2*sin(3*pi*x[0])

and boundary conditions given by

    u(0) = u(1) = 0
"""

from dolfin import *

# Create mesh
mesh = UnitInterval(50)

# save mesh
File("mesh.xml") << mesh
