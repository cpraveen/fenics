from dolfin import *
from param import *

mesh = UnitSquareMesh(n,n)
print 'Saving mesh to mesh.xml'
File('mesh.xml') << mesh
