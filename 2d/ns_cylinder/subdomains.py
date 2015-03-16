"""
Mark boundaries
"""

from dolfin import *
import math

set_log_level(1)

# These values must be same as those in cylinder_in_channel.geo file
xmin = -1.5
xmax = 2.2
H    =  0.2

# Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
class Cylinder(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Sub domain for inflow (right)
class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return math.fabs(x[0] - xmin) < DOLFIN_EPS and on_boundary

# Sub domain for outflow (left)
class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return math.fabs(x[0] - xmax) < DOLFIN_EPS and on_boundary

# Sub domain for cylinder
class Noslip(SubDomain):
    def inside(self, x, on_boundary):
        return (math.fabs(x[1] - H) < DOLFIN_EPS or \
                math.fabs(x[1] + H) < DOLFIN_EPS) and on_boundary

# Read mesh
mesh = Mesh("cylinder_in_channel.xml")

# Create mesh functions over the cell facets
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Mark all facets as sub domain 3
sub_domains.set_all(10)

# Mark cylinder as domain 4
cylinder = Cylinder()
cylinder.mark(sub_domains, 0)

# Mark inflow as sub domain 1, 01
inflow = Inflow()
inflow.mark(sub_domains, 1)

# Mark outflow as sub domain 2, 0.2, True
outflow = Outflow()
outflow.mark(sub_domains, 2)

# Mark no-slip facets as sub domain 0, 0.0
# Mark all faces as noslip
noslip = Noslip()
noslip.mark(sub_domains, 3)

# Save sub domains to file
File("subdomains.xml") << sub_domains

# Save sub domains to VTK files
File("subdomains.pvd") << sub_domains
