"""
Mark boundaries
"""

from dolfin import *

set_log_level(1)

# These values must be same as those in cylinder_in_channel.geo file
r1 = 1.0
r2 = 2.0
ht = 2.0

# Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
class Inner(SubDomain):
    def inside(self, x, on_boundary):
        return x[0]-r1<DOLFIN_EPS and on_boundary

# Sub domain for inflow (right)
class Outer(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] - r2 > -DOLFIN_EPS and on_boundary

# Sub domain for outflow (left)
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary

# Sub domain for cylinder
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] - ht > -DOLFIN_EPS and on_boundary

# Read mesh
mesh = Mesh("annulus.xml")

# Create mesh functions over the cell facets
sub_domains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)

# Mark all facets as sub domain 3
sub_domains.set_all(10)

# Mark cylinder as domain 4
inner = Inner()
inner.mark(sub_domains, 0)

# Mark inflow as sub domain 1, 01
outer = Outer()
outer.mark(sub_domains, 1)

# Mark outflow as sub domain 2, 0.2, True
bottom = Bottom()
bottom.mark(sub_domains, 2)

# Mark no-slip facets as sub domain 0, 0.0
# Mark all faces as noslip
top = Top()
top.mark(sub_domains, 3)

# Save sub domains to file
File("subdomains.xml") << sub_domains

# Save sub domains to VTK files
File("subdomains.pvd") << sub_domains
