"""
1) Flag a cell and refine it.
2) Flag and cell and its neighbors and refine
"""
from dolfin import *
import matplotlib.pyplot as plt

# Set-up mesh and connectivity
mesh = UnitSquareMesh(5,5)
tdim = mesh.topology().dim()
mesh.init(tdim - 1, tdim)

flag1 = CellFunction("bool", mesh)
flag2 = CellFunction("bool", mesh)

for cell in cells(mesh):
    flag1[cell] = False
    flag2[cell] = False

# Loop over cells
for cell in cells(mesh):
    if cell.index() == 15:
        flag1[cell] = True
        flag2[cell] = True
        for facet in facets(cell):
            # Filtered list of neighbor cells
            neighbor = filter(lambda ci: ci != cell.index(), 
                              facet.entities(tdim))
            ncell = Cell(mesh, neighbor[0])
            flag2[ncell] = True

mesh_new = adapt(mesh, flag2)

plt.figure()
plot(mesh)
plot(flag1)
plt.title("Flagged cell 15")

plt.figure()
plot(mesh_new)
plot(flag2)
plt.title("Flagged cell 15 and its neighbors")

plt.show()
