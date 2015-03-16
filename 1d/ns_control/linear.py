from dolfin import *
from param import *
from common import *

Uo = project(StationarySolution(), Qh)
U.assign(Uo)

B = - B2 - B3

dU = TrialFunction(Qh)
dB = derivative(B, U, dU)

# M du/dt = A u
M = PETScMatrix()
assemble( inner(dU, W)*dx, tensor=M )
A = PETScMatrix()
assemble(dB, tensor=A)
print M
print A

eigensolver = SLEPcEigenSolver(A, M)
eigensolver.solve()

Ndof = Qh.dim()

for i in range(Ndof):
   r, c = eigensolver.get_eigenvalue(i)
   print r, c
