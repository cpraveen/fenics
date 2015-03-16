from dolfin import *
from param  import *
from common import *
import scipy.io as sio
import numpy as np

Uf_stat = Function(Qh)
Uf_stat = project(StationarySolution(), Qh)

# Set initial condition
Uo = Function(Qh)
Uo = project(InitialCondition(), Qh)
U.assign(Uo)

# Needed by BDF2
Uoo = Function(Qh)

# Time and iteration counters
t, i = 0, 0

# First step : backward euler (BDF1)
B_BDF1 = (1.0/dt)*inner(U-Uo, W)*dx + B2 + B3

dU = TrialFunction(Qh)
dB_BDF1 = derivative(B_BDF1, U, dU)

problem = NonlinearVariationalProblem(B_BDF1, U, [], dB_BDF1)
solver = NonlinearVariationalSolver(problem)
solver.solve()
Uoo.assign(Uo)
Uo.assign(U)
t = t + dt
i = i + 1

# Second step onwards: BDF2
B_BDF2 = (1.0/dt)*inner(1.5*U-2.0*Uo+0.5*Uoo, W)*dx + B2 + B3
dB_BDF2 = derivative(B_BDF2, U, dU)

problem = NonlinearVariationalProblem(B_BDF2, U, [], dB_BDF2)
solver = NonlinearVariationalSolver(problem)

fo= File("sol.pvd")

energy = (U[0]**gamma/(gamma-1.0) + 0.5*U[1]**2/U[0])*dx
energy -= (Uf_stat[0]**gamma/(gamma-1.0) + 0.5*Uf_stat[1]**2/Uf_stat[0])*dx

# Load gain matrix
G = sio.loadmat('gain.mat')

while t < Tf:
   a = -np.dot(G['G'],Uo.vector().array()-Uf_stat.vector().array())
   ul.assign(a[0])
   ur.assign(a[1])
   solver.solve()
   Uoo.assign(Uo)
   Uo.assign(U)

   t = t + dt
   i = i + 1
   e = assemble(energy)
   print("Step = %g  %g  %18.12e %10.4e %10.4e" % (i, t, e, ul, ur))
   if i%10 == 0:
      fo << U
