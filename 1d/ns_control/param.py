from dolfin import *

degree = 2
nc   = 100
xmin = 0.0
xmax = 1.0

gamma  = 1.4
mu     = 0.01
Cip    = 2.0 * degree**2

# stationary solution
rho_stat = 1.0
u_stat = 0.0
U_stat = as_vector([rho_stat, rho_stat*u_stat])

# Shift to destabilize
omega = Constant(0.0);

# Boundary condition for velocity
ul     = Constant(0.0)
ur     = Constant(0.0)

# These two lines required for generating linear system.
# Otherwise comment them.
ul     = variable(ul)
ur     = variable(ur)

dt     = 1.0e-2   # time step
Tf     = 50.0     # final time
