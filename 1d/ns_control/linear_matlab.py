from dolfin import *
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as la
import scipy.io as sio

parameters.linear_algebra_backend = "uBLAS"

from param import *
from common import *

Uo = project(StationarySolution(), Qh)
U.assign(Uo)

B = - B2 - B3

dU = TrialFunction(Qh)
dB = derivative(B, U, dU)

M = assemble( inner(dU, W)*dx )
rows, cols, values = M.data()
M = sps.csc_matrix((values, cols, rows))
N = M.shape[0]
print "Size of M =",M.shape[0], M.shape[1]

A = assemble(dB)
rows, cols, values = A.data()
A = sps.csc_matrix((values, cols, rows))
N = A.shape[0]
print "Size of A =",A.shape[0],A.shape[1]

# Derivative wrt ul
B_ul = diff(B, ul)
B_ul = assemble(B_ul)
print B_ul
B_ul = B_ul.array()

# Derivative wrt ur
B_ur = diff(B, ur)
B_ur = assemble(B_ur)
print B_ur
B_ur = B_ur.array()

# Save matrices in matlab format
sio.savemat('linear.mat', mdict={'M':M, 'A':A, 'B_ul': B_ul, 'B_ur':B_ur})


# Now read it back
#M2 = sio.loadmat('M.mat')
#print M2
#a = M2['M']*Uo.vector().array()
#print a
