# import sympy as sym

from matplotlib import pyplot as plt
from math import sqrt,pi,cos,sin
import numpy as np


def cofactor_expansion(M):
    n = sym.shape(M)[0]
    if n == 1:
        return M[0,0]
    total = 0
    for i in range(n):
        minor = M.copy()
        minor.row_del(0)
        minor.col_del(i)
        total += (-1)**i * M[0,i] * cofactor_expansion(minor)
    return total



Avars = sym.symbols("a_1 a_2 a_3")
Bvars = [sym.symbols("b_{{{0}1}} b_{{{0}2}} b_{{{0}3}}".format(i)) for i in range(1,3+1)]
A = sym.Matrix(Avars)
Bs = [sym.Matrix(bvars) for bvars in Bvars]


M = sym.Matrix.vstack(
    *[
        sym.Matrix.hstack(
            (A - Bs[i]).T,
            sym.Matrix([sym.Rational(1,2)*(A + Bs[i]).dot(Bs[i] - A)])
        ) for i in range(3)
    ],
    sym.Matrix([0,0,0,1]).T
)
rhs = sym.Matrix([0,0,0,1])

lprint("", M)

Mdet_poly = cofactor_expansion(M)
lprint("", Mdet_poly)
for i in range(3):
    MM = M.copy()
    MM[:,i] = rhs
    MMdet_poly = cofactor_expansion(MM)
    lprint("$x_{}$".format(i+1), MMdet_poly/Mdet_poly)

# ldone()

def tetra(p,q,r,s):
    A = p
    Bs = [q,r,s]
    M = np.vstack((
        *[
            np.hstack((
                (A - Bs[i]).transpose(),
                np.array([0.5*np.dot(A + Bs[i], Bs[i] - A)])
            )) for i in range(3)
        ],
        np.array([0,0,0,1]).transpose()
    ))
    rhs = np.array([0,0,0,1])
    x = np.linalg.solve(M, rhs)[:3]
    for point in [p,q,r,s]:
        print(np.linalg.norm(x - p))
    
    
tetra(np.array([0,0,0]), np.array([2,0,0]), np.array([3,2,0]), np.array([1,1,3]))
