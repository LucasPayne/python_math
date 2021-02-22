from math import sqrt
import sympy as sym
import sympy.vector as symvec
import itertools
from math import floor, ceil
from printing import print_coeffs, print_matrix, matrix_row_strings
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols
Half = Rat(1,2)
Third = Rat(1,3)
Quarter = Rat(1,4)
def Rec(n):
    return Rat(1, n)

from sympy.diffgeom.rn import R2_r
from sympy.diffgeom import WedgeProduct
ex, ey = R2_r.base_vectors()
dx, dy = R2_r.base_oneforms()
print(WedgeProduct(dx, dy))
print(WedgeProduct(dx, dy)(ex, ey))

J = Mat([[Sym("J_{}{}".format(i, j)) for j in range(1,3)] for i in range(1,4)])
print_matrix(J)
print_matrix(J.T * J)



class k_vector:
    def __init__(self, n):
        self.n = n

def E(n, *vecs):
    assert(len(vecs) <= n)
    k = len(vecs)
    if len(set(vecs)) < k:
        return 0
    v = k_vector(n)
    


J = Mat([[1,0,1],[4,1,-1]]).T
C = symvec.CoordSys3D("C")
v = symvec.cross(C.i + C.k, 4*C.i + C.j - C.k) 
print(sqrt(v.dot(v)))
print(sqrt(sym.det(J.T * J)))

E,G,F,e,g,f = Sym("E G F e g f")
I = Mat([[E, F], [F, G]])
II = Mat([[e, f], [f, g]])
print_matrix(I)
print_matrix(II)
M = (I.inv() * II).expand()
print_matrix(M)
[V,D] = M.diagonalize()
print_matrix(V.expand())
print_matrix(D.expand())

