import sympy as sym
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



def cyclic_matrix(v):
    n = len(v)
    C = sym.Matrix([[0 for _ in range(n)] for __ in range(n)])
    for i in range(n):
        for j in range(n):
            C[(i+j)%n, i] = v[j]
    return C

def pointwise_mul(u, v):
    assert(u.rows == v.rows)
    return Mat([u[i, 0] * v[i, 0] for i in range(u.rows)])


vs = [
    [1,1,1,1],
    [1,-1,1,-1],
    [1,1j,-1,-1j],
    [-1,1j,1,-1j]
]
vsadj = [
    [1,1,1,1],
    [1,-1,1,-1],
    [1,-1j,-1,1j],
    [-1,-1j,1,1j]
]

# for v in vs:
#     vv = [x/2 for x in v]
#     print(cyclic_matrix(vv)*Mat(vv))

for v in vs:
    for vp in vs:
        vv = [x/2 for x in v]
        vvp = [x/2 for x in vp]
        print(cyclic_matrix(vv)*Mat(vvp))

P = Mat([[x/2 for x in row] for row in vs])
Padj = Mat([[x/2 for x in row] for row in vs])
print_matrix(P)
u = [9,1,3,8]
v = [1,2,3,4]
print(cyclic_matrix(u) * Mat(v))
print(cyclic_matrix(v) * Mat(u))

p = pointwise_mul(P*Mat(u), P*Mat(v))
print(p)
print((P * p).expand())
print_matrix(P)
print_matrix(P * P.T)
print_matrix(P.T * P)
