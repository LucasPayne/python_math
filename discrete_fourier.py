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



def circulant(v):
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
    [1,-1j,-1,1j]
]

# for v in vs:
#     vv = [x/2 for x in v]
#     print(circulant(vv)*Mat(vv))

for v in vs:
    for vp in vs:
        vv = [x/2 for x in v]
        vvp = [x/2 for x in vp]
        print(circulant(vv)*Mat(vvp))

# Why cuberoot 4 here? -----
# P appears 3 times in transformation. Vectors in P have square norm of 4.
P = Mat([[x/sym.cbrt(4) for x in row] for row in vs])
print_matrix(P)
u = [9,1,3,8]
v = [1,2,3,4]
print(circulant(u) * Mat(v))
print(circulant(v) * Mat(u))

print_matrix(P)
print_matrix(P * P.H)
print_matrix((P.H * pointwise_mul(P*Mat(u), P*Mat(v))).expand())
print_matrix(circulant(u) * Mat(v))



def circulant_eig(v):
    [V,D] = circulant(v).diagonalize()
    print("============================================================")
    print_matrix(V.expand())
    print_matrix(D.expand())
# for i in range(1,5):
#     circulant_eig(list(range(1,i+1)))


