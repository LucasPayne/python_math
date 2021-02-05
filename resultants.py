import sympy as sym
import itertools
from printing import print_coeffs, print_matrix
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols
Half = Rat(1,2)
Third = Rat(1,3)
Quarter = Rat(1,4)
def Rec(n):
    return Rat(1, n)


def sylvester_matrix(p, q):
    degp = len(p)-1
    degq = len(q)-1
    n = degp + degq
    M = [[0 for _ in range(n)] for __ in range(n)]
    for i in range(degq):
        for j in range(degp+1):
            M[i+j][i] = p[j]
    for i in range(degp):
        for j in range(degq+1):
            M[i+j][degq+i] = q[j]
    return sym.Matrix(M)

S = sylvester_matrix([1,-1,3,-3], [3,-1,-2])
print_matrix(S)
print(S.det())
