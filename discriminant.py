import sympy as sym
import numpy as np
import string
alphabet = string.ascii_lowercase



def substitution_general_discriminant(n):
    assert(n > 1)
    assert(n < 26)
    X = sym.symbols("X")
    syms = sym.symbols(" ".join(alphabet[:n+1]))
    P = sum(syms[i]*X**(n-i) for i in range(n+1))
    Pder = sym.diff(P, X)
    stationary = sym.solve(Pder, X)[0]
    Psubs = P.subs(X, stationary)
    
    print(80*"=")
    print("n = {}".format(n))
    print(80*"-")
    print("Polynomial:", P)
    print("Stationary location:", stationary)
    print("Discriminant:", Psubs)


def general_discriminant(degree):
    assert(degree > 1)
    assert(degree < 26)
    X = sym.symbols("X")
    syms = sym.symbols(" ".join(alphabet[:degree+1]))
    P = sum(syms[i]*X**(degree-i) for i in range(degree+1)).as_poly(X)
    Pder = sym.diff(P, X)
    P_vec = sym.Matrix(P.all_coeffs()).T
    Pder_vec = sym.Matrix(Pder.all_coeffs()).T
    sylvester_matrix = sym.zeros(2*degree-1,2*degree-1)
    print(P_vec)
    print(Pder_vec)

    n = degree
    m = degree-1
    for i in range(m):
        sylvester_matrix[i,i:i+n+1] = P_vec
    for i in range(n):
        sylvester_matrix[i+m,i:i+m+1] = Pder_vec
    discrim = sylvester_matrix.det()
    discrim = (-discrim/syms[0]).expand() # convention

    print("="*80)
    print("Degree {}".format(degree))
    print("-"*80)
    print("Sylvester matrix:", sylvester_matrix)
    print("Discriminant:", discrim)
    lprint("", sylvester_matrix)
    lprint("", discrim)

# for i in range(2, 5+1):
#     substitution_general_discriminant(i)

for i in range(2, 4+1):
    general_discriminant(i)

ldone()
