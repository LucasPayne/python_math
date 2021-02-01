import sympy as sym
from printing import print_coeffs, print_matrix
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols


def coefficient_matrix(polys):
    return sym.Matrix([poly.all_coeffs()[::-1] for poly in polys])


def bezier_basis(degree):
    n = degree + 1
    x = Sym('x')
    polys = [p.as_poly() for p in [sym.binomial(n, k) * x**k * (1-x)**(n-k)for k in range(n+1)]]
    print_matrix(coefficient_matrix(polys))


for i in range(5):
    bezier_basis(i)
