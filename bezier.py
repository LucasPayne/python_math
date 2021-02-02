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

# It is easier to work directly with the polynomials. The above method manually expands this power.
def bezier_segment_basis(degree):
    x,y = Sym("x y")
    p = ((x + y)**degree).expand()
    return p

def bezier_triangle_basis(degree):
    x,y,z = Sym("x y z")
    p = ((x + y + z)**degree).expand()
    return p

def bezier_tetrahedron_basis(degree):
    x,y,z,w = Sym("x y z w")
    p = ((x + y + z + w)**degree).expand()
    return p

for i in range(6):
    print(bezier_segment_basis(i))

for i in range(4):
    print(bezier_tetrahedron_basis(i))

for i in range(4):
    print(bezier_triangle_basis(i))
