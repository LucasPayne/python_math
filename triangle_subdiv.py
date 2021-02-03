import sympy as sym
from printing import print_coeffs, print_matrix
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols
Half = Rat(1,2)
Third = Rat(1,3)
Quarter = Rat(1,4)
def Rec(n):
    return Rat(1, n)


def trinomial(n, i,j,k):
    assert(i+j+k == n)
    return sym.binomial(n, i) * sym.binomial(n - i, j)


# ordered_sums(3, 5) will give all triples of sums 0+0+5, 0+1+4, ..., 3+1+1, etc., that add to 5.
def ordered_sums(terms, n):
    if terms == 1:
        yield tuple([n])
        return
    for i in range(0, n+1):
        for trailing in ordered_sums(terms-1, n-i):
            yield tuple([i]) + trailing

# Wrapper generator which gives a flat index as well.
def ordered_sums_indexed(terms, n):
    i = 0
    for t in ordered_sums(terms, n):
        yield tuple([i]) + t
        i += 1



def bezier_triangle_basis(degree):
    x,y,z = Sym("x y z")
    p = sym.poly((x + y + z)**degree)
    coeffs = []
    for i,j,k in ordered_sums(3, degree):
        coeffs.append(p.coeff_monomial(x**i * y**j * z**k) * x**i * y**j * z**k)
    return coeffs




# Given a triangular Bezier patch, there exist control points for the image of the middle triangle of the domain.
# These control points are convex combinations of the control points for the patch.
def triangular_bezier_subdivision_scheme(degree):
    print("============================================================")
    print("Triangular bezier patch interior triangle subdivision scheme")
    print("    degree = {}".format(degree))
    x,y,z = Sym("x y z")
    
    # An affine transformation of the parameter domain induces an affine transformation on the control points,
    # which can be represented as a matrix where each row sums to 1.
    # The weights in this matrix are computed by substituting the transformed x,y,z into the Bezier triangle polynomial and expanding, then
    # equating coefficients of each of the basis polynomials.
    ps = Sym(" ".join("p_{}{}{}".format(i,j,k) for i,j,k in ordered_sums(3, degree)))
    print(ps)
    q_poly = sym.poly(sum([(ps[index] * trinomial(degree, i,j,k) * ((x+y)/2)**i * ((y+z)/2)**j * ((z+x)/2)**k) for index,i,j,k in ordered_sums_indexed(3, degree)]))
    
    for i,j,k in ordered_sums(3, degree):
        c = sum(point * q_poly.coeff_monomial(x**i * y**j * z**k * point) / trinomial(degree, i,j,k) for point in ps)
        print("q_{}{}{} = {}".format(i,j,k, c))
    

for i in range(1,4):
    triangular_bezier_subdivision_scheme(i)



def quadratic_triangular_bspline_subdivision_scheme():
    X,Y,Z = Sym("X Y Z")
    x,y,z = Sym("x y z")

    basis = [
	X**2 / 2,
        Y**2 / 2,
        Z**2 / 2,
        X/2 + Y/2 + Z/2 + X*Y + X*Z + Y*Z
    ]

    ps = Sym("p_0 p_1 p_2 p_3")
    substituted_basis = [f.subs({
        X: (y+z)/2,
        Y: (x+z)/2,
        Z: (x+y)/2,
    }).expand() for f in basis]
    q_poly = sum(ps[i] * substituted_basis[i] for i in range(4)).expand()

    print(q_poly)
    print("--------------------------------------------------------------------------------")
    
    print(2 * q_poly.coeff(x, 2))
    print(2 * q_poly.coeff(y, 2))
    print(2 * q_poly.coeff(z, 2))
    f = x/2 + y/2 + z/2 + x*y + x*z + y*z

    g = q_poly
    g -= x**2 * q_poly.coeff(x, 2)
    print(sym.div(q_poly, f))
    g -= y**2 * q_poly.coeff(y, 2)
    g -= z**2 * q_poly.coeff(z, 2)
    g = g.expand()
    print(f)
    print(g)

    for f in substituted_basis:
        print(f)

    
    
quadratic_triangular_bspline_subdivision_scheme()

A = Mat([[0,0,0,Half,0,0,0,0,0],
         [0,0,0,0,Half,0,0,0,0],
         [0,0,0,0,0,Half,0,0,0],
         [Half,Half,Half,0,0,0,1,1,1]])
AM = Mat([[0,0,0,0,Rat(1,8),Rat(1,8),0,0,Quarter],
          [0,0,0,Rat(1,8),0,Rat(1,8),0,Quarter,0],
          [0,0,0,Rat(1,8),Rat(1,8),0,Quarter,0,0],
          [Half,Half,Half,Quarter,Quarter,Quarter,Rat(3,4),Rat(3,4),Rat(3,4)]])
print_matrix(AM * A.T * (A * A.T).inv())


A = Mat([[0,0,Half,0,0],
         [0,0,0,Half,0],
         [Half,Half,0,0,1]])
AM = Mat([[0,0,Half,0,0],
         [0,0,Rat(1,8),Rat(1,8),Quarter],
         [Rat(3,4), Rat(1,4), Rat(1,2), 0, Rat(1,2)]])

print_matrix(A)
print_matrix(AM)
print_matrix(AM * A.T * (A * A.T).inv())

x,y,z = Sym("x y z")
X,Y,Z = Sym("X Y Z")
print(((1-x)**2 / 2).expand())
print((x/2 + y/2 + x*y).subs(y, 1-x).expand())

A = Mat([[Half,-1,Half],
         [Half,1,-1],
         [0,0,Half]])
M = sym.diag(1, 2, 4).inv()
print_matrix(A * M * A.inv())


# def polynomial_transformation_matrix(num_vars, subst):
#     xs = sym.symbols(" ".join("x_{}".format(i) for i in range(1,num_vars+1)))


def fun(subst):
    x,y,z = Sym("x y z")
    X,Y,Z = Sym("X Y Z")
    polys = [X**2, Y**2, Z**2, X*Y, Y*Z, Z*X]
    M = []
    for p in polys:
        pp = p.subs(subst).expand().subs({
            x: X,
            y: Y,
            z: Z,
        }).as_poly()
        print(p, "|->", pp.as_expr())
        coeffs = []
        for coeff in polys:
            try:
                coeffs.append(pp.coeff_monomial(coeff))
            except:
                coeffs.append(0)
        M.append(coeffs)
    M = Mat(M)
    A = sym.diag(1,1,1,2,2,2)
    print_matrix(A * M * A.inv())
    

fun({
    X: (x+y)/2,
    Y: (y+z)/2,
    Z: (z+x)/2,
})
fun({
    X: x + y/2 + z/2,
    Y: y/2,
    Z: z/2,
})
