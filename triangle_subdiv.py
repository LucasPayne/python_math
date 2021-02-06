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

# Compute the affine combinations of control points that give a triangular subpatch.
# The subpatch is given by three barycentric points, in the space of the triangle.
# (For example, inputting 1,0,0, 0,1,0, 0,0,1 will give the identity weights.)
def bezier_triangle_subpatch(degree, a0, b0, c0, a1, b1, c1, a2, b2, c2):
    x,y,z = Sym("x y z")
    X,Y,Z = Sym("X Y Z")
    monomials = [X**i * Y**j * Z**k for i,j,k in ordered_sums(3, degree)]
    M = []
    for mono in monomials:
        p = mono.subs({
            X: a0*x + a1*y + a2*z,
            Y: b0*x + b1*y + b2*z,
            Z: c0*x + c1*y + c2*z,
        }).expand().subs({
            x: X,
            y: Y,
            z: Z,
        }).as_poly()
        # print(mono, "|->", p.as_expr())
        coeffs = []
        for mono2 in monomials:
            try:
                coeffs.append(p.coeff_monomial(mono2))
            except:
                coeffs.append(0)
        M.append(coeffs)
    M = Mat(M)
    A = sym.diag(*[trinomial(degree, i,j,k)  for i,j,k in ordered_sums(3, degree)])
    K = A * M * A.inv()
    print_matrix(K)

    for col, i,j,k in ordered_sums_indexed(3, degree):
        string = "Q_{}{}{} = ".format(i,j,k)
        string += " +\n        ".join("({})P_{}{}{}".format(K[row,col], ii,jj,kk) for row, ii,jj,kk in ordered_sums_indexed(3, degree) if K[row,col] != 0)
        print(string + "\n")


def bezier_triangle_middle_subpatch(degree):
    bezier_triangle_subpatch(degree, 
                             Rat(1,2), Rat(1,2), 0,
                             0,        Rat(1,2), Rat(1,2),
                             Rat(1,2), 0,        Rat(1,2))
def bezier_triangle_bottom_left_subpatch(degree):
    bezier_triangle_subpatch(degree, 
                             1,        0,        0,
                             Rat(1,2), Rat(1,2), 0,
                             Rat(1,2), 0,        Rat(1,2))
    

bezier_triangle_middle_subpatch(2)
bezier_triangle_bottom_left_subpatch(2)

bezier_triangle_subpatch(2,
    Half,Half,0,
    0,1,0,
    0,0,1,
)
bezier_triangle_subpatch(2,
    Half,Half,0,
    0,0,1,
    Half,0,Half
)

# bezier_triangle_middle_subpatch(3)
# bezier_triangle_bottom_left_subpatch(3)


from math import factorial

def rising_factorial(k, n):
    t = 1
    for i in range(n):
        t *= k + i
    return t

def stensor_size(k, rank):
    # The size is the k'th rank-figurate number.
    # (e.g. 2-figurate numbers are triangular, 3-figurate are pyramidal, etc.)
    return rising_factorial(k, rank) // factorial(rank)


def stensor_indices(k, rank):
    if k == 1:
        yield tuple([rank])
        return
    for i in range(0, rank+1):
        for trailing in stensor_indices(k-1, rank-i):
            yield tuple([i]) + trailing

def stensor_indices_size(k, rank):
    if k == 1:
        return 1
    return sum(stensor_indices_size(k-1, rank-i) for i in range(0, rank+1))


class stensor:
    def __init__(self, k, rank, init_array=None):
        if init_array == None:
            init_array = []
            self._vals = [0 for _ in range(stensor_size(k, rank))]
        else:
            assert(len(init_array) == stensor_size(k, rank))
            self._vals = init_array
        self.k = k
        self.rank = rank


    def to_tensor_index(self, index):
        assert(len(index) == self.k)
        assert(sum(index) == self.rank)
        tensor_index = [0 for _ in range(self.rank)]
        pos = 0
        for basis_index, multiplicity in enumerate(index):
            for i in range(multiplicity):
                tensor_index[pos + i] = basis_index
            pos += multiplicity
        return tensor_index

    def __getitem__(self, index):
        tensor_index = self.to_tensor_index(index)
        flat_index = 0
        for position, i in enumerate(tensor_index):
            flat_index += rising_factorial(self.k-1-i, self.rank - position) // factorial(self.rank - position)
        # return _vals[flat_index]
        return flat_index
    


# def bezier_triangle_subpatch_2(degree, M):
#     # degree: Degree of the the polynomial map. Determines the rank of the symmetric tensor and therefore the number of control points in the net.
#     k = 3 # Number of points in a basis simplex of the affine domain.
#     # M: (k+1)x(k+1) affine domain transformation matrix acting as Mp on points. Columns must sum to 1.
# 
# 
#     rank = degree # Rank of the symmetric tensor representing the polar form of the polynomial map.
# 
#     # Form the symmetric tensor of symbolic control points.
#     tensor_vars = sym.symbols(" ".join("N_" + "".join(str(i) for i in index)) for index in itertools.product(k, repeat=rank))
# 
# 
# 
#     control_point_vars = symmetric_tensor(k, rank,
#         sym.symbols(" ".join("P_" + "".join(i) for i in index) for index in ordered_sums(k, rank)
#     )
# 
#     # control_point_vars = symmetric_tensor(sym.symbols(" ".join("P_" + "".join(i) for i in index) for index in ordered_sums(k+1, rank)))
# 
# 
# 
# 
#     # tensor_weights = [[0 for _ in range(degree+1)] for __ in range(degree+1)]
#     # tensor_weights = [0 for _ in range((k+1)**rank)]
# 
# 
# 
#     for i,j in itertools.product(range(degree+1), repeat=2):
#         tensor_weights[i][j] = sum(tensor_vars[ii][jj] * M[ii, i] * M[jj, j] for ii,jj in itertools.product(range(degree+1), repeat=2))
# 
# 
#     mat = sym.Matrix(tensor_weights)
#     print_matrix(mat)
# 
#     for i,j in itertools.product(range(degree+1), repeat=2):
#         tri_index = [0,0,0]
#         tri_index[i] += 1
#         tri_index[j] += 1
#         net_var = sym.symbols("P_" + "".join(str(k) for k in tri_index))
#         mat = mat.subs(tensor_vars[i][j], net_var)
#     print_matrix(mat, lower_triangular=True)
        
        

# for degree in [2,3]:
#     print("Identity")
#     bezier_triangle_subpatch_2(degree, sym.eye(3))
#     print("Lower-left isometric")
#     bezier_triangle_subpatch_2(degree, Mat([
#         [1, Half, Half],
#         [0, Half, 0],
#         [0, 0,    Half]
#     ]))
#     print("Middle isometric")
#     bezier_triangle_subpatch_2(degree, Mat([
#         [Rat(1,2), Rat(1,2), 0],
#         [0,        Rat(1,2), Rat(1,2)],
#         [Rat(1,2), 0,        Rat(1,2)]
#     ]))



for index in stensor_indices(3, 2):
    print(index)

print("---")
for index in stensor_indices(3, 3):
    print(index)
print(stensor_indices_size(3,3))


k = 3
r = 4
T = stensor(k, r)
for index in stensor_indices(k, r):
    print(index, T.to_tensor_index(index))
for index in stensor_indices(k, r):
    print(T[index])
