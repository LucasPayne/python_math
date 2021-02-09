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


from functools import reduce
import operator
def prod(lis):
    return reduce(operator.mul, lis, 1)



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



def symmetric_tensor_indices(index):
    # Get all corresponding tensor indices to a symmetric tensor index.
    # The number of these is the multinomial coefficient (rank, index).
    any_to_give = False
    for position,i in enumerate(index):
        if i > 0:
            any_to_give = True
            new_index = [j for j in index]
            new_index[position] -= 1
            for trailing_index in symmetric_tensor_indices(new_index):
                yield [position] + trailing_index
    if not any_to_give:
        yield []

def stensor_to_tensor_index(k, rank, index):
    assert(len(index) == k)
    assert(sum(index) == rank)
    tensor_index = [0 for _ in range(rank)]
    pos = 0
    for basis_index, multiplicity in enumerate(index):
        for i in range(multiplicity):
            tensor_index[pos + i] = basis_index
        pos += multiplicity
    return tensor_index

def tensor_to_stensor_index(k, rank, index):
    stensor_index = [0 for _ in range(k)]
    for i in index:
        stensor_index[i] += 1
    return tuple(stensor_index)

def stensor_to_flat_index(k, rank, stensor_index):
    tensor_index = stensor_to_tensor_index(stensor_index)
    flat_index = 0
    for position, i in enumerate(tensor_index):
        flat_index += rising_factorial(k-1-i, rank - position) // factorial(rank - position)
    return flat_index

def tensor_to_flat_index(k, rank, tensor_index):
    # The tensor index is assumed to index into a symmetric tensor.
    # Brute force conversion.
    stensor_index = tensor_to_stensor_index(k, rank, tensor_index)
    for i,index in enumerate(stensor_indices(k, rank)):
        if index == stensor_index:
            return i
    assert(0)
    # #--------- this version doesn't work
    # flat_index = 0
    # for position, i in enumerate(tensor_index):
    #     flat_index += rising_factorial(k-1-i, rank - position) // factorial(rank - position)
    # return flat_index


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

        
    def indices(self):
        yield from stensor_indices(self.k, self.rank)

    def to_tensor_index(self, index):
        # print(len(index), self.k)
        assert(len(index) == self.k)
        # print(sum(index), self.rank)
        assert(sum(index) == self.rank)
        tensor_index = [0 for _ in range(self.rank)]
        pos = 0
        for basis_index, multiplicity in enumerate(index):
            for i in range(multiplicity):
                tensor_index[pos + i] = basis_index
            pos += multiplicity
        return tensor_index

    def to_stensor_index(self, tensor_index):
        stensor_index = [0 for _ in range(self.rank+1)]
        for i in tensor_index:
            stensor_index[i] += 1
        return tuple(stensor_index)


    def to_flat_index(self, index):
        # Brute force conversion to index into flat array.
        for i,other_index in enumerate(self.indices()):
            if index == other_index:
                return i
        assert(0)
        # ---this version doesn't work
        # tensor_index = self.to_tensor_index(index)
        # flat_index = 0
        # for position, i in enumerate(tensor_index):
        #     flat_index += rising_factorial(self.k-1-i, self.rank - position) // factorial(self.rank - position)
        # return flat_index

    def __getitem__(self, index):
        return self._vals[self.to_flat_index(index)]

    def set(self, index, value):
        self._vals[self.to_flat_index(index)] = value


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


def _multinomial_coefficient(n, values, prefix):
    if prefix == 2:
        return sym.binomial(n, values[1])
    return sym.binomial(n, values[prefix-1]) * _multinomial_coefficient(n - values[prefix-1], values, prefix-1)
def multinomial_coefficient(n, values):
    assert(sum(values) == n)
    return _multinomial_coefficient(n, values, len(values))


def transform_bezier_simplex(k, degree, M):
    print("------------------------------------------------------------")
    print("k={}, degree={}".format(k, degree))
    print("M =")
    print_matrix(M)
    assert(M.rows == k and M.cols == k)
    rank = degree
    control_points = stensor(k, rank, sym.symbols(" ".join("P_" + "".join(str(i) for i in index) for index in stensor_indices(k, rank))))
    masks = stensor(k, rank)
    for mask_index in masks.indices():
        mask_tensor_index = masks.to_tensor_index(mask_index)
        masks.set(mask_index, sum(
            sum(prod(M[tensor_index[j], mask_tensor_index[j]] for j in range(rank)) for tensor_index in symmetric_tensor_indices(index))
            * control_points[index]
            for index in control_points.indices()
        ))
        print(masks[mask_index])


m = Mat([
    [Rat(1,2), Rat(1,2), 0],
    [0,        Rat(1,2), Rat(1,2)],
    [Rat(1,2), 0,        Rat(1,2)]
])
for degree in range(1, 5):
    transform_bezier_simplex(3, degree, m)



transform_bezier_simplex(2, 4, Mat([[1,Half],[0,Half]]))

# print(multinomial_coefficient(6, [3,2,1]))
# c = 0
# test = [3,2,6,1]
# for index in symmetric_tensor_indices(test):
#     c += 1
#     print(index)
# print(c, multinomial_coefficient(sum(test), test))



def deboor_to_bezier_knot_mask(top_index, mask_extent):
    top_index = list(top_index)
    knot_mask = [tuple(top_index)]
    for depth in range(1, mask_extent+1):
        base_index = top_index[:]
        base_index[0] -= depth
        for comb in itertools.combinations_with_replacement(range(1, len(top_index)), depth):
            index = base_index[:]
            for i in comb:
                index[i] += 1
            knot_mask.append(tuple(index))
    return knot_mask
    

def deboor_to_bezier(domain_dim, degree, knots):
    assert(knots.rank == 2*degree-1)
    assert(knots.k == domain_dim+1)
    for deboor_net_index in stensor_indices(degree+1, domain_dim):
        top_knot_mask_index = list(deboor_net_index)
        top_knot_mask_index[0] += degree-1
        print("---")
        print(top_knot_mask_index)
        print("---")
        # knot_mask = [knot[knot_index
        for knot_index in deboor_to_bezier_knot_mask(top_knot_mask_index, degree-1):
            print(knot_index)
            print(knots[knot_index])

deboor_to_bezier_knot_mask((2,0,1), 2)
deboor_to_bezier(2, 2, stensor(3, 3, [index for index in stensor_indices(3,3)]))


def figurate_number(n, k):
    return rising_factorial(k, n) // factorial(n)


def deboor_to_bezier(domain_simplex_dim, continuity, knots=None):
    #----- This only works for domain_simplex_dim == 2. The attempt of generalization to higher simplex dimension did not work or really make sense.
    print("============================================================")
    print("deboor_to_bezier, domain_simplex_dim={}, continuity={}".format(domain_simplex_dim, continuity))
    print("------------------------------------------------------------")

    # Compute the de Boor net width as a figurate number, then check
    # that the knot tensor has the right shape.
    # ------------------------------------------------------------
    deboor_net_width = figurate_number(domain_simplex_dim-1, continuity+1)+1 #---
    num_deboor_points = figurate_number(domain_simplex_dim-1, deboor_net_width)
    degree = deboor_net_width-1
    knots_width = deboor_net_width + continuity
    num_knots = figurate_number(domain_simplex_dim-1, knots_width)
    print("deboor_net_width:", deboor_net_width)
    print("num_deboor_points:", num_deboor_points)
    print("degree:", degree)
    print("knots_width:", knots_width)
    print("num_knots:", num_knots)
    if knots == None:
        knots = stensor(domain_simplex_dim, knots_width-1)
        for index in knots.indices():
            knots.set(index, index) # just some distinct value for now
    elif type(knots) == list:
        knots = stensor(domain_simplex_dim, knots_width-1, knots)
    assert(knots.k == domain_simplex_dim)
    assert(knots.rank == knots_width-1)


    deboor_weights_matrix = []
    for deboor_index in stensor_indices(domain_simplex_dim, degree):
        # For each de Boor net index, get a list of corresponding knot
        # indices into the knot tensor.
        # ------------------------------------------------------------
        print("------------------------------------------------------------")
        print("deboor_index:", deboor_index)
        print("---")
        top_knot_mask_index = list(deboor_index)
        top_knot_mask_index[0] += continuity
        knot_mask_indices = deboor_to_bezier_knot_mask(top_knot_mask_index, continuity)
        print("knot_mask_indices:", knot_mask_indices)
        knot_mask = [knots[knot_index] for knot_index in knot_mask_indices]
        print("knot_mask:", knot_mask)
        
        # The values in the knot tensor are points in the affine domain. The de Boor net point p (of index deboor_index) is
        # p = f(...) where f is multi-affine and symmetric and the inputs are the knots in the knot mask of p, knot_mask.
        # Expanding p = f(...) gives p as an affine combination of Bezier points f(e_i,...), which are multiset inputs
        # of the affine domain basis simplex that the knots are expressed in.

        affine_weights = [0 for _ in range(num_deboor_points)] # one for each Bezier point.
        for multiindex in itertools.product(range(domain_simplex_dim), repeat=degree):
            flat_index = tensor_to_flat_index(domain_simplex_dim, degree, multiindex)
            affine_weights[flat_index] += prod(knot_mask[position][i] for position,i in enumerate(multiindex))
        print("weights:", affine_weights)

        deboor_weights_matrix.append(affine_weights)
        
        print(knot_mask)

    # Check that the rows sum to 1. If not, something definitely went wrong.
    for row in deboor_weights_matrix:
        assert(sum(row) == 1)

    deboor_weights_matrix = sym.Matrix(deboor_weights_matrix)
    print_matrix(deboor_weights_matrix)
    print_matrix(deboor_weights_matrix.inv())

    deboor_to_bezier_matrix = deboor_weights_matrix.inv()
    for i,index in enumerate(stensor_indices(domain_simplex_dim, degree)):
        print("{}: {}".format(index, ", ".join(str(c) for c in list(deboor_to_bezier_matrix.row(i)))))
    
    


def deboor_to_bezier_uniform_curve(continuity):
    deboor_to_bezier(2, continuity, [(continuity+1-i, -continuity+i) for i in range(2*(continuity+1))])

for i in range(6):
    deboor_to_bezier_uniform_curve(i)



def deboor_to_bezier_uniform_surface(continuity):
    knots = []
    degree = figurate_number(3-1, continuity+1)
    for index in stensor_indices(3, figurate_number(3-1, continuity+1)+continuity):
        knot = list(index)
        for i in range(len(knot)):
            knot[i] -= continuity
        knots.append(tuple(knot))
    deboor_to_bezier(3, continuity, knots)
    for k in knots:
        print(k)
    print(len(knots))

deboor_to_bezier_uniform_surface(0)
deboor_to_bezier_uniform_surface(1)



def bezier_curve_degree_elevation(m, n):
    # Multiply the Bezier curve polynomial of degree m by (u+v)^(n-m) and regroup terms
    # to get Bezier coefficients for the equivalent polynomial of degree n.
    assert(n > m)
    u,v = sym.symbols("u v")
    ps = sym.symbols(" ".join("p_{}".format(i) for i in range(m+1)))
    f = 0
    for i in range(m+1):
        f += ps[i] * sym.binomial(m, i) * u**i * v**(m-i)
    f = f.as_poly()
    print(f)
    g = (u + v)**(n - m) * f
    print(g)
    M = []
    for i in range(n+1):
        coeffs = []
        for j in range(m+1):
            coeffs.append(g.coeff_monomial(ps[j] * u**i * v**(n-i)) / sym.binomial(n, i))
        print(coeffs)
        M.append(coeffs)
    M = sym.Matrix(M)
    print_matrix(M)

bezier_curve_degree_elevation(2, 4)
bezier_curve_degree_elevation(1, 2)


M = Mat([[1,-2,1],[0,2,-2],[0,0,1]])
print_matrix(M)
print_matrix(M.inv())


def add_indices(index, add):
    assert(len(index) == len(add))
    return (index[i]+add[i] for i in range(len(index)))



def permutation_function(perm):
    return lambda lis: [lis[perm[i]] for i in range(len(lis))]

def cyclic_permutations(n):
    ns = list(range(n))
    for i in range(n-1):
        yield permutation_function(ns[i:]+ns[:i])


def print_triangle_stensor(s):
    assert(s.k == 3)
    lists = [[0 for _ in range(s.rank+1)] for i in range(s.rank+1)]
    for index in s.indices():
        lists[s.rank-index[0]][index[1]] = s[index]
    print_matrix(sym.Matrix(lists), lower_triangular=True)
        


def regular_triangular_bspline(continuity):
    assert(continuity % 2 == 0)

    patch_degree = 1
    patches_width = 3
    grid = stensor(3, 3)
    grid.set((1,1,1), 1)

    print_triangle_stensor(grid)
    
    # for iteration in range(continuity):
    #     for perm in cyclic_permutations(3):
    #         shift_subtract_grid = stensor(3, patch_degree*(patches_width+1))
    #         for index in grid:
    #             new_index = add_indices(index, perm((patch_degree, 0, 0)))
    #             shift_index = add_indices(new_index, perm((0, patch_degree, -patch_degree)))
    #             shift_subtract_grid.set(new_index, shift_subtract_grid[new_index] + grid[index])
    #             shift_subtract_grid.set(shift_index, shift_subtract_grid[shift_index] - grid[index])
            

            # new_patch_degree = patch_degree + 1
            # new_patches_width = new_patches_width + 1
            # new_grid = stensor(3, new_patch_degree * new_patches_width)
            # patch_degree = new_patch_degree
            # patches_width = new_patches_width
            # grid = new_grid

# lis = ["a", "b", "c"]
# for p in cyclic_permutations(3):
#     print(p(lis))

regular_triangular_bspline(2)

