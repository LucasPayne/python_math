from matplotlib import pyplot as plt
from math import factorial,sqrt,pi,cos,sin
import numpy as np
import sympy as sym
from sympy.functions.special.tensor_functions import LeviCivita as levi_civita
import itertools
from itertools import combinations_with_replacement as repeat_combinations
import random

plotmin = -5
plotmax = 5
def get_grid():
    xs = np.linspace(plotmin, plotmax, 500)
    ys = np.linspace(plotmin, plotmax, 500)
    X,Y = np.meshgrid(xs, ys)
    return X,Y

def multinomial(n, choose):
    # Compute a multinomial coefficient.
    # e.g. multinomial(5, (2,1,2)) gives 30 (5_choose_2 * (5-2)_choose_1 * (5-2-1)_choose_2).
    assert(sum(choose) == n)
    cur = n
    coefficient = 1
    for c in choose:
        coefficient *= sym.binomial(cur, c)
        cur -= c
    return coefficient

def prod(lis):
    if len(lis) == 0:
        return 1
    t = lis[0]
    for val in lis[1:]:
        t *= val
    return t


def plot_line(l):
    X,Y = get_grid()
    plt.contour(X, Y, l[0]*X + l[1]*Y + l[2], 0)

# ordered_sums(3, 5) will give all triples of sums 0+0+5, 0+1+4, ..., 3+1+1, etc., that add to 5.
def ordered_sums(terms, n):
    if terms == 1:
        yield tuple([n])
        return
    for i in range(0, n+1):
        for trailing in ordered_sums(terms-1, n-i):
            yield tuple([i]) + trailing




def point_point_join(p, q):
    # p: Homogeneous 3-vector.
    # q: Homogeneous 3-vector.
    # Returns: Line p^q represented by a 3x3 antisymmetric matrix.
    matrix = sym.zeros(3,3)
    for i,j in itertools.product(range(3), repeat=2):
        matrix[i,j] = p[i]*q[j] - p[j]*q[i]
    return matrix

# def line_line_intersect(l1, l2):
#     # l1: Line represented by a 3x3 Plucker matrix.
#     # l2: Line represented by a 3x3 Plucker matrix.
#     # Returns: Point 3-vector.
    

def plucker_to_line(matrix):
    # Convert 3x3 Plucker matrix to 3-vector.
    vector = np.zeros(3)
    for i,j in itertools.product(range(3), repeat=2):
        for k in range(3):
            vector[k] += levi_civita(i,j,k) * matrix[i, j]
    return vector

def five_point_conic(p,q,r,s,t):
    x,y,z = sym.symbols("x y z")
    monoms = lambda v: [v[i]*v[j] for i,j in repeat_combinations(range(3), 2)]
    M = sym.Matrix([
        *[monoms([*point, 1]) for point in (p,q,r,s,t)],
        monoms((x,y,z))
    ])
    conic_equation = M.det()
    f = sym.lambdify((x, y), conic_equation.subs(z, 1))
    X,Y = get_grid()
    Z = f(X, Y)
    plt.contour(X, Y, Z, 0)
    plt.scatter(*zip(p,q,r,s,t))
    # for var in [x,y]:
    #     fprime = sym.lambdify((x,y), sym.diff(conic_equation, var))
    #     plt.contour(X, Y, fprime(X,Y), 0)

    # Polarize the quadratic form.
    conic_poly = conic_equation.as_poly()
    conic_polar = sym.zeros(3,3)
    var = (x,y,z)
        
    for i,j in itertools.product(range(3), repeat=2):
        powers = [0,0,0]
        powers[i] += 1
        powers[j] += 1
        conic_polar[i,j] = conic_poly.coeff_monomial(var[i]*var[j])/multinomial(2,powers)
    return conic_equation, conic_polar


def plot_conic(conic_equation):
    x,y,z = sym.symbols("x y z")
    f = sym.lambdify((x, y), conic_equation.subs(z, 1))
    X,Y = get_grid()
    Z = f(X, Y)
    plt.contour(X, Y, Z, 0)


def conic():
    p = np.array([-1,0])
    q = np.array([2,-1])
    r = np.array([3,2])
    s = np.array([2,2.33])
    t = np.array([1,-1])
    conic_equation, conic_polar = five_point_conic(p,q,r,s,t)
    print(conic_equation)
    print(conic_polar)
    # l = point_point_join(np.array([*p,1]), np.array([*q,1]))
    l = point_point_join(np.array([0,0,1]), np.array([1,1,0]))
    R = l.T*conic_polar*l

    lamb = sym.symbols(r"\lambda")
    shifted_minor = R + lamb*l
    shifted_minor.col_del(0)
    shifted_minor.row_del(0)
    lamb_sol = sym.solve(shifted_minor.det(), lamb)[0]
    print(lamb_sol)
    shifted_R = R + lamb_sol*l
    intersects = shifted_R.col(0), shifted_R.row(0)
    intersects = [np.array([inter[0]/inter[2], inter[1]/inter[2]]) for inter in intersects]
    plt.scatter(*zip(*intersects), s=100, c="r")

    plot_line(plucker_to_line(l))
    a,b,c = sym.symbols("a b c")
    k = sym.Matrix([a,b,c])
    poly = (l*k).T * conic_polar * (l*k)
    poly = (-poly[0]).expand()
    # The result is a degenerate dual conic. This is a product of two lines,
    # whose dual points are the intersection of the original line with the conic.
    print(poly)
    

    


# def levi_civita_contraction(num_upper, num_lower, tensor):
#     N = tensor.shape[0]
#     assert(upper == N) #-also must be a square tensor.
#     contracted = np.zeros((N,)*num_lower)
#     for contraction_multiindex in itertools.product(range(N), repeat=num_lower):
#         for multiindex in itertools.product(range(N), repeat=num_upper):
#             contracted[*contraction_multiindex] += levi_civita(*contraction_multiindex)
#     levi_civita()
    


conic()
plt.show()
