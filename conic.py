from matplotlib import pyplot as plt
from math import factorial,sqrt,pi,cos,sin
import numpy as np
import sympy as sym
import itertools
from itertools import combinations_with_replacement as repeat_combinations

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
    conic_polar = np.zeros((3,3))
    var = (x,y,z)
        
    for i,j in itertools.product(range(3), repeat=2):
        powers = [0,0,0]
        powers[i] += 1
        powers[j] += 1
        conic_polar[i,j] = conic_poly.coeff_monomial(var[i]*var[j])/multinomial(2,powers)
    return conic_equation, conic_polar



def point_point_join(p, q):
    # p: Homogeneous 3-vector.
    # q: Homogeneous 3-vector.
    # Returns: Line p^q represented by a 3x3 antisymmetric matrix.
    matrix = np.zeros((3,3))
    for i,j in itertools.product(range(3), repeat=2):
        matrix[i,j] = p[i]*q[j] - p[j]*q[i]
    return matrix


def conic():
    p = np.array([0,0])
    q = np.array([2,-1])
    r = np.array([3,2])
    s = np.array([2,2])
    t = np.array([1,-1])
    conic_equation, conic_polar = five_point_conic(p,q,r,s,t)
    print(conic_equation)
    print(conic_polar)
    l = point_point_join(p, q)
    plot_line(l)


def points_determine_algebraic_form(n, points):
    assert(len(points) == (n+1)*(n+2)/2 - 1)
    x,y,z = sym.symbols("x y z")
    monoms = lambda v: [prod([v[i] for i in multiindex]) for multiindex in repeat_combinations(range(3), n)]
    M = sym.Matrix([
        *[monoms([*point, 1]) for point in points],
        monoms((x,y,z))
    ])
    equation = M.det()
    f = sym.lambdify((x, y), equation.subs(z, 1))
    X,Y = get_grid()
    Z = f(X, Y)
    plt.contour(X, Y, Z, 0)
    plt.scatter(*zip(*points))
    # for var in [x,y,z]:
    #     fprime = sym.lambdify((x,y,z), sym.diff(equation, var))
    #     plt.contour(X, Y, fprime(X,Y,1), 0)

    

    

# conic()

# points = [
#     np.array([0,0]),
#     np.array([2,-1.11]),
#     np.array([3.33,2]),
#     np.array([2,2.11]),
#     np.array([1,-1]),
#     np.array([0,-3.32]),
#     np.array([1.21321,-2.13213]),
#     np.array([0.22,-1.213]),
#     np.array([1,-1.0222213])
# ]
# points_determine_algebraic_form(3, points)

def algebraic_form(n):
    assert(n > 0)
    num_points = (n+1)*(n+2)//2 - 1
    points = [4*np.random.rand(2) for i in range(num_points)]
    points_determine_algebraic_form(n, points)
    nic_name = ["line", "conic", "cubic", "quartic", "quintic", "sextic", "septic"]
    plt.title("{} points determine a {}".format(num_points, nic_name[n-1] if n-1 < len(nic_name) else "{}-nic".format(n)))
    plt.show()

algebraic_form(1)
algebraic_form(2)
algebraic_form(3)
algebraic_form(4)
algebraic_form(5)
algebraic_form(6)
