from matplotlib import pyplot as plt
from math import sqrt,pi,cos,sin
import numpy as np
import sympy as sym
from itertools import combinations_with_replacement as repeat_combinations

plotmin = -5
plotmax = 5
def get_grid():
    xs = np.linspace(plotmin, plotmax, 100)
    ys = np.linspace(plotmin, plotmax, 100)
    X,Y = np.meshgrid(xs, ys)
    return X,Y


def plot_line(l):
    X,Y = get_grid()
    plt.contour(X, Y, l[0]*X + l[1]*Y + l[2], 0)


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
    
    
    return conic_equation


def conic():
    p = np.array([0,0])
    q = np.array([2,-1])
    r = np.array([3,2])
    s = np.array([2,2])
    t = np.array([1,-1])
    conic_equation = five_point_conic(p,q,r,s,t)
    l = np.array([3,-3.1222,1])
    plot_line(l)
    plt.show()
    

conic()
