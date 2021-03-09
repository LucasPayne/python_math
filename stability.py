import sympy as sym
import itertools
from matplotlib import pyplot as plt
import numpy as np
from scipy import linalg
from printing import print_coeffs, print_matrix
import string
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols
Half = Rat(1,2)
Third = Rat(1,3)
Quarter = Rat(1,4)
from math import factorial, sqrt, sin, cos, exp, floor, ceil
import time
def Rec(n):
    return Rat(1, n)
from finite_differences import finite_differences
import random
import cmath


def test_stability(k, y0, h, length_of_time, method="forward_euler"):
    assert(k.real < 0)
    num_nodes = ceil(length_of_time / h)
    nodes = np.linspace(0, length_of_time, num_nodes)
    h = nodes[1] - nodes[0]
    # Plot exact solution.
    exact = np.array([y0*cmath.exp(k*t) for t in nodes])
    plt.plot(nodes, np.real(exact), "r--")
    plt.plot(nodes, np.imag(exact), "b--")

    if method == "forward_euler":
        yi = lambda i: (1 + h*k)**i
    elif method == "trapezoidal_rule":
        yi = lambda i: ((1 + h*k/2)/(1 - h*k/2))**i
    y_approx = np.array([y0 * yi(i) for i in range(len(nodes))])

    plt.plot(nodes, np.real(y_approx), "r")
    plt.plot(nodes, np.imag(y_approx), "b")
    plt.title("{}, k={}, h={}".format(method, k, h))

    plt.show()


for method in ["forward_euler", "trapezoidal_rule"]:
    for k in [-1 + 0.5j, -1 + 2j, -1 + 1j]:
        # - A-stable. No matter the h, the solution should go to 0 as n->infinity.
        # - Very non-A-stable
        # - On barrier of A-stability region (not A-stable).
        n = 5
        for h in [0.1*i for i in range(1,n+1)]:
            test_stability(k, 1, h, 20, method)
