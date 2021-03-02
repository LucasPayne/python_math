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
from numpy_utility import diagonal_form

a = -10
b = 10
num_nodes = 300
nodes = np.linspace(a, b, num_nodes)
h = (b - a) / num_nodes
stencil_radius = 2
weights = [w.subs(sym.symbols('h'), h) for w in finite_differences(range(-stencil_radius, stencil_radius+1), 1)]

# initial_func = lambda x: sin(2*x) + 0.2*cos(5*x - 0.9) + 0.13*sin(10*x + 0.3) + 3*exp(-x**2 / 2)
initial_func = lambda x: 6*exp(-x**2/2)
u = np.array([initial_func(x) for x in nodes])
plt.plot(nodes, u, "k")
plt.ion()
dt = 0.01

# F = lambda x,u,i: u[i]
# F = lambda x,u,i: sin(x)

# scalar conservation law
# F = lambda x,u,i: u[i]*0.2 * (0.5 + 0.25*sin(x))
F = lambda x,u,i: u[i]**2/2

# inhomogeneous transport
transport_source = lambda x,u,i: 3*max(0, sin(x) - 0.1)

while True:
    D = np.zeros((num_nodes, num_nodes))
    for weight_index,shift in enumerate(range(-stencil_radius, stencil_radius+1)):
        if shift < 0:
            np.fill_diagonal(D[-shift:,:], weights[weight_index])
        else:
            np.fill_diagonal(D[:,shift:], weights[weight_index])
    # M = np.identity(num_nodes)
    # M -= dt * D
    # lower = stencil_radius
    # upper = stencil_radius
    # M_diagonal_form = diagonal_form(M, lower, upper)
    # u = linalg.solve_banded((lower, upper), M_diagonal_form, u)

    # scalar conservation
    F_v = np.array([F(x, u, i) for i,x in enumerate(nodes)])
    u = u + dt*D.dot(F_v)

    # transport
    # source = [transport_source(x, u, i) for i,x in enumerate(nodes)]
    # transport_speed = 1.5
    # u = u + dt*(transport_speed * D.dot(u) + source)

    plt.clf()
    plt.xlim(-10, 10)
    plt.ylim(-6, 6)
    plt.plot(nodes, u, "k")
    plt.plot(nodes, F_v, "b--")
    # plt.plot(nodes, source, "y--")
    plt.axhline(y=0)
    plt.pause(0.01)

