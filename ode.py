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
from math import factorial, sqrt, sin
def Rec(n):
    return Rat(1, n)
from plotting import plot_linspace_func

def Vec2(x, y):
    return np.array([x, y])

u = lambda X: Vec2(1, 2*X[0])

b = 6
num_points = 5
h = b / num_points
for case in range(2):
    x = Vec2(-3,0.5)
    plt.scatter([x[0]], [x[1]])
    def plot_tangent():
        dx = 1
        dy = 2*x[0]
        t = 2
        line_xs = [x[0]-dx*t, x[0]+dx*t]
        line_ys = [x[1]-dy*t, x[1]+dy*t]
        plt.plot(line_xs, line_ys)
    plot_tangent()
    for height in np.linspace(-15, 15, 40):
        plot_linspace_func(-4, 4, 150, lambda x: x**2 + height)
    for i in range(num_points):
        if case == 0:
            x = Vec2(x[0] + h, x[1] + 2*h*(x[0] + h))
        if case == 1:
            x = x + h*u(x)
        plot_tangent()
        plt.scatter([x[0]], [x[1]])
    plt.show()
