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

a = -5
b = 5
# f = sin
f = lambda x: 1 if 0 <= x and x <= 1 else 0
g = lambda x: abs(x)/2
nodes = 100
dt = (b - a)/(nodes - 1)
plot_linspace_func(a, b, nodes, f)
plot_linspace_func(a, b, nodes, g)
plot_linspace_func(a, b, nodes, lambda x: sum(f(y)*g(x-y)*dt for y in np.linspace(a, b, nodes)))

plt.show()

ns = [0,0,0, 1, 3, 1, -2, 8, 4, 0,0,0]
def get_ns(i):
    if i < 0 or i >= len(ns):
        return 0
    return ns[i]
diffs = [get_ns(i-1) - 2*get_ns(i) + get_ns(i+1) for i in range(len(ns))]
print(diffs)
print(sum(diffs))

