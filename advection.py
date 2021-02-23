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
def Rec(n):
    return Rat(1, n)



nodes = np.linspace(-10, 10, 200)
u = [exp(-x**2/2) if abs(x) < 4 else 0 for x in nodes]
plt.plot(nodes, u)

h = 0.5
for i in range(5):
    u = [u[i] - h*(u[i+1]-u[i]) for i in range(len(u)-1)] + [0]
    plt.plot(nodes, u)

plt.show()
