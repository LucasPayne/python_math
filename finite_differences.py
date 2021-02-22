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
from math import factorial
def Rec(n):
    return Rat(1, n)

def finite_differences(samples, order):
    assert(len(set(samples)) == len(samples))
    assert(order < len(samples))
    print("finite_differences {}".format(", ".join(str(x) for x in samples)))
    x = sym.symbols("x")
    h = sym.symbols("h")
    n = len(samples)
    power_matrix = []
    for i in range(n):
        power_matrix.append([(samples[i]*h)**j / factorial(j) for j in range(n)])
    power_matrix = sym.Matrix(power_matrix).T
    A = power_matrix.inv()
    weights = A.col(order).T
    return [weights[0, i] for i in range(weights.cols)]


intervals = 100
dt = 1 / intervals
h = Sym('h')
center = 2
samples = [-2, -1, 0, 1, 2]
weights = [w.subs(h, dt) for w in finite_differences(samples, 2)]
num_samples = len(samples)
system = np.zeros((intervals+1, intervals+1), dtype=float)

for shift in range(-1, 0):
    shifted_samples = [s - shift for s in samples]
    shifted_weights = finite_differences(shifted_samples, 2)
    for i in range(num_samples):
        system[shift+2, i] = shifted_weights[i].subs(h, dt)
for shift in range(1, 1+1):
    shifted_samples = [s - shift for s in samples]
    shifted_weights = finite_differences(shifted_samples, 2)
    for i in range(num_samples):
        system[intervals-2+shift, intervals+1-num_samples+i] = shifted_weights[i].subs(h, dt)
system[0,0] = 1
system[intervals,intervals] = 1

for i in range(2, intervals-1):
    for j in range(num_samples):
        system[i, j+i-2] = weights[j]

b = np.zeros(intervals+1, dtype=float)
b[0] = 5
b[intervals] = 10

x = linalg.solve(system, b)
plt.plot(np.linspace(0, 1, intervals+1), x)
plt.show()
