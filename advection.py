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


a = -10
b = 10
intervals = 500
nodes = np.linspace(a, b, intervals)
h = (b - a) / intervals
stencil_radius = 2
print(finite_differences(list(range(-stencil_radius,stencil_radius+1)), 2))
weights = [x.subs(Sym('h'), h) for x in finite_differences(list(range(-stencil_radius,stencil_radius+1)), 2)]
plt.xlim(-10, 10)
plt.ylim(-5, 5)
print(weights)
input()
def f(x):
    # return exp(-x**2/2)*sin(x)
    # return sqrt(abs(x))*exp(-x**2/2) + 0.2*sin(5*x)
    return 3*exp(-x**2/2)
u = [f(x) if abs(x) < 4 else 0 for x in nodes]
# vel = [sin(x)+1 for x in nodes]
vel = [1 for x in nodes]


# copied off a sympy post ----
def diagonal_form(a, upper = 1, lower= 1):
    """
    a is a numpy square matrix
    this function converts a square matrix to diagonal ordered form
    returned matrix in ab shape which can be used directly for scipy.linalg.solve_banded
    """
    n = a.shape[1]
    assert(np.all(a.shape ==(n,n)))
    
    ab = np.zeros((2*n-1, n))
    
    for i in range(n):
        ab[i,(n-1)-i:] = np.diagonal(a,(n-1)-i)
        
    for i in range(n-1): 
        ab[(2*n-2)-i,:i+1] = np.diagonal(a,i-(n-1))

    mid_row_inx = int(ab.shape[0]/2)
    upper_rows = [mid_row_inx - i for i in range(1, upper+1)]
    upper_rows.reverse()
    upper_rows.append(mid_row_inx)
    lower_rows = [mid_row_inx + i for i in range(1, lower+1)]
    keep_rows = upper_rows+lower_rows
    ab = ab[keep_rows,:]
    return ab


y = [exp(-x**2 / 2) for x in np.linspace(-10, 10, 50)]
h = 20/50
for i in range(4, 50-4+1):
    print((y[i-1] - 2*y[i] + y[i+1])/(h**2), (-y[i-2] + 16*y[i-1] - 30*y[i] + 16*y[i+1] - y[i+2])/(12*h**2))
input()

plt.ion()
dt = 0.1
while True:
    # u = [0] + [u[i] + vel[i]*dt*(u[i+1]-u[i-1])/(2*h) for i in range(1, len(u)-1)] + [0]

    # D = np.diag([-2 for _ in range(intervals)])
    # D += np.diagflat([1 for _ in range(intervals-1)], 1)
    # D += np.diagflat([1 for _ in range(intervals-1)], -1)

    D = np.zeros((intervals, intervals))
    for weight_index,shift in enumerate(range(-stencil_radius, stencil_radius+1)):
        if shift < 0:
            np.fill_diagonal(D[-shift:,:], weights[weight_index])
        else:
            np.fill_diagonal(D[:,shift:], weights[weight_index])
    
    M = np.identity(intervals)
    M -= dt * D

    banded = False
    start_time = time.time()
    if banded:
        lower = stencil_radius
        upper = stencil_radius
        M_diagonal_form = diagonal_form(M, lower, upper)
        u = linalg.solve_banded((1, 1), M_diagonal_form, u)
    else:
        u = linalg.solve(M, u)
    end_time = time.time()

    print("{:.9f}".format(end_time - start_time))

    plt.clf()
    plt.xlim(-10, 10)
    plt.ylim(-3, 3)
    plt.plot(nodes, u)
    plt.pause(0.01)
