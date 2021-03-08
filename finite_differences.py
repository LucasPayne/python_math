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
from math import factorial, sqrt, sin, cos, asin, acos
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
    print(sym.latex(power_matrix.T.subs(h, 1)))
    print(sym.latex(A.T.subs(h, 1)))
    weights = A.col(order).T
    return [weights[0, i] for i in range(weights.cols)]


# intervals = 100
# dt = 1 / intervals
# h = Sym('h')
# center = 2
# samples = [-2, -1, 0, 1, 2]
# weights = [w.subs(h, dt) for w in finite_differences(samples, 2)]
# num_samples = len(samples)
# system = np.zeros((intervals+1, intervals+1), dtype=float)
# 
# for shift in range(-1, 0):
#     shifted_samples = [s - shift for s in samples]
#     shifted_weights = finite_differences(shifted_samples, 2)
#     for i in range(num_samples):
#         system[shift+2, i] = shifted_weights[i].subs(h, dt)
# for shift in range(1, 1+1):
#     shifted_samples = [s - shift for s in samples]
#     shifted_weights = finite_differences(shifted_samples, 2)
#     for i in range(num_samples):
#         system[intervals-2+shift, intervals+1-num_samples+i] = shifted_weights[i].subs(h, dt)
# system[0,0] = 1
# system[intervals,intervals] = 1
# 
# for i in range(2, intervals-1):
#     for j in range(num_samples):
#         system[i, j+i-2] = weights[j]
# 
# b = np.zeros(intervals+1, dtype=float)
# b[0] = 5
# b[intervals] = 10
# 
# x = linalg.solve(system, b)
# plt.plot(np.linspace(0, 1, intervals+1), x)
# plt.show()

# print(finite_differences([-1, 0, 1], 2))
# print(finite_differences([-2, -1, 0, 1, 2], 2))
# print(finite_differences([-3, -2, -1, 0, 1, 2, 3], 2))

def plot_linspace_func(start, end, nodes, func):
    plt.plot(np.linspace(start, end, nodes), [func(x) for x in np.linspace(start, end, nodes)])


def solve_linear_dirichlet_bvp(coefficients, boundary, boundary_values, intervals, order_of_accuracy):
    datatype = float
    assert(len(boundary) == 2 and len(boundary_values) == 2 and boundary[0] < boundary[1])
    n = len(coefficients)
    # Central finite differences.
    # todo: Choose simplest central differences for the given accuracy.

    system = np.zeros((intervals+1,intervals+1), dtype=datatype)
    system[0,0] = 1
    system[intervals,intervals] = 1

    dt = 1/intervals
    b = np.zeros(intervals+1, dtype=datatype)
    b[0] = boundary_values[0]
    b[intervals] = boundary_values[1]
    for i in range(1, intervals):
        b[i] = -coefficients[0](boundary[0] + (boundary[1]-boundary[0])*i*dt)
    h = sym.symbols('h')

    for derivative_degree in range(1, n): # skip 0'th order term as it will be in b.
        samples = list(range(-1, 1+1)) # only up to order 2
        num_samples = len(samples)
        stencil = [x.subs(h, dt) for x in finite_differences(samples, derivative_degree)]
        stencil_radius = (len(stencil)-1)//2
        M = np.zeros((intervals+1, intervals+1), dtype=datatype)

        # Left boundary stencils.
        for index,shift in enumerate(range(-stencil_radius+1, 0)):
            shifted_samples = [s - shift for s in samples]
            shifted_weights = finite_differences(shifted_samples, derivative_degree)
            for i in range(num_samples):
                M[index+1, i] = shifted_weights[i].subs(h, dt)
        # Right boundary stencils.
        for index,shift in enumerate(range(1, stencil_radius)):
            shifted_samples = [s - shift for s in samples]
            shifted_weights = finite_differences(shifted_samples, derivative_degree)
            for i in range(num_samples):
                M[intervals-index-1, intervals+1-num_samples+i] = shifted_weights[i].subs(h, dt)
        # Middle stencils.
        for i in range(max(1, stencil_radius), min(intervals, intervals-stencil_radius+1)):
            for j in range(num_samples):
                M[i, j+i-stencil_radius] = stencil[j]
        # Coefficients are functions of x, so discretize this function and use it to weight the rows.
        row_multipliers = np.array([coefficients[derivative_degree](boundary[0] + (boundary[1]-boundary[0])*i*dt) for i in range(intervals+1)], dtype=datatype)
        M = np.multiply(M, row_multipliers[:, np.newaxis])
        print_matrix(Mat(M))
        # Add this matrix to the total system.
        system = np.add(system, M)
    print_matrix(Mat(system))
    print_matrix(Mat(b))

    x = linalg.solve(system, b)
    plt.plot(np.linspace(boundary[0], boundary[1], intervals+1), x, "k")

def main():
    # solve_linear_dirichlet_bvp([lambda x: 10 if x > 3/10 and x < 7/10 else 0, lambda x: 0, lambda x: 1], (0,1), (0, 0), 100, 2)
    # solve_linear_dirichlet_bvp([lambda x: 10 if x > 3/10 and x < 7/10 else 0, lambda x: 0, lambda x: 1], (0,1), (0, 0), 10, 2)
    
    # solve_linear_dirichlet_bvp([lambda x: -sin(x), lambda x: 1, lambda x: 1], (0,1), (1, 0), 10, 2)

    solve_linear_dirichlet_bvp([lambda x: 0, lambda x: 1, lambda x: 1], (0, 1), (2, 3), 100, 2)
    exact = lambda x: sin(x)/sin(1) + 2
    plt.plot(np.linspace(0, 1, 101), [exact(x) for x in np.linspace(0, 1, 101)], "b--")

    u = sym.symbols("u", cls=sym.Function)
    x = sym.symbols("x")
    diffeq = sym.Eq(u(x).diff(x, x) + u(x), 0)
    sol = sym.dsolve(diffeq)
    print(sol)

    plt.legend(["Solution with second order finite differences", "Exact solution"], loc="upper left")
    plt.title("")
    plt.xlabel("$x$")
    # plt.ylabel("$y$")
    plt.savefig("3a.pdf")
    
    
    # solve_linear_dirichlet_bvp([lambda x: 10*x, lambda x: 0, lambda x: 1], (0,1), (0, 1), 10, 4)
    
    # solve_linear_dirichlet_bvp([lambda x: x, lambda x: 0, lambda x: 1], (0,1), (0, 0), 100, 4)
    # plot_linspace_func(boundary[0], boundary[1], intervals+1, lambda x: (1/6)*(x - x**3))
    plt.show()

if __name__ == "__main__":
    main()
