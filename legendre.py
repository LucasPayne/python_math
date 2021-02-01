# Least squares polynomial regression using Legendre polynomials.
# Inner products are done using two-point Gaussian quadrature.
#
# Possibly higher degree Gaussian quadrature should be used when fitting higher degree polynomials,
# so that the sampled polynomial vectors are numerically orthogonal.

import sympy as sym
import numpy as np
import itertools
import random
from matplotlib import pyplot as plt

degree = 8

x = sym.symbols('x')
# Normalized Legendre polynomials (integral from -1 to 1 of P[n]^2 is 1)
P = [1 / sym.sqrt(2)]
for i in range(1, degree+1):
    f = x**i
    for p in P:
        f -= p * sym.integrate(f * p, (x, -1, 1))
    norm = sym.sqrt(sym.integrate(f * f, (x, -1, 1)))
    f /= norm
    P.append(f)
    print(f.expand())

# P = [
#     1 / sym.sqrt(2),
#     (sym.sqrt(3) / sym.sqrt(2)) * x,
#     (5/(2*sym.sqrt(10)))*(3*x**2 - 1),
#     sym.sqrt(sym.Rational(7, 2)) * sym.Rational(1, 2) * (5*x**3 - 3*x)
# ]


from math import sin, cos, sqrt, exp
def f(x):
    random.seed(x)
    return sin(5*x**2) + 0.1*random.gauss(0, 1)

dx = 0.01
xs = []
for s in np.arange(-1, 1, dx):
    xs += [s + dx*(0.5 - 1/(2*sqrt(3))),
           s + dx*(0.5 + 1/(2*sqrt(3)))]

Pd = [[float(P[n].subs(x, xs[i])) for i in range(len(xs))] for n in range(len(P))]
for a,b in itertools.product(Pd, repeat=2):
    t = 0
    for i in range(len(a)):
        t += a[i] * b[i] * dx * 0.5
    print(t)


# Project f to normalized Legendre polynomials by inner products.
coeffs = []
for n in range(len(P)):
    integral = 0
    for i in range(len(xs)):
        integral += Pd[n][i] * f(xs[i]) * dx * 0.5
    coeffs.append(integral)
print(coeffs)

N = 200
plot_xs = [-1 + 2*i/N for i in range(N)]
plt.plot(plot_xs, [sum(float(coeffs[n]*P[n].subs(x, plot_xs[i])) for n in range(len(coeffs))) for i in range(N)])
plt.plot(plot_xs, [f(plot_xs[i]) for i in range(N)])
plt.show()
