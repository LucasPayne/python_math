#...
import sympy as sym
from cmath import *


def pade_approximant(coefficients, numerator_degree):
    m = numerator_degree
    n = len(coefficients) - m

    x = sym.symbols('x')
    poly = sum([coefficients[i] * x**i for i in range(len(coefficients))])
    
    factor = poly
    divided = x**(m+n+1)

    while 1:
        q, r = sym.div(divided, factor, domain = 'QQ')
        print("q:", q)
        print("r:", r)
        if sym.degree(r) <= m:
            break
        divided = factor
        factor = r




pade_approximant([1,0,3,0,-6], 1)

x = sym.symbols('x')
p = 1 + (sym.Rational(252, 25) - sym.Rational(72, 5)*x**2)*(-x**2/6 - sym.Rational(1, 12))
print(p)
print(p.expand())

print((4/(25*p)).expand())

l1 = 1/10 + 0.5*sqrt(-17/75)
l2 = 1/10 - 0.5*sqrt(-17/75)
for n in range(10):
    print(sqrt(-75/17)/15 * (1/l1**(n+1) - 1/l2**(n+1)))


print(sqrt(l1.real*l1.real + l1.imag*l1.imag))

for n in range(10):
    print(1j*sqrt(2)*(1j/sqrt(2))**n * (1 - (-1)**n))


def cont(f, n=1):
    if n == 500:
        return 1
    return f(n) + 1/cont(f, n+1)

print(cont(lambda n: n))
print(1/cont(lambda n: n))
