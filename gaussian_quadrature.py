import sympy as sym

a = -1/sym.sqrt(3)
b = 1/sym.sqrt(3)

x, fa, fb = sym.symbols('x fa fb')

f = fa * (x - b)/(a -b) + fb * (x - a)/(b - a)
print(sym.integrate(f, (x, -1, 1)))

