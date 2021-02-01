import sympy as sym

def chebyshev_polynomial(n):
    x = sym.symbols('x')
    return (((x + sym.I*sym.sqrt(1-x**2))**n + (x - sym.I*sym.sqrt(1-x**2))**n)/2).expand()

for n in range(10):
    print(chebyshev_polynomial(n))
