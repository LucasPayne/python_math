import sympy as sym

def depress_cubic():
    a,b,c,d = sym.symbols("a b c d")
    x = sym.symbols("x")
    P = a*x**3 + b*x**2 + c*x + d
    Pdepressed = P.subs(x, x - b/(3*a)).as_poly(x)
    print(P)
    print(Pdepressed)

def discriminant():
    p,q = sym.symbols("p q")
    x = sym.symbols("x")
    R = x**3 + p*x + q
    v = sym.sqrt(-p/3)
    expr = p**2 - (R.subs(x, v) - p)**2
    print(expr.simplify())

print("depress_cubic")
depress_cubic()
print("discriminant")
discriminant()
