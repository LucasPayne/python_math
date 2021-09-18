# import sympy as sym

x, theta, lamb = sym.symbols(r"x \theta \lambda")
y = sym.Function("y")(x)
yp = sym.diff(y, x)

J = sym.Matrix([[1, yp*sym.cos(theta), yp*sym.sin(theta)],
                [0, -y*sym.sin(theta), y*sym.cos(theta)]])
M = J * J.T
# print(sym.sqrt(M.det()).simplify())

L = y*sym.sqrt(1+yp**2)

el = sym.diff(L, y) - sym.diff(sym.diff(L, yp), x)
el = el.simplify()
con = -sym.diff(sym.sqrt(1 + yp**2), x)

lprint("", el.subs(y, sym.cosh(x)).doit().simplify())
lprint("", con.subs(y, sym.cosh(x)).doit().simplify())
lprint("", (el + lamb*con).subs(y, sym.cosh(x)).doit().simplify())


lprint("", el)
lprint("", con)
ldone()
