Rat = sym.Rational

A = sym.Matrix([[1, 0, Rat(1,3)], [0, Rat(1,3), Rat(1,5)], [Rat(1,3), Rat(1,5), 0]])
lprint("A:", A)
x = A.inv() * sym.Matrix([-Rat(1,5), 0, 0])
lprint("x:", x)

x = sym.symbols("x")
f = Rat(6,8) * sym.sqrt(10) * (x**2 - Rat(1,3))
lprint("f:", f)
g = sym.integrate(f**2, (x, -1, 1))
lprint("integral -1 to 1:", g)


y,r,R = sym.symbols("y, r, R")
f = ((x - R)**2 + y**2)*((x + R)**2 + y**2)
lprint("f:", f)
lprint("f:", sym.expand(f))



ldone()
