Rat = sym.Rational

h = sym.symbols("h")
A = sym.Matrix([[1,1,1], [0, h, -h], [0, h*h/2, h*h/2]])

lprint("A:", A)
lprint(":", A.inv()*sym.Matrix([0,0,1]))


A = sym.Matrix([[1,1,1,1,1], [0, h, -h, 2*h, -2*h],
                         [0, h*h/2, h*h/2, h*h, h*h],
                         [0, h*h*h/6, -h*h*h/6, 4*h*h*h/3, -4*h*h*h/3],
                         [0, h**4/24, h**4/24, 2*h**4/3, 2*h**4/3]])
lprint("A:", A)
lprint(":", A.inv()*sym.Matrix([0,0,0,0,1]))
lprint(":", A.inv()*sym.Matrix([0,0,1,0,0]))

A = sym.Matrix([[1, -h, h*h/2, -h*h*h/6], [1, 0,0,0], [1,h,h*h/2, h*h*h/6], [1, 1.5*h, (9/4) * h*h/2, (27/8) * h*h*h/2]])
lprint("A", A)
lprint("A.inv()", A.inv())

A = sym.Matrix([[1,-h,h*h/2],[1,0,0],[1,h,h*h/2]])
lprint("A", A)
lprint("A.inv()", A.inv())


ldone()
