Rat = sym.Rational

# 
# A = sym.Matrix([[0.5, 0.5, 0], [-2, 2, 0], [2, -4, 2]])
# B = sym.Matrix([[0.5, 0.5, 0], [-1, 1, 0], [0.5, -1, 0.5]])
# 
# lprint("ok ", A.inv() * B)
# 
# 
# t = sym.symbols('t')
# n13 = sym.expand((0.5 + t - t**2)*((2 - t)/3) + 0.5*((1-t)**2) * ((2+t)/3))
# lprint("n13: ", n13)
# n23 = sym.expand(2/3 - (t - 1)**2 - (t-1)**3 / 2)
# lprint("n23: ", n23)
# 
# 
# A3 = sym.Matrix([[1/6, 2/3, 1/6, 0], [-1/2, 0, 1/2, 0], [1/2, -1, 1/2, 0], [-1/6, 1/2, -1/2, 1/6]])
# B3 = sym.diag(1,2,4,8) * A3 
# 
# lprint("A3 ", A3)
# lprint("B3 ", B3)
# lprint("$A3^-1 * B3$", A3.inv() * B3)

# A = sym.Matrix([[0.5, -1, 0.5], [0.5, 1, -1], [0, 0, 0.5]])
# lprint("ok ", A * sym.diag(1, 0.5, 0.25) * A.inv())

A = Rat(1,6) * sym.Matrix([[1, -3, 3, -1], [4, 0, -6, 3], [1, 3, 3, -3], [0,0,0,1]])
lprint("A: ", A)
M = sym.diag(1, 2, 4, 8).inv()
K = A * M * A.inv()
lprint("K: ", K)

G = []
for i in range(4):
    row = []
    for j in range(4):
        row += [sym.symbols("P_{}{}".format(i, j))]
    G += [row]
G = sym.Matrix(G)
lprint("G: ", G)
Q = K.T * G * K
lprint("$K^T G K$:", Q)
# lprint("$K^T G K$:", (A * M * A.inv()).T * G * (A * M * A.inv()))


A = sym.Matrix([[Rat(1,2), -1, Rat(1,2)], [Rat(1,2), 1, -1], [0, 0, Rat(1/2)]])
lprint("A: ", A)
M = sym.diag(1, 2, 4).inv()
K = A * M * A.inv()
lprint("K: ", K)

G = []
for i in range(3):
    row = []
    for j in range(3):
        row += [sym.symbols("P_{}{}".format(i, j))]
    G += [row]
G = sym.Matrix(G)
lprint("G: ", G)
Q = K.T * G * K
lprint("$K^T G K$:", Q)


ldone()
