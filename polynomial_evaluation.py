# Trying to derive a factoring scheme giving a quartic multiplication with two multiplies.
# (e.g. maybe ((x+A)(x+B)+x+C)*(x+D) gives a general quartic.)
import sympy as sym
import itertools

# a1,b1, a2,b2 = sym.symbols("a_1 b_1 a_2 b_2")
# x = sym.symbols("x")
# alpha1,beta1, alpha2,beta2 = sym.symbols(r"\alpha_1 \beta_1 \alpha_2 \beta_2")
# A,B = sym.symbols(r"A B")
# 
# P = (x**2 + a1*x + b1)*(x**2 + a2*x + b2)
# P = P.expand()
# Q = ((x + A)*(x + B) + alpha1*x + beta1) * ((x + A)*(x + B) + alpha2*x + beta2)
# Q = Q.expand()
# 
# for i in range(4):
#     # lprint("${}$".format(i), P.coeff(x, i))
#     lprint("${}$".format(i), Q.coeff(x, i))
# 
# ldone()

# a,b,c,d = sym.symbols("a b c d")
# x = sym.symbols("x")
# P = ((x + a)*(x + b) + c)**2 + d
# P = P.expand()
# for i in range(4+1):
#     lprint("${}$".format(i), P.coeff(x, i))
# 
# 
# M = sym.Matrix([[0,0,0,1,0,0,1,0,0,0,2,0,0,1],
#                 [0,0,0,0,0,0,0,0,2,2,0,2,2,0],
#                 [0,0,2,0,1,1,0,4,0,0,0,0,0,0],
#                 [2,2,0,0,0,0,0,0,0,0,0,0,0,0]])
# 
# lprint("", M * M.transpose())
# # lprint("", M)
# # lprint("", M * M.transpose())
# S = M.transpose() * (M * M.transpose()).inv()
# lprint("", S)
# c0,c1,c2,c3 = sym.symbols("c_0 c_1 c_2 c_3")
# b = sym.Matrix([c0, c1, c2, c3])
# v = S*b
# lprint("", v)
# 
# lprint("", M*v)
# 
# 
# ldone()

a,b,c,d = sym.symbols("a b c d")
x = sym.symbols("x")
P = sym.Poly(((x + a)*(x + b) + c + x)**2 + d,  x)
M = []

monoms = []
for na,nb,nc,nd in itertools.product(range(0,2+1), repeat=4):
    monoms.append(a**na * b**nb * c**nc * d**nd)
    
for i in range(4+1):
    row = []
    coeff = P.coeff_monomial(x**i)
    p = sym.Poly(coeff, a,b,c,d)
    for na,nb,nc,nd in itertools.product(range(0,2+1), repeat=4):
        monom = a**na * b**nb * c**nc * d**nd
        row.append(p.coeff_monomial(monom))
    M.append(row)
M = sym.Matrix(M)
        
# lprint("", M * M.transpose())
# lprint("", M)
# lprint("", M * M.transpose())
S = M.transpose() * (M * M.transpose()).inv()
# lprint("", S)
c0,c1,c2,c3,c4 = sym.symbols("c_0 c_1 c_2 c_3 c_4")
B = sym.Matrix([c0, c1, c2, c3, c4])
v = S*B
# lprint("", v)
# lprint("", M*v)

A = v[monoms.index(a)]
B = v[monoms.index(b)]
C = v[monoms.index(c)]
D = v[monoms.index(d)]
print(A)
print(B)
print(C)
print(D)
print(P)
PP = P.subs({a:A, b:B, c:C, d:D}).as_expr().simplify()
print(PP)


# ldone()
