import sympy as sym
from printing import print_matrix
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols


m = Mat([[Rat(1,2), Rat(1,2), 0],
         [Rat(1,6), Rat(2,3), Rat(1,6)],
         [0, Rat(1,2), Rat(1,2)]])

[P,D] = m.diagonalize()
print_matrix(P)
print_matrix(D)


K = P * Mat([[0,0,0],[0,0,0],[0,0,1]]) * P.inv()
print_matrix(K)

C = [1, 5, 9]
c = [1, 5, 9]
for i in range(100):
    c = [(c[0]+c[1])/2, c[0]/6 + 2*c[1]/3 + c[2]/6, (c[1]+c[2])/2]
print(c[1])
print(C[0]/5 + 3*C[1]/5 + C[2]/5)





