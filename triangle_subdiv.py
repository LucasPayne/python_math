import sympy as sym
from printing import print_coeffs, print_matrix
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols
Half = Rat(1,2)
Third = Rat(1,3)
Quarter = Rat(1,4)
def Rec(n):
    return Rat(1, n)


def trinomial(n, i,j,k):
    assert(i+j+k == n)
    return sym.binomial(n, i) * sym.binomial(n - i, j)


# ordered_sums(3, 5) will give all triples of sums 0+0+5, 0+1+4, ..., 3+1+1, etc., that add to 5.
def ordered_sums(terms, n):
    if terms == 1:
        yield tuple([n])
        return
    for i in range(0, n+1):
        for trailing in ordered_sums(terms-1, n-i):
            yield tuple([i]) + trailing


def bezier_triangle_basis(degree):
    x,y,z = Sym("x y z")
    p = sym.poly((x + y + z)**degree)
    coeffs = []
    for (i,j,k) in ordered_sums(3, degree):
        coeffs.append(p.coeff_monomial(x**i * y**j * z**k) * x**i * y**j * z**k)
    return coeffs




degree = 2
coeffs = bezier_triangle_basis(degree)
x,y,z = Sym("x y z")



coeffs_q = [(Rat(1,2**degree) * trinomial(degree, i,j,k) * (x + z)**i * (x+y)**j * (x+z)**k).expand() for i,j,k in ordered_sums(3, degree)]

ps = Sym(" ".join("p_{}{}{}".format(i,j,k) for i,j,k in ordered_sums(3, degree)))
q_poly = sum(ps[i]*coeffs_q[i] for i in range(len(coeffs_q))).as_poly()

for i,j,k in ordered_sums(3, degree):
    c = sum(point * q_poly.coeff_monomial(x**i * y**j * z**k * point**1) for point in ps)
    print("q_{}{}{} = {}".format(i,j,k, c))

# for i in range(len(coeffs_q)):
#     print(q_poly.coeff(ps[i], 1))
# 
# for c in coeffs:
#     print(c)


# M = Mat([[1,2,1,-2,-2,1],
#          [-2,0,2,4,0,-2],
#          [1,-2,1,-2,2,1],
#          [2,0,-2,0,4,-2],
#          [-2,4,-2,0,0,2],
#          [1,-2,1,2,-2,1]])
# print_matrix(M.inv())

print_matrix(Mat([[Half, Half, 0], [0, Half, Half], [Half, 0, Half]]).inv())
