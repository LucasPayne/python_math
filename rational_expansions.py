import sympy as sym
Rat = sym.Rational
Sym = sym.symbols

# # 1/(3 + 2x - x^2)
# 
# # Recurrence relation from linear system.
# A = [Rat(1, 3), -Rat(2, 9)]
# for i in range(2, 10):
#     A.append(Rat(1, 3)*A[i-2] - Rat(2, 3)*A[i-1])
# print(A)
# 
# # Closed form from partial fractions.
# for n in range(1, 10):
#     print(Rat(1, 4*3**(n+1)) + (-1)**n * Rat(1, 4))


def toeplitz(c, r):
    assert len(c)==len(r)
    assert(c[0] == r[0])
    n = len(c)
    M = [[0 for __ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n-i):
            M[j][i+j] = r[i]
    for i in range(1,n):
        for j in range(n-i):
            M[i+j][j] = c[i]
    return sym.Matrix(M)

def print_matrix(M):
    strings = [[str(sym.simplify(c)) for c in M.row(i)] for i in range(M.rows)]
    for col in range(M.cols):
        max_len = 0
        for row in range(M.rows):
            if len(strings[row][col]) > max_len:
                max_len = len(strings[row][col])
        for row in range(M.rows):
            l = len(strings[row][col])
            strings[row][col] += ",   " + " "*max(0, max_len-l)
    for row in strings:
        print("".join(row))


def f(n):
    M = toeplitz([1,-2,1]+[0]*n, [1,0,0]+[0]*n)
    y = sym.symbols("y")
    M -= y*sym.eye(3+n)
    print(sym.factor(M.det()))
for i in range(5):
    f(i)


A = [1, 4]
for i in range(2, 10):
    A.append(4*A[i-1] - 3*A[i-2])
print(A)


x = Sym('x')
print(sym.solve(x**2 - 2*x - 2))

l1 = 1 + sym.sqrt(3)
l2 = 1 - sym.sqrt(3)
M = sym.Matrix([[l1, l2], [-1, -1]])
print_matrix(M.inv() * sym.Matrix([1,1]))
a = M[0]
b = M[1]
def g(n):
    return 1/l1 * a * (1 / l1)**n + 1/l2 * b * (1 / l2)**n

for i in range(10):
    print(sym.simplify(g(i)))

for n in range(10):
    print(((n+1)*(n+2))//2 * 2**n)
