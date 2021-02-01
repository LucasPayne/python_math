from math import factorial

def nCk(n, k):
    if k == 0:
        return 1
    numerator = n
    x = n
    for i in range(k-1):
        x -= 1
        numerator *= x
    return numerator / factorial(k)
        


# recursive version
def S(n, a):
    if a == 0:
        return n
    return -1/(a+1) * (sum([nCk(a+1, k) * S(n, k) for k in range(a)]) + 1 - (n+1)**(a+1))


# iterative version (can handle larger powers)
def Si(n, a):
    table = [0 for _ in range(a+1)]
    table[0] = n
    for alpha in range(2, a+1):
        s = -sum([nCk(alpha + 1, k) * table[k] for k in range(alpha)]) + (n+1)**(alpha + 1) - 1
        table[alpha] = (1 / (alpha + 1)) * s
    return table[a]


def Stest(n, a):
    t = 0
    for i in range(1, n+1):
        t += i**a
    print("{}, {}, {}".format(t, S(n, a), Si(n, a)))


import sympy as sym

m = sym.Matrix([[1,0,0,0,0],
                [1,2,0,0,0],
                [1,3,3,0,0],
                [1,4,6,4,0],
                [1,5,10,10,5]])
n = sym.symbols('n')
v = sym.Matrix([(n+1)**(i+1) - 1 for i in range(5)])
print(sym.Matrix([e.expand() for e in m.inv() * v]))



