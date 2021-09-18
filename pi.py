from math import factorial, sqrt

def term(n, x):
    return factorial(2*n)/((2**n * factorial(n))**2) * x**(2*n+1) / (2*n + 1)

t = 0
x = 0.5 * sqrt(2 - sqrt((5 + sqrt(5))/2))
for n in range(100):
    t += term(n, x)
    print(t*20)
