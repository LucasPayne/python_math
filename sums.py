from math import factorial

def harmonic(n):
    x = 0
    for i in range(1,n+1):
        x += 1/i
    return x

def T_closed(n):
    s = 1
    for k in range(1,n+1):
        s += 1/factorial(k-1)
    return s * (factorial(n)**2)/factorial(n+1)

def T(n):
    if n==0:
        return 1
    return (n**2 / (n+1)) * T(n-1) + n/(n+1)

def alpha(n):
    x = 1
    for k in range(1,n+1):
        x *= (k+1)/k**2
    print(x)
    print(factorial(n+1)/(factorial(n)**2))
    return x



def sum_squares(n):
    t = 0
    for k in range(1, n+1):
        t += k*k
    return t

def sum_squares_2(n):
    t = 0
    for k in range(1, n+1):
        t += (n-k+1)*(2*k-1)
    return t

def sum_squares_3(n):
    return ((n+1)**3 - 1 - n - ((3*n*(n+1))/2))/3


def sum_powers(power, n):
    if power == 0:
        return n
    if power == 1:
        return n*(n+1)/2
    s = 0
    for k in range(power):
        c = factorial(power+1) / (factorial(k) * factorial(power+1-k))
        s += c * sum_powers(k, n)
    return ((n+1)**(power+1) - 1 - s)/(power + 1)

def sum_cubes(n):
    t = 0
    for i in range(1, n+1):
        t += i ** 3
    print(t)
    return ((n+1)**4 - 1 - 6*sum_squares_3(n) - 4*(n*(n+1))/2 - n)/4

def S1(n, x):
    # Sum k=0..n of (k+1)x^3k
    t = 0
    for k in range(n+1):
        t += (k+1)*x**(3*k)
    print(t)
    # Rational closed form derived from perturbation method.
    return (1 + x**3 * ((1-x**(3*(n+1)))/(1-x**3)) - (n+1)*x**(3*(n+1))) / (1 - x**3)

def S2(n):
    # Sum k=0..n of k2^k
    t = 0
    for k in range(n+1):
        t += k * 2**k
    print(t)
    # Closed form from perturbation method.
    return (n+1)*2**(n+1) + 2 - 2**(n+2)
    


for i in range(2,11):
    print("{}   :   {}".format(T(i), T_closed(i)))
