from math import gcd

def pythagorean_triple(x, y):
    a = 2*x*y
    b = abs(x*x - y*y)
    c = x*x + y*y
    g = gcd(a,gcd(b,c))
    a //= g
    b //= g
    c //= g
    # print("a:", a)
    # print("b:", b)
    # print("c:", c)
    print("a^2 + b^2 = {}^2 + {}^2 = {}".format(a, b, a*a + b*b))
    print("c^2 = {}^2 = {}".format(c, c*c))
    

for i in range(10):
    pythagorean_triple(i,10)


