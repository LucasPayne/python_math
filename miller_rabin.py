from random import randrange

def miller_rabin(n):
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    a = randrange(1, n)
    c = a**d % n
    if c != 1 and c != n - 1:
        return False
    for r in range(1, n):
        if a**(2**r * d) % n != n - 1:
            return False
    return True

for i in range(2,100):
    print(i, miller_rabin(i))
            

