
def prime_sieve(n):
    vals = [1 for i in range(n)]
    for i in range(2, floor(sqrt(n)+1)):
        j = 2
        while True:
            if j*i >= n:
                break
            vals[j*i] = 0
            j += 1
    return vals

