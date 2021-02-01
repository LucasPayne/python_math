from math import sqrt, floor

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

s = prime_sieve(100000000)
best_p = 0
best_c = 0
for i,v in enumerate(s):
    if v != 1:
        continue
    p = i
    c = 1
    while p > 2:
        print(p)
        p = (p - 1)//2
        if s[p] == 1:
            c += 1
        else:
            break
    if c > best_c:
        best_c = c
        best_p = i

print(best_c, best_p)
        
            
