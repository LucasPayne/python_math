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
        

seq = [1]

for n in range(25):
    b = (-1)**n * sum([nCk(2*(n+1), 2*k) * (-1)**k * seq[k] for k in range(n+1)])
    seq.append(b)

for b in seq:
    print(b)
    

# binomial convolution of this sequence with sinx to get tanx coefficients.
tseq = []

for n in range(25):
    t = (-1)**n * sum([nCk(2*n+1, 2*k) * (-1)**k * seq[k] for k in range(n+1)])
    tseq.append(t)

print("tangent")
for t in tseq:
    print(t)
