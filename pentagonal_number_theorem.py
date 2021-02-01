import itertools
from math import ceil, floor, sqrt

# Euler's pentagonal number theorem for a non-local recurrence of the number of partitions.
def partitions(num):
    # Pre-compute required pentagonal numbers.
    pentagonals = [(k*(3*k - 1))//2 for k in range(1, floor(0.5 + sqrt(0.25 + 6*num)))]
    seq = [1]
    for n in range(2, num+1):
        p = 0
        sgn = 1
        for k in itertools.count():
            pent = pentagonals[k]
            anti_pent = pent + k + 1
            if n-pent-1 >= 0:
                p += sgn * seq[n-pent-1] 
            if n-anti_pent-1 >= 0:
                p += sgn * seq[n-anti_pent-1]
            sgn *= -1
            if n-anti_pent-1 <= 0:
                break
        seq.append(p)
    return seq

for i in range(1, 7):
    print(partitions(i))
print(partitions(100))
