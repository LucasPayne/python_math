

t1 = 1
t2 = 0.25
t3 = 0.25
for n in range(1, 15):
    t1 += 1 / (n+1)**3
    t2 += (1 - n**3/(n*(n+1)*(n+2)))/n**3
    t3 += 1 / n**3 - 1/(n*(n+1)*(n+2))
    print(t1, "   ", t2, "   ", t3)


