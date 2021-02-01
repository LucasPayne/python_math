import random
from math import sqrt
min_dis = 0.2
ns = []
counter = 0
while counter < 12:
    (a, b) = (random.uniform(-1,1), random.uniform(-1,1))
    if a**2 + b**2 > 1:
        continue
    bad = False
    for (c,d) in ns:
        if sqrt((c-a)**2 + (d-b)**2) < min_dis:
            bad = True
    if bad:
        continue
    print("    vec2({}, {}),".format(a, b))
    counter += 1
    ns += [(a,b)]

