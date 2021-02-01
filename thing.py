
for i in range(1, 27):
    j = i
    c = 0
    while j != 1:
        print(j)
        j *= i
        j %= 27
        c += 1
        if c > 100:
            break
    print("{}: {}", i, c)
    input()
        

