

def sat(nums):
    n = 1
    table = [n for n in nums]
    print(nums)
    while n < len(nums):
        for i in range(len(nums)-1, n-1, -1):
            table[i] += table[i - n]
        print(table)
        n *= 2


sat(list(range(8)))

