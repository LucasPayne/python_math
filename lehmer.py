import itertools


def permutation_to_lehmer_code(perm):
    n = len(perm)
    lehm = [v for v in perm]
    for i in range(n-1):
        for j in range(i+1,n):
            if lehm[j] > lehm[i]:
                lehm[j] -= 1
    return lehm
    

def lehmer_code_to_permutation(lehm):
    n = len(lehm)
    perm = [0 for i in range(n)]
    select = [i for i in range(n)]
    for i in range(n):
        perm[i] = select[lehm[i]]
        for j in range(lehm[i],n):
            select[j] = select[(j+1)%n]
    return perm
        
        

def test(perm):
    lehm = permutation_to_lehmer_code(perm)
    perm2 = lehmer_code_to_permutation(lehm)
    print(perm, lehm, perm2)

test([0,1,2,3])
test([1,0,2,3]) # [1,0,0,0]
test([2,3,0,1]) # [2,2,0,0]




