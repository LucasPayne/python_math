import sympy as sym
from sympy.combinatorics.permutations import Permutation
import numpy as np
import itertools
from math import factorial

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
    return Permutation(perm)

def signed_permutations(vals):
    N = len(vals)
    lehmer_code = [0 for i in range(N)]
    while True:
        lehmer_pos = N-2
        yield (-2*(sum(lehmer_code)%2)+1, lehmer_code_to_permutation(lehmer_code)(vals))
        lehmer_code[lehmer_pos] += 1
        while lehmer_code[lehmer_pos] > N-lehmer_pos-1:
            if lehmer_pos == 0:
                return
            lehmer_pos -= 1
            lehmer_code[lehmer_pos] += 1
            for i in range(lehmer_pos+1,N):
                lehmer_code[i] = 0
    


# A flat is represented by an antisymmetric tensor. This tensor has a small number of degrees of freedom,
# so if the dual flat (e.g. plane to point in R^3) has less entries, that tensor can be stored explicitly,
# while the dual tensor is contracted with implicitly, by taking a tensor index and computing the corresponding value
# with the correct sign.

# Note:
# In R^3, a line is represented by the 4x4 Plucker matrix, while its dual is also a line. So the above method
# does not save any space. However, there are 6 values non-repeated in the matrix, and four degrees of freedom.
# The four degrees of freedom might be harder to encode. What is a nice way to think of a related 6-vector in computations?

# Note:
# s = L.dot(k3) # Join line with plane
# This is simple since L is already in tensor form.
# Even if, for example, a plane k stores 4 values only, an interface would
# allows all usual tensor operations w/r/t the 4x4x4 antisymmetric tensor.
# Ideally this tensor interface would be facilitated by some feature of numpy etc.

def plane_coefficient(k1, i,j,k):
    # A plane in R^3 is represented by a 4x4x4 antisymmetric tensor.)
    # Given its dual point k1 = (a,b,c,d), this returns the entry of the corresponding tensor at (i,j,k).
    indices = (i,j,k)
    if len(set(indices)) < 3:
        return 0
    n = list(set([0,1,2,3]) - set(indices))[0]
    pi = Permutation([(index - n)%4 - 1 for index in indices])

    sign = 1 if pi.is_even else -1
    if n % 2 == 0:
        sign *= -1
    return sign * k1[n]
    
        
def plane_intersection(k1, k2):
    # k1: Plane represented by a 4-vector.
    # k2: Plane represented by a 4-vector.
    # Compute the regressive product k1 v k2 = k1 k2*.
    # This is the contraction k1^(ijk)k2_(i)\partial_(jk).
    # The result is a 4x4 antisymmetric rank-two Plucker matrix.
    plucker_matrix = np.zeros((4,4))
    for j,k in itertools.product(range(4), repeat=2):
        plucker_matrix[j,k] = sum(plane_coefficient(k1, i,j,k)*k2[i] for i in range(4))
    return plucker_matrix

def line_point_join(L, p):
    # L: 4x4 Plucker matrix.
    # p: 4-vector.
    # Compute the exterior product L ^ p.
    # The result is a plane represented by a 4x4x4 antisymmetric tensor.
    # This is returned as a 4-vector over the basis 123,234,341,412.
    plane = np.zeros(4)
    for n in range(4):
        t = 0
        for sign,perm in signed_permutations([(n+1+i)%4 for i in range(3)]):
            if n%2 == 0:
                sign *= -1 #---Sign?
            t += sign * L[perm[0],perm[1]]*p[perm[2]]
        plane[n] = t
    return plane


def point_point_join(p, q):
    # p: 4-vector.
    # q: 4-vector.
    # Compute the exterior product p ^ q.
    # This is a line represented by a 4x4 Plucker matrix.
    plucker_matrix = np.zeros((4,4))
    for i,j in itertools.product(range(4), repeat=2):
        plucker_matrix[i,j] = p[i]*q[j] - p[j]*q[i]
    return plucker_matrix
        
        
    
    

# k1 = np.array([1,2,3,4])
# k2 = np.array([2,3,4,5])
# L = plane_intersection(k1, k2)
# print(L)
# print(L.dot(k1))
# print(L.dot(k2))
# 
# k3 = np.array([2,3,0,-1])
# s = L.dot(k3) # Join line with plane. The regressive product here is simple since L is already in the tensor form.
# print(s)
# print(s.dot(k3))
# 
# print(line_point_join(L, s))
# print(L)
# plane = line_point_join(L, np.array([0,0,0,1]))
# print(plane)
# print(L.dot(plane))
# 
# 
# # perms = [perm for perm in alternating_permutations([0,1,2])]
# # print(perms)
# # print(len(perms))
# # print(len(set(perms)))
# 
# # l = plane_intersection(np.array([1,0,0,0]),np.array([0,1,0,0]))
# # print(l)
# # print(l.dot(np.array([0,0,1,-1]))) # z = 1 plane intersection
# # p = np.array([1,0,0,1])
# # plane = line_point_join(l, p)
# # print(plane)
# 

p = np.array([0,0,0,1])
q = np.array([1,0,0,0])
L = point_point_join(p, q)
print(L)
print(line_point_join(L, p))
print(line_point_join(L, q))
print(line_point_join(L, np.array([0,1,0,1])))
# 

# for sign,perm in signed_permutations([0,1,2,3]):
#     print(sign, perm)
