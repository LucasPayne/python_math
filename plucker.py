import sympy as sym
from sympy.combinatorics.permutations import Permutation
import numpy as np
import itertools

def alternating_permutations(lis):
    # Generate permutations with alternating sign (+even, -odd, +even, -odd, ...).
    N = len(lis)
    indices = list(range(N))
    for n in range(N):
    


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
    # A plane in R^3 is represented by a 4x4x4 antisymmetric tensor.
    # Given its dual point k1 = (a,b,c,d), this returns the entry of the corresponding tensor at (i,j,k).
    indices = (i,j,k)
    if (len(set(indices)) < 3):
        return 0
    n = list(set([0,1,2,3]) - set(indices))[0]
    pi = Permutation([(index - n)%4 - 1 for index in indices])

    sign = 1 if pi.is_even else -1
    if n % 2 == 0:
        sign *= -1
    return sign * k1[n]
    
        
def plane_intersection(k1, k2):
    # Contraction k1^(ijk)k2_(i)\partial_(jk).
    # This is one form of the regressive product k1 v k2.
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
        for i,j,k in itertools.permutations([(n+i)%4 for i in range(3)]):
            sign = 1
            t += sign * L[i,j]*p[k]
        plane[n] = t
    return plane
        
        
        
    
    

k1 = np.array([1,2,3,4])
k2 = np.array([2,3,4,5])
L = plane_intersection(k1, k2)
print(L)
print(L.dot(k1))
print(L.dot(k2))

k3 = np.array([2,3,0,-1])
s = L.dot(k3) # Join line with plane. The regressive product here is simple since L is already in the tensor form.
print(s)
print(s.dot(k3))

line_point_join(L, s)
