import numpy as np
from scipy import linalg

# Exponential of a skew-symmetric matrix is orthogonal.

def matrix_exp(M):
    Mpow = np.eye(M.shape[0])
    Mexp = Mpow
    for i in range(20):
        Mpow = Mpow.dot((1/(i+1)) * M)
        Mexp += Mpow
    return Mexp
        
        
S = np.array([[0, -1, 2],
              [1, 0, 1.22],
              [-2, -1.22, 0]])
print(S)
M = matrix_exp(S)
print(M)
print(M.T.dot(M))
M = linalg.expm(S)
print(M)
print(M.T.dot(M))
