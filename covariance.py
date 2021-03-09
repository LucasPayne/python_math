import sympy as sym
import itertools
from matplotlib import pyplot as plt
import numpy as np
from scipy import linalg
from printing import print_coeffs, print_matrix
import string
Rat = sym.Rational
Mat = sym.Matrix
Sym = sym.symbols
Half = Rat(1,2)
Third = Rat(1,3)
Quarter = Rat(1,4)
from math import factorial, sqrt, sin, cos, exp, floor, ceil
import time
def Rec(n):
    return Rat(1, n)
from finite_differences import finite_differences
import random

from mpl_toolkits.mplot3d import Axes3D


n = 100
scales = [random.uniform(0, 1) for _ in range(3)]
scale_dir = np.random.rand(3)
scale_dir *= 1/sqrt(sum(x**2 for x in scale_dir))
scale = random.uniform(0,0.5)
points = np.array([[scales[i] * np.random.normal(0, 1) for i in range(3)] for __ in range(n)])
for row in range(points.shape[0]):
    points[row,:] = points[row,:] + (scale-1) * scale_dir * points[row,:].dot(scale_dir)

# print_matrix(sym.Matrix(points))
mean = 1/n * sum(points[row,:] for row in range(points.shape[0]))
# C = 1/n * sum(np.outer(points[row,:] - mean, points[row,:] - mean) for row in range(points.shape[0]))
centralized_points = np.array([points[row,:] - mean for row in range(points.shape[0])])
C = (1/n) * centralized_points.T.dot(centralized_points) # Easier way to write the above.
print_matrix(Mat(C))
print("det(C) =", linalg.det(C))


U,s,Vh = linalg.svd(C)
# Compute principal plane.
C2 = sum(s[i] * np.outer(U[:,i], Vh[i,:].T) for i in range(2))
print_matrix(Mat(C2))
print("det(C2) =", linalg.det(C2)) # Should be 0.



fig = plt.figure()
ax = plt.gca(projection="3d")
X = U[:,0]
Y = U[:,1]
plane_points = np.array([x*X + y*Y + mean for x,y in np.mgrid[-3:3:25j, -3:3:25j].reshape(2,-1).T])

# ax.plot_surface(xx, yy, z, alpha=0.3)
ax.scatter(points[:,0], points[:,1], points[:,2])
ax.scatter(plane_points[:,0], plane_points[:,1], plane_points[:,2], alpha=0.2)
plt.show()
