from math import sin,cos
import numpy as np


def rotation_matrix(theta):
    return np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])

def line_transform(z1, z2, M):
    N = 50
    xs = np.linspace(z1[0], z2[0], N)
    ys = np.linspace(z1[1], z2[1], N)
    points = [np.array(p) for p in list(zip(xs, ys))]
    transformed_points = M.dot(p) for p in 


line_transform((0,0), (1,1))
