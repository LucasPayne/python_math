from matplotlib import pyplot as plt
from math import sqrt,pi,cos,sin
import numpy as np

def side_length(p, q):
    return sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2)

def draw_circle(c, radius):
    ps = [(c[0] + radius*cos(theta), c[1] + radius*sin(theta)) for theta in np.linspace(0,2*pi,100)]
    plt.plot(*zip(*ps))


def draw_tri(p, q, r):
    plt.plot(*zip(p,q,r,p))
    a = np.linalg.norm(q - r)
    b = np.linalg.norm(r - p)
    c = np.linalg.norm(p - q)

    center = (a*p + b*q + c*r)/(a + b + c)
    plt.scatter([center[0]], [center[1]])
    l = np.dot(center - p, (q-p)/np.linalg.norm(q-p))
    d = sqrt(np.dot(center - p, center - p) - l**2)
    draw_circle(center, d)

    r1 = np.array([p[1] - q[1], q[0] - p[0]])
    h1 = np.array([*((p+q)/2), 1])

    r2 = np.array([p[1] - r[1], r[0] - p[0]])
    h2 = np.array([*((p+r)/2), 1])

    point = np.cross(np.cross(r1, h1), np.cross(r2, h2))
    point = point[:2]/point[2]
    plt.scatter(*point)
    radius = np.linalg.norm(point - p)
    draw_circle(point, radius)



draw_tri(np.array([0,0]), np.array([2,0]), np.array([3,2]))
plt.gca().set_aspect(1)
plt.show()



