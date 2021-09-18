from matplotlib import pyplot as plt
from math import sqrt,pi,cos,sin
import numpy as np

def side_length(p, q):
    return sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2)

def draw_circle(c, radius):
    ps = [(c[0] + radius*cos(theta), c[1] + radius*sin(theta)) for theta in np.linspace(0,2*pi,100)]
    print(ps)
    plt.plot(*zip(*ps))


def draw_tri(p, q, r):
    plt.plot(*zip(p,q,r,p))
    a = np.linalg.norm(q - r)
    b = np.linalg.norm(r - p)
    c = np.linalg.norm(p - q)
    c = (a*p + b*q + c*r)/(a + b + c)

    plt.scatter([c[0]], [c[1]])

    l = np.dot(c - p, (q-p)/np.linalg.norm(q-p))
    d = sqrt(np.dot(c - p, c - p) - l**2)
    print(d)

    draw_circle(c, d)


draw_tri(np.array([0,0]), np.array([2,0]), np.array([3,2]))
plt.gca().set_aspect(1)
plt.show()



