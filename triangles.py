from matplotlib import pyplot as plt
from math import sqrt,pi,cos,sin
import numpy as np
import sympy as sym


def cofactor_expansion(M):
    n = sym.shape(M)[0]
    if n == 1:
        return M[0,0]
    total = 0
    for i in range(n):
        minor = M.copy()
        minor.row_del(0)
        minor.col_del(i)
        total += (-1)**i * M[0,i] * cofactor_expansion(minor)
    return total


def side_length(p, q):
    return sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2)

def draw_circle(c, radius):
    ps = [(c[0] + radius*cos(theta), c[1] + radius*sin(theta)) for theta in np.linspace(0,2*pi,100)]
    plt.plot(*zip(*ps))


def circumcircle_center(p, q, r):
    x,y = sym.symbols("x y")
    var_p = sym.symbols("p_x p_y")
    var_q = sym.symbols("q_x q_y")
    var_r = sym.symbols("r_x r_y")
    M = sym.Matrix([
        *[[point[0], point[1], np.dot(point, point), 1] for point in (p,q,r)],
        # *[[point[0], point[1], np.dot(point, point), 1] for point in (var_p,var_q,var_r)],
        [x, y, x**2 + y**2, 1]
    ])
    circle_equation = M.det()
    grad = [sym.diff(circle_equation, var) for var in (x,y)]
    sol = sym.solve(grad)
    return np.array([sol[x], sol[y]])

def circumcircle_linear_equation():
    x,y = sym.symbols("x y")
    var_p = sym.symbols("p_x p_y")
    var_q = sym.symbols("q_x q_y")
    var_r = sym.symbols("r_x r_y")
    M = sym.Matrix([
        *[[point[0], point[1], 1, np.dot(point, point)][::-1] for point in (var_p,var_q,var_r)],
        [x, y, 1, x**2 + y**2][::-1]
    ]).T
    poly = cofactor_expansion(M)
    print(poly)

    # circle_equation = M.det()
    # grad = [sym.diff(circle_equation, var) for var in (x,y)]
    # sol = sym.solve(grad)
    # return np.array([sol[x], sol[y]])


def circumcircle_center_vector_algebra(p, q, r):
    # This formula is derived by
    #     - Forming the determinant for the circumcircle equation.
    #     - Writing this expanded determinant in terms of vector algebra (dot, cross, triple product).
    #     - Differentiating with respect to x and y.
    #     - Equating this gradient to zero, then solving the linear equations.
    p_h = np.array([*p, 1])
    q_h = np.array([*q, 1])
    r_h = np.array([*r, 1])
    s = np.dot(p, p)*np.cross(q_h, r_h) \
        + np.dot(q, q)*np.cross(r_h, p_h) \
        + np.dot(r, r)*np.cross(p_h, q_h)
    s = s/(2*(np.dot(p_h, np.cross(q_h, r_h))))
    return s[:2]
    



def draw_tri(p, q, r):
    plt.plot(*zip(p,q,r,p))
    a = np.linalg.norm(q - r)
    b = np.linalg.norm(r - p)
    c = np.linalg.norm(p - q)
    
    # Incircle
    center = (a*p + b*q + c*r)/(a + b + c)
    plt.scatter([center[0]], [center[1]])
    l = np.dot(center - p, (q-p)/np.linalg.norm(q-p))
    d = sqrt(np.dot(center - p, center - p) - l**2)
    draw_circle(center, d)

    r1 = np.array([p[1] - q[1], q[0] - p[0]])
    h1 = np.array([*((p+q)/2), 1])

    r2 = np.array([p[1] - r[1], r[0] - p[0]])
    h2 = np.array([*((p+r)/2), 1])

    # Circumcircle
    point = np.cross(np.cross(r1, h1), np.cross(r2, h2))
    point = point[:2]/point[2]
    plt.scatter(*point)
    radius = np.linalg.norm(point - p)
    draw_circle(point, radius)

    point_2 = circumcircle_center(p, q, r)
    print(point)
    print(point_2)

    circumcircle_linear_equation()
    
    point_3 = circumcircle_center_vector_algebra(p,q,r)
    print(point_3)
    

    




draw_tri(np.array([0,0]), np.array([2,0]), np.array([3,2]))
plt.gca().set_aspect(1)
plt.show()


