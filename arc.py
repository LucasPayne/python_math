import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, asin, acos, sin, cos, atan, tan
import scipy as sci
from scipy import integrate

x_nodes = np.linspace(-1,1,100)

def arc(x):
    return x*sqrt(1 - x**2) + 2*sci.integrate.quad(lambda x: sqrt(1-x**2), x, 1)[0]
plt.plot(x_nodes, [arc(x) for x in x_nodes])
plt.plot(x_nodes, [acos(x) for x in x_nodes])
plt.show()

def arc(x):
    return sci.integrate.quad(lambda x: 1/(1 + x**2), 0, x)[0]
plt.plot(x_nodes, [arc(x) for x in x_nodes])
plt.plot(x_nodes, [atan(x) for x in x_nodes])
plt.show()




