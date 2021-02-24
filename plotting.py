from matplotlib import pyplot as plt
import numpy as np
import sympy as sym

def plot_linspace_func(start, end, nodes, func):
    plt.plot(np.linspace(start, end, nodes), [func(x) for x in np.linspace(start, end, nodes)])
