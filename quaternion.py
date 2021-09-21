import sympy as sym
from sympy.algebras.quaternion import Quaternion

c0 = sym.sqrt(3)
q = Quaternion(2, c0, c0, c0)
p = (1,2,-3)
pp = Quaternion.rotate_point((1,2,-3), q)
print(q)
print(p)
print(pp)

print(q * Quaternion(0,1,2,-3))
