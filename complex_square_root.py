import cmath
from math import sqrt

def complex_sqrt(a, b):
    modulus = sqrt(a**2 + b**2)
    sqrt_modulus = sqrt(modulus)
    x = a / modulus
    b_sgn = 1 if b >= 0 else -1
    return (sqrt_modulus * sqrt((1 + x)/2),
            sqrt_modulus * b_sgn * sqrt((1 - x)/2))

def complex_cube_root(a, b):
    cube_root_modulus = pow(a**2 + b**2, 1/6)
    y = 

    return (cube_root_modulus * y,
            cube_root_modulus * sqrt(1 - y**2))

def test(a, b):
    print(cmath.sqrt(a + b*1j))
    print("   ", complex_sqrt(a, b))

test(1, 2)
test(1, 3)
test(0, -1)
test(2, -3)
test(-3, -3)
test(-2, 3)
