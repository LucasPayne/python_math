from math import pi, sqrt, cos, sin

def radius_from_subtend(degrees):
    theta = degrees * pi / 180;
    return sqrt((cos(theta)-1)**2 + sin(theta)**2)/sqrt((1+cos(theta))**2 + sin(theta)**2)
