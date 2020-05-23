# various utility functions

import numpy as np
import matplotlib.pyplot as plt

def orient(p1, p2, p3):
    sign, _ =  np.linalg.slogdet([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])
    return sign

def psuedo_angle(p1, p2):
    if p1 == p2:
        return (0,0)

    length = ((p2[0] - p1[0])**2 + (p2[1] - p1[1]) ** 2) ** (1/2)
    height = (p2[1] - p1[1]) / length

    if p2[0] >= p1[0]:
        return (height, length)
    else:
        return (2 - height, length)

def draw_hull(hull):
    h = len(hull)
    for i in range(-1, h-1):
        plt.plot([hull[i][0], hull[i+1][0]], [hull[i][1], hull[i+1][1]], 'r')

def psuedo_distance(a, b, c):
    return abs(np.linalg.det([[a[0], b[0], c[0]], [a[1], b[1], c[1]], [1, 1, 1]]))
