# various utility functions

import numpy as np
import matplotlib.pyplot as plt

def orient(p1, p2, p3):
    sign, _ =  np.linalg.slogdet([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])
    return sign

def psuedo_angle(point):
    point = point / np.sqrt(point[0]**2 + point[1]**2)
    if point[0] >= 0:
        return point[1]
    else:
        return 2 - point[1]

def draw_hull(hull):
    h = len(hull)
    for i in range(-1, h-1):
        plt.plot([hull[i][0], hull[i+1][0]], [hull[i][1], hull[i+1][1]], 'r')
