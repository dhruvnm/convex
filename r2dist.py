# Generates points in R2 in various distributions

import numpy as np
import pdb

def uniform_disk(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.sqrt(np.random.uniform(0, 1))
        x[i] = r * np.cos(a)
        y[i] = r * np.sin(a)
    return (x, y)

def uniform_circle(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        x[i] = radius * np.cos(a)
        y[i] = radius * np.sin(a)
    return (x, y)

def uniform_kgon(n, k):
    x = [None] * n
    y = [None] * n


    inner_angle = 2*np.pi / k
    isoc_angle = (np.pi - inner_angle) / 2
    isoc_side = np.sin(isoc_angle) / np.sin(inner_angle)

    # p1 = (0, 0)
    p2 = (isoc_side, 0)
    p3 = (isoc_side * np.cos(inner_angle), isoc_side * np.sin(inner_angle))

    for i in range(n):
        u1 = np.random.uniform(0, 1)
        p4 = (u1*p2[0], u1*p2[1])
        u2 = np.sqrt(np.random.uniform(0, 1))
        p5 = ((1 - u2)*p3[0] + u2*p4[0], (1 - u2)*p3[1] + u2*p4[1])

        r = np.sqrt(p5[0]**2 + p5[1]**2)
        theta = np.arctan(p5[1]/p5[0])

        u3 = np.floor(k*np.random.uniform(0, 1))
        theta += inner_angle * u3

        x[i] = r * np.cos(theta)
        y[i] = r * np.sin(theta)

    return (x, y)

def gaussian(n, sigma):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        x[i] = np.random.normal(0, sigma)
        y[i] = np.random.normal(0, sigma)
    return (x, y)
