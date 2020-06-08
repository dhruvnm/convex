# Generates points in R2 in various distributions

import numpy as np

def uniform_disk(n, radius=1, center=(0,0)):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.sqrt(np.random.uniform(0, 1))
        x[i] = r * np.cos(a) + center[0]
        y[i] = r * np.sin(a) + center[1]
    return (x, y)

def uniform_circle(n, radius=1, center=(0,0)):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        x[i] = radius * np.cos(a) + center[0]
        y[i] = radius * np.sin(a) + center[1]
    return (x, y)

def uniform_kgon(n, k, center=(0,0)):
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

        x[i] = r * np.cos(theta) + center[0]
        y[i] = r * np.sin(theta) + center[1]

    return (x, y)

def gaussian(n, mu_x=0, mu_y=0, sigma_x=1, sigma_y=1):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        x[i] = np.random.normal(mu_x, sigma_x)
        y[i] = np.random.normal(mu_y, sigma_y)
    return (x, y)

def uniform_annulus(n, inner_radius=0.5, outer_radius=1, center=(0,0)):
    x = [None] * n
    y = [None] * n
    C = 2/(outer_radius**2 - inner_radius**2)
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = np.sqrt(2*np.random.uniform(0,1)/C + inner_radius**2)
        x[i] = r * np.cos(a) + center[0]
        y[i] = r * np.sin(a) + center[1]
    return (x, y)

