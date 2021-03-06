# Generates points in R2 in various distributions

import numpy as np
from scipy import stats

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
        u1 = np.sqrt(np.random.uniform(0, 1))
        u2 = np.random.uniform(0, 1)
        px = u1*(1-u2)*p2[0] + u2*u1*p3[0]
        py = u1*(1-u2)*p2[1] + u2*u1*p3[1]

        r = np.sqrt(px**2 + py**2)
        theta = np.arctan(px/py)

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

def exponential(n, beta_x=1, beta_y=1):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        x[i] = np.random.exponential(beta_x)
        y[i] = np.random.exponential(beta_y)
    return (x, y)

def lognormal(n, mu_x=0, mu_y=0, sigma_x=1, sigma_y=1):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        x[i] = np.random.lognormal(mu_x, sigma_x)
        y[i] = np.random.lognormal(mu_y, sigma_y)
    return (x, y)

def johnsonsu(n, a_x=1, a_y=1, b_x=1, b_y=1):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        x[i] = stats.johnsonsu.rvs(a_x, b_x)
        y[i] = stats.johnsonsu.rvs(a_y, b_y)
    return (x, y)

def extreme_value(n, c_x=0, c_y=0):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        x[i] = stats.genextreme.rvs(c_x)
        y[i] = stats.genextreme.rvs(c_y)
    return (x, y)