# Generates points in R2 in various distributions

import numpy as np

def v_sparse_hull_circle(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.random.uniform(0, 1) ** 4
        x[i] = r * np.cos(a)
        y[i] = r * np.sin(a)
    return (x, y)

def sparse_hull_circle(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.random.uniform(0, 1)
        x[i] = r * np.cos(a)
        y[i] = r * np.sin(a)
    return (x, y)

def uniform_circle(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.sqrt(np.random.uniform(0, 1))
        x[i] = r * np.cos(a)
        y[i] = r * np.sin(a)
    return (x, y)

def dense_hull_circle(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.random.uniform(0, 1) ** 0.25
        x[i] = r * np.cos(a)
        y[i] = r * np.sin(a)
    return (x, y)

def v_dense_hull_circle(n, radius):
    x = [None] * n
    y = [None] * n
    for i in range(n):
        a = np.random.uniform(0, 2*np.pi)
        r = radius * np.random.uniform(0, 1) ** (1/16)
        x[i] = r * np.cos(a)
        y[i] = r * np.sin(a)
    return (x, y)
