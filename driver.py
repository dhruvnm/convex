from r2dist import *
from convex import grahams_scan, andrews_monotone_chain, divide_and_conquer, jarvis_march, chans_algorithm, quickhull, chans_algorithm_mod
from utility import draw_hull
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

def generate_pointsets(n, r, distribution, num_dists):
    pointsets = []
    for _ in range(num_dists):
        x, y = distribution(n, r)
        pointsets.append(list(zip(x,y)))

    return pointsets

def performance(pointsets, algorithm):
    times = []
    for points in pointsets:
        tic = perf_counter()
        algorithm(points)
        toc = perf_counter()
        times.append(toc - tic)

    return min(times), np.median(times), np.mean(times), max(times)

(x, y) = uniform_annulus(1000)
hull = quickhull(list(zip(x,y)))
plt.figure()
plt.scatter(x, y)
plt.show()
plt.figure()
plt.scatter(x, y)
draw_hull(hull)
plt.show()

str = "{}\n    Min: {} Median: {} Average: {} Max: {}"

# mn, med, avg, mx = performance(pointsets, grahams_scan)
# print(str.format("Graham's Scan", mn, med, avg, mx))

#mn, med, avg, mx = performance(pointsets, andrews_monotone_chain)
#print(str.format("Andrew's Monotone Chain", mn, med, avg, mx))

#mn, med, avg, mx = performance(pointsets, divide_and_conquer)
#print(str.format("Divide and Conquer", mn, med, avg, mx))

#mn, med, avg, mx = performance(pointsets, quickhull)
#print(str.format("QuickHull", mn, med, avg, mx))

#mn, med, avg, mx = performance(pointsets, jarvis_march)
#print(str.format("Jarvis March", mn, med, avg, mx))

#mn, med, avg, mx = performance(pointsets, chans_algorithm)
#print(str.format("Chan's Algorithm", mn, med, avg, mx))

# mn, med, avg, mx = performance(pointsets, chans_algorithm_mod)
# print(str.format("Chan's Algorithm Modified", mn, med, avg, mx))
