from r2dist import *
from convex import grahams_scan, andrews_monotone_chain, divide_and_conquer, jarvis_march, chans_algorithm
from utility import draw_hull
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

#point_sets = [None] * 100
#for i in range(100):
#    point_sets[i] = uniform_circle(1000, 20)
#
#total_time_graham = 0
#total_time_andrew = 0
#for point_set in point_sets:
#    points = list(zip(point_set[0], point_set[1]))
#    tic = perf_counter()
#    grahams_scan(points)
#    toc = perf_counter()
#    total_time_graham += toc - tic
#
#for point_set in point_sets:
#    points = list(zip(point_set[0], point_set[1]))
#    tic = perf_counter()
#    andrews_monotone_chain(points)
#    toc = perf_counter()
#    total_time_andrew += toc - tic
#
#avg_graham = total_time_graham / 100
#avg_andrew = total_time_andrew / 100
#
#print("Average Graham:", avg_graham)
#print("Average Andrew:", avg_andrew)

(x, y) = uniform_circle(1000, 5)
points = list(zip(x,y))
plt.scatter(x, y)
hull = chans_algorithm(points)
#p = jarvis_binary_search((-4, -4), hull)
draw_hull(hull)
#plt.plot([-4, p[0]],[-4, p[1]], 'g')
plt.show()
