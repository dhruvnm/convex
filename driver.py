# Used as a testing sandbox

from r2dist import *
from convex import grahams_scan, andrews_monotone_chain, divide_and_conquer, jarvis_march, chans_algorithm, quickhull, chans_algorithm_mod, symmetric_hull
from utility import draw_hull
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

(x, y) = uniform_kgon(1000, 6)
plt.figure()
plt.scatter(x, y)
plt.show()
hull = grahams_scan(list(zip(x,y)), True)
plt.figure()
plt.scatter(x, y)
draw_hull(hull)
plt.show()
