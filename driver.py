from r2dist import *
from convex import *
from utility import *
import numpy as np
import matplotlib.pyplot as plt

(x, y) = uniform_circle(1000, 5)
points = list(zip(x,y))
hull = divide_and_conquer(points)
plt.scatter(x, y)
draw_hull(hull)
plt.show()
