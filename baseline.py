# Run all algorithms on all distributions and tabulate the statistics

from r2dist import uniform_disk, uniform_circle, uniform_kgon, gaussian, \
    uniform_annulus, exponential, lognormal, johnsonsu, extreme_value
from convex import grahams_scan, andrews_monotone_chain, divide_and_conquer, \
    jarvis_march, chans_algorithm, quickhull, chans_algorithm_mod, symmetric_hull
from statistics import median, mean
from time import perf_counter_ns
import csv

# 1000 point distributions
N = 1000

# Number of pointsets for each distribution
M = 5

# Create M distributions of each kind with N points
distributions = [uniform_disk, \
                 uniform_circle, \
                 lambda n : uniform_kgon(n, 3), \
                 lambda n : uniform_kgon(n, 4), \
                 lambda n : uniform_kgon(n, 6), \
                 gaussian, \
                 uniform_annulus, \
                 exponential, \
                 lognormal, \
                 johnsonsu, \
                 extreme_value]

pointsets = []

for distribution in distributions:
    pointsets.append([])
    for _ in range(M):
        (x, y) = distribution(N)
        pointsets[-1].append((x, y))

# Get hulls and store statistical data
algorithms = [grahams_scan, \
              andrews_monotone_chain, \
              divide_and_conquer, \
              jarvis_march, \
              chans_algorithm, \
              chans_algorithm_mod, \
              quickhull, \
              symmetric_hull]

min_list = []
max_list = []
mean_list = []
median_list = []

for algorithm in algorithms:
    min_list.append([])
    max_list.append([])
    mean_list.append([])
    median_list.append([])
    for points in pointsets:
        times = []
        for (x, y) in points:
            points = list(zip(x, y))
            t_start = perf_counter_ns()
            algorithm(points)
            t_stop = perf_counter_ns()
            times.append(t_stop - t_start)
        min_list[-1].append(min(times))
        max_list[-1].append(max(times))
        mean_list[-1].append(mean(times))
        median_list[-1].append(median(times))

# Create a CSV file with output
alg_labels = ["Graham's Scan", \
              "Andrew's Monotone Chain", \
              "Divide and Conquer", \
              "Jarvis March", \
              "Chan's Algorithm", \
              "Chan's Algorithm Modified", \
              "Quickhull", \
              "SymmetricHull"]

dist_labels = ["", \
               "Uniform Disk", \
               "Uniform Circle", \
               "Uniform Triangle", \
               "Uniform Square", \
               "Uniform Hexagon", \
               "Gaussian", \
               "Uniform Annulus", \
               "Exponential", \
               "Lognormal", \
               "Johnson's SU", \
               "Extreme Value"]

with open ("min.csv", mode="w", newline="") as min_file:
    writer = csv.writer(min_file, delimiter=",")
    
    writer.writerow(dist_labels)
    for i in range(len(alg_labels)):
        row = [alg_labels[i]] + min_list[i]
        writer.writerow(row)

with open ("max.csv", mode="w", newline="") as max_file:
    writer = csv.writer(max_file, delimiter=",")
    
    writer.writerow(dist_labels)
    for i in range(len(alg_labels)):
        row = [alg_labels[i]] + max_list[i]
        writer.writerow(row)

with open ("mean.csv", mode="w", newline="") as mean_file:
    writer = csv.writer(mean_file, delimiter=",")
    
    writer.writerow(dist_labels)
    for i in range(len(alg_labels)):
        row = [alg_labels[i]] + mean_list[i]
        writer.writerow(row)

with open ("median.csv", mode="w", newline="") as median_file:
    writer = csv.writer(median_file, delimiter=",")
    
    writer.writerow(dist_labels)
    for i in range(len(alg_labels)):
        row = [alg_labels[i]] + median_list[i]
        writer.writerow(row)
