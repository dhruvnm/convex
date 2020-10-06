# Run all algorithms on all distributions and tabulate the statistics

from r2dist import uniform_disk, uniform_circle, uniform_kgon, gaussian, \
    uniform_annulus, exponential, lognormal, johnsonsu, extreme_value
from convex import grahams_scan, andrews_monotone_chain, divide_and_conquer, \
    jarvis_march, chans_algorithm, quickhull, chans_algorithm_mod, symmetric_hull
from statistics import median, mean
from time import perf_counter
import csv

# Define the number of points in each distributioon
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

# Each algorithm with and without akl-toussant heuristic
algorithms = [grahams_scan, \
              lambda p : grahams_scan(p, True), \
              andrews_monotone_chain, \
              lambda p : andrews_monotone_chain(p, True), \
              divide_and_conquer, \
              lambda p : divide_and_conquer(p, True), \
              jarvis_march, \
              lambda p : jarvis_march(p, True), \
              chans_algorithm, \
              lambda p : chans_algorithm(p, True), \
              lambda p : chans_algorithm(p, False, symmetric_hull), \
              lambda p : chans_algorithm(p, True, symmetric_hull), \
              chans_algorithm_mod, \
              lambda p : chans_algorithm_mod(p, True), \
              lambda p : chans_algorithm_mod(p, False, symmetric_hull), \
              lambda p : chans_algorithm_mod(p, True, symmetric_hull), \
              quickhull, \
              lambda p : quickhull(p, True), \
              symmetric_hull, \
              lambda p : symmetric_hull(p, True)]

min_list = []
max_list = []
mean_list = []
median_list = []

alg_labels = ["Graham's Scan", \
              "Graham's Scan w/A-T", \
              "Andrew's Monotone Chain", \
              "Andrew's Monotone Chain w/A-T", \
              "Divide and Conquer", \
              "Divide and Conquer w/A-T", \
              "Jarvis March", \
              "Jarvis March w/A-T", \
              "Chan's Algorithm", \
              "Chan's Algorithm w/A-T", \
              "Chan's Algorithm w/Symmetric Hull", \
              "Chan's Algorithm w/A-T and Symmetric Hull", \
              "Chan's Algorithm Modified", \
              "Chan's Algorithm Modified w/A-T", \
              "Chan's Algorithm Modified w/Symmetric Hull", \
              "Chan's Algorithm Modified w/A-T and Symmetric Hull", \
              "Quickhull", \
              "Quickhull w/A-T", \
              "SymmetricHull", \
              "SymmetricHull w/A-T"]

# Create a CSV file with output
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

for i, algorithm in enumerate(algorithms):
    print("Starting Algorithm: ", alg_labels[i])
    min_list.append([])
    max_list.append([])
    mean_list.append([])
    median_list.append([])
    for j, points in enumerate(pointsets):
        times = []
        for (x, y) in points:
            points = list(zip(x, y))
            t_start = perf_counter()
            algorithm(points)
            t_stop = perf_counter()
            times.append(t_stop - t_start)
        min_list[-1].append(min(times))
        max_list[-1].append(max(times))
        mean_list[-1].append(mean(times))
        median_list[-1].append(median(times))
    print("Finishing Algorithm: ", alg_labels[i])

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
