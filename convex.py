from utility import orient, psuedo_angle, psuedo_distance
import numpy as np
import pdb
import matplotlib.pyplot as plt

def grahams_scan(points):
    if len(points) <= 3:
        return points

    p = (0, float('Inf'))
    for point in points:
        if point[1] < p[1]:
            p = point
        elif point[1] == p[1]:
            if point[0] < p[0]:
                p = point

    points.sort(key = lambda x : psuedo_angle(p, x))
    stack = [points[0], points[1], points[2]]

    for point in points[3:]:
        while len(stack) > 1 and orient(stack[-2], stack[-1], point) < 0:
            stack.pop()
        stack.append(point)

    while len(stack) > 2 and orient(stack[-2], stack[-1], stack[0]) < 0:
        stack.pop()

    return stack

def andrews_monotone_chain(points):
    if len(points) <= 3:
        return points

    points.sort(key = lambda x : x[0])

    upper = []
    for point in points:
        while len(upper) > 1 and orient(upper[-2], upper[-1], point) > 0:
            upper.pop()
        upper.append(point)

    upper = upper[::-1]

    lower = []
    for point in points:
        while len(lower) > 1 and orient(lower[-2], lower[-1], point) < 0:
            lower.pop()
        lower.append(point)

    return upper[1:] + lower[1:]

def divide_and_conquer(points):
    if len(points) <= 3:
        return points

    points.sort(key = lambda x : x[0])
    upper = dac_helper(points, 1)
    lower = dac_helper(points, -1)

    upper = upper[::-1]
    return upper[1:] + lower[1:]

def dac_helper(points, bad_dir):
    if len(points) < 3:
        return points

    if len(points) == 3:
        if orient(points[0], points[1], points[2]) == bad_dir:
            return [points[0], points[2]]
        else:
            return points

    m = len(points) // 2
    left = dac_helper(points[:m], bad_dir)
    right = dac_helper(points[m:], bad_dir)

    lp, rp = len(left) - 1, 0
    found_upper_tangent = False
    while not found_upper_tangent:
        if lp > 0 and orient(left[lp - 1], left[lp], right[rp]) == bad_dir:
            lp -= 1
        elif rp < len(right) - 1 and orient(left[lp], right[rp], right[rp + 1]) == bad_dir:
            rp += 1
        else:
            found_upper_tangent = True

    return left[:lp+1] + right[rp:]

def jarvis_march(points):
    if len(points) <= 3:
        return points

    p = (0, float('Inf'))
    x = float('Inf')
    for point in points:
        if point[1] < p[1]:
            p = point
        if point[0] < x:
            x = point[0]

    hull = [(x - 100, 0), p]
    hull_computed = False

    while not hull_computed:
        next_point = points[0]
        for point in points[1:]:
            if orient(hull[-2], hull[-1], point) > 0 and orient(hull[-1], next_point, point) <= 0:
                next_point = point
        if next_point == hull[1]:
            hull_computed = True
        else:
            hull.append(next_point)

    return hull[1:]

def chans_algorithm(points):
    if len(points) <= 3:
        return points

    p = (0, float('Inf'))
    x = float('Inf')
    for point in points:
        if point[1] < p[1]:
            p = point
        if point[0] < x:
            x = point[0]

    n = len(points)
    m = int(np.ceil(np.log2(np.log2(n))))
    i = 1
    success = False
    for i in range(m):
        h = 2 ** (2 ** (i+1))
        final_hull = [(x - 100, 0), p]

        r = int(np.ceil(n/h))
        mini_hulls = [None] * r
        start = 0
        end = h
        for j in range(r-1):
            mini_hulls[j] = grahams_scan(points[start:end])
            start = end
            end += h
        mini_hulls[-1] = grahams_scan(points[start:])

        for j in range(h):
            next_point = jarvis_binary_search(final_hull[-1], mini_hulls[0])
            for hull in mini_hulls[1:]:
                point = jarvis_binary_search(final_hull[-1], hull)
                if orient(final_hull[-2], final_hull[-1], point) > 0 and orient(final_hull[-1], next_point, point) <= 0:
                    next_point = point
            if next_point == final_hull[1]:
                success = True
                break
            else:
                final_hull.append(next_point)

        if success:
            break

    return final_hull[1:]

def chans_algorithm_mod(points):
    if len(points) <= 3:
        return points

    p = (0, float('Inf'))
    x = float('Inf')
    for point in points:
        if point[1] < p[1]:
            p = point
        if point[0] < x:
            x = point[0]

    n = len(points)
    m = int(np.ceil(np.log2(np.log2(n))))
    i = 1
    success = False
    for i in range(m):
        h = 2 ** (2 ** (i+1))
        final_hull = [(x - 100, 0), p]

        r = int(np.ceil(n/h))
        mini_hulls = [None] * r
        start = 0
        end = h
        for j in range(r-1):
            mini_hulls[j] = grahams_scan(points[start:end])
            start = end
            end += h
        mini_hulls[-1] = grahams_scan(points[start:])

        for j in range(h):
            next_point = jarvis_binary_search(final_hull[-1], mini_hulls[0])
            for hull in mini_hulls[1:]:
                point = jarvis_binary_search(final_hull[-1], hull)
                if orient(final_hull[-2], final_hull[-1], point) > 0 and orient(final_hull[-1], next_point, point) <= 0:
                    next_point = point
            if next_point == final_hull[1]:
                success = True
                break
            else:
                final_hull.append(next_point)

        if success:
            break

        points = []
        for hull in mini_hulls:
            points += hull
        n = len(points)

    return final_hull[1:]


def jarvis_binary_search(q, hull):
    if len(hull) == 1:
        return hull[0]

    n = len(hull)
    i = 0
    k = n - 1
    found_point = False

    while not found_point:
        m = (i + k) // 2
        if i == k:
            found_point = True
        elif q == hull[i]:
            i += 1
            found_point = True
        elif q == hull[k]:
            i = (k + 1) % n
            found_point = True
        elif i == m:
            if orient(q, hull[i], hull[i+1]) < 0:
                i = k
            found_point = True
        elif orient(q, hull[i], hull[i+1]) >= 0 and orient(q, hull[k], hull[(k+1)%n]) < 0:
            found_point = True
        elif orient(q, hull[i], hull[i+1]) < 0 and orient(q, hull[k], hull[(k+1)%n]) >= 0:
            if orient(q, hull[m], hull[m+1]) >= 0:
                k = m
            else:
                i = m
        elif orient(q, hull[i], hull[i+1]) < 0 and orient(q, hull[k], hull[(k+1)%n]) < 0:
            if orient(q, hull[m], hull[m+1]) >= 0:
                k = m
            elif orient(q, hull[i], hull[m]) >= 0:
                k = m
            else:
                i = m
        else:
            if orient(q, hull[m], hull[m+1]) < 0:
                i = m
            elif orient(q, hull[i], hull[m]) < 0:
                k = m
            else:
                i = m

    return hull[i]

class Node:
    def __init__(self, data):
        self.data = data
        self.next = None

def quickhull(points):
    if len(points) <= 3:
        return points

    leftmost = (float('Inf'), 0)
    rightmost = (float('-Inf'), 0)
    for point in points:
        if point[0] < leftmost[0]:
            leftmost = point
        if point[0] > rightmost[0]:
            rightmost = point

    L = Node(leftmost)
    R = Node(rightmost)
    L.next = R
    R.next = L

    top = []
    bot = []

    for point in points:
        side = orient(L.data, R.data, point)
        if side > 0:
            top.append(point)
        elif side < 0:
            bot.append(point)

    qh_helper(top, L, R)
    qh_helper(bot, R, L)

    hull = [L.data]
    start = L
    L = L.next
    while L is not start:
        hull.append(L.data)
        L = L.next

    return hull


def qh_helper(points, P, Q):
    if not points:
        return

    max_dist = float('-Inf')
    max_dist_point = None

    for point in points:
        dist = psuedo_distance(P.data, Q.data, point)
        if  dist > max_dist:
            max_dist = dist
            max_dist_point = point

    C = Node(max_dist_point)
    Q.next = C
    C.next = P

    left = []
    right = []

    for point in points:
        if orient(P.data, C.data, point) > 0:
            left.append(point)
        elif orient(C.data, Q.data, point) > 0:
            right.append(point)

    qh_helper(left, P, C)
    qh_helper(right, C, Q)
    return
