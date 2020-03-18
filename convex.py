from utility import orient, psuedo_angle
import numpy as np

def graham_scan(points):
    points.sort(key = psuedo_angle)

    p = (0, float('Inf'))
    for (i, point) in enumerate(points):
        if point[1] < p[1]:
            p = point
            start = i

    stack = []
    n = len(points)

    for i in range(n):
        while len(stack) > 1 and orient(stack[-2], stack[-1], points[(i + start) % n]) < 0:
            stack.pop()
        stack.append(points[(i + start) % n])

    while len(stack) > 2 and orient(stack[-2], stack[-1], stack[0]) < 0:
        stack.pop()

    return stack

def andrew_monotone_chain(points):
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
