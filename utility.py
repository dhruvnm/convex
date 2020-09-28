# various utility functions

import numpy as np
import matplotlib.pyplot as plt

def orient(p1, p2, p3):
    """
    This computes the orientation of p1, p2, and p3. In 2D, this check to see if
    the points form a left turn, right turn, or no turn at all when travelling
    from p1, to p2, to p3.
    Returns 1 if left turn
            -1 if right turn
            0 if no turn.
    """
    sign, _ =  np.linalg.slogdet([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])
    return sign

def pseudo_angle(p1, p2):
    """
    This computes a pseudo angle between the the vector from p1 to p2 against
    the horizontal through p1. When calculating multiple pseudo angles for a
    pointset, we use a constant p1. It should be the lowest point in the set. If
    a pointset is sorted by this pseudo angle, it will also be sorted by the
    actual angle described above.

    The return value will be a tuple containing the pseudo angle and the
    distance between p1 and and p2. This allows sorting by angle primarily and
    by distance from p1 secondarily.
    """

    # if the points are the same, the angle is 0 and the distnace is 0.
    if p1 == p2:
        return (0,0)

    # Calculate the distance from p1 to p2
    length = ((p2[0] - p1[0])**2 + (p2[1] - p1[1]) ** 2) ** (1/2)
    # Scale the height down such that the distance from p1 to p2 is 1
    # This is to make sure the distance has no impact on the pseudo angle.
    height = (p2[1] - p1[1]) / length

    # If the point is the the right of p1, then increasing vertical distance
    # implies increasing angle. If it is to the left, decreasing vertical
    # implied increasing angle. Points to the left have a larger angle than
    # points to the right. We use this information to return a pseudo angle
    # so that the points will be sored correctly.
    if p2[0] >= p1[0]:
        return (height, length)
    else:
        return (2 - height, length)

def draw_hull(hull):
    """Draw the hull on the current plot. For visualization purposes."""
    h = len(hull)
    for i in range(-1, h-1):
        plt.plot([hull[i][0], hull[i+1][0]], [hull[i][1], hull[i+1][1]], 'r')

def pseudo_distance(a, b, c):
    """
    Computes the pseudo distance between the point c and the line formed by a
    and b. The determinant computed below corresponds to the distance. 
    """
    return abs(np.linalg.det([[a[0], b[0], c[0]], [a[1], b[1], c[1]], [1, 1, 1]]))

def akl_toussaint(points):
    """
    Runs the Akl-Toussaint heuristic to get rid of points inside the convex
    quadrilateral formed by the 4 most extreme points. 
    """
    # Search for endpoints
    p_top = (0, float('-Inf'))
    p_bot = (0, float('Inf'))
    p_right = (float('-Inf'), 0)
    p_left = (float('Inf'), 0)
    for point in points:
        if p_top[1] < point[1]:
            p_top = point
        if p_bot[1] > point[1]:
            p_bot = point
        if p_right[0] < point[0]:
            p_right = point
        if p_left[0] > point[0]:
            p_left = point

    # If quadrilateral is not well formed do not reduce
    quad = [p_top, p_left, p_bot, p_right]
    if len(set(quad)) < 4:
        return points

    new_points = []
    for point in points:
        if not inside_quad(point, quad):
            new_points.append(point)

    return new_points, p_top, p_bot, p_right, p_left

def inside_quad(p, quad):
    p1 = quad[-1]
    if p == p1:
        return False

    for p2 in quad:
        if p == p2:
            return False

        if orient(p1, p2, p) <= 0:
            return False

        p1 = p2

    return True
    