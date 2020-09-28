from utility import orient, pseudo_angle, pseudo_distance, akl_toussaint
import numpy as np
import matplotlib.pyplot as plt
from llist import dllist, dllistnode

def grahams_scan(points, ak=False):
    """
    Compute the convex hull of points using Graham's Scan.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, _, _, _ = akl_toussaint(points)
         

    # First we determine the lowest point in the pointset
    # Break ties by taking the point with smaller x-coord
    p = (0, float('Inf'))
    for point in points:
        if point[1] < p[1]:
            p = point
        elif point[1] == p[1]:
            if point[0] < p[0]:
                p = point

    # Shoot a ray horizontally to the right from p
    # Sort the points in order of increasing angle from the ray
    points.sort(key = lambda x : pseudo_angle(p, x))

    stack = [points[0], points[1], points[2]]
    for point in points[3:]:
        # Remove points on the current hull until there is a left turn
        # with respect to the last two points on the hull and the current point
        while len(stack) > 1 and orient(stack[-2], stack[-1], point) < 0:
            stack.pop()
        stack.append(point)

    return stack

def andrews_monotone_chain(points, ak=False):
    """
    Compute the convex hull of points using Andrew's Monotone Chain Algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, _, _, _ = akl_toussaint(points)

    # Sort the points based on increasing x coordinate.
    points.sort(key = lambda x : x[0])

    # Compute the upper hull
    upper = []
    for point in points:
        # Remove points on the current hull until there is a right turn
        # with respect to the last two points on the hull and the current point
        while len(upper) > 1 and orient(upper[-2], upper[-1], point) > 0:
            upper.pop()
        upper.append(point)

    upper = upper[::-1] # Reverse upper hull since we want left turns

    # Computer the lower hull
    lower = []
    for point in points:
        # Remove points on the current hull until there is a left turn
        # with respect to the last two points on the hull and the current point
        while len(lower) > 1 and orient(lower[-2], lower[-1], point) < 0:
            lower.pop()
        lower.append(point)

    # Concatenate the two hulls to form the complete hull
    # Note they will both have the leftmost and rightmost points
    # Thus each of them ignores one of the points.
    return upper[1:] + lower[1:]

def divide_and_conquer(points, ak=False):
    """
    Compute the convex hull of points using the Divide and Conquer algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, _, _, _ = akl_toussaint(points)

    # Sort the points by increasing x-coordinate
    points.sort(key = lambda x : x[0])
    upper = dac_helper(points, 1) # Compute upper hull
    lower = dac_helper(points, -1) # Compute lower hull

    upper = upper[::-1] # Reverse to get left turns
    return upper[1:] + lower[1:] # Concatenate and return complete hull

def dac_helper(points, bad_dir):
    """
    Recursive helper method to compute a partial hull of points
    bad_dir determines which direction we DON'T want when computing the hull
        1 means we don't want left turns
        -1 means we don't want right turns
    Returns a list of points on the partial hull.
    """
    if len(points) < 3:
        return points

    if len(points) == 3:
        # Check to see if the middle point is inside the hull
        if orient(points[0], points[1], points[2]) == bad_dir:
            return [points[0], points[2]]
        else:
            return points

    m = len(points) // 2
    # Recursively compute partial hulls on the left and right half of the points
    left = dac_helper(points[:m], bad_dir)
    right = dac_helper(points[m:], bad_dir)

    lp, rp = len(left) - 1, 0
    found_tangent = False
    # Look for the tangent line that has both partial hulls on the same side.
    # Walk a test line between the two hulls until it is found
    while not found_tangent:
        if lp > 0 and orient(left[lp - 1], left[lp], right[rp]) == bad_dir:
            lp -= 1
        elif rp < len(right) - 1 and orient(left[lp], right[rp], right[rp + 1]) == bad_dir:
            rp += 1
        else:
            found_tangent = True

    # Remove points below the tangent line, concatenate, and return the hull
    return left[:lp+1] + right[rp:]

def jarvis_march(points, ak=False):
    """
    Compute the convex hull of points using the Jarvis March algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, p, _, x = akl_toussaint(points)
        x = x[0]
    else:
        # Find the lowest point in the pointset.
        p = (0, float('Inf'))
        x = float('Inf')
        for point in points:
            if point[1] < p[1]:
                p = point
            if point[0] < x:
                x = point[0]

    # The hull is initialized with a sentinel point that is effectively very
    # far to the left. The first actual point is the lowest point that was
    # computed above.
    hull = [(x - 100, 0), p]
    hull_computed = False

    # Look for the next point on the hull. Essentially looking for the point
    # that cause the least amount of "turning" when compared to the previous
    # edge. This ensures, all other points will be within the hull. The loop
    # ends when we get back to p.
    while not hull_computed:
        next_point = points[0]
        for point in points[1:]:
            if orient(hull[-2], hull[-1], point) > 0 and orient(hull[-1], next_point, point) <= 0:
                next_point = point
        if next_point == hull[1]:
            hull_computed = True
        else:
            hull.append(next_point)

    return hull[1:] # Return the hull minus the sentinel point.

def chans_algorithm(points, ak=False):
    """
    Compute the convex hull using Chan's Algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, p, _, x = akl_toussaint(points)
        x = x[0]
    else:
        # Find the lowest point in the pointset.
        p = (0, float('Inf'))
        x = float('Inf')
        for point in points:
            if point[1] < p[1]:
                p = point
            if point[0] < x:
                x = point[0]

    n = len(points)
    # Max number of iterations needed for the algorithm to succeed
    m = int(np.ceil(np.log2(np.log2(n))))
    i = 1
    success = False
    for i in range(m):
        h = 2 ** (2 ** (i+1)) # the "guess" for how many points on the hull.
        final_hull = [(x - 100, 0), p] # hull has sentinel and p to start

        r = int(np.ceil(n/h)) # number of mini-hulls to compute
        mini_hulls = [None] * r
        start = 0
        end = h
        # Compute all mini hulls using graham's scan
        for j in range(r-1):
            mini_hulls[j] = grahams_scan(points[start:end])
            start = end
            end += h
        mini_hulls[-1] = grahams_scan(points[start:])

        # For each mini-hull, find the correct representative point such that
        # the line through the last point on the hull and the representative is
        # tangent to the mini hull. This is done through a binary search.
        # Then we perform a jarvis march on all of these possible points to
        # find the correct one.
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

    return final_hull[1:] # return hull minus the sentinel

def chans_algorithm_mod(points, ak=False):
    """
    Compute the convex hull using a modified version of Chan's Algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, p, _, x = akl_toussaint(points)
        x = x[0]
    else:
        # Find the lowest point in the pointset.
        p = (0, float('Inf'))
        x = float('Inf')
        for point in points:
            if point[1] < p[1]:
                p = point
            if point[0] < x:
                x = point[0]

    n = len(points)
    # Max number of iterations needed for the algorithm to succeed
    m = int(np.ceil(np.log2(np.log2(n))))
    i = 1
    success = False
    for i in range(m):
        h = 2 ** (2 ** (i+1)) # the "guess" for how many points on the hull.
        final_hull = [(x - 100, 0), p] # hull has sentinel and p to start

        r = int(np.ceil(n/h)) # number of mini-hulls to compute
        mini_hulls = [None] * r
        start = 0
        end = h
        # Compute all mini hulls using graham's scan
        for j in range(r-1):
            mini_hulls[j] = grahams_scan(points[start:end])
            start = end
            end += h
        mini_hulls[-1] = grahams_scan(points[start:])

        # For each mini-hull, find the correct representative point such that
        # the line through the last point on the hull and the representative is
        # tangent to the mini hull. This is done through a binary search.
        # Then we perform a jarvis march on all of these possible points to
        # find the correct one.
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

        # This is an optimization. Clearly points iwthin each mini hull will
        # never be on the final hull. Thus for the next iteration, we only
        # consider the points on each mini hull.
        points = []
        for hull in mini_hulls:
            points += hull
        n = len(points)

    return final_hull[1:] # return hull minus the sentinel


def jarvis_binary_search(q, hull):
    """
    This is a helper method for Chan's algorithm.
    It utilizes binary search to find a point on the hull such that the support
    line through q and this point will have the points on the hull to the left.
    This point is returned.
    """
    if len(hull) == 1:
        return hull[0]

    n = len(hull)
    i = 0
    k = n - 1
    found_point = False

    # A giant case analysis to find the correct point
    # Will possibly detail how this all works at a later date.
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

def quickhull(points, ak=False):
    """
    Compute the convex hull of points using the Quickhull algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, _, _, rightmost, leftmost = akl_toussaint(points)
    else:
        # Find the leftmost and rightmost points in the set.
        leftmost = (float('Inf'), 0)
        rightmost = (float('-Inf'), 0)
        for point in points:
            if point[0] < leftmost[0]:
                leftmost = point
            if point[0] > rightmost[0]:
                rightmost = point

    # The leftmost and rightmost points are guaranteed to be on the hull.
    # We form a circular linked list that will represent this "initial" hull.
    hull = dllist([leftmost, rightmost])
    L = hull.first
    R = hull.last

    # Separate the points into two groups. Those that fall above the current
    # hull and those that fall below.
    top = []
    bot = []
    for point in points:
        side = orient(L.value, R.value, point)
        if side > 0:
            top.append(point)
        elif side < 0:
            bot.append(point)

    # Recursively compute the upper and lower hulls
    qh_helper(top, hull, L, R)
    qh_helper(bot, hull, R, L)

    return list(hull)


def qh_helper(points, hull, P, Q):
    """
    Recursive helper method to find next point on the hull.
    P and Q represent an edge on the current hull.
    points is all points in between P and Q with respect to x-coordinate.
    This method updates the convex hull and returns.
    """
    if not points:
        return

    max_dist = float('-Inf')
    max_dist_point = None

    # Find the point that is the farthest distance from the line formed P and Q.
    for point in points:
        dist = pseudo_distance(P.value, Q.value, point)
        if  dist > max_dist:
            max_dist = dist
            max_dist_point = point

    # The point that is the farther distance is clearly on the convex hull.
    # We add it between P and Q.
    C = dllistnode(max_dist_point)
    hull.insertnode(C, Q)

    # Separate the points into those between P and C and those between C and Q.
    # Discard all points below either of these edges since they obviously
    # will not be on the hull.
    left = []
    right = []
    for point in points:
        if orient(P.value, C.value, point) > 0:
            left.append(point)
        elif orient(C.value, Q.value, point) > 0:
            right.append(point)

    # Recursively compute the rest of the hull and return.
    qh_helper(left, hull, P, C)
    qh_helper(right, hull, C, Q)
    return

def symmetric_hull(points, ak=False):
    """
    Compute the convex hull of points using the SymmetricHull algorithm.
    Returns a list of points on the hull.
    """
    if len(points) <= 3:
        return points

    if ak:
        points, p_top, p_bot, p_right, p_left = akl_toussaint(points)
    else:
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

    # Lexicographic sort by y-coordinate
    points.sort(key = lambda x : (x[1], x[0]))

    # Intialize stacks
    t1= [(p_top, float('Inf'))]
    t2= [(p_top, float('-Inf'))]
    t3= [(p_bot, float('Inf'))]
    t4= [(p_bot, float('-Inf'))]

    # Q3 and Q4
    stop_right = False
    stop_left = False
    for point in points:
        if not stop_right or not stop_left:
            # Q4 case
            if point[0] > p_bot[0] and not stop_right:
                t = t4
                x_compare = lambda p, t: True if p[0] > t[-1][0][0] else False
                s_compare = lambda s, t: True if s < t[-1][1] else False
                if point[0] == p_right[0]:
                    stop_right = True
            
            # Q3 case
            elif point[0] <= p_bot[0] and not stop_left:
                t = t3
                x_compare = lambda p, t: True if p[0] < t[-1][0][0] else False
                s_compare = lambda s, t: True if s > t[-1][1] else False
                if point[0] == p_left[0]:
                    stop_left = True
            else:
                continue
            
            # Ensure the new point adheres to monotonic slope increase/decrease
            if x_compare(point, t):
                slope = (point[1] - t[-1][0][1]) / (point[0] - t[-1][0][0])
                while s_compare(slope, t):
                    t.pop()
                    slope = (point[1] - t[-1][0][1]) / (point[0] - t[-1][0][0])
                t.append((point, slope))
        else:
            break

    # Q1 and Q2
    stop_right = False
    stop_left = False
    for point in reversed(points):
        if not stop_right or not stop_left:
            # Q1 case
            if point[0] > p_top[0] and not stop_right:
                t = t1
                x_compare = lambda p, t: True if p[0] > t[-1][0][0] else False
                s_compare = lambda s, t: True if s > t[-1][1] else False
                if point[0] == p_right[0]:
                    stop_right = True

            # Q2 case
            elif point[0] <= p_bot[0] and not stop_left:
                t = t2
                x_compare = lambda p, t: True if p[0] < t[-1][0][0] else False
                s_compare = lambda s, t: True if s < t[-1][1] else False
                if point[0] == p_left[0]:
                    stop_left = True
            else:
                continue
            
            # Ensure the new point adheres to monotonic slope increase/decrease
            if x_compare(point, t):
                slope = (point[1] - t[-1][0][1]) / (point[0] - t[-1][0][0])
                while s_compare(slope, t):
                    t.pop()
                    slope = (point[1] - t[-1][0][1]) / (point[0] - t[-1][0][0])
                t.append((point, slope))
        else:
            break

    # Concatenate hulls, remove slopes, and return
    hull = t1[::-1] + t2[1:-1] + t3[::-1] + t4[1:-1]
    ret = [x[0] for x in hull]
    return ret