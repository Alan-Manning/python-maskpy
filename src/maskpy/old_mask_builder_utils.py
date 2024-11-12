from operator import itemgetter

import numpy as np
from numpy import cos, pi, sin, tan
from numpy.typing import NDArray

# @dataclass(init=False)
# class Layer:
#     name: str
#     number: int
#     datatype: int
#
#     def __init__(self, name: str, number: int, datatype: int):
#         if not isinstance(name, str):
#             raise TypeError(f"Layer's name should be of type str, not {type(name)}.")
#
#         if not isinstance(number, int):
#             raise TypeError(f"Layer's number should be of type int, not {type(number)}.")
#
#         if not isinstance(datatype, int):
#             raise TypeError(f"Layer's datatype should be of type int, not {type(datatype)}.")
#
#         if (number < 0) or (number > (2**16)):
#             raise ValueError(f"Layer's number should be within range 0 -> {int(2**16)}.")
#
#         self.name = name
#         self.number = number
#         self.datatype = datatype


# @dataclass
# class ConfigOverride:
#     name: str
#     value: float | int
#


def determine_polygon_orientation(points: list[list[float | int]]) -> int:
    """Determine the orientation of a polygon, that is whether the points are
    defined in a clockwise or anti-clockwise direction.

    Parameters
    ----------
    points : list[list[float | int]]
        list of [x,y] lists that define the coordinates of all the points the
        polygon to be sized.

    Returns
    -------
    orientation : int
        This is +1 for clockwise orientation, and -1 for anti-clockwise.
    """

    # answer_for_n = (x_n+1 - x_n) * (y_n+1 + y_n)
    # sign of summing all answer_for_n is orientation

    x_cs = np.array(list(map(itemgetter(0), points)))
    y_cs = np.array(list(map(itemgetter(1), points)))

    x_ns = np.roll(x_cs, 1)
    y_ns = np.roll(y_cs, 1)

    sums_for_each_vertex = (x_ns - x_cs) * (y_ns + y_cs)
    final_sum = np.sum(sums_for_each_vertex)

    if final_sum > 0:
        return +1
    else:
        return -1


def round_polygon(
    points: list[list[float | int]],
    bend_radius: float | int,
    bend_points_gap: float | int = 10.0,
) -> list[list[float | int]]:
    """.

    Parameters
    ----------
    points : list[list[float | int]]
        list of [x,y] lists that define the coordinates of all the points the
        polygon.

    bend_radius : float, int
        radius for the bends to be created smoothing the feedline out.

    KwArgs
    ------
    bend_points_gap = 10
        This is the distance (in microns) between points in the rounded section
        of the line. This value can be made smaller or larger to make the
        rounded section smoother or more polygonal respectively. The default
        10um value however is smooth enough without creating a huge number of
        points.

    Returns
    -------
    new_points : list[list[float | int]]
        list of [x,y] lists that define coordinates that the rounded feedline
        passes through. Similar in form to the input list
        ("feedline_pass_through_points") however this new list will contain a
        lot more points. The
    """

    if bend_radius < 0:
        raise ValueError(f"bend_radius should be a positive value.")

    R = bend_radius
    epsilon = 0.001

    rounded_polygon_points: list[list[float]] = []
    rounded_polygon_points.append([points[0][0], points[0][1]])

    for i in range(len(points) - 2):
        x1 = points[i][0]
        y1 = points[i][1]

        x2 = points[i + 1][0]
        y2 = points[i + 1][1]

        x3 = points[i + 2][0]
        y3 = points[i + 2][1]

        len1 = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        len2 = np.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)

        if len1 == len2 == 0:
            continue

        v1x = (x2 - x1) / len1
        v1y = (y2 - y1) / len1

        v2x = (x3 - x2) / len2
        v2y = (y3 - y2) / len2

        cross_product = (v1x * v2y) - (v1y * v2x)
        if abs(cross_product) < epsilon:
            continue
        elif cross_product > 0.0:
            bend = "left"
        else:
            bend = "right"

        if bend == "left":
            px1 = x1 - v1y * R
            py1 = y1 + v1x * R
            px2 = x2 - v2y * R
            py2 = y2 + v2x * R
        else:
            # if bend == "right":
            px1 = x1 + v1y * R
            py1 = y1 - v1x * R
            px2 = x2 + v2y * R
            py2 = y2 - v2x * R

        den = (v1x * v2y) - (v2x * v1y)

        k1 = (v2y * (px2 - px1) - v2x * (py2 - py1)) / den
        k2 = (v1y * (px2 - px1) - v1x * (py2 - py1)) / den

        cx = px1 + k1 * v1x
        cy = py1 + k1 * v1y

        tx1 = x1 + k1 * v1x
        ty1 = y1 + k1 * v1y
        tx2 = x2 + k2 * v2x
        ty2 = y2 + k2 * v2y

        ang_inc = bend_points_gap / R

        init_angle = np.arctan2((ty1 - cy), (tx1 - cx))

        init_ang_unit_vec = [(tx1 - cx), (ty1 - cy)] / np.linalg.norm([(tx1 - cx), (ty1 - cy)])
        final_ang_unit_vec = [(tx2 - cx), (ty2 - cy)] / np.linalg.norm([(tx2 - cx), (ty2 - cy)])

        ang_diff = np.arccos(np.dot(init_ang_unit_vec, final_ang_unit_vec))

        init_angle = np.arctan2(init_ang_unit_vec[1], init_ang_unit_vec[0])

        loop_ang = init_angle

        try:
            int(ang_diff / ang_inc)
        except Exception as e:
            locs = locals()
            print(f"locals = {locs}")
            print(f"locals")
            for key, val in locs.items():
                print(f"{key} = {val}")
            print(f"======")
            print(f"ang_diff = {ang_diff}")
            print(f"ang_inc = {ang_inc}")
            print(f"")
            print(f"init_ang_unit_vec = {init_ang_unit_vec}")
            print(f"final_ang_unit_vec = {final_ang_unit_vec}")
            raise e
        for _ in range(int(ang_diff / ang_inc)):
            bx = cx + R * cos(loop_ang)
            by = cy + R * sin(loop_ang)
            rounded_polygon_points.append([bx, by])

            if bend == "left":
                loop_ang += ang_inc
            else:
                loop_ang -= ang_inc

    rounded_polygon_points.append([points[-1][0], points[-1][1]])

    return rounded_polygon_points


def size_polygon(
    points: list[list[float | int]],
    size_diff: float | tuple[float, float],
    cutoff_angle: float = 130,
    # cutoff_angle: float = 70,
) -> list[list[float | int]]:
    """Increase the size of a polygon (analogous to KLayout's Edit -> Selection -> Size Shapes).

    Parameters
    ----------
    points : list[list[float | int]]
        list of [x,y] lists that define the coordinates of all the points the
        polygon to be sized.

    size_diff: float | tuple[float, float]
        The diff in size for the new polygon. If a single value is passed then
        scaling dx and dy is the same. If a tuple of two values is passed then
        scaling for dx is the first value and dy is the second.

    KwArgs
    ------
    cutoff_angle: float
        The angle **in degrees** between the normals of the edges at which the
        corner will be cutoff to avoid large sharp points.

    Returns
    -------
    new_points: list[list[float | int]]
        The new points of the sized polygon.
    """
    cutoff_angle = cutoff_angle * (pi / 180)

    if isinstance(size_diff, tuple):
        if (size_diff[0] * size_diff[1]) < 0:
            raise ValueError(f"elements of size_diff tuple should be of the same sign. Current values are ({size_diff[0]},{size_diff[1]}).")
        dx = size_diff[0]
        dy = size_diff[1]
    elif isinstance(size_diff, (float, int)):
        dx = size_diff
        dy = size_diff
    else:
        raise TypeError(f"size_diff should be of type float or tuple[float, float]. Not of type {type(size_diff)}.")

    if dx == 0 and dy == 0:
        return points

    size_sign = +1 if dx > 0 else -1

    direction = determine_polygon_orientation(points)

    sign = size_sign * direction

    new_points: list[list[float]] = []

    for i in range(len(points)):

        # p = previous point
        # c = current point
        # n = next point

        if i == 0:
            p = np.array(points[-1])
        else:
            p = np.array(points[i - 1])

        c = np.array(points[i])
        if i == (len(points) - 1):
            n = np.array(points[0])
        else:
            n = np.array(points[i + 1])

        p_c = c - p
        c_n = n - c

        c_edge_angle = np.arctan2(p_c[1], p_c[0]) % (2 * pi)
        n_edge_angle = np.arctan2(c_n[1], c_n[0]) % (2 * pi)

        c_edge_normal = (c_edge_angle - (pi / 2)) % (2 * pi)
        n_edge_normal = (n_edge_angle - (pi / 2)) % (2 * pi)

        angle_between_normals = (n_edge_normal - c_edge_normal) % (2 * pi)

        dL_in_c_norm = sign * np.sqrt((dx * cos(c_edge_normal)) ** 2 + (dy * sin(c_edge_normal)) ** 2)
        dL_in_n_norm = sign * np.sqrt((dx * cos(n_edge_normal)) ** 2 + (dy * sin(n_edge_normal)) ** 2)

        test_x = -34473.85200
        test_y = -15735.90600
        if (test_x - 0.2 < c[0] < test_x + 0.2) and (test_y - 0.2 < c[1] < test_y + 0.2):
            print(angle_between_normals)

        if angle_between_normals > cutoff_angle:
            # make the mitre corner
            mitre_length_in_c_edge_direction = np.sqrt((dx * cos(c_edge_angle)) ** 2 + (dy * sin(c_edge_angle)) ** 2)
            p_c_unit_vec = p_c / np.linalg.norm(p_c)

            point_moved_in_c_edge_normal = np.array(
                [c[0] + (dL_in_c_norm * cos(c_edge_normal)), c[1] + (dL_in_c_norm * sin(c_edge_normal))]
            )

            new_point_1 = point_moved_in_c_edge_normal + (mitre_length_in_c_edge_direction * p_c_unit_vec)

            dL_in_n_edge_direction = np.sqrt((dx * cos(n_edge_angle)) ** 2 + (dy * sin(n_edge_angle)) ** 2)
            n_c_unit_vec = -c_n / np.linalg.norm(c_n)

            point_moved_in_n_edge_normal = np.array(
                [c[0] + (dL_in_n_norm * cos(n_edge_normal)), c[1] + (dL_in_n_norm * sin(n_edge_normal))]
            )

            new_point_2 = point_moved_in_n_edge_normal + (dL_in_n_edge_direction * n_c_unit_vec)

            new_points.append([new_point_1[0], new_point_1[1]])
            new_points.append([new_point_2[0], new_point_2[1]])

        else:
            # make the sharp corner point
            w = dL_in_c_norm / cos((pi / 2) - angle_between_normals)
            m = dL_in_n_norm * tan((pi / 2) - angle_between_normals)

            angle_t = np.arctan((w - m) / dL_in_n_norm)
            angle_T = n_edge_normal - angle_t

            angle_for_new_point = (angle_T) % (2 * pi)

            dL_in_point = sign * np.sqrt(
                ((dL_in_c_norm / cos((pi / 2) - angle_between_normals)) - (dL_in_n_norm * tan((pi / 2) - angle_between_normals))) ** 2
                + (dL_in_n_norm) ** 2
            )

            point_dx = dL_in_point * cos(angle_for_new_point)
            point_dy = dL_in_point * sin(angle_for_new_point)

            new_points.append([c[0] + point_dx, c[1] + point_dy])

    # round = True
    # if round:
    #     return round_polygon(new_points, max(dx, dy))

    return new_points


def get_points_bounding_box(*points: list[list[float | int]]) -> NDArray:
    """Get the bounding box points for a set of points or many sets of points.
    Parameters
    ----------
    *points : list[list[float]]
        This is any number of lists of [x,y] lists defining the x,y for each
        point.

    Returns
    -------
    bbox: NDArray | None
        This is a numpy array of shape 2x2 or None if points have no
        polygons. The array will be [[x_min, y_min], [x_max, y_max]].
    """

    x_mins: list[float] = []
    x_maxs: list[float] = []
    y_mins: list[float] = []
    y_maxs: list[float] = []
    for point_set in points:
        if len(point_set) == 0:  # empty points check
            continue

        points_array = np.array(point_set)

        x_min = np.min(points_array[:, 0])
        y_min = np.min(points_array[:, 1])
        x_max = np.max(points_array[:, 0])
        y_max = np.max(points_array[:, 1])

        x_mins.append(x_min)
        y_mins.append(y_min)
        x_maxs.append(x_max)
        y_maxs.append(y_max)

    try:
        return np.array([[min(x_mins), min(y_mins)], [max(x_maxs), max(y_maxs)]])
    except ValueError as ve:
        print(f"Unable to find bbox for points provided:\n{points}")
        raise ve


def mirror_points_around_yaxis(points: list[list[float] | tuple[float, float]]) -> list[list[float]]:
    """Mirrors a set of points around the y axis. Just flipping the x sign.

    Parameters
    ----------
    points : list[list[float]]
        list of [x,y] lists defining the x,y for each point to mirror.

    Returns
    -------
    mirrored_points: list[list[float]]
        The new mirrored points around the y axis.
    """
    mirrored_points = points.copy()
    for i in range(len(points)):
        old = list(points[i])
        px = -old[0]
        py = old[1]
        mirrored_points[i] = [px, py]

    return mirrored_points


def mirror_points_around_xaxis(points: list[list[float]]) -> list[list[float]]:
    """Mirrors a set of points around the x axis. Just flipping the y sign.

    points : list[list[float]]
        list of [x,y] lists defining the x,y for each point to mirror.

    Returns
    -------
    mirrored_points: list[list[float]]
        The new mirrored points around the x axis.
    """
    mirrored_points = points.copy()
    for i in range(len(points)):
        old = list(points[i])
        px = old[0]
        py = -old[1]
        mirrored_points[i] = [px, py]

    return mirrored_points


def inside_hexagon(
    xy: list[float],
    hex_rad: float,
    hex_center: list[float] = [0, 0],
):
    """Checks if an xy point sits within a hexagon with given radius.

    Parameters
    ----------
    xy : list
        list of [x,y] coordinates of a point to check if its inside the
        hexagon.

    hex_rad : float, int
        The radius of the hexagon to check if point sits inside.

    KwArgs
    ------
    hex_center = [0,0]
        list of the [x,y] coordinates for the center of the haxagon.

    Returns
    -------
    boolean :
        True if the xy point given sits within the hexagon.
        False if the xy point sits on the boundary or outside the hexagon.
    """

    R = hex_rad
    x, y = map(abs, [xy[0] - hex_center[0], xy[1] - hex_center[1]])
    cx, cy = hex_center
    return y < np.sqrt(3) * min((R - x), (R / 2))


def deg_to_rad(deg: int | float) -> float:
    """Converts an angle given in degrees to radians."""
    return (deg / 180.0) * pi


def rad_to_deg(rad: int | float) -> float:
    """Converts an angle given in radians to degrees."""
    return rad * (180.0 / pi)


def rotate(
    ox: float,
    oy: float,
    px: float,
    py: float,
    angle: float,
) -> tuple[float, float]:
    """Rotate a point counterclockwise by a given angle around a given
    origin.

    Parameters
    ----------
    ox, oy: float
        The origin x, y to rotate around

    px, py: float,
        The point x, y to rotate.

    angle: float,
        The angle in radians to rotate.

    """
    qx = ox + cos(angle) * (px - ox) - sin(angle) * (py - oy)
    qy = oy + sin(angle) * (px - ox) + cos(angle) * (py - oy)
    return qx, qy


def rotate_points_list(
    points: list[list[float]],
    rot_angle: float,
    ox: float = 0.0,
    oy: float = 0.0,
) -> list[list[float]]:
    """Rotates a list of points lists all by angle. Rotates counterclockwise by
    a given angle around a given origin. The angle should be given in radians.

    Parameters
    ----------
    points : list[list[float]]
        list of many [x,y] lists to be moved.

    rot_angle: float,
        The angle (**in radians**) to rotate counterclockwise by.

    KwArgs
    ------
    ox: float = 0.0,
        The x origin of rotation. Default is 0.0.

    oy: float = 0.0,
        The y origin of rotation. Default is 0.0.

    Returns
    -------
    new_points : list[list[float]]
        Returns a new list of points that have been roated. List is of the
        same form as input points list.
    """
    points_copy = points.copy()
    new_points = []
    for i in range(len(points_copy)):  # rotating each element
        old = list(points_copy[i])
        px = old[0]
        py = old[1]
        qx = ox + cos(rot_angle) * (px - ox) - sin(rot_angle) * (py - oy)
        qy = oy + sin(rot_angle) * (px - ox) + cos(rot_angle) * (py - oy)
        new_points.append([qx, qy])

    return new_points


def rotate_and_move_single_point(
    point: list[float | int],
    rot_angle: float,
    dx: float | int,
    dy: float | int,
    ox: float | int = 0.0,
    oy: float | int = 0.0,
) -> list[float]:
    """Rotates and moves a single point [x,y] by angle and dx,dy.

    Rotates counterclockwise by a given angle around a given origin.
    The angle should be given in radians. Rotation is done first.

    Parameters
    ----------
    points : list[float]
        list of [x,y] coordinates to be moved and rotated.

    rot_angle: float,
        The angle (**in radians**) to rotate counterclockwise by.

    dx: float | int,
        The delta x to move the points in the points list.

    dy: float | int,
        The delta y to move the points in the points list.

    KwArgs
    ------
    ox: float | int = 0.0
        The x origin of rotation. Default is 0.0.

    oy: float | int = 0.0
        The y origin of rotation. Default is 0.0.

    Returns
    -------
    new_points : list[float]
        Returns a new list of point's coordinates that have been roated and
        then moved. List is of the same form as input points list.
    """
    old = list(point.copy())
    px = old[0]
    py = old[1]
    qx = ox + cos(rot_angle) * (px - ox) - sin(rot_angle) * (py - oy)
    qy = oy + sin(rot_angle) * (px - ox) + cos(rot_angle) * (py - oy)
    new_point = [qx + dx, qy + dy]

    return new_point


def rotate_and_move_points_list(
    points: list[list[float]],
    rot_angle: float,
    dx: float | int,
    dy: float | int,
    ox: float | int = 0.0,
    oy: float | int = 0.0,
) -> list[list[float]]:
    """Rotates and moves a list of points lists all by angle and dx,dy. The
    rotation is done before the translation. Rotates counterclockwise by a
    given angle around a given origin. The angle should be given in radians.

    Parameters
    ----------
    points : list[list[float]]
        list of many [x,y] lists to be moved.

    rot_angle: float | int
        The angle (**in radians**) to rotate counterclockwise by.

    dx: float | int
        The delta x to move the points in the points list.

    dy: float | int
        The delta y to move the points in the points list.

    KwArgs
    ------
    ox: float | int = 0.0
        The x origin of rotation. Default is 0.0.

    oy: float | int = 0.0
        The y origin of rotation. Default is 0.0.

    Returns
    -------
    new_points : list[list[float]]
        Returns a new list of points that have been roated and then moved.
        List is of the same form as input points list.
    """
    points_copy = points.copy()
    new_points = []
    for i in range(len(points_copy)):  # shifting each element
        old_point = list(points_copy[i])
        # px = old[0]
        # py = old[1]
        # qx = ox + cos(rot_angle) * (px - ox) - sin(rot_angle) * (py - oy)
        # qy = oy + sin(rot_angle) * (px - ox) + cos(rot_angle) * (py - oy)
        new_point = rotate_and_move_single_point(old_point, rot_angle, dx, dy, ox=ox, oy=oy)
        new_points.append(new_point)

    return new_points


def move_points_list(
    points: list[list[float]],
    dx: float,
    dy: float,
) -> list[list[float]]:
    """Moves a list of points lists all by dx,dy.

    Parameters
    ----------
    points : list[list[float]]
        list of many [x,y] lists to be moved.

    dx : float
        The delta x to move the points in the points list.

    dy : float
        The delta y to move the points in the points list.

    Returns
    -------
    new_points : list
        Returns a new list of points that have been moved. List is of the same
        form as input points list.
    """
    points_copy = points.copy()
    new_points = []
    for i in range(len(points_copy)):  # shifting each element
        old = list(points_copy[i])
        px = old[0]
        py = old[1]
        new_points.append([px + dx, py + dy])

    return new_points


def move_single_point(
    point: list[float],
    dx: float,
    dy: float,
) -> list[float]:
    """Moves an [x,y] point by dx,dy.

    Parameters
    ----------
    point : list[float]
        Point [x,y] list to be moved.

    dx : float
        The delta x to move the point by.

    dy : float
        The delta y to move the point by.

    Returns
    -------
    new_point : list[float]
        Returns a new point that has been moved. List is of the same form as
        input point list.
    """
    point_copy = point.copy()
    px = point_copy[0]
    py = point_copy[1]
    new_point = [px + dx, py + dy]

    return new_point
