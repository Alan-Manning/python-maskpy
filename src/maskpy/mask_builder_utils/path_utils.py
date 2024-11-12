from typing import Sequence

import numpy as np
from numpy import cos, sin


def get_rounded_path_points(
    path_points: Sequence[Sequence[float | int]],
    bend_radius: float | int,
    no_points_in_bend: int = 10,
) -> list[list[float]]:
    """Generates a list of [x,y] points of a path with rounded corners from
    an input list of [x,y] points for a path. This rounds corners with given
    bend radius, the bend is put before the corner points. e.g. bending the
    path [[0,1], [1,1], [1,0]] with bend radius of 1 will form a
    circular arc from [0,1] to [1,0] like the first quadrent of a circle with
    radius 1.

    Parameters
    ----------
    path_points: Sequence[Sequence[float | int]]
        Sequence of [x,y] Sequences which are the coordinates that the path
        passes through. **Note**, these points will not be neccessarily be
        included in the final returned rounded feedline.

    bend_radius: float, int
        The radius for the bends to be created in the path.

    KwArgs
    ------
    no_points_in_bend: int = 10.0
        This is the number of points to add for each corners bend. This value
        can be made smaller or larger to make the rounded corners smoother or
        more polygonal respectively. The default value of 10 however most of
        the time creating smooth enough corners without creating a huge number
        of points.

    returns
    -------
    rounded_feedline_points: list[list[float]]
        list of [x,y] lists that define coordinates for the rounded path.
    """
    rounded_path_points: list[list[float]] = []
    rounded_path_points.append(
        [
            path_points[0][0],
            path_points[0][1],
        ]
    )

    R = bend_radius

    for i in range(len(path_points) - 2):
        x1 = path_points[i][0]
        y1 = path_points[i][1]

        x2 = path_points[i + 1][0]
        y2 = path_points[i + 1][1]

        x3 = path_points[i + 2][0]
        y3 = path_points[i + 2][1]

        len1 = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        len2 = np.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)

        v1x = (x2 - x1) / len1
        v1y = (y2 - y1) / len1

        v2x = (x3 - x2) / len2
        v2y = (y3 - y2) / len2

        if v1x * v2y - v1y * v2x > 0.0:
            bend_sign = +1
        else:
            bend_sign = -1

        px1 = x1 - bend_sign * (v1y * R)
        py1 = y1 + bend_sign * (v1x * R)
        px2 = x2 - bend_sign * (v2y * R)
        py2 = y2 + bend_sign * (v2x * R)

        den = v1x * v2y - v2x * v1y

        k1 = (v2y * (px2 - px1) - v2x * (py2 - py1)) / den
        k2 = (v1y * (px2 - px1) - v1x * (py2 - py1)) / den

        cx = px1 + k1 * v1x
        cy = py1 + k1 * v1y

        tx1 = x1 + k1 * v1x
        ty1 = y1 + k1 * v1y
        tx2 = x2 + k2 * v2x
        ty2 = y2 + k2 * v2y

        init_angle = np.arctan2((ty1 - cy), (tx1 - cx))

        init_ang_unit_vec = [(tx1 - cx), (ty1 - cy)] / np.linalg.norm([(tx1 - cx), (ty1 - cy)])
        final_ang_unit_vec = [(tx2 - cx), (ty2 - cy)] / np.linalg.norm([(tx2 - cx), (ty2 - cy)])

        ang_diff = np.arccos(np.dot(init_ang_unit_vec, final_ang_unit_vec))

        init_angle = np.arctan2(init_ang_unit_vec[1], init_ang_unit_vec[0])

        loop_angles = np.linspace(init_angle, init_angle + (bend_sign * ang_diff), no_points_in_bend)

        for loop_angle in loop_angles:
            bx = cx + R * cos(loop_angle)
            by = cy + R * sin(loop_angle)
            rounded_path_points.append([bx, by])

    rounded_path_points.append([path_points[-1][0], path_points[-1][1]])

    return rounded_path_points
