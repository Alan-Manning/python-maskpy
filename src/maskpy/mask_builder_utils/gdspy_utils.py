import numpy as np
from gdspy import FlexPath
from numpy import cos, pi, sin
from numpy.lib import math
from numpy.typing import NDArray


def create_miter_join(p0, v0, p1, v1, p2, w) -> list[list[float]]:
    """Creates a miter on the corner of the bends in a gdspy FlexPath
    object. takes 6 arguments (vertex and direction vector from both
    segments being joined, the center and width of the path) and return a
    list of vertices that make the join.

    Parameters
    ----------
    p0 : array-like[2]
        Vertex [x, y] describing the end of the first path segment.

    v0 : rray-like[2]
        Vector [dx, dy] describing the first path section.

    p1 : array-like[2]
        Vertex [x, y] describing the start of the second path segment.

    v1 : array-like[2]
        Vector [dx, dy] describing the second path section.

    p2 : array-like[2]
        Vertex [x, y] describing center of the joins between the path segments.

    w : float | int
        The width of the path.

    Returns
    ----------
    corner_join_vertices : list[list[float]]
        List of array-like[2] describing the vertices that make the corner join.
    """

    # Calculate intersection point p between lines defined by
    # p0 + u0 * v0 (for all u0) and p1 + u1 * v1 (for all u1)
    den = v1[1] * v0[0] - v1[0] * v0[1]
    lim = 1e-12 * (v0[0] ** 2 + v0[1] ** 2) * (v1[0] ** 2 + v1[1] ** 2)
    if den**2 < lim:
        # Lines are parallel: use mid-point
        u0 = u1 = 0
        p = 0.5 * (p0 + p1)
    else:
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        u0 = (v1[1] * dx - v1[0] * dy) / den
        u1 = (v0[1] * dx - v0[0] * dy) / den
        p = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
    if u0 <= 0 and u1 >= 0:
        # Inner corner
        return [p]
    # Outer corner
    angle0 = np.arctan2(v0[1], v0[0])
    angle1 = np.arctan2(v1[1], v1[0])

    corner_join_vertices = [
        [p0[0] + (w / 2) * cos(angle0 + pi), p0[1] + (w / 2) * sin(angle0 + pi)],
        [p1[0] + (w / 2) * cos(angle1), p1[1] + (w / 2) * sin(angle1)],
    ]

    return corner_join_vertices


def get_flexpath_bounding_box(flexpath: FlexPath) -> NDArray | None:
    """Get the bounding box points for a gdspy FlexPath object.
    Parameters
    ----------
    flexpath: FlexPath
        The flexpath to get the bbox for.

    Returns
    -------
    bbox: NDArray | None
        This is a numpy array of shape 2x2 or None if flexpath has no
        polygons. The array will be [[x_min, y_min], [x_max, y_max]].
    """

    polygon_set = flexpath.to_polygonset()
    if polygon_set is None:
        print(f"get_flexpath_bounding_box got flexpath `{flexpath}` with None type PolygonSet. Returning None.")
        return None

    bbox = polygon_set.get_bounding_box()
    return bbox


def get_polys_from_flexpath(flexpath: FlexPath):
    """Gets the polygon points for the shape created by a gdspy FlexPath
    object.

    Parameters
    ----------
    flexpath : FlexPath
        An instanciated FlexPath object.

    Returns
    -------
    polys : list
        list of [x,y] lists defining the points in the polygon.
    """
    polys = []
    flex_path_ploys = flexpath.get_polygons()
    for poly in flex_path_ploys:
        poly_set = []
        for point in poly:
            poly_set.append([point[0], point[1]])
            if math.isnan(point[0]):
                raise ValueError(f"Encountered NaN x value in flexpath polygon points. point:{point}, in FlexPath:{flexpath}")
            if math.isnan(point[1]):
                raise ValueError(f"Encountered NaN y value in flexpath polygon points. point:{point}, in FlexPath:{flexpath}")

        polys.append(poly_set)

    return polys
