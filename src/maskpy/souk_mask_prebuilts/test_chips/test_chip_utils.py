from collections.abc import Sequence


def get_6_inch_test_chip_centers(
    center_point: Sequence[float | int] = (0, 0),
) -> list[list[float]]:
    """Get the x,y coordinates for the test chip centers on an 6 inch souk test
    chip mask. There are 14 chips with 4 rows and 4 columns with the upper left
    most and upper right most chips excluded. The chip centers are returned in
    the order from top left going right, then the next row down going right and
    so on. See notes for diagram.

    KwArgs
    ------
    center_point: Sequence[float | int] = (0, 0)
        The placement of the test chip centers is relative to (0, 0) but this
        can be altered here.

    Returns
    -------
    chip_centers: list[list[float]]
        The chip center coords [x,y] for the top left chip reading right.

    Notes
    -----
    The chip placements. the ordering is shown in each cell and x shows the
    rough position of the center point.
    >>>      +----+----+
    >>>      | 0  | 1  |
    >>> +----+----+----+----+
    >>> | 2  | 3  | 4  | 5  |
    >>> +----+----+----+----+
    >>> | 6  | 7  x 8  | 9  |
    >>> +----+----+----+----+
    >>> | 10 | 11 | 12 | 13 |
    >>> +----+----+----+----+
    """

    OFFSET_BETWEEN_CHIP_CENTERS = 27300.0

    chip_0_x = (-0.5 * OFFSET_BETWEEN_CHIP_CENTERS) + center_point[0]
    chip_0_y = 47875.0 + center_point[1]

    # Top 2 chips
    chip_centers: list[list[float]] = []

    chip_centers.append(
        [
            chip_0_x,
            chip_0_y,
        ]
    )
    chip_centers.append(
        [
            chip_0_x + OFFSET_BETWEEN_CHIP_CENTERS,
            chip_0_y,
        ]
    )

    # All chips in the grid below
    for v in range(1, 4):
        for h in range(4):
            chip_centers.append(
                [
                    chip_0_x + ((-1 + h) * OFFSET_BETWEEN_CHIP_CENTERS),
                    chip_0_y + (-v * OFFSET_BETWEEN_CHIP_CENTERS),
                ]
            )
    return chip_centers
