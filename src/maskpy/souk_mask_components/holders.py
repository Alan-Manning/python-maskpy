import gdspy
import numpy as np
from numpy import cos, pi, sin

from .. import mask_builder_utils as mbu


def make_Toms_6_inch_holder_and_get_ports(
    mask_builder,
    chip_center_xy: list[float | int] | tuple[float | int, float | int],
    wafer_flat_length: float | int,
    middle_flat_angle: float | int = (-pi / 2),
) -> dict[int, list[float]]:
    """Generates a dictionary of port x, y, and rotaions for the 6inch
    wafer holder design. This contains 18 port locations where the three on
    each edge of the center point and vertexes of a hexagon edge. This
    dictionary has keys "0" through "18" where "0" is botttom most middle
    right of the wafer where the flat is located at the bottom, the
    numbering continues going counter-clockwise. Each dict key contains the
    x,y and roation of the port conection points.

    Parameters
    ----------
    center_xy: list[float | int] | tuple[float | int, float | int]
        list or tuple containing the [x,y] coordinate for the center point
        which to place the octagonal holder around.

    wafer_flat_length: float | int
        The length of the flat on the silicon wafer.

    KwArgs
    ------
    middle_flat_angle: float | int
        Angle (**in radians**) the falt is located. Default is (-pi/2),
        ie the bottom middle side of the hexagon wafer. This angle is the
        angle the middle of the flat makes with the chip center.

    Returns
    -------
    ports_dict: dict[int, list[float]]
        dictionary containing the x,y,rot(**in radians**) for each port in
        the 6inch holder.
    """

    six_inch = 150000
    three_inch = six_inch / 2

    wafer_diameter = six_inch
    wafer_radius = three_inch

    holder_outer_hex_radius = six_inch * 4 / 6
    holder_outer_port_width = 5000

    port_edge_pitch = 18000

    ports_dict: dict[int, list[float]] = {}
    holder_outer_hex_points = []

    holder_positives = gdspy.Cell("HOLDER_POSITIVES")
    holder_cutouts = gdspy.Cell("HOLDER_CUTOUTS")

    # loop making the hexagonal outer and getting port locations and rotations
    for i in range(6):
        ang_cent_to_edge = i * (pi / 3)
        ang_edge_to_cent = ang_cent_to_edge - pi

        holder_outer_hex_points.append(
            [
                chip_center_xy[0] + holder_outer_hex_radius * cos(ang_cent_to_edge - (pi / 6)),
                chip_center_xy[1] + holder_outer_hex_radius * sin(ang_cent_to_edge - (pi / 6)),
            ]
        )

        offset_from_cent_angle = np.arcsin(port_edge_pitch / wafer_radius)

        top_edge_port_xy_rot = [
            chip_center_xy[0] + 1.0025 * wafer_radius * cos(ang_cent_to_edge + offset_from_cent_angle),
            chip_center_xy[1] + 1.0025 * wafer_radius * sin(ang_cent_to_edge + offset_from_cent_angle),
            ang_edge_to_cent,
        ]

        mid_edge_port_xy_rot = [
            chip_center_xy[0] + 1.0025 * wafer_radius * cos(ang_cent_to_edge),
            chip_center_xy[1] + 1.0025 * wafer_radius * sin(ang_cent_to_edge),
            ang_edge_to_cent,
        ]

        bot_edge_port_xy_rot = [
            chip_center_xy[0] + 1.0025 * wafer_radius * cos(ang_cent_to_edge - offset_from_cent_angle),
            chip_center_xy[1] + 1.0025 * wafer_radius * sin(ang_cent_to_edge - offset_from_cent_angle),
            ang_edge_to_cent,
        ]

        ports_dict[(3 * i) + 0] = bot_edge_port_xy_rot
        ports_dict[(3 * i) + 1] = mid_edge_port_xy_rot
        ports_dict[(3 * i) + 2] = top_edge_port_xy_rot

        top_port_rect = gdspy.Rectangle(
            [top_edge_port_xy_rot[0] - wafer_radius, top_edge_port_xy_rot[1] + (holder_outer_port_width / 2)],
            [top_edge_port_xy_rot[0] + wafer_radius, top_edge_port_xy_rot[1] - (holder_outer_port_width / 2)],
            layer=mask_builder.layers.Chip_holder.number,
            datatype=mask_builder.layers.Chip_holder.datatype,
        )
        top_port_rect.rotate(top_edge_port_xy_rot[2], [top_edge_port_xy_rot[0], top_edge_port_xy_rot[1]])
        holder_cutouts.add(top_port_rect)

        mid_port_rect = gdspy.Rectangle(
            [mid_edge_port_xy_rot[0] - wafer_radius, mid_edge_port_xy_rot[1] + (holder_outer_port_width / 2)],
            [mid_edge_port_xy_rot[0] + wafer_radius, mid_edge_port_xy_rot[1] - (holder_outer_port_width / 2)],
            layer=mask_builder.layers.Chip_holder.number,
            datatype=mask_builder.layers.Chip_holder.datatype,
        )
        mid_port_rect.rotate(mid_edge_port_xy_rot[2], [mid_edge_port_xy_rot[0], mid_edge_port_xy_rot[1]])
        holder_cutouts.add(mid_port_rect)

        bot_port_rect = gdspy.Rectangle(
            [bot_edge_port_xy_rot[0] - wafer_radius, bot_edge_port_xy_rot[1] + (holder_outer_port_width / 2)],
            [bot_edge_port_xy_rot[0] + wafer_radius, bot_edge_port_xy_rot[1] - (holder_outer_port_width / 2)],
            layer=mask_builder.layers.Chip_holder.number,
            datatype=mask_builder.layers.Chip_holder.datatype,
        )
        bot_port_rect.rotate(bot_edge_port_xy_rot[2], [bot_edge_port_xy_rot[0], bot_edge_port_xy_rot[1]])
        holder_cutouts.add(bot_port_rect)

    holder_outer_hex = gdspy.Polygon(
        holder_outer_hex_points,
        layer=mask_builder.layers.Chip_holder.number,
        datatype=mask_builder.layers.Chip_holder.datatype,
    )
    holder_positives.add(holder_outer_hex)

    wafer_circle_cutout = mask_builder.make_wafer_shape(
        mask_builder.Chip_holder,
        wafer_radius=wafer_diameter / 2,
        wafer_flat_length=wafer_flat_length,
        chip_center=chip_center_xy,
        flat_angle=middle_flat_angle,
    )
    holder_cutouts.add(wafer_circle_cutout)

    holder = gdspy.boolean(
        holder_positives,
        holder_cutouts,
        "not",
        layer=mask_builder.layers.Chip_holder.number,
        datatype=mask_builder.layers.Chip_holder.datatype,
    )
    mask_builder.Main.add(holder)

    return ports_dict


def generate_octagonal_holder_port_dict(
    mask_builder,
    center_xy: list[float | int] | tuple[float | int, float | int],
    octagon_rotation: float | int,
    wafer_radius: float | int = 75000.0,
) -> dict[str, list[float]]:
    """Generates a dictionary of port x, y, and rotaions for the octagonal
    holder design. This contains 16 port locations where the two on each
    edge of the octagon are 10000um from the center point of that edge.
    This dictionary has keys "0" through "16" where "0" is the left most
    face upper, the numbering continues going counter-clockwise. each key
    contains the x,y and roation of the octagon port conection points.

    Parameters
    ----------
    center_xy: list
        list or tuple containing the [x,y] coordinate for the center point
        which to place the octagonal holder around.

    octagon_rotation: float | int
        the angle (**in Degrees**) for the rotaion of the holder around the
        center_xy argument.

    KwArgs
    ------
    wafer_radius: float | int
        This is the radius of the wafer where the holder positions will be.
        Default is 75000.0 which is a 3inch radius for a six inch wafer.

    Returns
    -------
    oct_ports_dict: dict[str, list[float]]
        dictionary containing the x,y,rot(degrees) for each port in the
        octagonal holder.
    """

    center_wafer_x, center_wafer_y = center_xy
    # the roation of the the octagonal holder for the ports to connect into
    # setting up a dictionary of [x, y, rotation] for each possible conection in the octagonal holder
    oct_ports_dict: dict[str, list[float]] = {}

    oct_ports_dict["0"] = list(
        mbu.rotate(
            center_wafer_x,
            center_wafer_y,
            (center_wafer_x - np.sqrt(wafer_radius**2 - 10000**2)),
            (center_wafer_y + 10000),
            mbu.deg_to_rad(octagon_rotation),
        )
    ) + [
        octagon_rotation
    ]  # sets the x,y and roation of the first two octagon port conection points
    oct_ports_dict["1"] = list(
        mbu.rotate(
            center_wafer_x,
            center_wafer_y,
            (center_wafer_x - np.sqrt(wafer_radius**2 - 10000**2)),
            (center_wafer_y - 10000),
            mbu.deg_to_rad(octagon_rotation),
        )
    ) + [octagon_rotation]

    for i in range(1, 8):  # makes the rest of the ports around the octagonal holder
        oct_ports_dict[f"{2 * i}"] = list(
            mbu.rotate(
                center_wafer_x,
                center_wafer_y,
                oct_ports_dict[f"{2 * (i - 1)}"][0],
                oct_ports_dict[f"{2 * (i - 1)}"][1],
                pi / 4,
            )
        ) + [(oct_ports_dict[f"{2 * (i - 1)}"][2] + 45) % 360]
        oct_ports_dict[f"{(2 * i) + 1}"] = list(
            mbu.rotate(
                center_wafer_x,
                center_wafer_y,
                oct_ports_dict[f"{2 * i - 1}"][0],
                oct_ports_dict[f"{2 * i - 1}"][1],
                pi / 4,
            )
        ) + [(oct_ports_dict[f"{2 * i - 1}"][2] + 45) % 360]

    for i in range(8):  # adds a rectangle around the outside of the wafer shoing where the ports would be+
        mask_builder.Main.add(
            gdspy.Rectangle(
                (oct_ports_dict[f"{2 * i}"][0] - 1000, oct_ports_dict[f"{2 * i}"][1] - 250),
                (oct_ports_dict[f"{2 * i}"][0], oct_ports_dict[f"{2 * i}"][1] + 250),
            ).rotate(
                mbu.deg_to_rad(oct_ports_dict[f"{2 * i}"][2]),
                (oct_ports_dict[f"{2 * i}"][0], oct_ports_dict[f"{2 * i}"][1]),
            )
        )
        mask_builder.Main.add(
            gdspy.Rectangle(
                (oct_ports_dict[f"{(2 * i) + 1}"][0] - 1000, oct_ports_dict[f"{(2 * i) + 1}"][1] - 250),
                (oct_ports_dict[f"{(2 * i) + 1}"][0], oct_ports_dict[f"{(2 * i) + 1}"][1] + 250),
            ).rotate(
                mbu.deg_to_rad(oct_ports_dict[f"{(2 * i) + 1}"][2]),
                (oct_ports_dict[f"{(2 * i) + 1}"][0], oct_ports_dict[f"{(2 * i) + 1}"][1]),
            )
        )

    return oct_ports_dict
