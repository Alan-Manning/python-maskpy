from typing import Literal

import gdspy

from .. import mask_builder_utils as mbu
from ..souk_mask_configs import SoukMaskConfig, get_mask_default_config


def add_combiner_section_and_get_conect_point(
    mask_builder,
    x: float | int,
    y: float | int,
    rot: float | int,
    outer_ring_conection_gap: float | int,
    outer_ring_linewidth: float | int,
    combiner_section_90ghz_config_override: dict[str, float | int] | None = None,
    combiner_section_150ghz_config_override: dict[str, float | int] | None = None,
    combiner_type: Literal["90GHZ", "150GHZ"] = "90GHZ",
    mirror_combiner: bool = False,
    return_configurator_points: bool = False,
) -> list[list[float | int]] | tuple[list[list[float | int]], dict | None]:
    """Adds the phase combiner to the outer ring of the filter bank. This
    is by default the 90GHz or optionally the 150GHz. The specific
    geometries dimensions is determined by the parameters in the config.
    This function will return the coordinate of the conection point.

    Parameters
    ----------
    x, y: float, int
        x, y coordinate for the base of the combiner structure, this is the
        very bottom middle where it connects to the outer ring of the filter
        bank. i.e. the middle of the gap in the outer ring.

    rot: float, int
        The angle (**in radians**) of the rotation of the whole combiner
        geometry.

    outer_ring_conection_gap: float, int
        The gap distance in the outer ring so the combiner can connect to
        either side.

    outer_ring_linewidth: float, int
        The line width of the outer ring of the filter bank this combiner
        attaches to.

    KwArgs
    ------
    combiner_section_90ghz_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    combiner_section_150ghz_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    combiner_type: Literal["90GHZ", "150GHZ"] = "90GHZ"
        This defines what type of phase combiner to add. This can take string
        values of either "90GHZ" or "150GHZ".

    mirror_combiner: bool = False,
        This defines if the combiner should be mirrored around the line
        perpedicular to the outer ring of the filter bank. The default is
        False or can be True.

    return_configurator_points: bool = False
        return a the points for use in the configurator.

    Returns
    -------
    points_to_conect_kids_to: list[list[float | int]]
        list of [x,y] coordinate for the conection point which connects to a
        KID.
    """

    if combiner_type == "90GHZ":
        config = get_mask_default_config(
            SoukMaskConfig.COMBINER_SECTION_90GHZ,
            config_override=combiner_section_90ghz_config_override,
        )
    elif combiner_type == "150GHZ":
        config = get_mask_default_config(
            SoukMaskConfig.COMBINER_SECTION_150GHZ,
            config_override=combiner_section_150ghz_config_override,
        )
    else:
        raise ValueError("'combiner_type' argument should be either '90GHZ' or '150GHZ'")

    # Making First section
    first_linewidth = config["first_linewidth"]  # 5
    distance_from_outer_ring = config["distance_from_outer_ring"]  # 50
    distance_from_outer_ring += (first_linewidth / 2) + (outer_ring_linewidth / 2)

    first_section_height = config["first_section_height"]  # 560
    first_section_width = config["first_section_width"]  # 155
    first_section_back_height = config["first_section_back_height"]  # 105

    first_section_top_points = [
        (-outer_ring_linewidth / 2, outer_ring_conection_gap / 2),
        (distance_from_outer_ring, outer_ring_conection_gap / 2),
        (distance_from_outer_ring, first_section_height / 2),
        (distance_from_outer_ring + first_section_width, first_section_height / 2),
        (distance_from_outer_ring + first_section_width, first_section_height / 2 - first_section_back_height),
    ]

    first_section_bot_points = [
        (-outer_ring_linewidth / 2, -outer_ring_conection_gap / 2),
        (distance_from_outer_ring, -outer_ring_conection_gap / 2),
        (distance_from_outer_ring, -first_section_height / 2),
        (distance_from_outer_ring + first_section_width, -first_section_height / 2),
        (distance_from_outer_ring + first_section_width, -first_section_height / 2 + first_section_back_height),
    ]

    if mirror_combiner:
        mirrored_first_section_top_points = mbu.mirror_points_around_xaxis(first_section_top_points)
        mirrored_first_section_bot_points = mbu.mirror_points_around_xaxis(first_section_bot_points)

        first_section_top = gdspy.FlexPath(
            mirrored_first_section_top_points,
            first_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
        first_section_bot = gdspy.FlexPath(
            mirrored_first_section_bot_points,
            first_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        first_section_top = gdspy.FlexPath(
            first_section_top_points,
            first_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
        first_section_bot = gdspy.FlexPath(
            first_section_bot_points,
            first_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    first_section_top.rotate(rot, (0, 0))
    first_section_top.translate(x, y)
    first_section_bot.rotate(rot, (0, 0))
    first_section_bot.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(first_section_top, mask_builder.layers.Nb_Antenna)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(first_section_bot, mask_builder.layers.Nb_Antenna)
    # self.Main.add(first_section_top)
    # self.Main.add(first_section_bot)

    # Making second section
    second_section_linewidth = config["second_section_linewidth"]  # 4.5

    second_section_top_width1 = config["second_section_top_width1"]  # 25.25
    second_section_top_height1 = config["second_section_top_height1"]  # 16.5
    second_section_top_width2 = config["second_section_top_width2"]  # 20.5
    second_section_top_height2 = config["second_section_top_height2"]  # 16.5
    second_section_top_width3 = config["second_section_top_width3"]  # 29
    second_section_top_height3 = config["second_section_top_height3"]  # 16.5
    second_section_top_width4 = config["second_section_top_width4"]  # 20.5
    second_section_top_height4 = config["second_section_top_height4"]  # 11.5
    second_section_top_width5 = config["second_section_top_width5"]  # 26

    second_section_top_first_point = [first_section_top_points[-1][0], first_section_top_points[-1][1] + second_section_linewidth / 2]

    second_section_top_points = [
        second_section_top_first_point,
        [second_section_top_first_point[0] + second_section_top_width1, second_section_top_first_point[1]],
        [second_section_top_first_point[0] + second_section_top_width1, second_section_top_first_point[1] - second_section_top_height1],
        [
            second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2,
            second_section_top_first_point[1] - second_section_top_height1,
        ],
        [
            second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2,
            second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2,
        ],
        [
            second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2 + second_section_top_width3,
            second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2,
        ],
        [
            second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2 + second_section_top_width3,
            second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2 - second_section_top_height3,
        ],
        [
            second_section_top_first_point[0]
            + second_section_top_width1
            + second_section_top_width2
            + second_section_top_width3
            + second_section_top_width4,
            second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2 - second_section_top_height3,
        ],
        [
            second_section_top_first_point[0]
            + second_section_top_width1
            + second_section_top_width2
            + second_section_top_width3
            + second_section_top_width4,
            second_section_top_first_point[1]
            - second_section_top_height1
            + second_section_top_height2
            - second_section_top_height3
            + second_section_top_height4,
        ],
        [
            second_section_top_first_point[0]
            + second_section_top_width1
            + second_section_top_width2
            + second_section_top_width3
            + second_section_top_width4
            + second_section_top_width5,
            second_section_top_first_point[1]
            - second_section_top_height1
            + second_section_top_height2
            - second_section_top_height3
            + second_section_top_height4,
        ],
    ]

    if mirror_combiner:
        mirrored_second_section_top_points = mbu.mirror_points_around_xaxis(second_section_top_points)
        second_section_top = gdspy.FlexPath(
            mirrored_second_section_top_points,
            second_section_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        second_section_top = gdspy.FlexPath(
            second_section_top_points,
            second_section_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    second_section_top.rotate(rot, (0, 0))
    second_section_top.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(second_section_top, mask_builder.layers.Nb_Antenna)
    # self.Main.add(second_section_top)

    second_section_bot_width1 = config["second_section_bot_width1"]  # 24.25
    second_section_bot_height1 = config["second_section_bot_height1"]  # 15.5
    second_section_bot_width2 = config["second_section_bot_width2"]  # 20.5
    second_section_bot_height2 = config["second_section_bot_height2"]  # 15.5
    second_section_bot_width3 = config["second_section_bot_width3"]  # 24
    second_section_bot_height3 = config["second_section_bot_height3"]  # 15.5
    second_section_bot_width4 = config["second_section_bot_width4"]  # 20.5
    second_section_bot_height4 = config["second_section_bot_height4"]  # 15.5
    second_section_bot_width5 = config["second_section_bot_width5"]  # 32

    second_section_bot_first_point = [first_section_bot_points[-1][0], first_section_bot_points[-1][1] - second_section_linewidth / 2]

    second_section_bot_points = [
        second_section_bot_first_point,
        [second_section_bot_first_point[0] + second_section_bot_width1, second_section_bot_first_point[1]],
        [second_section_bot_first_point[0] + second_section_bot_width1, second_section_bot_first_point[1] + second_section_bot_height1],
        [
            second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2,
            second_section_bot_first_point[1] + second_section_bot_height1,
        ],
        [
            second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2,
            second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2,
        ],
        [
            second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2 + second_section_bot_width3,
            second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2,
        ],
        [
            second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2 + second_section_bot_width3,
            second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2 + second_section_bot_height3,
        ],
        [
            second_section_bot_first_point[0]
            + second_section_bot_width1
            + second_section_bot_width2
            + second_section_bot_width3
            + second_section_bot_width4,
            second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2 + second_section_bot_height3,
        ],
        [
            second_section_bot_first_point[0]
            + second_section_bot_width1
            + second_section_bot_width2
            + second_section_bot_width3
            + second_section_bot_width4,
            second_section_bot_first_point[1]
            + second_section_bot_height1
            - second_section_bot_height2
            + second_section_bot_height3
            - second_section_bot_height4,
        ],
        [
            second_section_bot_first_point[0]
            + second_section_bot_width1
            + second_section_bot_width2
            + second_section_bot_width3
            + second_section_bot_width4
            + second_section_bot_width5,
            second_section_bot_first_point[1]
            + second_section_bot_height1
            - second_section_bot_height2
            + second_section_bot_height3
            - second_section_bot_height4,
        ],
    ]

    if mirror_combiner:
        mirrored_second_section_bot_points = mbu.mirror_points_around_xaxis(second_section_bot_points)
        second_section_bot = gdspy.FlexPath(
            mirrored_second_section_bot_points,
            second_section_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        second_section_bot = gdspy.FlexPath(
            second_section_bot_points,
            second_section_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    second_section_bot.rotate(rot, (0, 0))
    second_section_bot.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(second_section_bot, mask_builder.layers.Nb_Antenna)
    # self.Main.add(second_section_bot)

    start_of_second_section_vertical_linewidth = config["start_of_second_section_vertical_linewidth"]  # 3.5
    start_of_second_section_vertical_points = [
        [
            first_section_top_points[-1][0] - first_linewidth / 2 + start_of_second_section_vertical_linewidth / 2,
            first_section_top_points[-1][1],
        ],
        [
            first_section_bot_points[-1][0] - first_linewidth / 2 + start_of_second_section_vertical_linewidth / 2,
            first_section_bot_points[-1][1],
        ],
    ]
    if mirror_combiner:
        mirrored_start_of_second_section_vertical_points = mbu.mirror_points_around_xaxis(start_of_second_section_vertical_points)
        start_of_second_section_vertical = gdspy.FlexPath(
            mirrored_start_of_second_section_vertical_points,
            start_of_second_section_vertical_linewidth,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        start_of_second_section_vertical = gdspy.FlexPath(
            start_of_second_section_vertical_points,
            start_of_second_section_vertical_linewidth,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    start_of_second_section_vertical.rotate(rot, (0, 0))
    start_of_second_section_vertical.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(start_of_second_section_vertical, mask_builder.layers.Nb_Antenna)
    # self.Main.add(start_of_second_section_vertical)

    third_section_linewidth = config["third_section_linewidth"]  # 3

    end_of_second_section_vertical_linewidth = config["end_of_second_section_vertical_linewidth"]  # 6

    end_of_second_section_vertical_points = [
        [second_section_top_points[-1][0], second_section_top_points[-1][1] + second_section_linewidth / 2],
        [
            second_section_top_points[-1][0] + end_of_second_section_vertical_linewidth,
            second_section_top_points[-1][1] - second_section_linewidth / 2 + third_section_linewidth,
        ],
        [
            second_section_bot_points[-1][0] + end_of_second_section_vertical_linewidth,
            second_section_bot_points[-1][1] + second_section_linewidth / 2 - third_section_linewidth,
        ],
        [second_section_bot_points[-1][0], second_section_bot_points[-1][1] - second_section_linewidth / 2],
    ]
    if mirror_combiner:
        mirrored_end_of_second_section_vertical_points = mbu.mirror_points_around_xaxis(end_of_second_section_vertical_points)
        end_of_second_section_vertical = gdspy.Polygon(
            mirrored_end_of_second_section_vertical_points,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        end_of_second_section_vertical = gdspy.Polygon(
            end_of_second_section_vertical_points,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    end_of_second_section_vertical.rotate(rot, (0, 0))
    end_of_second_section_vertical.translate(x, y)
    mask_builder.Main.add(end_of_second_section_vertical)

    # Making third section
    third_section_width1 = config["third_section_width1"]  # 21.5
    third_section_height1 = config["third_section_height1"]  # 56
    third_section_width2 = config["third_section_width2"]  # 23.5
    third_section_height2 = config["third_section_height2"]  # 56
    third_section_width3 = config["third_section_width3"]  # 23.5

    third_section_start_points = [
        end_of_second_section_vertical_points[1][0],
        end_of_second_section_vertical_points[1][1] - third_section_linewidth / 2,
    ]
    third_section_end_points = [
        end_of_second_section_vertical_points[2][0],
        end_of_second_section_vertical_points[2][1] + third_section_linewidth / 2,
    ]

    third_section_points = [
        third_section_start_points,
        [third_section_start_points[0] + third_section_width1, third_section_start_points[1]],
        [third_section_start_points[0] + third_section_width1, third_section_start_points[1] - third_section_height1],
        [
            third_section_start_points[0] + third_section_width1 + third_section_width2,
            third_section_start_points[1] - third_section_height1,
        ],
        [
            third_section_start_points[0] + third_section_width1 + third_section_width2,
            third_section_start_points[1] - third_section_height1 + third_section_height2,
        ],
        [
            third_section_start_points[0] + third_section_width1 + third_section_width2 + third_section_width3,
            third_section_start_points[1] - third_section_height1 + third_section_height2,
        ],
        [
            third_section_end_points[0] + third_section_width1 + third_section_width2 + third_section_width3,
            third_section_end_points[1] + third_section_height1 - third_section_height2,
        ],
        [
            third_section_end_points[0] + third_section_width1 + third_section_width2,
            third_section_end_points[1] + third_section_height1 - third_section_height2,
        ],
        [
            third_section_end_points[0] + third_section_width1 + third_section_width2,
            third_section_end_points[1] + third_section_height1,
        ],
        [third_section_end_points[0] + third_section_width1, third_section_end_points[1] + third_section_height1],
        [third_section_end_points[0] + third_section_width1, third_section_end_points[1]],
        third_section_end_points,
    ]
    if mirror_combiner:
        mirrored_third_section_points = mbu.mirror_points_around_xaxis(third_section_points)
        third_section = gdspy.FlexPath(
            mirrored_third_section_points,
            third_section_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        third_section = gdspy.FlexPath(
            third_section_points,
            third_section_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    third_section.rotate(rot, (0, 0))
    third_section.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(third_section, mask_builder.layers.Nb_Antenna)
    # self.Main.add(third_section)

    conection_to_kid_path_linewidth = config["conection_to_kid_path_linewidth"]  # 5
    conection_to_kid_path_start_piece_length = config["conection_to_kid_path_start_piece_length"]  # 10

    conection_to_kid_path_start_path_points = [
        [
            third_section_points[5][0] + third_section_linewidth / 2 - conection_to_kid_path_linewidth / 2,
            third_section_points[5][1] - third_section_linewidth / 2,
        ],
        [
            third_section_points[5][0] + third_section_linewidth / 2 - conection_to_kid_path_linewidth / 2,
            third_section_points[5][1] - third_section_linewidth / 2 + conection_to_kid_path_start_piece_length,
        ],
    ]

    if mirror_combiner:
        mirrored_conection_to_kid_path_start_path_points = mbu.mirror_points_around_xaxis(conection_to_kid_path_start_path_points)
        conection_to_kid_path_start = gdspy.FlexPath(
            mirrored_conection_to_kid_path_start_path_points,
            conection_to_kid_path_linewidth,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        conection_to_kid_path_start = gdspy.FlexPath(
            conection_to_kid_path_start_path_points,
            conection_to_kid_path_linewidth,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    conection_to_kid_path_start.rotate(rot, (0, 0))
    conection_to_kid_path_start.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(conection_to_kid_path_start, mask_builder.layers.Nb_Antenna)

    if mirror_combiner:
        points_to_conect_kids_to = mbu.rotate_and_move_single_point(mirrored_conection_to_kid_path_start_path_points[-1], rot, x, y)
    else:
        points_to_conect_kids_to = mbu.rotate_and_move_single_point(conection_to_kid_path_start_path_points[-1], rot, x, y)

    # Making the meander section and connect to frame

    meander_height = config["meander_height"]  # 525
    meander_linewidth = config["meander_linewidth"]  # 5

    meander_initial_gap = config["meander_initial_gap"]  # 10

    meander_final_gap = config["meander_final_gap"]
    meander_gap_spacing = config["meander_gap_spacing"]

    meander_conect_bot_width = config["meander_conect_bot_width"]  # 30
    meander_conect_right_height = config["meander_conect_right_height"]  # 265
    meander_conect_final_width = config["meander_conect_final_width"]  # 7.5

    meander_conect_linewidth = config["meander_conect_linewidth"]  # 5

    no_of_full_straights = int(config["no_of_full_straights"])  # 10
    direction = 1
    meander_points = []
    meander_offset_x = distance_from_outer_ring + meander_initial_gap + (meander_linewidth / 2) + (first_linewidth / 2)

    for i in range(no_of_full_straights - 1):
        direction *= -1
        meander_points.append(
            [meander_offset_x + i * meander_gap_spacing, direction * meander_height / 2 + direction * meander_linewidth / 2]
        )
        meander_points.append(
            [meander_offset_x + i * meander_gap_spacing, -direction * meander_height / 2 - direction * meander_linewidth / 2]
        )

    meander_points.append([meander_points[-1][0] + meander_final_gap, meander_points[-1][1]])
    meander_points.append([meander_points[-1][0], meander_points[-1][1] - meander_height - (meander_linewidth * 1.5)])

    if mirror_combiner:
        mirrored_meander_points = mbu.mirror_points_around_xaxis(meander_points)
        meander = gdspy.FlexPath(
            mirrored_meander_points,
            meander_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Aluminium_Direct.number,
            datatype=mask_builder.layers.Aluminium_Direct.datatype,
        )
    else:
        meander = gdspy.FlexPath(
            meander_points,
            meander_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Aluminium_Direct.number,
            datatype=mask_builder.layers.Aluminium_Direct.datatype,
        )

    meander_bbox_points = mbu.get_flexpath_bounding_box(meander)
    if meander_bbox_points is not None:
        x_min, y_min = meander_bbox_points[0]
        x_max, y_max = meander_bbox_points[1]

        nb_patch_box_oversize_x = config.get("nb_patch_box_oversize_x", 6)
        nb_patch_box_oversize_y = config.get("nb_patch_box_oversize_y", 6)

        # x_min -= config.get("nb_patch_box_oversize_x", 6)
        # x_max += config.get("nb_patch_box_oversize_x", 6)
        # y_min -= config.get("nb_patch_box_oversize_y", 6)
        # y_max += config.get("nb_patch_box_oversize_y", 6)

        x_min -= nb_patch_box_oversize_x
        x_max += nb_patch_box_oversize_x
        y_min -= nb_patch_box_oversize_y
        y_max += nb_patch_box_oversize_y

        meander_bbox_patch = gdspy.Rectangle(
            [x_min, y_min],
            [x_max, y_max],
            layer=mask_builder.layers.Nb_Patch.number,
            datatype=mask_builder.layers.Nb_Patch.datatype,
        )

        meander_bbox_patch.rotate(rot, (0, 0))
        meander_bbox_patch.translate(x, y)

        # self.nb_patch_positives.add(meander_bbox_patch)
        mask_builder.Main.add(meander_bbox_patch)
        # self.layers.Nb_Patch

    meander.rotate(rot, (0, 0))
    meander.translate(x, y)

    mask_builder.make_flexpath_into_polygons_and_add_to_main(meander, mask_builder.layers.Aluminium_Direct)

    meander_to_frame_box_size = config["meander_to_frame_box_size"]

    meander_to_frame_box = gdspy.Rectangle(
        [
            meander_points[-1][0] - (meander_to_frame_box_size / 2),
            meander_points[-1][1] - (meander_to_frame_box_size / 2) + (meander_linewidth / 2),
        ],
        [
            meander_points[-1][0] + (meander_to_frame_box_size / 2),
            meander_points[-1][1] + (meander_to_frame_box_size / 2) + (meander_linewidth / 2),
        ],
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    if mirror_combiner:
        meander_to_frame_box.mirror([0, 0], [1, 0])

    meander_to_frame_box.rotate(rot, (0, 0))
    meander_to_frame_box.translate(x, y)
    mask_builder.Main.add(meander_to_frame_box)

    connect_meander_to_frame_points = [
        [meander_points[-1][0], meander_points[-1][1]],
        [meander_points[-1][0], -meander_height / 2 - meander_linewidth / 2],
        [meander_points[-1][0] + meander_conect_bot_width, -meander_height / 2 - meander_linewidth / 2],
        [meander_points[-1][0] + meander_conect_bot_width, -meander_height / 2 - meander_linewidth / 2 + meander_conect_right_height],
        [
            meander_points[-1][0] + meander_conect_bot_width + meander_conect_final_width,
            -meander_height / 2 - meander_linewidth / 2 + meander_conect_right_height,
        ],
    ]
    if mirror_combiner:
        mirrored_connect_meander_to_frame_points = mbu.mirror_points_around_xaxis(connect_meander_to_frame_points)
        connect_meander_to_frame = gdspy.FlexPath(
            mirrored_connect_meander_to_frame_points,
            meander_conect_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        connect_meander_to_frame = gdspy.FlexPath(
            connect_meander_to_frame_points,
            meander_conect_linewidth,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    connect_meander_to_frame.rotate(rot, (0, 0))
    connect_meander_to_frame.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(connect_meander_to_frame, mask_builder.layers.Nb_Antenna)
    # self.Main.add(connect_meander_to_frame)

    # meander_fork
    meander_last_fork_wdith = config["meander_last_fork_wdith"]  # 5
    meander_last_fork_height = config["meander_last_fork_height"]  # 172
    meander_last_fork_linewdith = config["meander_last_fork_linewdith"]  # 2.5
    meander_last_fork_start_height = config["meander_last_fork_start_height"]  # 95

    meander_fork_start_xy = [
        connect_meander_to_frame_points[2][0] - meander_linewidth / 2,
        connect_meander_to_frame_points[2][1] + meander_last_fork_start_height + meander_linewidth / 2 + meander_last_fork_linewdith / 2,
    ]

    meander_last_fork_points = [
        meander_fork_start_xy,
        [meander_fork_start_xy[0] - meander_last_fork_wdith - meander_last_fork_linewdith / 2, meander_fork_start_xy[1]],
        [
            meander_fork_start_xy[0] - meander_last_fork_wdith - meander_last_fork_linewdith / 2,
            meander_fork_start_xy[1] + meander_last_fork_height + meander_last_fork_linewdith / 2,
        ],
    ]
    if mirror_combiner:
        mirrored_meander_last_fork_points = mbu.mirror_points_around_xaxis(meander_last_fork_points)
        meander_last_fork = gdspy.FlexPath(
            mirrored_meander_last_fork_points,
            meander_last_fork_linewdith,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
    else:
        meander_last_fork = gdspy.FlexPath(
            meander_last_fork_points,
            meander_last_fork_linewdith,
            corners=mbu.create_miter_join,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )

    meander_last_fork.rotate(rot, (0, 0))
    meander_last_fork.translate(x, y)
    mask_builder.make_flexpath_into_polygons_and_add_to_main(meander_last_fork, mask_builder.layers.Nb_Antenna)
    # self.Main.add(meander_last_fork)

    meander_last_fork_top_box_size = config["meander_last_fork_top_box_size"]

    meander_last_fork_top_box = gdspy.Rectangle(
        [meander_last_fork_points[-1][0] - meander_last_fork_top_box_size / 2, meander_last_fork_points[-1][1]],
        [
            meander_last_fork_points[-1][0] + meander_last_fork_top_box_size / 2,
            meander_last_fork_points[-1][1] + meander_last_fork_top_box_size,
        ],
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    if mirror_combiner:
        meander_last_fork_top_box.mirror([0, 0], [1, 0])

    meander_last_fork_top_box.rotate(rot, (0, 0))
    meander_last_fork_top_box.translate(x, y)
    mask_builder.Main.add(meander_last_fork_top_box)

    meander_last_fork_via_box_size_x = config["meander_last_fork_via_box_size_x"]
    meander_last_fork_via_box_size_y = config["meander_last_fork_via_box_size_y"]
    meander_last_fork_via_box_y_offset = config["meander_last_fork_via_box_y_offset"]

    meander_last_fork_via_box = gdspy.Rectangle(
        [
            meander_last_fork_points[-1][0] - meander_last_fork_via_box_size_x / 2,
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset,
        ],
        [
            meander_last_fork_points[-1][0] + meander_last_fork_via_box_size_x / 2,
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset + meander_last_fork_via_box_size_y,
        ],
        layer=mask_builder.layers.SiN_dep.number,
        datatype=mask_builder.layers.SiN_dep.datatype,
    )
    if mirror_combiner:
        meander_last_fork_via_box.mirror([0, 0], [1, 0])

    meander_last_fork_via_box.rotate(rot, (0, 0))
    meander_last_fork_via_box.translate(x, y)
    mask_builder.silicon_nitride_cutouts.add(meander_last_fork_via_box)

    if not return_configurator_points:
        return points_to_conect_kids_to, None

    configurator_points = {}

    # --------------------------------------------------------------------- First Section

    start = mbu.rotate_and_move_single_point(
        [first_section_top_points[0][0], first_section_top_points[0][1] - (first_linewidth / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [first_section_top_points[0][0], first_section_top_points[0][1] + (first_linewidth / 2)], rot, x, y
    )
    configurator_points["first_linewidth"] = {
        "text": "first_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [first_section_top_points[1][0] - distance_from_outer_ring, first_section_top_points[1][1]], rot, x, y
    )
    end = mbu.rotate_and_move_single_point([first_section_top_points[1][0], first_section_top_points[1][1]], rot, x, y)
    configurator_points["distance_from_outer_ring"] = {
        "text": "distance_from_outer_ring",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([first_section_top_points[2][0], first_section_top_points[2][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([first_section_bot_points[2][0], first_section_bot_points[2][1]], rot, x, y)
    configurator_points["first_section_height"] = {
        "text": "first_section_height",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([first_section_top_points[2][0], first_section_top_points[2][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([first_section_top_points[3][0], first_section_top_points[3][1]], rot, x, y)
    configurator_points["first_section_width"] = {
        "text": "first_section_width",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([first_section_top_points[3][0], first_section_top_points[3][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([first_section_top_points[4][0], first_section_top_points[4][1]], rot, x, y)
    configurator_points["first_section_back_height"] = {
        "text": "first_section_back_height",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Second Section Top
    start = mbu.rotate_and_move_single_point(
        [
            second_section_top_points[0][0] + (second_section_linewidth / 2),
            second_section_top_points[0][1] - (second_section_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            second_section_top_points[0][0] + (second_section_linewidth / 2),
            second_section_top_points[0][1] + (second_section_linewidth / 2),
        ],
        rot,
        x,
        y,
    )

    configurator_points["second_section_linewidth_top"] = {
        "text": "second_section_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[0][0], second_section_top_points[0][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[1][0], second_section_top_points[1][1]], rot, x, y)
    configurator_points["second_section_top_width1"] = {
        "text": "second_section_top_width1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[1][0], second_section_top_points[1][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[2][0], second_section_top_points[2][1]], rot, x, y)
    configurator_points["second_section_top_height1"] = {
        "text": "second_section_top_height1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[2][0], second_section_top_points[2][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[3][0], second_section_top_points[3][1]], rot, x, y)
    configurator_points["second_section_top_width2"] = {
        "text": "second_section_top_width2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[3][0], second_section_top_points[3][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[4][0], second_section_top_points[4][1]], rot, x, y)
    configurator_points["second_section_top_height2"] = {
        "text": "second_section_top_height2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[4][0], second_section_top_points[4][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[5][0], second_section_top_points[5][1]], rot, x, y)
    configurator_points["second_section_top_width3"] = {
        "text": "second_section_top_width3",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[5][0], second_section_top_points[5][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[6][0], second_section_top_points[6][1]], rot, x, y)
    configurator_points["second_section_top_height3"] = {
        "text": "second_section_top_height3",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[6][0], second_section_top_points[6][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[7][0], second_section_top_points[7][1]], rot, x, y)
    configurator_points["second_section_top_width4"] = {
        "text": "second_section_top_width4",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[7][0], second_section_top_points[7][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[8][0], second_section_top_points[8][1]], rot, x, y)
    configurator_points["second_section_top_height4"] = {
        "text": "second_section_top_height4",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_top_points[8][0], second_section_top_points[8][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_top_points[9][0], second_section_top_points[9][1]], rot, x, y)
    configurator_points["second_section_top_width5"] = {
        "text": "second_section_top_width5",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Second Section bot
    start = mbu.rotate_and_move_single_point(
        [
            second_section_bot_points[0][0] + (second_section_linewidth / 2),
            second_section_bot_points[0][1] - (second_section_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            second_section_bot_points[0][0] + (second_section_linewidth / 2),
            second_section_bot_points[0][1] + (second_section_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    configurator_points["second_section_linewidth_bot"] = {
        "text": "second_section_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[0][0], second_section_bot_points[0][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[1][0], second_section_bot_points[1][1]], rot, x, y)
    configurator_points["second_section_bot_width1"] = {
        "text": "second_section_bot_width1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[1][0], second_section_bot_points[1][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[2][0], second_section_bot_points[2][1]], rot, x, y)
    configurator_points["second_section_bot_height1"] = {
        "text": "second_section_bot_height1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[2][0], second_section_bot_points[2][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[3][0], second_section_bot_points[3][1]], rot, x, y)
    configurator_points["second_section_bot_width2"] = {
        "text": "second_section_bot_width2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[3][0], second_section_bot_points[3][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[4][0], second_section_bot_points[4][1]], rot, x, y)
    configurator_points["second_section_bot_height2"] = {
        "text": "second_section_bot_height2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[4][0], second_section_bot_points[4][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[5][0], second_section_bot_points[5][1]], rot, x, y)
    configurator_points["second_section_bot_width3"] = {
        "text": "second_section_bot_width3",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[5][0], second_section_bot_points[5][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[6][0], second_section_bot_points[6][1]], rot, x, y)
    configurator_points["second_section_bot_height3"] = {
        "text": "second_section_bot_height3",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[6][0], second_section_bot_points[6][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[7][0], second_section_bot_points[7][1]], rot, x, y)
    configurator_points["second_section_bot_width4"] = {
        "text": "second_section_bot_width4",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[7][0], second_section_bot_points[7][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[8][0], second_section_bot_points[8][1]], rot, x, y)
    configurator_points["second_section_bot_height4"] = {
        "text": "second_section_bot_height4",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([second_section_bot_points[8][0], second_section_bot_points[8][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([second_section_bot_points[9][0], second_section_bot_points[9][1]], rot, x, y)
    configurator_points["second_section_bot_width5"] = {
        "text": "second_section_bot_width5",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Second Section Vertical Line
    start = mbu.rotate_and_move_single_point(
        [
            start_of_second_section_vertical_points[0][0] + (start_of_second_section_vertical_linewidth / 2),
            start_of_second_section_vertical_points[0][1],
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            start_of_second_section_vertical_points[0][0] - (start_of_second_section_vertical_linewidth / 2),
            start_of_second_section_vertical_points[0][1],
        ],
        rot,
        x,
        y,
    )
    configurator_points["start_of_second_section_vertical_linewidth"] = {
        "text": "start_of_second_section_vertical_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [
            end_of_second_section_vertical_points[1][0] - end_of_second_section_vertical_linewidth,
            end_of_second_section_vertical_points[1][1],
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [end_of_second_section_vertical_points[1][0], end_of_second_section_vertical_points[1][1]], rot, x, y
    )
    configurator_points["end_of_second_section_vertical_linewidth"] = {
        "text": "end_of_second_section_vertical_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Third Section Top
    start = mbu.rotate_and_move_single_point(
        [third_section_points[0][0], third_section_points[0][1] - (third_section_linewidth / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [third_section_points[0][0], third_section_points[0][1] + (third_section_linewidth / 2)], rot, x, y
    )
    configurator_points["third_section_linewidth_top"] = {
        "text": "third_section_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[0][0], third_section_points[0][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[1][0], third_section_points[1][1]], rot, x, y)
    configurator_points["third_section_width1_top"] = {
        "text": "third_section_width1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[1][0], third_section_points[1][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[2][0], third_section_points[2][1]], rot, x, y)
    configurator_points["third_section_height1_top"] = {
        "text": "third_section_height1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[2][0], third_section_points[2][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[3][0], third_section_points[3][1]], rot, x, y)
    configurator_points["third_section_width2_top"] = {
        "text": "third_section_width2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[3][0], third_section_points[3][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[4][0], third_section_points[4][1]], rot, x, y)
    configurator_points["third_section_height2_top"] = {
        "text": "third_section_height2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[4][0], third_section_points[4][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[5][0], third_section_points[5][1]], rot, x, y)
    configurator_points["third_section_width3_top"] = {
        "text": "third_section_width3",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Third Section Bot
    start = mbu.rotate_and_move_single_point(
        [third_section_points[11][0], third_section_points[11][1] - (third_section_linewidth / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [third_section_points[11][0], third_section_points[11][1] + (third_section_linewidth / 2)], rot, x, y
    )
    configurator_points["third_section_linewidth_bot"] = {
        "text": "third_section_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[11][0], third_section_points[11][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[10][0], third_section_points[10][1]], rot, x, y)
    configurator_points["third_section_width1_bot"] = {
        "text": "third_section_width1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[10][0], third_section_points[10][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[9][0], third_section_points[9][1]], rot, x, y)
    configurator_points["third_section_height1_bot"] = {
        "text": "third_section_height1",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[9][0], third_section_points[9][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[8][0], third_section_points[8][1]], rot, x, y)
    configurator_points["third_section_width2_bot"] = {
        "text": "third_section_width2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[8][0], third_section_points[8][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[7][0], third_section_points[7][1]], rot, x, y)
    configurator_points["third_section_height2_bot"] = {
        "text": "third_section_height2",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([third_section_points[7][0], third_section_points[7][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([third_section_points[6][0], third_section_points[6][1]], rot, x, y)
    configurator_points["third_section_width3_bot"] = {
        "text": "third_section_width3",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Connection To Kid Path
    start = mbu.rotate_and_move_single_point(
        [
            conection_to_kid_path_start_path_points[0][0] - (conection_to_kid_path_linewidth / 2),
            conection_to_kid_path_start_path_points[0][1],
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            conection_to_kid_path_start_path_points[0][0] + (conection_to_kid_path_linewidth / 2),
            conection_to_kid_path_start_path_points[0][1],
        ],
        rot,
        x,
        y,
    )
    configurator_points["conection_to_kid_path_linewidth"] = {
        "text": "conection_to_kid_path_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [conection_to_kid_path_start_path_points[0][0], conection_to_kid_path_start_path_points[0][1]], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [conection_to_kid_path_start_path_points[1][0], conection_to_kid_path_start_path_points[1][1]], rot, x, y
    )
    configurator_points["conection_to_kid_path_start_piece_length"] = {
        "text": "conection_to_kid_path_start_piece_length",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Meander Section
    start = mbu.rotate_and_move_single_point(
        [meander_points[3][0] - (meander_linewidth / 2), meander_points[3][1] + (meander_linewidth / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [meander_points[3][0] + (meander_linewidth / 2), meander_points[3][1] + (meander_linewidth / 2)], rot, x, y
    )
    configurator_points["meander_linewidth"] = {
        "text": "meander_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([meander_points[2][0], meander_points[2][1] - (meander_linewidth / 2)], rot, x, y)
    end = mbu.rotate_and_move_single_point([meander_points[3][0], meander_points[3][1] + (meander_linewidth / 2)], rot, x, y)
    configurator_points["meander_height"] = {
        "text": "meander_height",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([meander_points[0][0] - (meander_linewidth / 2), meander_points[0][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point(
        [meander_points[0][0] - (meander_linewidth / 2) - meander_initial_gap, meander_points[0][1]], rot, x, y
    )
    configurator_points["meander_initial_gap"] = {
        "text": "meander_initial_gap",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([meander_points[0][0], meander_points[0][1] + (meander_linewidth)], rot, x, y)
    end = mbu.rotate_and_move_single_point([meander_points[3][0], meander_points[3][1] + (meander_linewidth)], rot, x, y)
    configurator_points["meander_gap_spacing"] = {
        "text": "meander_gap_spacing",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([meander_points[1][0], meander_points[1][1] - 4 * (meander_linewidth)], rot, x, y)
    end = mbu.rotate_and_move_single_point([meander_points[-2][0], meander_points[-2][1] - 4 * (meander_linewidth)], rot, x, y)
    configurator_points["no_of_full_straights"] = {
        "text": "no_of_full_straights",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([meander_points[-2][0], meander_points[-2][1] - (meander_linewidth / 2)], rot, x, y)
    end = mbu.rotate_and_move_single_point([meander_points[-3][0], meander_points[-3][1] - (meander_linewidth / 2)], rot, x, y)
    configurator_points["meander_final_gap"] = {
        "text": "meander_final_gap",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Meander To Frame Conenct Box
    start = mbu.rotate_and_move_single_point(
        [meander_points[-1][0] - (meander_to_frame_box_size / 2), meander_points[-1][1] - (meander_linewidth / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [meander_points[-1][0] + (meander_to_frame_box_size / 2), meander_points[-1][1] - (meander_linewidth / 2)], rot, x, y
    )
    configurator_points["meander_to_frame_box_size_bot"] = {
        "text": "meander_to_frame_box_size",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [
            meander_points[-1][0] - (meander_to_frame_box_size / 2),
            meander_points[-1][1] - (meander_to_frame_box_size / 2) + (meander_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            meander_points[-1][0] - (meander_to_frame_box_size / 2),
            meander_points[-1][1] + (meander_to_frame_box_size / 2) + (meander_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_to_frame_box_size_side"] = {
        "text": "meander_to_frame_box_size",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Meander To Frame Conenct Path
    start = mbu.rotate_and_move_single_point(
        [
            connect_meander_to_frame_points[1][0] + (meander_conect_linewidth / 2),
            connect_meander_to_frame_points[1][1] - (meander_conect_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            connect_meander_to_frame_points[1][0] + (meander_conect_linewidth / 2),
            connect_meander_to_frame_points[1][1] + (meander_conect_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_conect_linewidth"] = {
        "text": "meander_conect_linewidth",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([connect_meander_to_frame_points[1][0], connect_meander_to_frame_points[1][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([connect_meander_to_frame_points[2][0], connect_meander_to_frame_points[2][1]], rot, x, y)
    configurator_points["meander_conect_bot_width"] = {
        "text": "meander_conect_bot_width",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([connect_meander_to_frame_points[2][0], connect_meander_to_frame_points[2][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([connect_meander_to_frame_points[3][0], connect_meander_to_frame_points[3][1]], rot, x, y)
    configurator_points["meander_conect_right_height"] = {
        "text": "meander_conect_right_height",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point([connect_meander_to_frame_points[3][0], connect_meander_to_frame_points[3][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point([connect_meander_to_frame_points[4][0], connect_meander_to_frame_points[4][1]], rot, x, y)
    configurator_points["meander_conect_final_width"] = {
        "text": "meander_conect_final_width",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Meander To Frame Conenct Fork
    start = mbu.rotate_and_move_single_point([meander_last_fork_points[0][0], meander_last_fork_points[0][1]], rot, x, y)
    end = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[1][0] + (meander_last_fork_linewdith / 2), meander_last_fork_points[1][1]], rot, x, y
    )
    configurator_points["meander_last_fork_wdith"] = {
        "text": "meander_last_fork_wdith",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[1][0], meander_last_fork_points[1][1] + (meander_last_fork_linewdith / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point([meander_last_fork_points[2][0], meander_last_fork_points[2][1]], rot, x, y)
    configurator_points["meander_last_fork_height"] = {
        "text": "meander_last_fork_height",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[0][0], meander_last_fork_points[0][1] - (meander_last_fork_linewdith / 2)], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[0][0], meander_last_fork_points[0][1] + (meander_last_fork_linewdith / 2)], rot, x, y
    )
    configurator_points["meander_last_fork_linewdith"] = {
        "text": "meander_last_fork_linewdith",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [
            connect_meander_to_frame_points[2][0] - (meander_conect_linewidth / 2) - (meander_last_fork_wdith / 2),
            connect_meander_to_frame_points[2][1] + (meander_conect_linewidth / 2),
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[0][0] - (meander_last_fork_wdith / 2),
            meander_last_fork_points[0][1] - (meander_last_fork_linewdith / 2),
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_last_fork_start_height"] = {
        "text": "meander_last_fork_start_height",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    # --------------------------------------------------------------------- Meander To Frame Conenct Fork Via Box
    start = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[-1][0] - (meander_last_fork_top_box_size / 2), meander_last_fork_points[-1][1]], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[-1][0] + (meander_last_fork_top_box_size / 2), meander_last_fork_points[-1][1]], rot, x, y
    )
    configurator_points["meander_last_fork_top_box_size_bot"] = {
        "text": "meander_last_fork_top_box_size",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[-1][0] - (meander_last_fork_top_box_size / 2), meander_last_fork_points[-1][1]], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[-1][0] - (meander_last_fork_top_box_size / 2),
            meander_last_fork_points[-1][1] + (meander_last_fork_top_box_size),
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_last_fork_top_box_size_side"] = {
        "text": "meander_last_fork_top_box_size",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[-1][0] - (meander_last_fork_via_box_size_x / 2),
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset + meander_last_fork_via_box_size_y,
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[-1][0] + (meander_last_fork_via_box_size_x / 2),
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset + meander_last_fork_via_box_size_y,
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_last_fork_via_box_size_x"] = {
        "text": "meander_last_fork_via_box_size_x",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[-1][0] + (meander_last_fork_via_box_size_x / 2),
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset,
        ],
        rot,
        x,
        y,
    )
    end = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[-1][0] + (meander_last_fork_via_box_size_x / 2),
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset + meander_last_fork_via_box_size_y,
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_last_fork_via_box_size_y"] = {
        "text": "meander_last_fork_via_box_size_y",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    start = mbu.rotate_and_move_single_point(
        [meander_last_fork_points[-1][0] + (meander_last_fork_top_box_size / 2), meander_last_fork_points[-1][1]], rot, x, y
    )
    end = mbu.rotate_and_move_single_point(
        [
            meander_last_fork_points[-1][0] + (meander_last_fork_top_box_size / 2),
            meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset,
        ],
        rot,
        x,
        y,
    )
    configurator_points["meander_last_fork_via_box_y_offset"] = {
        "text": "meander_last_fork_via_box_y_offset",
        "start": [start[0], start[1]],
        "end": [end[0], end[1]],
    }

    return points_to_conect_kids_to, configurator_points
