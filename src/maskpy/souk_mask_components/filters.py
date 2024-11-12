import gdspy
import numpy as np
from numpy import cos, pi, sin

from .. import mask_builder_utils as mbu
from ..souk_mask_configs import SoukMaskConfig, get_mask_default_config


def add_Lo_pass_filters(
    mask_builder,
    x: float | int,
    y: float | int,
    inner_ring_line_width: float | int,
    inner_ring_radius: float | int,
    init_angle: float | int,
    direction: str,
    Lo_pass_filters_config_override: dict[str, float | int] | None = None,
    return_configurator_points: bool = False,
):
    """Adds the low pass filter arms around the inner ring of the filter
    bank.

    Parameters
    ----------
    x, y: float, int
        x, y coordinate of the center of the antenna to be places around.

    inner_ring_line_width: float, int
        Line width of the inner ring that the filters conect to.

    inner_ring_radius: float, int
        Radius of the inner ring that the filters conect to.

    init_angle: float, int
        Initial angle (*in radians*) to start at when making the filter arms.

    direction: str
        The direction around the inner ring in which to make the filters point.
        Allowed values are `clockwise` or `anti-clockwise`.

    KwArgs
    ------
    Lo_pass_filters_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    return_configurator_points=False
        return a the points for use in the configurator.
    """
    config = get_mask_default_config(SoukMaskConfig.LO_PASS_FILTERS, config_override=Lo_pass_filters_config_override)

    # LoPass arm details
    Lo_pass_arm1_linewidth = config["Lo_pass_arm1_linewidth"]  # 9.5
    Lo_pass_arm1_height = config["Lo_pass_arm1_height"]  # 20
    Lo_pass_arm1_length = config["Lo_pass_arm1_length"]  # 245
    Lo_pass_arm1_offset = config["Lo_pass_arm1_offset"]  # 152
    Lo_pass_arm1_offset_arc_angle = Lo_pass_arm1_offset / inner_ring_radius

    Lo_pass_arm2_linewidth = config["Lo_pass_arm2_linewidth"]  # 20.5
    Lo_pass_arm2_height = config["Lo_pass_arm2_height"]  # 41.5
    Lo_pass_arm2_length = config["Lo_pass_arm2_length"]  # 218
    Lo_pass_arm2_offset = config["Lo_pass_arm2_offset"]  # 278
    Lo_pass_arm2_offset_arc_angle = Lo_pass_arm2_offset / inner_ring_radius

    Lo_pass_arm3_linewidth = config["Lo_pass_arm3_linewidth"]  # 20
    Lo_pass_arm3_height = config["Lo_pass_arm3_height"]  # 40
    Lo_pass_arm3_length = config["Lo_pass_arm3_length"]  # 219
    Lo_pass_arm3_offset = config["Lo_pass_arm3_offset"]  # 278
    Lo_pass_arm3_offset_arc_angle = Lo_pass_arm3_offset / inner_ring_radius

    Lo_pass_arm4_linewidth = config["Lo_pass_arm4_linewidth"]  # 20.5
    Lo_pass_arm4_height = config["Lo_pass_arm4_height"]  # 41.5
    Lo_pass_arm4_length = config["Lo_pass_arm4_length"]  # 218
    Lo_pass_arm4_offset = config["Lo_pass_arm4_offset"]  # 278
    Lo_pass_arm4_offset_arc_angle = Lo_pass_arm4_offset / inner_ring_radius

    Lo_pass_arm5_linewidth = config["Lo_pass_arm5_linewidth"]  # 9.5
    Lo_pass_arm5_height = config["Lo_pass_arm5_height"]  # 120.5
    Lo_pass_arm5_length = config["Lo_pass_arm5_length"]  # 145
    Lo_pass_arm5_offset = config["Lo_pass_arm5_offset"]  # 278
    Lo_pass_arm5_offset_arc_angle = Lo_pass_arm5_offset / inner_ring_radius

    via_overlap_length = config["via_overlap_length"]  # 5
    via_extend_length = config["via_extend_length"]  # 5

    if direction == "clockwise":
        sign = -1
    elif direction == "anti-clockwise":
        sign = +1

    arm_linewidths = [
        Lo_pass_arm1_linewidth,
        Lo_pass_arm2_linewidth,
        Lo_pass_arm3_linewidth,
        Lo_pass_arm4_linewidth,
        Lo_pass_arm5_linewidth,
    ]
    arm_heights = [Lo_pass_arm1_height, Lo_pass_arm2_height, Lo_pass_arm3_height, Lo_pass_arm4_height, Lo_pass_arm5_height]
    arm_lengths = [Lo_pass_arm1_length, Lo_pass_arm2_length, Lo_pass_arm3_length, Lo_pass_arm4_length, Lo_pass_arm5_length]
    offset_arc_angles = [
        Lo_pass_arm1_offset_arc_angle,
        Lo_pass_arm2_offset_arc_angle,
        Lo_pass_arm3_offset_arc_angle,
        Lo_pass_arm4_offset_arc_angle,
        Lo_pass_arm5_offset_arc_angle,
    ]
    offset_arc_cumsum_angles = np.cumsum(offset_arc_angles)

    configurator_points = {}

    for i in range(5):
        angle = init_angle + sign * offset_arc_cumsum_angles[i]
        x_pos = x + (inner_ring_radius * cos(angle))
        y_pos = y + (inner_ring_radius * sin(angle))

        arm_geom_points = [
            [(sign) * arm_linewidths[i] / 2, 0],
            [(sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
            [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2)],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i],
                arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
            ],
            [(-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i], arm_heights[i] + (inner_ring_line_width / 2)],
            [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
            [(-sign) * arm_linewidths[i] / 2, 0],
        ]

        arm_geom = mbu.rotate_and_move_points_list(arm_geom_points, angle - pi / 2, x_pos, y_pos)

        configurator_points[f"Lo_pass_arm{i+1}_linewidth"] = {
            "text": f"Lo_pass_arm{i+1}_linewidth",
            "start": [arm_geom[2][0], arm_geom[2][1]],
            "end": [arm_geom[5][0], arm_geom[5][1]],
        }

        configurator_points[f"Lo_pass_arm{i+1}_length"] = {
            "text": f"Lo_pass_arm{i+1}_length",
            "start": [arm_geom[4][0], arm_geom[4][1]],
            "end": [arm_geom[5][0], arm_geom[5][1]],
        }

        configurator_points[f"Lo_pass_arm{i+1}_height"] = {
            "text": f"Lo_pass_arm{i+1}_height",
            "start": [arm_geom[5][0], arm_geom[5][1]],
            "end": [arm_geom[6][0], arm_geom[6][1]],
        }

        arm = gdspy.Polygon(
            arm_geom,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
        mask_builder.Main.add(arm)

        via_points = [
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
            ],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                arm_heights[i] + (inner_ring_line_width / 2),
            ],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                arm_heights[i] + (inner_ring_line_width / 2),
            ],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
            ],
        ]
        via_geom = mbu.rotate_and_move_points_list(via_points, angle - pi / 2, x_pos, y_pos)

        configurator_points[f"via_extend_length_{i+1}"] = {
            "text": "via_extend_length",
            "start": [via_geom[0][0] + (via_geom[3][0] - via_geom[0][0]) / 2, via_geom[0][1] + (via_geom[3][1] - via_geom[0][1]) / 2],
            "end": [via_geom[3][0], via_geom[3][1]],
        }

        configurator_points[f"via_overlap_length_{i+1}"] = {
            "text": "via_overlap_length",
            "start": [via_geom[0][0] + (via_geom[3][0] - via_geom[0][0]) / 2, via_geom[0][1] + (via_geom[3][1] - via_geom[0][1]) / 2],
            "end": [via_geom[0][0], via_geom[0][1]],
        }

        via_box = gdspy.Polygon(
            via_geom,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(via_box)

    if not return_configurator_points:
        return

    configurator_points[f"Lo_pass_arm1_offset"] = {
        "text": f"Lo_pass_arm1_offset",
        "start": [x + inner_ring_radius * cos(init_angle), y + inner_ring_radius * sin(init_angle)],
        "end": [
            x + inner_ring_radius * cos(init_angle + (sign * offset_arc_cumsum_angles[0])),
            y + inner_ring_radius * sin(init_angle + (sign * offset_arc_cumsum_angles[0])),
        ],
    }
    for i in range(1, 5):
        configurator_points[f"Lo_pass_arm{i+1}_offset"] = {
            "text": f"Lo_pass_arm{i+1}_offset",
            "start": [
                x + inner_ring_radius * cos(init_angle + (sign * offset_arc_cumsum_angles[i - 1])),
                y + inner_ring_radius * sin(init_angle + (sign * offset_arc_cumsum_angles[i - 1])),
            ],
            "end": [
                x + inner_ring_radius * cos(init_angle + (sign * offset_arc_cumsum_angles[i])),
                y + inner_ring_radius * sin(init_angle + (sign * offset_arc_cumsum_angles[i])),
            ],
        }

    return configurator_points


def add_Hi_pass_filters(
    mask_builder,
    x: float | int,
    y: float | int,
    inner_ring_line_width: float | int,
    inner_ring_radius: float | int,
    init_angle: float | int,
    direction: str,
    Hi_pass_filters_config_override: dict[str, float | int] | None = None,
    return_configurator_points: bool = False,
):
    """Adds the high pass filter arms around the inner ring of the filter
    bank.

    Parameters
    ----------
    x, y: float, int
        x, y coordinate of the center of the antenna to be places around.

    inner_ring_line_width: float, int
        Line width of the inner ring that the filters conect to.

    inner_ring_radius: float, int
        Radius of the inner ring that the filters conect to.

    init_angle: float, int
        Initial angle (*in radians*) to start at when making the filter arms.

    direction: str
        The direction around the inner ring in which to make the filters point.
        Allowed values are `clockwise` or `anti-clockwise`.

    KwArgs
    ------
    Hi_pass_filters_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    return_configurator_points=False
        return a the points for use in the configurator.
    """
    config = get_mask_default_config(SoukMaskConfig.HI_PASS_FILTERS, config_override=Hi_pass_filters_config_override)

    # HiPass arm details
    Hi_pass_arm1_linewidth = config["Hi_pass_arm1_linewidth"]  # 22
    Hi_pass_arm1_height = config["Hi_pass_arm1_height"]  # 23
    Hi_pass_arm1_length = config["Hi_pass_arm1_length"]  # 120.5
    Hi_pass_arm1_offset = config["Hi_pass_arm1_offset"]  # 185
    Hi_pass_arm1_offset_arc_angle = Hi_pass_arm1_offset / inner_ring_radius

    Hi_pass_arm2_linewidth = config["Hi_pass_arm2_linewidth"]  # 41.5
    Hi_pass_arm2_height = config["Hi_pass_arm2_height"]  # 34.5
    Hi_pass_arm2_length = config["Hi_pass_arm2_length"]  # 97
    Hi_pass_arm2_offset = config["Hi_pass_arm2_offset"]  # 160
    Hi_pass_arm2_offset_arc_angle = Hi_pass_arm2_offset / inner_ring_radius

    Hi_pass_arm3_linewidth = config["Hi_pass_arm3_linewidth"]  # 40.5
    Hi_pass_arm3_height = config["Hi_pass_arm3_height"]  # 34.5
    Hi_pass_arm3_length = config["Hi_pass_arm3_length"]  # 99.5
    Hi_pass_arm3_offset = config["Hi_pass_arm3_offset"]  # 160
    Hi_pass_arm3_offset_arc_angle = Hi_pass_arm3_offset / inner_ring_radius

    Hi_pass_arm4_linewidth = config["Hi_pass_arm4_linewidth"]  # 41.5
    Hi_pass_arm4_height = config["Hi_pass_arm4_height"]  # 34.5
    Hi_pass_arm4_length = config["Hi_pass_arm4_length"]  # 97
    Hi_pass_arm4_offset = config["Hi_pass_arm4_offset"]  # 160
    Hi_pass_arm4_offset_arc_angle = Hi_pass_arm4_offset / inner_ring_radius

    Hi_pass_arm5_linewidth = config["Hi_pass_arm5_linewidth"]  # 22
    Hi_pass_arm5_height = config["Hi_pass_arm5_height"]  # 23
    Hi_pass_arm5_length = config["Hi_pass_arm5_length"]  # 120.5
    Hi_pass_arm5_offset = config["Hi_pass_arm5_offset"]  # 160
    Hi_pass_arm5_offset_arc_angle = Hi_pass_arm5_offset / inner_ring_radius

    via_overlap_length = config["via_overlap_length"]  # 5
    via_extend_length = config["via_extend_length"]  # 5

    if direction == "clockwise":
        sign = -1
    elif direction == "anti-clockwise":
        sign = +1

    arm_linewidths = [
        Hi_pass_arm1_linewidth,
        Hi_pass_arm2_linewidth,
        Hi_pass_arm3_linewidth,
        Hi_pass_arm4_linewidth,
        Hi_pass_arm5_linewidth,
    ]
    arm_heights = [Hi_pass_arm1_height, Hi_pass_arm2_height, Hi_pass_arm3_height, Hi_pass_arm4_height, Hi_pass_arm5_height]
    arm_lengths = [Hi_pass_arm1_length, Hi_pass_arm2_length, Hi_pass_arm3_length, Hi_pass_arm4_length, Hi_pass_arm5_length]
    offset_arc_angles = [
        Hi_pass_arm1_offset_arc_angle,
        Hi_pass_arm2_offset_arc_angle,
        Hi_pass_arm3_offset_arc_angle,
        Hi_pass_arm4_offset_arc_angle,
        Hi_pass_arm5_offset_arc_angle,
    ]
    offset_arc_cumsum_angles = np.cumsum(offset_arc_angles)

    configurator_points = {}

    for i in range(5):
        angle = init_angle + sign * offset_arc_cumsum_angles[i]
        x_pos = x + (inner_ring_radius * cos(angle))
        y_pos = y + (inner_ring_radius * sin(angle))

        arm_geom_points = [
            [(sign) * arm_linewidths[i] / 2, 0],
            [(sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
            [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2)],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i],
                arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
            ],
            [(-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i], arm_heights[i] + (inner_ring_line_width / 2)],
            [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
            [(-sign) * arm_linewidths[i] / 2, 0],
        ]

        arm_geom = mbu.rotate_and_move_points_list(arm_geom_points, angle - pi / 2, x_pos, y_pos)

        configurator_points[f"Hi_pass_arm{i+1}_linewidth"] = {
            "text": f"Hi_pass_arm{i+1}_linewidth",
            "start": [arm_geom[2][0], arm_geom[2][1]],
            "end": [arm_geom[5][0], arm_geom[5][1]],
        }

        configurator_points[f"Hi_pass_arm{i+1}_length"] = {
            "text": f"Hi_pass_arm{i+1}_length",
            "start": [arm_geom[4][0], arm_geom[4][1]],
            "end": [arm_geom[5][0], arm_geom[5][1]],
        }

        configurator_points[f"Hi_pass_arm{i+1}_height"] = {
            "text": f"Hi_pass_arm{i+1}_height",
            "start": [arm_geom[5][0], arm_geom[5][1]],
            "end": [arm_geom[6][0], arm_geom[6][1]],
        }

        arm = gdspy.Polygon(
            arm_geom,
            layer=mask_builder.layers.Nb_Antenna.number,
            datatype=mask_builder.layers.Nb_Antenna.datatype,
        )
        mask_builder.Main.add(arm)

        via_points = [
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
            ],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                arm_heights[i] + (inner_ring_line_width / 2),
            ],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                arm_heights[i] + (inner_ring_line_width / 2),
            ],
            [
                (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
            ],
        ]

        via_geom = mbu.rotate_and_move_points_list(via_points, angle - pi / 2, x_pos, y_pos)

        configurator_points[f"via_extend_length_{i+1}"] = {
            "text": "via_extend_length",
            "start": [via_geom[0][0] + (via_geom[3][0] - via_geom[0][0]) / 2, via_geom[0][1] + (via_geom[3][1] - via_geom[0][1]) / 2],
            "end": [via_geom[3][0], via_geom[3][1]],
        }

        configurator_points[f"via_overlap_length_{i+1}"] = {
            "text": "via_overlap_length",
            "start": [via_geom[0][0] + (via_geom[3][0] - via_geom[0][0]) / 2, via_geom[0][1] + (via_geom[3][1] - via_geom[0][1]) / 2],
            "end": [via_geom[0][0], via_geom[0][1]],
        }

        via_box = gdspy.Polygon(
            via_geom,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(via_box)

    if not return_configurator_points:
        return

    configurator_points[f"Hi_pass_arm1_offset"] = {
        "text": f"Hi_pass_arm1_offset",
        "start": [x + inner_ring_radius * cos(init_angle), y + inner_ring_radius * sin(init_angle)],
        "end": [
            x + inner_ring_radius * cos(init_angle + (sign * offset_arc_cumsum_angles[0])),
            y + inner_ring_radius * sin(init_angle + (sign * offset_arc_cumsum_angles[0])),
        ],
    }
    for i in range(1, 5):
        configurator_points[f"Hi_pass_arm{i+1}_offset"] = {
            "text": f"Hi_pass_arm{i+1}_offset",
            "start": [
                x + inner_ring_radius * cos(init_angle + (sign * offset_arc_cumsum_angles[i - 1])),
                y + inner_ring_radius * sin(init_angle + (sign * offset_arc_cumsum_angles[i - 1])),
            ],
            "end": [
                x + inner_ring_radius * cos(init_angle + (sign * offset_arc_cumsum_angles[i])),
                y + inner_ring_radius * sin(init_angle + (sign * offset_arc_cumsum_angles[i])),
            ],
        }

    return configurator_points


def add_filter_bank_ring_overlap_and_get_conections(
    mask_builder,
    x: float | int,
    y: float | int,
    ant_center_x: float | int,
    ant_center_y: float | int,
    overlap_no: int,
    rot: float | int,
    filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None,
    return_configurator_points: bool = False,
):
    """Adds a ring overlap conection to bridge the inner and outer rings
    over one another. This function will return the conection points where
    the inner and outer rings should connect to.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate of the center of the ring overlap.

    ant_center_x, ant_center_y: float, int
        The x, y coordinate of the center of the antenna structure,
        i.e. the center of the horn.

    overlap_no: int
        The number of the ring overlap. This determines where to draw it around
        the antenna. Starting at 0 for left middle and +1 for each subsequent
        overlap going anti-clockwise. Should not be more than 3, values more
        than this wrap back to 0 (left middle placement) because the overlap_no
        operates like it is modulo 4.

    rot: float | int
        The angle (**in radians**) which the overlap geometry should be rotated
        at. This rot angle defined as the anti-clockwise angle made with the
        positive x-axis.

    KwArgs
    ------
    filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    return_configurator_points: bool = False
        return a the points for use in the configurator.

    Returns
    -------
    conections_dict: dict
        This dictionary contains keys that map to an [x,y] list which are the
        coordinates defining the conection points where the inner and outer
        ring should connect to this overlap structure.

        This dict has keys: **'inner_conect_0', 'inner_conect_1',
        'outer_conect_0', 'outer_conect_1'**.
    """
    config = get_mask_default_config(SoukMaskConfig.FILTER_BANK_RING_OVERLAP, config_override=filter_bank_ring_overlap_config_override)

    outer_box_width = config["outer_box_width"]  # 100
    outer_box_height = config["outer_box_height"]  # 100
    outer_box_inner_cutout_width = config["outer_box_inner_cutout_width"]  # 16
    outer_box_inner_cutout_height = config["outer_box_inner_cutout_height"]  # 6
    linewidth = config["linewidth"]  # 5

    half_conect_taper_extra_height = config["half_conect_taper_extra_height"]  # 5
    half_conect_taper_stright_width = config["half_conect_taper_stright_width"]  # 5
    half_conect_taper_diag_width = config["half_conect_taper_diag_width"]  # 5.5
    half_conect_taper_inner_height = config["half_conect_taper_inner_height"]  # 4
    half_conect_distacne_from_center = config["half_conect_distacne_from_center"]  # 3

    half_conect_bridge_rect_width = config["half_conect_bridge_rect_width"]  # 14
    half_conect_bridge_rect_height = config["half_conect_bridge_rect_height"]  # 4

    half_conect_bridge_pad_width = config["half_conect_bridge_pad_width"]  # 2
    half_conect_bridge_pad_height = config["half_conect_bridge_pad_height"]  # 2
    half_conect_bridge_pad_offset_from_center = config["half_conect_bridge_pad_offset_from_center"]  # 5

    full_conect_taper_extra_width = config["full_conect_taper_extra_width"]  # 5
    full_conect_taper_stright_height = config["full_conect_taper_stright_height"]  # 5
    full_conect_taper_diag_height = config["full_conect_taper_diag_height"]  # 5.5
    full_conect_taper_start_from_center = config["full_conect_taper_start_from_center"]  # 3
    full_conect_center_width = config["full_conect_center_width"]  # 4

    conect_half_conect_left = mbu.rotate_and_move_single_point([-outer_box_width / 2, 0], rot, x, y)
    conect_half_conect_right = mbu.rotate_and_move_single_point([outer_box_width / 2, 0], rot, x, y)
    conect_full_conect_bot = mbu.rotate_and_move_single_point([0, -outer_box_height / 2], rot, x, y)
    conect_full_conect_top = mbu.rotate_and_move_single_point([0, outer_box_height / 2], rot, x, y)

    conections_rotated_to_quadrant_1 = [
        mbu.rotate(ant_center_x, ant_center_y, conect_half_conect_left[0], conect_half_conect_left[1], (overlap_no * -pi / 2)),
        mbu.rotate(ant_center_x, ant_center_y, conect_half_conect_right[0], conect_half_conect_right[1], (overlap_no * -pi / 2)),
        mbu.rotate(ant_center_x, ant_center_y, conect_full_conect_bot[0], conect_full_conect_bot[1], (overlap_no * -pi / 2)),
        mbu.rotate(ant_center_x, ant_center_y, conect_full_conect_top[0], conect_full_conect_top[1], (overlap_no * -pi / 2)),
    ]

    conections_sorted_list = sorted(conections_rotated_to_quadrant_1.copy(), key=lambda k: [k[1], k[0]])

    conections_dict = {}
    conections_dict["inner_conect_0"] = mbu.rotate(
        ant_center_x, ant_center_y, conections_sorted_list[0][0], conections_sorted_list[0][1], (overlap_no * pi / 2)
    )
    conections_dict["outer_conect_0"] = mbu.rotate(
        ant_center_x, ant_center_y, conections_sorted_list[1][0], conections_sorted_list[1][1], (overlap_no * pi / 2)
    )
    conections_dict["inner_conect_1"] = mbu.rotate(
        ant_center_x, ant_center_y, conections_sorted_list[2][0], conections_sorted_list[2][1], (overlap_no * pi / 2)
    )
    conections_dict["outer_conect_1"] = mbu.rotate(
        ant_center_x, ant_center_y, conections_sorted_list[3][0], conections_sorted_list[3][1], (overlap_no * pi / 2)
    )

    outer_box_poly_points = [
        [-outer_box_width / 2, -outer_box_height / 2],
        [-outer_box_width / 2, outer_box_height / 2],
        [outer_box_width / 2, outer_box_height / 2],
        [outer_box_width / 2, outer_box_inner_cutout_height / 2],
        [-outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
        [-outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
        [outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
        [outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
        [outer_box_width / 2, outer_box_inner_cutout_height / 2],
        [outer_box_width / 2, -outer_box_height / 2],
    ]

    outer_box_poly_points = mbu.rotate_and_move_points_list(outer_box_poly_points, rot, x, y)
    outer_box = gdspy.Polygon(
        outer_box_poly_points, layer=mask_builder.layers.Aluminium.number, datatype=mask_builder.layers.Aluminium.datatype
    )
    # self.Main.add(outer_box)

    outer_box_inner_cutout_points = [
        [-outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
        [-outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
        [outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
        [outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
    ]

    outer_box_inner_cutout_points = mbu.rotate_and_move_points_list(outer_box_inner_cutout_points, rot, x, y)

    outer_box_inner_cutout = gdspy.Polygon(
        outer_box_inner_cutout_points,
        layer=mask_builder.layers.Nb_Groundplane.number,
        datatype=mask_builder.layers.Nb_Groundplane.datatype,
    )
    mask_builder.ground_plane_cutouts.add(outer_box_inner_cutout)

    half_conect_poly_points_right = [
        [half_conect_distacne_from_center, half_conect_taper_inner_height / 2],
        [half_conect_distacne_from_center + half_conect_taper_diag_width, linewidth / 2 + half_conect_taper_extra_height],
        [
            half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width,
            linewidth / 2 + half_conect_taper_extra_height,
        ],
        [half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width, linewidth / 2],
        [outer_box_width / 2, linewidth / 2],
        [outer_box_width / 2, -linewidth / 2],
        [half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width, -linewidth / 2],
        [
            half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width,
            -linewidth / 2 - half_conect_taper_extra_height,
        ],
        [half_conect_distacne_from_center + half_conect_taper_diag_width, -linewidth / 2 - half_conect_taper_extra_height],
        [half_conect_distacne_from_center, -half_conect_taper_inner_height / 2],
    ]

    half_conect_poly_points_right = mbu.rotate_and_move_points_list(half_conect_poly_points_right, rot, x, y)
    half_conect_right = gdspy.Polygon(
        half_conect_poly_points_right,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(half_conect_right)

    half_conect_poly_points_left = [
        [-half_conect_distacne_from_center, half_conect_taper_inner_height / 2],
        [-half_conect_distacne_from_center - half_conect_taper_diag_width, linewidth / 2 + half_conect_taper_extra_height],
        [
            -half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width,
            linewidth / 2 + half_conect_taper_extra_height,
        ],
        [-half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width, linewidth / 2],
        [-outer_box_width / 2, linewidth / 2],
        [-outer_box_width / 2, -linewidth / 2],
        [-half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width, -linewidth / 2],
        [
            -half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width,
            -linewidth / 2 - half_conect_taper_extra_height,
        ],
        [-half_conect_distacne_from_center - half_conect_taper_diag_width, -linewidth / 2 - half_conect_taper_extra_height],
        [-half_conect_distacne_from_center, -half_conect_taper_inner_height / 2],
    ]

    half_conect_poly_points_left = mbu.rotate_and_move_points_list(half_conect_poly_points_left, rot, x, y)
    half_conect_left = gdspy.Polygon(
        half_conect_poly_points_left,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(half_conect_left)

    half_conect_bridge_rect = gdspy.Rectangle(
        [-half_conect_bridge_rect_width / 2, -half_conect_bridge_rect_height / 2],
        [half_conect_bridge_rect_width / 2, half_conect_bridge_rect_height / 2],
        layer=mask_builder.layers.Nb_Groundplane.number,
        datatype=mask_builder.layers.Nb_Groundplane.datatype,
    )
    half_conect_bridge_rect.rotate(rot, (0, 0))
    half_conect_bridge_rect.translate(x, y)
    mask_builder.Main.add(half_conect_bridge_rect)

    half_conect_bridge_pad_right = gdspy.Rectangle(
        [half_conect_bridge_pad_offset_from_center, -half_conect_bridge_pad_height / 2],
        [half_conect_bridge_pad_offset_from_center + half_conect_bridge_pad_width, half_conect_bridge_pad_height / 2],
        layer=mask_builder.layers.SiN_dep.number,
        datatype=mask_builder.layers.SiN_dep.datatype,
    )
    half_conect_bridge_pad_right.rotate(rot, (0, 0))
    half_conect_bridge_pad_right.translate(x, y)
    # self.Main.add(half_conect_bridge_pad_right)
    mask_builder.silicon_nitride_cutouts.add(half_conect_bridge_pad_right)

    half_conect_bridge_pad_left = gdspy.Rectangle(
        [half_conect_bridge_pad_offset_from_center, -half_conect_bridge_pad_height / 2],
        [half_conect_bridge_pad_offset_from_center + half_conect_bridge_pad_width, half_conect_bridge_pad_height / 2],
        layer=mask_builder.layers.SiN_dep.number,
        datatype=mask_builder.layers.SiN_dep.datatype,
    )
    half_conect_bridge_pad_left.rotate(pi, (0, 0))
    half_conect_bridge_pad_left.rotate(rot, (0, 0))
    half_conect_bridge_pad_left.translate(x, y)
    # self.Main.add(half_conect_bridge_pad_left)
    mask_builder.silicon_nitride_cutouts.add(half_conect_bridge_pad_left)

    full_conect_poly_points = [
        [full_conect_center_width / 2, full_conect_taper_start_from_center],
        [linewidth / 2 + full_conect_taper_extra_width, full_conect_taper_start_from_center + full_conect_taper_diag_height],
        [
            linewidth / 2 + full_conect_taper_extra_width,
            full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height,
        ],
        [linewidth / 2, full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height],
        [linewidth / 2, outer_box_height / 2],
        [-linewidth / 2, outer_box_height / 2],
        [-linewidth / 2, full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height],
        [
            -linewidth / 2 - full_conect_taper_extra_width,
            full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height,
        ],
        [-linewidth / 2 - full_conect_taper_extra_width, full_conect_taper_start_from_center + full_conect_taper_diag_height],
        [-full_conect_center_width / 2, full_conect_taper_start_from_center],
        [-full_conect_center_width / 2, -full_conect_taper_start_from_center],
        [-linewidth / 2 - full_conect_taper_extra_width, -full_conect_taper_start_from_center - full_conect_taper_diag_height],
        [
            -linewidth / 2 - full_conect_taper_extra_width,
            -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height,
        ],
        [-linewidth / 2, -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height],
        [-linewidth / 2, -outer_box_height / 2],
        [linewidth / 2, -outer_box_height / 2],
        [linewidth / 2, -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height],
        [
            linewidth / 2 + full_conect_taper_extra_width,
            -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height,
        ],
        [linewidth / 2 + full_conect_taper_extra_width, -full_conect_taper_start_from_center - full_conect_taper_diag_height],
        [full_conect_center_width / 2, -full_conect_taper_start_from_center],
    ]

    full_conect_poly_points = mbu.rotate_and_move_points_list(full_conect_poly_points, rot, x, y)

    full_conect = gdspy.Polygon(
        full_conect_poly_points,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(full_conect)

    if not return_configurator_points:
        return conections_dict

    configurator_points = {}

    configurator_points["outer_box_width"] = {
        "text": "outer_box_width",
        "start": [x - (outer_box_width / 2), y + (outer_box_height / 2)],
        "end": [x + (outer_box_width / 2), y + (outer_box_height / 2)],
    }

    configurator_points["outer_box_height"] = {
        "text": "outer_box_height",
        "start": [x - (outer_box_width / 2), y + (outer_box_height / 2)],
        "end": [x - (outer_box_width / 2), y - (outer_box_height / 2)],
    }

    configurator_points["outer_box_inner_cutout_width"] = {
        "text": "outer_box_inner_cutout_width",
        "start": [outer_box_inner_cutout_points[2][0], outer_box_inner_cutout_points[2][1]],
        "end": [outer_box_inner_cutout_points[1][0], outer_box_inner_cutout_points[1][1]],
    }

    configurator_points["outer_box_inner_cutout_height"] = {
        "text": "outer_box_inner_cutout_height",
        "start": [outer_box_inner_cutout_points[0][0], outer_box_inner_cutout_points[0][1]],
        "end": [outer_box_inner_cutout_points[1][0], outer_box_inner_cutout_points[1][1]],
    }

    configurator_points["linewidth"] = {
        "text": "linewidth",
        "start": [full_conect_poly_points[3][0], (full_conect_poly_points[3][1] + full_conect_poly_points[4][1]) / 2],
        "end": [full_conect_poly_points[6][0], (full_conect_poly_points[6][1] + full_conect_poly_points[5][1]) / 2],
    }

    configurator_points["linewidth_B"] = {
        "text": "linewidth",
        "start": [(half_conect_poly_points_left[4][0] + half_conect_poly_points_left[3][0]) / 2, half_conect_poly_points_left[3][1]],
        "end": [(half_conect_poly_points_left[5][0] + half_conect_poly_points_left[6][0]) / 2, half_conect_poly_points_left[5][1]],
    }

    configurator_points["half_conect_taper_extra_height"] = {
        "text": "half_conect_taper_extra_height",
        "start": [half_conect_poly_points_left[2][0], half_conect_poly_points_left[2][1]],
        "end": [half_conect_poly_points_left[3][0], half_conect_poly_points_left[3][1]],
    }

    configurator_points["half_conect_taper_stright_width"] = {
        "text": "half_conect_taper_stright_width",
        "start": [half_conect_poly_points_left[1][0], half_conect_poly_points_left[1][1]],
        "end": [half_conect_poly_points_left[2][0], half_conect_poly_points_left[2][1]],
    }

    configurator_points["half_conect_taper_diag_width"] = {
        "text": "half_conect_taper_diag_width",
        "start": [half_conect_poly_points_left[0][0], half_conect_poly_points_left[1][1]],
        "end": [half_conect_poly_points_left[1][0], half_conect_poly_points_left[1][1]],
    }

    configurator_points["half_conect_taper_inner_height"] = {
        "text": "half_conect_taper_inner_height",
        "start": [half_conect_poly_points_left[0][0], half_conect_poly_points_left[0][1]],
        "end": [half_conect_poly_points_left[-1][0], half_conect_poly_points_left[-1][1]],
    }

    configurator_points["half_conect_distacne_from_center"] = {
        "text": "half_conect_distacne_from_center",
        "start": [x, (half_conect_poly_points_left[0][1] + half_conect_poly_points_left[-1][1]) / 2],
        "end": [half_conect_poly_points_left[0][0], (half_conect_poly_points_left[0][1] + half_conect_poly_points_left[-1][1]) / 2],
    }

    configurator_points["half_conect_bridge_rect_width"] = {
        "text": "half_conect_bridge_rect_width",
        "start": [-half_conect_bridge_rect_width / 2, -half_conect_bridge_rect_height / 2],
        "end": [half_conect_bridge_rect_width / 2, -half_conect_bridge_rect_height / 2],
    }

    configurator_points["half_conect_bridge_rect_height"] = {
        "text": "half_conect_bridge_rect_height",
        "start": [half_conect_bridge_rect_width / 2, -half_conect_bridge_rect_height / 2],
        "end": [half_conect_bridge_rect_width / 2, +half_conect_bridge_rect_height / 2],
    }

    configurator_points["half_conect_bridge_pad_width"] = {
        "text": "half_conect_bridge_pad_width",
        "start": [half_conect_bridge_pad_offset_from_center, -half_conect_bridge_pad_height / 2],
        "end": [half_conect_bridge_pad_offset_from_center + half_conect_bridge_pad_width, -half_conect_bridge_pad_height / 2],
    }

    configurator_points["half_conect_bridge_pad_height"] = {
        "text": "half_conect_bridge_pad_height",
        "start": [half_conect_bridge_pad_offset_from_center, -half_conect_bridge_pad_height / 2],
        "end": [half_conect_bridge_pad_offset_from_center, +half_conect_bridge_pad_height / 2],
    }

    configurator_points["half_conect_bridge_pad_offset_from_center"] = {
        "text": "half_conect_bridge_pad_offset_from_center",
        "start": [x, y],
        "end": [half_conect_bridge_pad_offset_from_center, y],
    }

    configurator_points["full_conect_taper_extra_width"] = {
        "text": "full_conect_taper_extra_width",
        "start": [full_conect_poly_points[3][0], full_conect_poly_points[3][1]],
        "end": [full_conect_poly_points[2][0], full_conect_poly_points[2][1]],
    }

    configurator_points["full_conect_taper_stright_height"] = {
        "text": "full_conect_taper_stright_height",
        "start": [full_conect_poly_points[1][0], full_conect_poly_points[1][1]],
        "end": [full_conect_poly_points[2][0], full_conect_poly_points[2][1]],
    }

    configurator_points["full_conect_taper_diag_height"] = {
        "text": "full_conect_taper_diag_height",
        "start": [full_conect_poly_points[1][0], full_conect_poly_points[0][1]],
        "end": [full_conect_poly_points[1][0], full_conect_poly_points[1][1]],
    }

    configurator_points["full_conect_taper_start_from_center"] = {
        "text": "full_conect_taper_start_from_center",
        "start": [x, full_conect_poly_points[0][1]],
        "end": [x, y],
    }

    configurator_points["full_conect_center_width"] = {
        "text": "full_conect_center_width",
        "start": [full_conect_poly_points[0][0], full_conect_poly_points[0][1]],
        "end": [full_conect_poly_points[9][0], full_conect_poly_points[9][1]],
    }

    return conections_dict, configurator_points


def add_filter_bank_and_get_conection_points(
    mask_builder,
    x: float | int,
    y: float | int,
    filter_bank_config_override: dict[str, float | int] | None = None,
    filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None,
    Hi_pass_filters_config_override: dict[str, float | int] | None = None,
    Lo_pass_filters_config_override: dict[str, float | int] | None = None,
    combiner_section_90ghz_config_override: dict[str, float | int] | None = None,
    combiner_section_150ghz_config_override: dict[str, float | int] | None = None,
    with_combiner: bool = True,
    with_crossover: bool = True,
    only_1_pol: bool = False,
    return_configurator_points: bool = False,
    return_configurator_points_for_Lo_pass: bool = False,
    return_configurator_points_for_Hi_pass: bool = False,
    return_configurator_points_for_combiner_150ghz: bool = False,
    return_configurator_points_for_combiner_90ghz: bool = False,
    **kwargs,
):
    """Adds the filter bank structure to the chip centered at the x,y
    coordinate given. By default this will be drawn with phase combiners
    and ring overlap crossovers and with both polorizations filtered.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center filter bank structure around.

    KwArgs
    ------
    filter_bank_config_override,
    Hi_pass_filters_config_override,
    Lo_pass_filters_config_override,
    combiner_section_90ghz_config_override,
    combiner_section_150ghz_config_override,
    filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    with_combiner = True
        Adds the combiner along with the filter bank rings by default. When
        Fasle the outer ring has no gap in it and the conection points will be
        the base of where the combiners would have gone, i.e. the middle of the
        outer ring arc sections.

    with_crossover = True
        Adds the crossover between inner and outer rings of the filter bank by
        default. When Fasle, the inner is only added for the first and 3rd
        antenna leaving the other two to conect to nothing. The outer ring is
        removed for the first and second antenna. The Hugh filters remain.

    only_1_pol = False
        If set True this will only draw one polorization of the filter bank.
        This will not draw the relevant inner and outer rings for one
        polorization and adds an inductive meander to absorb energy from the
        unused antenna polorization.

        **Note**:
        **There will still be some part of the structure that wont be used left
        on the mask**, namely the Hugh filters. These extra parts will have to
        be **removed manually**.

    return_configurator_points=False
        return a the points for use in the configurator.

    return_configurator_points_for_Lo_pass=False
        return a the points for use in the configurator.

    return_configurator_points_for_Hi_pass=False
        return a the points for use in the configurator.

    return_configurator_points_for_combiner_150ghz=False
        return a the points for use in the configurator.

    return_configurator_points_for_combiner_90ghz=False
        return a the points for use in the configurator.

    Returns
    -------
    points_to_conect_kids_dict: dict
        This dictionary contains keys that map to an [x,y] list which are the
        coordinates defining the conection points where the KIDs should connect
        to the filter bank.

        This dict has keys: **'TR', 'TL', 'BL', 'BR'**.
    """
    filter_bank_config = get_mask_default_config(SoukMaskConfig.FILTER_BANK, config_override=filter_bank_config_override)

    # ring dimension and properties
    inner_ring_radius = filter_bank_config["inner_ring_radius"]  # 3920/2
    inner_ring_line_width = filter_bank_config["inner_ring_line_width"]  # 5
    outer_ring_radius = filter_bank_config["outer_ring_radius"]  # 4315/2
    outer_ring_line_width = filter_bank_config["outer_ring_line_width"]  # 5

    ring_overlap_distance_from_center = filter_bank_config["ring_overlap_distance_from_center"]  # inner_ring_radius + 150

    inner_ring_overlap_gap = filter_bank_config["inner_ring_overlap_gap"]  # 325
    outer_ring_overlap_gap = filter_bank_config["outer_ring_overlap_gap"]  # 325
    if with_combiner:
        outer_ring_conector_gap = filter_bank_config["outer_ring_conector_gap"]  # 20
    else:
        outer_ring_conector_gap = 0

    inner_ring_arc_length = inner_ring_overlap_gap
    inner_ring_arc_angle = inner_ring_arc_length / inner_ring_radius

    outer_ring_arc_length = outer_ring_overlap_gap
    outer_ring_arc_angle = outer_ring_arc_length / outer_ring_radius

    outer_combiner_arc_length = outer_ring_conector_gap
    outer_combiner_arc_angle = outer_combiner_arc_length / outer_ring_radius

    inner_ring_overlap_theta = inner_ring_arc_angle / 2
    outer_ring_overlap_theta = outer_ring_arc_angle / 2
    outer_combiner_theta = outer_combiner_arc_angle / 2

    inner_straight_extend_distance = filter_bank_config["inner_straight_extend_distance"]  # 20
    inner_bend_radius = filter_bank_config["inner_bend_radius"]  # 15
    outer_straight_extend_distance = filter_bank_config["outer_straight_extend_distance"]  # 10
    outer_bend_radius = filter_bank_config["outer_bend_radius"]  # 9

    Hugh_filter_thin_length1 = filter_bank_config["Hugh_filter_thin_length1"]  # 51
    Hugh_filter_thin_length2 = filter_bank_config["Hugh_filter_thin_length2"]  # 106

    Hugh_filter_thick_length1 = filter_bank_config["Hugh_filter_thick_length1"]  # 48
    Hugh_filter_thick_length2 = filter_bank_config["Hugh_filter_thick_length2"]  # 58.5

    Hugh_filter_thin_lengths = np.array([Hugh_filter_thin_length1, Hugh_filter_thin_length2])
    Hugh_filter_thin_arc_angles = Hugh_filter_thin_lengths / outer_ring_radius
    Hugh_filter_thick_lengths = np.array([Hugh_filter_thick_length1, Hugh_filter_thick_length2])
    Hugh_filter_thick_arc_angles = Hugh_filter_thick_lengths / outer_ring_radius
    Hugh_filter_thin_width = filter_bank_config["Hugh_filter_thin_width"]  # 2
    Hugh_filter_thick_width = filter_bank_config["Hugh_filter_thick_width"]  # 20
    half_Hugh_filter_total_arc_angle = (
        np.sum(Hugh_filter_thin_arc_angles) + Hugh_filter_thick_arc_angles[0] + (Hugh_filter_thick_arc_angles[1] / 2)
    )

    ring_overlap_0_rot = filter_bank_config["ring_overlap_0_rot"]
    ring_overlap_1_rot = filter_bank_config["ring_overlap_1_rot"]
    ring_overlap_2_rot = filter_bank_config["ring_overlap_2_rot"]
    ring_overlap_3_rot = filter_bank_config["ring_overlap_3_rot"]

    ring_overlap_rotations = [
        mbu.deg_to_rad(ring_overlap_0_rot),
        mbu.deg_to_rad(ring_overlap_1_rot),
        mbu.deg_to_rad(ring_overlap_2_rot),
        mbu.deg_to_rad(ring_overlap_3_rot),
    ]

    points_to_conect_kids_dict = {}
    conect_kids_dict_key_names = ["TR", "TL", "BL", "BR"]

    for i in range(4):
        if not only_1_pol or i in [0, 2]:
            inner_arc_section = gdspy.Round(
                [x, y],
                inner_ring_radius + (inner_ring_line_width / 2),
                inner_ring_radius - (inner_ring_line_width / 2),
                initial_angle=inner_ring_overlap_theta,
                final_angle=inner_ring_overlap_theta + (pi / 2) - inner_ring_arc_angle,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.Main.add(inner_arc_section)

        if i in [0, 1] and with_crossover:
            outer_arc_section1 = gdspy.Round(
                [x, y],
                outer_ring_radius + (outer_ring_line_width / 2),
                outer_ring_radius - (outer_ring_line_width / 2),
                initial_angle=outer_ring_overlap_theta,
                final_angle=(i * pi / 2) + (pi / 4) - outer_combiner_theta,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            outer_arc_section2 = gdspy.Round(
                [x, y],
                outer_ring_radius + (outer_ring_line_width / 2),
                outer_ring_radius - (outer_ring_line_width / 2),
                initial_angle=(i * pi / 2) + (pi / 4) + outer_combiner_theta,
                final_angle=outer_ring_overlap_theta + (pi / 2) - outer_ring_arc_angle,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.Main.add(outer_arc_section1)
            mask_builder.Main.add(outer_arc_section2)

        if i in [2, 3]:
            if (not only_1_pol) or (with_crossover):
                outer_arc_section11 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (outer_ring_line_width / 2),
                    outer_ring_radius - (outer_ring_line_width / 2),
                    initial_angle=outer_ring_overlap_theta,
                    final_angle=(i * pi / 2) + (pi / 8) - (half_Hugh_filter_total_arc_angle),
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                outer_arc_section12 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (outer_ring_line_width / 2),
                    outer_ring_radius - (outer_ring_line_width / 2),
                    initial_angle=(i * pi / 2) + (pi / 8) + (half_Hugh_filter_total_arc_angle),
                    final_angle=(i * pi / 2) + (pi / 4) - outer_combiner_theta,
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                mask_builder.Main.add(outer_arc_section11)
                mask_builder.Main.add(outer_arc_section12)

            outer_arc_section21 = gdspy.Round(
                [x, y],
                outer_ring_radius + (outer_ring_line_width / 2),
                outer_ring_radius - (outer_ring_line_width / 2),
                initial_angle=(i * pi / 2) + (pi / 4) + outer_combiner_theta,
                final_angle=(i * pi / 2) + (pi / 8) + (pi / 4) - (half_Hugh_filter_total_arc_angle),
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            outer_arc_section22 = gdspy.Round(
                [x, y],
                outer_ring_radius + (outer_ring_line_width / 2),
                outer_ring_radius - (outer_ring_line_width / 2),
                initial_angle=(i * pi / 2) + (pi / 8) + (pi / 4) + (half_Hugh_filter_total_arc_angle),
                final_angle=outer_ring_overlap_theta + (pi / 2) - outer_ring_arc_angle,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.Main.add(outer_arc_section21)
            mask_builder.Main.add(outer_arc_section22)

            for j in range(2):
                mid_Hugh_angle = pi / 8 if j == 0 else pi / 8 + pi / 4
                angle_inc = 0

                if (only_1_pol and not with_crossover) and j == 0:
                    if i != 3:
                        continue

                    # add arc to hugh from lopass bot left
                    inner_conect_ang = inner_ring_overlap_theta - inner_ring_arc_angle
                    outer_conect_ang = outer_ring_overlap_theta - outer_ring_arc_angle
                    extra_length_connect_arc = 200
                    conection_low_pass_to_hugh_points = [
                        [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                        [
                            x + inner_ring_radius * cos(inner_conect_ang) + extra_length_connect_arc * cos(inner_conect_ang + pi / 2),
                            y + inner_ring_radius * sin(inner_conect_ang) + extra_length_connect_arc * sin(inner_conect_ang + pi / 2),
                        ],
                        [
                            x + outer_ring_radius * cos(inner_conect_ang) + extra_length_connect_arc * cos(outer_conect_ang + pi / 2),
                            y + outer_ring_radius * sin(inner_conect_ang) + extra_length_connect_arc * sin(outer_conect_ang + pi / 2),
                        ],
                        [x + outer_ring_radius * cos(outer_conect_ang), y + outer_ring_radius * sin(outer_conect_ang)],
                    ]
                    conection_low_pass_to_hugh = gdspy.FlexPath(
                        conection_low_pass_to_hugh_points,
                        inner_ring_line_width,
                        corners="circular bend",
                        bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                        layer=mask_builder.layers.Nb_Antenna.number,
                        datatype=mask_builder.layers.Nb_Antenna.datatype,
                    )
                    mask_builder.make_flexpath_into_polygons_and_add_to_main(conection_low_pass_to_hugh, mask_builder.layers.Nb_Antenna)
                    BL_KID_connect_1_pol = [x + outer_ring_radius * cos(-3 * pi / 4), y + outer_ring_radius * sin(-3 * pi / 4)]

                    # add from hugh bot right
                    inner_conect_ang = -(pi / 4)
                    outer_conect_ang = -(pi / 4)
                    extra_length_connect_arc = 200
                    outer_arc_conection_radius = outer_ring_radius + (outer_ring_radius - inner_ring_radius)
                    conection_from_hugh_points = [
                        [x + outer_ring_radius * cos(inner_conect_ang), y + outer_ring_radius * sin(inner_conect_ang)],
                        [
                            x + outer_ring_radius * cos(inner_conect_ang) + extra_length_connect_arc * cos(inner_conect_ang - pi / 2),
                            y + outer_ring_radius * sin(inner_conect_ang) + extra_length_connect_arc * sin(inner_conect_ang - pi / 2),
                        ],
                        [
                            x
                            + outer_arc_conection_radius * cos(inner_conect_ang)
                            + extra_length_connect_arc * cos(outer_conect_ang - pi / 2),
                            y
                            + outer_arc_conection_radius * sin(inner_conect_ang)
                            + extra_length_connect_arc * sin(outer_conect_ang - pi / 2),
                        ],
                        [
                            x + outer_arc_conection_radius * cos(outer_conect_ang),
                            y + outer_arc_conection_radius * sin(outer_conect_ang),
                        ],
                    ]
                    conection_from_hugh = gdspy.FlexPath(
                        conection_from_hugh_points,
                        inner_ring_line_width,
                        corners="circular bend",
                        bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                        layer=mask_builder.layers.Nb_Antenna.number,
                        datatype=mask_builder.layers.Nb_Antenna.datatype,
                    )
                    mask_builder.make_flexpath_into_polygons_and_add_to_main(conection_from_hugh, mask_builder.layers.Nb_Antenna)
                    BR_KID_connect_1_pol = [
                        x + outer_arc_conection_radius * cos(-pi / 4),
                        y + outer_arc_conection_radius * sin(-pi / 4),
                    ]

                    # add path conecting right mid inner to outer
                    inner_conect_ang = (pi / 2) + inner_ring_overlap_theta
                    outer_conect_ang = (pi / 2) + outer_ring_overlap_theta - outer_ring_arc_angle
                    extra_length_connect_path = 100
                    inner_to_outer_connect_path_points = [
                        [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                        [
                            x + inner_ring_radius * cos(inner_conect_ang) + extra_length_connect_path * cos(inner_conect_ang - pi / 2),
                            y + inner_ring_radius * sin(inner_conect_ang) + extra_length_connect_path * sin(inner_conect_ang - pi / 2),
                        ],
                        [
                            x + outer_ring_radius * cos(outer_conect_ang) + extra_length_connect_path * cos(outer_conect_ang + pi / 2),
                            y + outer_ring_radius * sin(outer_conect_ang) + extra_length_connect_path * sin(outer_conect_ang + pi / 2),
                        ],
                        [x + outer_ring_radius * cos(outer_conect_ang), y + outer_ring_radius * sin(outer_conect_ang)],
                    ]

                    inner_to_outer_connect_path = gdspy.FlexPath(
                        inner_to_outer_connect_path_points,
                        inner_ring_line_width,
                        corners="circular bend",
                        bend_radius=(outer_ring_radius - inner_ring_radius),
                        layer=mask_builder.layers.Nb_Antenna.number,
                        datatype=mask_builder.layers.Nb_Antenna.datatype,
                    )
                    mask_builder.make_flexpath_into_polygons_and_add_to_main(inner_to_outer_connect_path, mask_builder.layers.Nb_Antenna)

                    # add arc from hipass top
                    inner_conect_ang = pi + inner_ring_overlap_theta - inner_ring_arc_angle
                    outer_conect_ang = pi + outer_ring_overlap_theta - outer_ring_arc_angle
                    extra_inner_length_connect_arc = 400
                    extra_outer_length_connect_arc_after = 800
                    conection_high_pass_points = [
                        [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                        [
                            x + inner_ring_radius * cos(inner_conect_ang) + extra_inner_length_connect_arc * cos(inner_conect_ang + pi / 2),
                            y + inner_ring_radius * sin(inner_conect_ang) + extra_inner_length_connect_arc * sin(inner_conect_ang + pi / 2),
                        ],
                        [x + outer_ring_radius * cos(outer_conect_ang), y + outer_ring_radius * sin(outer_conect_ang)],
                        [
                            x + outer_ring_radius * cos(outer_conect_ang) + extra_outer_length_connect_arc_after,
                            y + outer_ring_radius * sin(outer_conect_ang),
                        ],
                    ]
                    conection_high_pass = gdspy.FlexPath(
                        conection_high_pass_points,
                        inner_ring_line_width,
                        corners="circular bend",
                        bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                        layer=mask_builder.layers.Nb_Antenna.number,
                        datatype=mask_builder.layers.Nb_Antenna.datatype,
                    )
                    mask_builder.make_flexpath_into_polygons_and_add_to_main(conection_high_pass, mask_builder.layers.Nb_Antenna)
                    TR_KID_connect_1_pol = conection_high_pass_points[-1]

                    # add small section on the left mid toward top left
                    inner_conect_ang = pi + inner_ring_arc_angle / 2
                    extra_length_straight = 100

                    small_connection_left_mid_points = [
                        [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                        [
                            x + inner_ring_radius * cos(inner_conect_ang) + extra_length_straight * cos(inner_conect_ang - pi / 2),
                            y + inner_ring_radius * sin(inner_conect_ang) + extra_length_straight * sin(inner_conect_ang - pi / 2),
                        ],
                        [
                            x + inner_ring_radius * cos(inner_conect_ang) + extra_length_straight * cos(inner_conect_ang - pi / 2),
                            y
                            + inner_ring_radius * sin(inner_conect_ang)
                            + extra_length_straight * sin(inner_conect_ang - pi / 2)
                            + extra_length_straight,
                        ],
                    ]

                    small_connection_left_mid = gdspy.FlexPath(
                        small_connection_left_mid_points,
                        inner_ring_line_width,
                        corners="circular bend",
                        bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                        layer=mask_builder.layers.Nb_Antenna.number,
                        datatype=mask_builder.layers.Nb_Antenna.datatype,
                    )
                    mask_builder.make_flexpath_into_polygons_and_add_to_main(small_connection_left_mid, mask_builder.layers.Nb_Antenna)
                    TL_KID_connect_1_pol = small_connection_left_mid_points[-1]

                    mask_builder.Main.add(
                        gdspy.Round(
                            TL_KID_connect_1_pol,
                            10,
                            layer=mask_builder.layers.General_labeling.number,
                            datatype=mask_builder.layers.General_labeling.datatype,
                        )
                    )
                    mask_builder.Main.add(
                        gdspy.Round(
                            TR_KID_connect_1_pol,
                            10,
                            layer=mask_builder.layers.General_labeling.number,
                            datatype=mask_builder.layers.General_labeling.datatype,
                        )
                    )
                    mask_builder.Main.add(
                        gdspy.Round(
                            BL_KID_connect_1_pol,
                            10,
                            layer=mask_builder.layers.General_labeling.number,
                            datatype=mask_builder.layers.General_labeling.datatype,
                        )
                    )
                    mask_builder.Main.add(
                        gdspy.Round(
                            BR_KID_connect_1_pol,
                            10,
                            layer=mask_builder.layers.General_labeling.number,
                            datatype=mask_builder.layers.General_labeling.datatype,
                        )
                    )
                    one_pol_connect_dict = {
                        "TL": TL_KID_connect_1_pol,
                        "TR": TR_KID_connect_1_pol,
                        "BL": BL_KID_connect_1_pol,
                        "BR": BR_KID_connect_1_pol,
                    }

                    continue

                middle_Hugh_thick_arc = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thick_width / 2),
                    outer_ring_radius - (Hugh_filter_thick_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle - (Hugh_filter_thick_arc_angles[1] / 2),
                    final_angle=(i * pi / 2) + mid_Hugh_angle + (Hugh_filter_thick_arc_angles[1] / 2),
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                mask_builder.Main.add(middle_Hugh_thick_arc)
                angle_inc += Hugh_filter_thick_arc_angles[1] / 2

                middle_Hugh_thin_arc1 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thin_width / 2),
                    outer_ring_radius - (Hugh_filter_thin_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[1],
                    final_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc,
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                middle_Hugh_thin_arc2 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thin_width / 2),
                    outer_ring_radius - (Hugh_filter_thin_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc,
                    final_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc + Hugh_filter_thin_arc_angles[1],
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                mask_builder.Main.add(middle_Hugh_thin_arc1)
                mask_builder.Main.add(middle_Hugh_thin_arc2)
                angle_inc += Hugh_filter_thin_arc_angles[1]

                outer_Hugh_thick_arc1 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thick_width / 2),
                    outer_ring_radius - (Hugh_filter_thick_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0],
                    final_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc,
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                outer_Hugh_thick_arc2 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thick_width / 2),
                    outer_ring_radius - (Hugh_filter_thick_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc,
                    final_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc + Hugh_filter_thick_arc_angles[0],
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                mask_builder.Main.add(outer_Hugh_thick_arc1)
                mask_builder.Main.add(outer_Hugh_thick_arc2)
                angle_inc += Hugh_filter_thick_arc_angles[0]

                outer_Hugh_thin_arc1 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thin_width / 2),
                    outer_ring_radius - (Hugh_filter_thin_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0],
                    final_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc,
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                outer_Hugh_thin_arc2 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (Hugh_filter_thin_width / 2),
                    outer_ring_radius - (Hugh_filter_thin_width / 2),
                    initial_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc,
                    final_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc + Hugh_filter_thin_arc_angles[0],
                    layer=mask_builder.layers.Nb_Antenna.number,
                    datatype=mask_builder.layers.Nb_Antenna.datatype,
                )
                mask_builder.Main.add(outer_Hugh_thin_arc1)
                mask_builder.Main.add(outer_Hugh_thin_arc2)
                angle_inc += Hugh_filter_thin_arc_angles[0]

        if i == 0:
            mask_builder.add_Lo_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "clockwise",
                Lo_pass_filters_config_override=Lo_pass_filters_config_override,
            )
            configurator_points_for_Hi_pass = mask_builder.add_Hi_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "anti-clockwise",
                Hi_pass_filters_config_override=Hi_pass_filters_config_override,
                return_configurator_points=return_configurator_points_for_Hi_pass,
            )
        if i == 1 and not only_1_pol:
            mask_builder.add_Hi_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "clockwise",
                Hi_pass_filters_config_override=Hi_pass_filters_config_override,
            )
            configurator_points_for_Lo_pass = mask_builder.add_Lo_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "anti-clockwise",
                Lo_pass_filters_config_override=Lo_pass_filters_config_override,
                return_configurator_points=return_configurator_points_for_Lo_pass,
            )
        if i == 2:
            mask_builder.add_Hi_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "clockwise",
                Hi_pass_filters_config_override=Hi_pass_filters_config_override,
            )
            mask_builder.add_Lo_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "anti-clockwise",
                Lo_pass_filters_config_override=Lo_pass_filters_config_override,
            )
        if i == 3 and not only_1_pol:
            mask_builder.add_Lo_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "clockwise",
                Lo_pass_filters_config_override=Lo_pass_filters_config_override,
            )
            mask_builder.add_Hi_pass_filters(
                x,
                y,
                inner_ring_line_width,
                inner_ring_radius,
                (i * pi / 2) + pi / 4,
                "anti-clockwise",
                Hi_pass_filters_config_override=Hi_pass_filters_config_override,
            )

        combiner_xpos = x + (outer_ring_radius * cos(pi / 4 + i * (pi / 2)))
        combiner_ypos = y + (outer_ring_radius * sin(pi / 4 + i * (pi / 2)))

        conections_dict = {}

        if with_combiner:
            if i in [1, 3]:
                mirror_combiner = False
            else:
                mirror_combiner = True

            if i in [0, 1]:
                conection_point, configurator_points_for_combiner_150ghz = mask_builder.add_combiner_section_and_get_conect_point(
                    combiner_xpos,
                    combiner_ypos,
                    (pi / 4 + i * pi / 2),
                    outer_ring_conector_gap,
                    outer_ring_line_width,
                    combiner_section_150ghz_config_override=combiner_section_150ghz_config_override,
                    combiner_type="150GHZ",
                    mirror_combiner=mirror_combiner,
                    return_configurator_points=return_configurator_points_for_combiner_150ghz,
                )
            else:
                conection_point, configurator_points_for_combiner_90ghz = mask_builder.add_combiner_section_and_get_conect_point(
                    combiner_xpos,
                    combiner_ypos,
                    (pi / 4 + i * pi / 2),
                    outer_ring_conector_gap,
                    outer_ring_line_width,
                    combiner_section_90ghz_config_override=combiner_section_90ghz_config_override,
                    combiner_type="90GHZ",
                    mirror_combiner=mirror_combiner,
                    return_configurator_points=return_configurator_points_for_combiner_90ghz,
                )
        else:
            conection_point = [
                x + outer_ring_radius * cos((i * pi / 2) + (pi / 4) - outer_combiner_theta),
                y + outer_ring_radius * sin((i * pi / 2) + (pi / 4) - outer_combiner_theta),
            ]

        points_to_conect_kids_dict[conect_kids_dict_key_names[i]] = conection_point

        if with_crossover:
            ring_overlap_xpos = x + (ring_overlap_distance_from_center * cos(i * (pi / 2)))
            ring_overlap_ypos = y + (ring_overlap_distance_from_center * sin(i * (pi / 2)))
            ring_overlap_number = i

            conections_dict = mask_builder.add_filter_bank_ring_overlap_and_get_conections(
                ring_overlap_xpos,
                ring_overlap_ypos,
                x,
                y,
                ring_overlap_number,
                ring_overlap_rotations[i],
                filter_bank_ring_overlap_config_override=filter_bank_ring_overlap_config_override,
            )

            # inner conections to ring overlap
            inner_0_conection_points = [
                conections_dict["inner_conect_0"],
                [
                    conections_dict["inner_conect_0"][0] + inner_straight_extend_distance * cos(-3 * pi / 4 + i * pi / 2),
                    conections_dict["inner_conect_0"][1] + inner_straight_extend_distance * sin(-3 * pi / 4 + i * pi / 2),
                ],
                [
                    x
                    + (
                        inner_ring_radius
                        * cos((i * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)
                    ),
                    y
                    + (
                        inner_ring_radius
                        * sin((i * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)
                    ),
                ],
                [
                    x + (inner_ring_radius * cos((i * pi / 2) - inner_ring_arc_angle / 2)),
                    y + (inner_ring_radius * sin((i * pi / 2) - inner_ring_arc_angle / 2)),
                ],
            ]

            inner_0_path = gdspy.FlexPath(
                inner_0_conection_points,
                inner_ring_line_width,
                corners="circular bend",
                bend_radius=inner_bend_radius,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.make_flexpath_into_polygons_and_add_to_main(inner_0_path, mask_builder.layers.Nb_Antenna)
            # self.Main.add(inner_0_path)

            inner_1_conection_points = [
                conections_dict["inner_conect_1"],
                [
                    conections_dict["inner_conect_1"][0] + inner_straight_extend_distance * cos(3 * pi / 4 + i * pi / 2),
                    conections_dict["inner_conect_1"][1] + inner_straight_extend_distance * sin(3 * pi / 4 + i * pi / 2),
                ],
                [
                    x
                    + (
                        inner_ring_radius
                        * cos((i * pi / 2) + inner_ring_arc_angle / 2 - inner_straight_extend_distance / inner_ring_radius)
                    ),
                    y
                    + (
                        inner_ring_radius
                        * sin((i * pi / 2) + inner_ring_arc_angle / 2 - inner_straight_extend_distance / inner_ring_radius)
                    ),
                ],
                [
                    x + (inner_ring_radius * cos((i * pi / 2) + inner_ring_arc_angle / 2)),
                    y + (inner_ring_radius * sin((i * pi / 2) + inner_ring_arc_angle / 2)),
                ],
            ]

            inner_1_path = gdspy.FlexPath(
                inner_1_conection_points,
                inner_ring_line_width,
                corners="circular bend",
                bend_radius=inner_bend_radius,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.make_flexpath_into_polygons_and_add_to_main(inner_1_path, mask_builder.layers.Nb_Antenna)
            # self.Main.add(inner_1_path)

            # outer conections to ring overlap
            outer_0_conection_points = [
                conections_dict["outer_conect_0"],
                [
                    conections_dict["outer_conect_0"][0] + outer_straight_extend_distance * cos(-pi / 4 + i * pi / 2),
                    conections_dict["outer_conect_0"][1] + outer_straight_extend_distance * sin(-pi / 4 + i * pi / 2),
                ],
                [
                    x
                    + (
                        outer_ring_radius
                        * cos((i * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)
                    ),
                    y
                    + (
                        outer_ring_radius
                        * sin((i * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)
                    ),
                ],
                [
                    x + (outer_ring_radius * cos((i * pi / 2) - outer_ring_arc_angle / 2)),
                    y + (outer_ring_radius * sin((i * pi / 2) - outer_ring_arc_angle / 2)),
                ],
            ]

            outer_0_path = gdspy.FlexPath(
                outer_0_conection_points,
                inner_ring_line_width,
                corners="circular bend",
                bend_radius=outer_bend_radius,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.make_flexpath_into_polygons_and_add_to_main(outer_0_path, mask_builder.layers.Nb_Antenna)
            # self.Main.add(outer_0_path)

            outer_1_conection_points = [
                conections_dict["outer_conect_1"],
                [
                    conections_dict["outer_conect_1"][0] + outer_straight_extend_distance * cos(pi / 4 + i * pi / 2),
                    conections_dict["outer_conect_1"][1] + outer_straight_extend_distance * sin(pi / 4 + i * pi / 2),
                ],
                [
                    x
                    + (
                        outer_ring_radius
                        * cos((i * pi / 2) + outer_ring_arc_angle / 2 - outer_straight_extend_distance / outer_ring_radius)
                    ),
                    y
                    + (
                        outer_ring_radius
                        * sin((i * pi / 2) + outer_ring_arc_angle / 2 - outer_straight_extend_distance / outer_ring_radius)
                    ),
                ],
                [
                    x + (outer_ring_radius * cos((i * pi / 2) + outer_ring_arc_angle / 2)),
                    y + (outer_ring_radius * sin((i * pi / 2) + outer_ring_arc_angle / 2)),
                ],
            ]

            outer_1_path = gdspy.FlexPath(
                outer_1_conection_points,
                inner_ring_line_width,
                corners="circular bend",
                bend_radius=outer_bend_radius,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.make_flexpath_into_polygons_and_add_to_main(outer_1_path, mask_builder.layers.Nb_Antenna)
            # self.Main.add(outer_1_path)

        inner_ring_overlap_theta += pi / 2
        outer_ring_overlap_theta += pi / 2

    if only_1_pol and not with_crossover:
        return one_pol_connect_dict

    if not return_configurator_points:
        return points_to_conect_kids_dict

    if return_configurator_points_for_Lo_pass:
        return points_to_conect_kids_dict, configurator_points_for_Lo_pass

    if return_configurator_points_for_Hi_pass:
        return points_to_conect_kids_dict, configurator_points_for_Hi_pass

    if return_configurator_points_for_combiner_150ghz:
        return points_to_conect_kids_dict, configurator_points_for_combiner_150ghz

    if return_configurator_points_for_combiner_90ghz:
        return points_to_conect_kids_dict, configurator_points_for_combiner_90ghz

    configurator_points = {}

    annotate_ang = pi / 2 - 0.15
    configurator_points["outer_ring_radius"] = {
        "text": "outer_ring_radius",
        "start": [x, y],
        "end": [
            x + outer_ring_radius * cos(annotate_ang),
            y + outer_ring_radius * sin(annotate_ang),
        ],
    }

    annotate_ang -= 0.3
    configurator_points["inner_ring_radius"] = {
        "text": "inner_ring_radius",
        "start": [x, y],
        "end": [
            x + inner_ring_radius * cos(annotate_ang),
            y + inner_ring_radius * sin(annotate_ang),
        ],
    }

    configurator_points["ring_overlap_distance_from_center"] = {
        "text": "ring_overlap_distance_from_center",
        "start": [x, y],
        "end": [
            x + ring_overlap_distance_from_center * cos(0.0),
            y + ring_overlap_distance_from_center * sin(0.0),
        ],
    }

    annotate_ang = pi / 2 - 0.3

    configurator_points["inner_ring_line_width"] = {
        "text": "inner_ring_line_width",
        "start": [
            x + (inner_ring_radius - inner_ring_line_width / 2) * cos(annotate_ang),
            y + (inner_ring_radius - inner_ring_line_width / 2) * sin(annotate_ang),
        ],
        "end": [
            x + (inner_ring_radius + inner_ring_line_width / 2) * cos(annotate_ang),
            y + (inner_ring_radius + inner_ring_line_width / 2) * sin(annotate_ang),
        ],
    }

    configurator_points["outer_ring_line_width"] = {
        "text": "outer_ring_line_width",
        "start": [
            x + (outer_ring_radius - outer_ring_line_width / 2) * cos(annotate_ang),
            y + (outer_ring_radius - outer_ring_line_width / 2) * sin(annotate_ang),
        ],
        "end": [
            x + (outer_ring_radius + outer_ring_line_width / 2) * cos(annotate_ang),
            y + (outer_ring_radius + outer_ring_line_width / 2) * sin(annotate_ang),
        ],
    }

    configurator_points["inner_ring_overlap_gap"] = {
        "text": "inner_ring_overlap_gap",
        "start": [
            x + inner_ring_radius * cos(inner_ring_overlap_theta),
            y + inner_ring_radius * sin(inner_ring_overlap_theta),
        ],
        "end": [
            x + inner_ring_radius * cos(-inner_ring_overlap_theta),
            y + inner_ring_radius * sin(-inner_ring_overlap_theta),
        ],
    }

    configurator_points["outer_ring_overlap_gap"] = {
        "text": "outer_ring_overlap_gap",
        "start": [
            x + outer_ring_radius * cos(outer_ring_overlap_theta),
            y + outer_ring_radius * sin(outer_ring_overlap_theta),
        ],
        "end": [
            x + outer_ring_radius * cos(-outer_ring_overlap_theta),
            y + outer_ring_radius * sin(-outer_ring_overlap_theta),
        ],
    }

    configurator_points["outer_ring_conector_gap"] = {
        "text": "outer_ring_conector_gap",
        "start": [
            x + outer_ring_radius * cos((pi / 4) + outer_combiner_theta),
            y + outer_ring_radius * sin((pi / 4) + outer_combiner_theta),
        ],
        "end": [
            x + outer_ring_radius * cos((pi / 4) - outer_combiner_theta),
            y + outer_ring_radius * sin((pi / 4) - outer_combiner_theta),
        ],
    }

    configurator_points["inner_straight_extend_distance"] = {
        "text": "inner_straight_extend_distance",
        "start": [
            x + (inner_ring_radius * cos((3 * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)),
            y + (inner_ring_radius * sin((3 * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)),
        ],
        "end": [
            x + (inner_ring_radius * cos((3 * pi / 2) - inner_ring_arc_angle / 2)),
            y + (inner_ring_radius * sin((3 * pi / 2) - inner_ring_arc_angle / 2)),
        ],
    }

    configurator_points["inner_straight_extend_distance_B"] = {
        "text": "inner_straight_extend_distance",
        "start": [
            conections_dict["inner_conect_0"][0],
            conections_dict["inner_conect_0"][1],
        ],
        "end": [
            conections_dict["inner_conect_0"][0] + inner_straight_extend_distance * cos(-3 * pi / 4 + 3 * pi / 2),
            conections_dict["inner_conect_0"][1] + inner_straight_extend_distance * sin(-3 * pi / 4 + 3 * pi / 2),
        ],
    }

    dy = (
        -configurator_points["inner_straight_extend_distance_B"]["end"][1] + configurator_points["inner_straight_extend_distance"]["end"][1]
    )
    dx = (
        -configurator_points["inner_straight_extend_distance_B"]["end"][0] + configurator_points["inner_straight_extend_distance"]["end"][0]
    )
    angle_of_middle_section = np.arctan2(dy, dx)

    dy = (
        -configurator_points["inner_straight_extend_distance_B"]["start"][1]
        + configurator_points["inner_straight_extend_distance_B"]["end"][1]
    )
    dx = (
        -configurator_points["inner_straight_extend_distance_B"]["start"][0]
        + configurator_points["inner_straight_extend_distance_B"]["end"][0]
    )
    angle_of_overlap_ext_section = np.arctan2(dy, dx)

    dy = (
        -configurator_points["inner_straight_extend_distance"]["start"][1] + configurator_points["inner_straight_extend_distance"]["end"][1]
    )
    dx = (
        -configurator_points["inner_straight_extend_distance"]["start"][0] + configurator_points["inner_straight_extend_distance"]["end"][0]
    )
    angle_of_arc_ext_section = np.arctan2(dy, dx)

    configurator_points["inner_bend_radius"] = {
        "text": "inner_bend_radius",
        "start": [
            conections_dict["inner_conect_0"][0]
            + inner_straight_extend_distance * cos(-3 * pi / 4 + 3 * pi / 2)
            + inner_bend_radius * cos((pi / 2) + (angle_of_overlap_ext_section + angle_of_middle_section) / 2),
            conections_dict["inner_conect_0"][1]
            + inner_straight_extend_distance * sin(-3 * pi / 4 + 3 * pi / 2)
            + inner_bend_radius * sin((pi / 2) + (angle_of_overlap_ext_section + angle_of_middle_section) / 2),
        ],
        "end": [
            conections_dict["inner_conect_0"][0] + inner_straight_extend_distance * cos(-3 * pi / 4 + 3 * pi / 2),
            conections_dict["inner_conect_0"][1] + inner_straight_extend_distance * sin(-3 * pi / 4 + 3 * pi / 2),
        ],
    }

    configurator_points["inner_bend_radius_B"] = {
        "text": "inner_bend_radius",
        "start": [
            x
            + (inner_ring_radius * cos((3 * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius))
            + inner_bend_radius * cos((pi / 2) + (angle_of_arc_ext_section + angle_of_middle_section) / 2),
            y
            + (inner_ring_radius * sin((3 * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius))
            + inner_bend_radius * sin((pi / 2) + (angle_of_arc_ext_section + angle_of_middle_section) / 2),
        ],
        "end": [
            x + (inner_ring_radius * cos((3 * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)),
            y + (inner_ring_radius * sin((3 * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)),
        ],
    }

    configurator_points["outer_straight_extend_distance"] = {
        "text": "outer_straight_extend_distance",
        "start": [
            x + (outer_ring_radius * cos((3 * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)),
            y + (outer_ring_radius * sin((3 * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)),
        ],
        "end": [
            x + (outer_ring_radius * cos((3 * pi / 2) - outer_ring_arc_angle / 2)),
            y + (outer_ring_radius * sin((3 * pi / 2) - outer_ring_arc_angle / 2)),
        ],
    }

    configurator_points["outer_straight_extend_distance_B"] = {
        "text": "outer_straight_extend_distance",
        "start": [
            conections_dict["outer_conect_0"][0],
            conections_dict["outer_conect_0"][1],
        ],
        "end": [
            conections_dict["outer_conect_0"][0] + outer_straight_extend_distance * cos(-3 * pi / 2 + 3 * pi / 4),
            conections_dict["outer_conect_0"][1] + outer_straight_extend_distance * sin(-3 * pi / 2 + 3 * pi / 4),
        ],
    }

    dy = (
        -configurator_points["outer_straight_extend_distance_B"]["end"][1] + configurator_points["outer_straight_extend_distance"]["end"][1]
    )
    dx = (
        -configurator_points["outer_straight_extend_distance_B"]["end"][0] + configurator_points["outer_straight_extend_distance"]["end"][0]
    )
    angle_of_middle_section = np.arctan2(dy, dx)

    dy = (
        -configurator_points["outer_straight_extend_distance_B"]["start"][1]
        + configurator_points["outer_straight_extend_distance_B"]["end"][1]
    )
    dx = (
        -configurator_points["outer_straight_extend_distance_B"]["start"][0]
        + configurator_points["outer_straight_extend_distance_B"]["end"][0]
    )
    angle_of_overlap_ext_section = np.arctan2(dy, dx)

    dy = (
        -configurator_points["outer_straight_extend_distance"]["start"][1] + configurator_points["outer_straight_extend_distance"]["end"][1]
    )
    dx = (
        -configurator_points["outer_straight_extend_distance"]["start"][0] + configurator_points["outer_straight_extend_distance"]["end"][0]
    )
    angle_of_arc_ext_section = np.arctan2(dy, dx)

    configurator_points["outer_bend_radius"] = {
        "text": "outer_bend_radius",
        "start": [
            conections_dict["outer_conect_0"][0]
            + outer_straight_extend_distance * cos(-3 * pi / 2 + 3 * pi / 4)
            + outer_bend_radius * cos((pi / 2) + (angle_of_overlap_ext_section + angle_of_middle_section) / 2),
            conections_dict["outer_conect_0"][1]
            + outer_straight_extend_distance * sin(-3 * pi / 2 + 3 * pi / 4)
            + outer_bend_radius * sin((pi / 2) + (angle_of_overlap_ext_section + angle_of_middle_section) / 2),
        ],
        "end": [
            conections_dict["outer_conect_0"][0] + outer_straight_extend_distance * cos(-3 * pi / 2 + 3 * pi / 4),
            conections_dict["outer_conect_0"][1] + outer_straight_extend_distance * sin(-3 * pi / 2 + 3 * pi / 4),
        ],
    }

    configurator_points["outer_bend_radius_B"] = {
        "text": "outer_bend_radius",
        "start": [
            x
            + (outer_ring_radius * cos((3 * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius))
            + outer_bend_radius * cos((-pi / 2) + (angle_of_arc_ext_section + angle_of_middle_section) / 2),
            y
            + (outer_ring_radius * sin((3 * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius))
            + outer_bend_radius * sin((-pi / 2) + (angle_of_arc_ext_section + angle_of_middle_section) / 2),
        ],
        "end": [
            x + (outer_ring_radius * cos((3 * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)),
            y + (outer_ring_radius * sin((3 * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)),
        ],
    }

    angle_inc -= Hugh_filter_thin_arc_angles[0]
    configurator_points["Hugh_filter_thin_length1"] = {
        "text": "Hugh_filter_thin_length1",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    configurator_points["Hugh_filter_thin_width"] = {
        "text": "Hugh_filter_thin_width",
        "start": [
            x
            + (outer_ring_radius - (Hugh_filter_thin_width / 2)) * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
            y
            + (outer_ring_radius - (Hugh_filter_thin_width / 2)) * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
        ],
        "end": [
            x
            + (outer_ring_radius + (Hugh_filter_thin_width / 2)) * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
            y
            + (outer_ring_radius + (Hugh_filter_thin_width / 2)) * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
        ],
    }

    angle_inc -= Hugh_filter_thick_arc_angles[0]
    configurator_points["Hugh_filter_thick_length1"] = {
        "text": "Hugh_filter_thick_length1",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    configurator_points["Hugh_filter_thick_width"] = {
        "text": "Hugh_filter_thick_width",
        "start": [
            x
            + (outer_ring_radius - (Hugh_filter_thick_width / 2))
            * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
            y
            + (outer_ring_radius - (Hugh_filter_thick_width / 2))
            * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
        ],
        "end": [
            x
            + (outer_ring_radius + (Hugh_filter_thick_width / 2))
            * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
            y
            + (outer_ring_radius + (Hugh_filter_thick_width / 2))
            * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
        ],
    }

    angle_inc -= Hugh_filter_thin_arc_angles[1]
    configurator_points["Hugh_filter_thin_length2"] = {
        "text": "Hugh_filter_thin_length2",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[1]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[1]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    angle_inc -= Hugh_filter_thick_arc_angles[1]
    configurator_points["Hugh_filter_thick_length2"] = {
        "text": "Hugh_filter_thick_length2",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[1]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[1]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    angle_inc -= Hugh_filter_thin_arc_angles[1]
    configurator_points["Hugh_filter_thin_length2_B"] = {
        "text": "Hugh_filter_thin_length2",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[1]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[1]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    angle_inc -= Hugh_filter_thick_arc_angles[0]
    configurator_points["Hugh_filter_thick_length1_B"] = {
        "text": "Hugh_filter_thick_length1",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    angle_inc -= Hugh_filter_thin_arc_angles[0]
    configurator_points["Hugh_filter_thin_length1_B"] = {
        "text": "Hugh_filter_thin_length1",
        "start": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0]),
        ],
        "end": [
            x + outer_ring_radius * cos((pi) + mid_Hugh_angle - angle_inc),
            y + outer_ring_radius * sin((pi) + mid_Hugh_angle - angle_inc),
        ],
    }

    ang_mult = 0
    configurator_points["ring_overlap_0_rot"] = {
        "text": "ring_overlap_0_rot",
        "start": [
            x + ring_overlap_distance_from_center * cos((ang_mult * pi)),
            y + ring_overlap_distance_from_center * sin((ang_mult * pi)),
        ],
        "end": [
            x
            + ring_overlap_distance_from_center * cos((ang_mult * pi))
            + (outer_straight_extend_distance * cos(mbu.deg_to_rad(ring_overlap_0_rot))),
            y
            + ring_overlap_distance_from_center * sin((ang_mult * pi))
            + (outer_straight_extend_distance * sin(mbu.deg_to_rad(ring_overlap_0_rot))),
        ],
    }

    ang_mult = 1 / 2
    configurator_points["ring_overlap_1_rot"] = {
        "text": "ring_overlap_1_rot",
        "start": [
            x + ring_overlap_distance_from_center * cos((ang_mult * pi)),
            y + ring_overlap_distance_from_center * sin((ang_mult * pi)),
        ],
        "end": [
            x
            + ring_overlap_distance_from_center * cos((ang_mult * pi))
            + (outer_straight_extend_distance * cos(mbu.deg_to_rad(ring_overlap_1_rot))),
            y
            + ring_overlap_distance_from_center * sin((ang_mult * pi))
            + (outer_straight_extend_distance * sin(mbu.deg_to_rad(ring_overlap_1_rot))),
        ],
    }

    ang_mult = 1
    configurator_points["ring_overlap_2_rot"] = {
        "text": "ring_overlap_2_rot",
        "start": [
            x + ring_overlap_distance_from_center * cos((ang_mult * pi)),
            y + ring_overlap_distance_from_center * sin((ang_mult * pi)),
        ],
        "end": [
            x
            + ring_overlap_distance_from_center * cos((ang_mult * pi))
            + (outer_straight_extend_distance * cos(mbu.deg_to_rad(ring_overlap_2_rot))),
            y
            + ring_overlap_distance_from_center * sin((ang_mult * pi))
            + (outer_straight_extend_distance * sin(mbu.deg_to_rad(ring_overlap_2_rot))),
        ],
    }

    ang_mult = 3 / 2
    configurator_points["ring_overlap_3_rot"] = {
        "text": "ring_overlap_3_rot",
        "start": [
            x + ring_overlap_distance_from_center * cos((ang_mult * pi)),
            y + ring_overlap_distance_from_center * sin((ang_mult * pi)),
        ],
        "end": [
            x
            + ring_overlap_distance_from_center * cos((ang_mult * pi))
            + (outer_straight_extend_distance * cos(mbu.deg_to_rad(ring_overlap_3_rot))),
            y
            + ring_overlap_distance_from_center * sin((ang_mult * pi))
            + (outer_straight_extend_distance * sin(mbu.deg_to_rad(ring_overlap_3_rot))),
        ],
    }

    return points_to_conect_kids_dict, configurator_points
