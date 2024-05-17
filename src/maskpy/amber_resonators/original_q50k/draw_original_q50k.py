from typing import Callable

import gdspy
import numpy as np

# from ...souk_mask_builder import SoukMaskBuilder
from .utils_original_q50k import _get_config_checking_override


def draw(
    mask_builder,
    x: float,
    y: float,
    rot_angle: float,
    f0: float,
    config_override: dict | None = None,
    mux_func_override: Callable | None = None,
    mirror=False,
    IDC_and_frame_material: str = "IDC_Nb",
    meander_material: str = "IDC_Nb",
):
    """Draw the original q50k amber resonator to the Main cell at the x,y
    cooardinate given.

    ?
    The resonator is placed where the base middle of the inductive meander
    is at this x,y. The resonator geometry is defined by the dimensions within
    the Main_config_file_dict. By default it will, but optionally can
    choose not to, add all the neccessay cutouts for the structure.
    ?

    KwArgs
    ------
    config_override: dict | None = None

    """

    IDC_and_frame_material_lookup = {"IDC_Nb": mask_builder.IDC_Nb, "Nb": mask_builder.Nb_Antenna, "Al": mask_builder.Aluminium}
    material_idc_and_frame = IDC_and_frame_material_lookup[IDC_and_frame_material]

    meander_material_lookup = {"Al": mask_builder.Aluminium, "IDC_Nb": mask_builder.IDC_Nb, "Nb": mask_builder.Nb_Antenna}
    material_meander = meander_material_lookup[meander_material]

    config = _get_config_checking_override(config_override)

    # Making the meander section
    meander_lw = config["meander_lw"]
    meander_corner_bend_radius = config["meander_corner_bend_radius"]

    meander_bot_width = config["meander_bot_width"]

    number_of_right_meander_sections = 9
    number_of_left_meander_sections = 9

    # Making the meander path corner points
    right_side_meander_path_points = [
        [(meander_bot_width / 2), -(meander_lw / 2)],
    ]

    for section_no in range(1, number_of_right_meander_sections + 1):
        # Adding the vertical section
        previous_corner_xy = right_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0]
        corner_y = previous_corner_xy[1] + config[f"meander_right_height_{section_no}"]
        right_side_meander_path_points.append([corner_x, corner_y])

        # Adding the horizontal section
        previous_corner_xy = right_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0] + config[f"meander_right_width_{section_no}"]
        corner_y = previous_corner_xy[1]
        right_side_meander_path_points.append([corner_x, corner_y])

    left_side_meander_path_points = [
        [-(meander_bot_width / 2), -(meander_lw / 2)],
    ]

    for section_no in range(1, number_of_left_meander_sections + 1):
        # Adding the vertical section
        previous_corner_xy = left_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0]
        corner_y = previous_corner_xy[1] + config[f"meander_left_height_{section_no}"]
        left_side_meander_path_points.append([corner_x, corner_y])

        # Adding the horizontal section
        previous_corner_xy = left_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0] + config[f"meander_left_width_{section_no}"]
        corner_y = previous_corner_xy[1]
        left_side_meander_path_points.append([corner_x, corner_y])

    # Adding the final horizontal section to the left side. It has a last height but not a last width
    left_previous_corner_xy = left_side_meander_path_points[-1]
    corner_x = left_previous_corner_xy[0]
    corner_y = left_previous_corner_xy[1] + config[f"meander_left_height_10"]
    left_side_meander_path_points.append([corner_x, corner_y])

    meander_path_points = right_side_meander_path_points[::-1] + left_side_meander_path_points

    # Making the left and right side of the frame
    right_side_end_of_meander = meander_path_points[0]
    left_side_end_of_meander = meander_path_points[-1]

    coupler_arm_lw = config["coupler_arm_lw"]
    idc_arms_offset_from_meander = config["idc_arms_offset_from_meander"]
    top_idc_arm_to_coupler_arm = config["top_idc_arm_to_coupler_arm"]
    idc_arm_lw = config["idc_arm_lw"]
    idc_arm_spacing = config["idc_arm_spacing"]

    idc_arm_lengths = np.ones(21) * 3960
    no_of_idc_arms = len(idc_arm_lengths)

    coupler_len = 2056

    left_frame_height_overrun = config["left_frame_height_overrun"]
    left_frame_width = config["left_frame_width"]

    left_frame_height = (
        idc_arms_offset_from_meander + (no_of_idc_arms * idc_arm_lw) + ((no_of_idc_arms - 1) * idc_arm_spacing) + left_frame_height_overrun
    )
    left_frame_points = [
        [left_side_end_of_meander[0] - (left_frame_width / 2), left_side_end_of_meander[1]],  # Bot Left Corner
        [left_side_end_of_meander[0] - (left_frame_width / 2), left_side_end_of_meander[1] + left_frame_height],  # Top Left Corner
        [left_side_end_of_meander[0] + (left_frame_width / 2), left_side_end_of_meander[1] + left_frame_height],  # Top Right Corner
        [left_side_end_of_meander[0] + (left_frame_width / 2), left_side_end_of_meander[1]],  # Bot Right Corner
    ]

    right_frame_height = (
        idc_arms_offset_from_meander + (no_of_idc_arms * idc_arm_lw) + ((no_of_idc_arms - 1) * idc_arm_spacing) + top_idc_arm_to_coupler_arm
    )
    right_frame_width = config["right_frame_width"]

    # this is the right side frame and the coupler arm defined from the bot right corner going anti-clockwise.
    right_frame_points = [
        [right_side_end_of_meander[0], right_side_end_of_meander[1] + (meander_lw / 2)],  # Bot right corner
        [
            right_side_end_of_meander[0],
            right_side_end_of_meander[1] + (meander_lw / 2) + right_frame_height + (coupler_arm_lw),
        ],  # Top right corner
        [
            right_side_end_of_meander[0] - coupler_len,
            right_side_end_of_meander[1] + (meander_lw / 2) + right_frame_height + (coupler_arm_lw),
        ],  # Top left corner inside coupler fork
        [
            right_side_end_of_meander[0] - coupler_len,
            right_side_end_of_meander[1] + (meander_lw / 2) + right_frame_height,
        ],  # Bot left corner inside coupler fork
        [
            right_side_end_of_meander[0] - right_frame_width,
            right_side_end_of_meander[1] + (meander_lw / 2) + right_frame_height,
        ],  # Bot side of coupler next to frame
        [
            right_side_end_of_meander[0] - right_frame_width,
            right_side_end_of_meander[1] + (meander_lw / 2),
        ],  # Bot left corner
    ]

    # Making the idc arm polygons
    initial_idc_arm_y = left_side_end_of_meander[1] + idc_arms_offset_from_meander
    left_side_idc_arm_x = left_frame_points[-1][0]
    right_side_idc_arm_x = right_frame_points[-1][0]

    idc_arm_poly_points_dict: dict[str, list[list[float]]] = {}
    for i, arm_length in enumerate(idc_arm_lengths):

        bot_side_y_for_arm = initial_idc_arm_y + (i * (idc_arm_lw + idc_arm_spacing))

        # idc arms always start on the left hand side
        if i % 2 == 0:
            frame_side_x_for_arm = left_side_idc_arm_x
            sign = +1
        else:
            frame_side_x_for_arm = right_side_idc_arm_x
            sign = -1

        # Making the IDC arm points. Defined from the bot (R | L) corner of the
        # arm going (Clockwise | AntiClockwise). Bracketed expresions here is
        # if the arm is attached to the (RightFrame | LeftFame).
        poly_points = [
            [frame_side_x_for_arm, bot_side_y_for_arm],
            [frame_side_x_for_arm + (sign * arm_length), bot_side_y_for_arm],
            [frame_side_x_for_arm + (sign * arm_length), bot_side_y_for_arm + idc_arm_lw],
            [frame_side_x_for_arm, bot_side_y_for_arm + idc_arm_lw],
        ]
        idc_arm_poly_points_dict[f"arm_{i+1}"] = poly_points

    # Making the coupler fork polygon
    coupler_fork_bot_arm_lw = config["coupler_fork_bot_arm_lw"]
    coupler_fork_top_arm_lw = config["coupler_fork_top_arm_lw"]
    coupler_fork_stub_lw = config["coupler_fork_stub_lw"]
    coupler_fork_x_offset_from_coupler_arm = config["coupler_fork_x_offset_from_coupler_arm"]
    coupler_fork_bot_arm_y_offset_from_coupler_arm = config["coupler_fork_bot_arm_y_offset_from_coupler_arm"]
    coupler_fork_top_arm_y_offset_from_coupler_arm = config["coupler_fork_top_arm_y_offset_from_coupler_arm"]
    coupler_fork_top_arm_offset_to_ground = config["coupler_fork_top_arm_offset_to_ground"]

    print(right_frame_points[2])
    coupler_fork_TR = [
        right_frame_points[2][0] - coupler_fork_x_offset_from_coupler_arm,
        right_frame_points[2][1]
        + coupler_fork_top_arm_y_offset_from_coupler_arm
        + coupler_fork_top_arm_lw
        + coupler_fork_top_arm_offset_to_ground,
    ]

    coupler_fork_top_arm_len = 1000
    coupler_fork_bot_arm_len = 1000

    # Defined from the top right corner of the coupler fork going Clockwise
    coupler_fork_poly_points = [
        [coupler_fork_TR[0], coupler_fork_TR[1]],
        [coupler_fork_TR[0], coupler_fork_TR[1] - coupler_fork_top_arm_offset_to_ground],
        [
            coupler_fork_TR[0] + coupler_fork_top_arm_len,
            coupler_fork_TR[1] - coupler_fork_top_arm_offset_to_ground,
        ],  # Top left corner of coupler fork bot arm
        [
            coupler_fork_TR[0] + coupler_fork_top_arm_len,
            coupler_fork_TR[1] - coupler_fork_top_arm_offset_to_ground - coupler_fork_top_arm_lw,
        ],
        [
            coupler_fork_TR[0],
            coupler_fork_TR[1] - coupler_fork_top_arm_offset_to_ground - coupler_fork_top_arm_lw,
        ],  # Bot left corner of coupler fork bot arm
        [
            coupler_fork_TR[0],
            coupler_fork_TR[1]
            - coupler_fork_top_arm_offset_to_ground
            - coupler_fork_top_arm_lw
            - coupler_fork_top_arm_y_offset_from_coupler_arm
            - coupler_arm_lw
            - coupler_fork_bot_arm_y_offset_from_coupler_arm,
        ],  # Top left corner of coupler fork bot arm
        [
            coupler_fork_TR[0] + coupler_fork_bot_arm_len,
            coupler_fork_TR[1]
            - coupler_fork_top_arm_offset_to_ground
            - coupler_fork_top_arm_lw
            - coupler_fork_top_arm_y_offset_from_coupler_arm
            - coupler_arm_lw
            - coupler_fork_bot_arm_y_offset_from_coupler_arm,
        ],
        [
            coupler_fork_TR[0] + coupler_fork_bot_arm_len,
            coupler_fork_TR[1]
            - coupler_fork_top_arm_offset_to_ground
            - coupler_fork_top_arm_lw
            - coupler_fork_top_arm_y_offset_from_coupler_arm
            - coupler_arm_lw
            - coupler_fork_bot_arm_y_offset_from_coupler_arm
            - coupler_fork_bot_arm_lw,
        ],
        [
            coupler_fork_TR[0] - coupler_fork_stub_lw,
            coupler_fork_TR[1]
            - coupler_fork_top_arm_offset_to_ground
            - coupler_fork_top_arm_lw
            - coupler_fork_top_arm_y_offset_from_coupler_arm
            - coupler_arm_lw
            - coupler_fork_bot_arm_y_offset_from_coupler_arm
            - coupler_fork_bot_arm_lw,
        ],  # Bot left corner of coupler fork
        [coupler_fork_TR[0] - coupler_fork_stub_lw, coupler_fork_TR[1]],
    ]

    print(len(coupler_fork_poly_points))
    print(coupler_fork_poly_points)

    #################################
    # Adding everything to the mask #
    #################################
    # Adding the meander
    if mirror:
        new_meander_path_points = mask_builder.mirror_points_around_yaxis(meander_path_points)
        new_meander_path_points = mask_builder.rotate_and_move_points_list(new_meander_path_points, rot_angle, x, y)
    else:
        new_meander_path_points = mask_builder.rotate_and_move_points_list(meander_path_points, rot_angle, x, y)

    meander_path = gdspy.FlexPath(
        new_meander_path_points, meander_lw, corners="circular bend", bend_radius=meander_corner_bend_radius, **material_meander
    )

    mask_builder.make_flexpath_into_polygons_and_add_to_main(
        meander_path, layer=material_meander["layer"], datatype=material_meander["datatype"]
    )

    # Adding the left and right frame
    if mirror:
        new_left_frame_points = mask_builder.mirror_points_around_yaxis(left_frame_points)
        new_left_frame_points = mask_builder.rotate_and_move_points_list(new_left_frame_points, rot_angle, x, y)
    else:
        new_left_frame_points = mask_builder.rotate_and_move_points_list(left_frame_points, rot_angle, x, y)

    left_frame_polygon = gdspy.Polygon(
        new_left_frame_points, layer=material_idc_and_frame["layer"], datatype=material_idc_and_frame["datatype"]
    )
    mask_builder.Main.add(left_frame_polygon)

    if mirror:
        new_right_frame_points = mask_builder.mirror_points_around_yaxis(right_frame_points)
        new_right_frame_points = mask_builder.rotate_and_move_points_list(new_right_frame_points, rot_angle, x, y)
    else:
        new_right_frame_points = mask_builder.rotate_and_move_points_list(right_frame_points, rot_angle, x, y)

    right_frame_polygon = gdspy.Polygon(
        new_right_frame_points, layer=material_idc_and_frame["layer"], datatype=material_idc_and_frame["datatype"]
    )
    mask_builder.Main.add(right_frame_polygon)

    # Adding the IDC arms
    for arm_poly_points in idc_arm_poly_points_dict.values():
        if mirror:
            new_arm_poly_points = mask_builder.mirror_points_around_yaxis(arm_poly_points)
            new_arm_poly_points = mask_builder.rotate_and_move_points_list(new_arm_poly_points, rot_angle, x, y)
        else:
            new_arm_poly_points = mask_builder.rotate_and_move_points_list(arm_poly_points, rot_angle, x, y)

        arm_polygon = gdspy.Polygon(new_arm_poly_points, layer=material_idc_and_frame["layer"], datatype=material_idc_and_frame["datatype"])
        mask_builder.Main.add(arm_polygon)

    # Adding the coupler fork
    if mirror:
        new_coupler_fork_poly_points = mask_builder.mirror_points_around_yaxis(coupler_fork_poly_points)
        new_coupler_fork_poly_points = mask_builder.rotate_and_move_points_list(new_coupler_fork_poly_points, rot_angle, x, y)
    else:
        new_coupler_fork_poly_points = mask_builder.rotate_and_move_points_list(coupler_fork_poly_points, rot_angle, x, y)

    coupler_fork_polygon = gdspy.Polygon(
        new_coupler_fork_poly_points, layer=material_idc_and_frame["layer"], datatype=material_idc_and_frame["datatype"]
    )
    mask_builder.Main.add(coupler_fork_polygon)

    return
