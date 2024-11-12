from typing import Callable

import gdspy
import numpy as np

from ... import mask_builder_utils as mbu
from ...amber_muxing import get_mux_func_for_resonator_type
from ...layers import Layer
from ...logging import styled_type_error
from ..resonator_types import AmberResonatorType
from ..utils.get_config import get_resonator_config

# from ...souk_mask_builder import SoukMaskBuilder  # circular import.


def draw(
    # mask_builder: SoukMaskBuilder,
    mask_builder,
    resonator_type: AmberResonatorType,
    x: float,
    y: float,
    rot_angle: float,
    f0: float,
    mux_func_override: Callable | None = None,
    resonator_config_override: dict | None = None,
    mirror=False,
    IDC_and_frame_material: Layer | str | None = None,
    meander_material: Layer | str | None = None,
    coupler_fork_material: Layer | str | None = None,
    add_grnd_cutout: bool = True,
    add_SiN_dep_dielectric_cutout: bool = True,
    add_SiO_cutout: bool = True,
    add_SiN_membrane_cutout: bool = True,
    add_backside_check: bool = False,
    add_inductor_cover: bool = False,
):
    """Draw an Amber resonator to the Main cell at the x,y cooardinate given.

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
    if not isinstance(resonator_type, AmberResonatorType):
        styled_type_error(resonator_type, "resonator_type", AmberResonatorType)

    accepted_resonator_types = [
        AmberResonatorType.ORIGINAL_NB_IND2,
        AmberResonatorType.ORIGINAL_NB_IND4,
        AmberResonatorType.ORIGINAL_AL_IND2,
        AmberResonatorType.ORIGINAL_AL_IND4,
    ]

    if resonator_type not in accepted_resonator_types:
        raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

    resonator_config = get_resonator_config(resonator_type, resonator_config_override=resonator_config_override)

    if mux_func_override is None:
        mux_func = get_mux_func_for_resonator_type(resonator_type)
    else:
        mux_func = mux_func_override

    if isinstance(IDC_and_frame_material, str):
        IDC_and_frame_material = mask_builder.layers.get_layer_from_layer_name(IDC_and_frame_material)

    if isinstance(meander_material, str):
        meander_material = mask_builder.layers.get_layer_from_layer_name(meander_material)

    if isinstance(coupler_fork_material, str):
        coupler_fork_material = mask_builder.layers.get_layer_from_layer_name(coupler_fork_material)

    if IDC_and_frame_material is None:
        match resonator_type:
            case AmberResonatorType.ORIGINAL_AL_IND2:
                IDC_and_frame_material = mask_builder.layers.Aluminium
            case AmberResonatorType.ORIGINAL_AL_IND4:
                IDC_and_frame_material = mask_builder.layers.Aluminium
            case AmberResonatorType.ORIGINAL_NB_IND2:
                IDC_and_frame_material = mask_builder.layers.IDC_Nb
            case AmberResonatorType.ORIGINAL_NB_IND4:
                IDC_and_frame_material = mask_builder.layers.IDC_Nb

    if meander_material is None:
        match resonator_type:
            case AmberResonatorType.ORIGINAL_AL_IND2:
                meander_material = mask_builder.layers.Aluminium
            case AmberResonatorType.ORIGINAL_AL_IND4:
                meander_material = mask_builder.layers.Aluminium
            case AmberResonatorType.ORIGINAL_NB_IND2:
                meander_material = mask_builder.layers.IDC_Nb
            case AmberResonatorType.ORIGINAL_NB_IND4:
                meander_material = mask_builder.layers.IDC_Nb

    if coupler_fork_material is None:
        coupler_fork_material = mask_builder.layers.Nb_Groundplane

    # Making the meander section
    meander_lw = resonator_config["meander_lw"]
    meander_corner_bend_radius = resonator_config["meander_corner_bend_radius"]

    meander_bot_width = resonator_config["meander_bot_width"]

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
        corner_y = previous_corner_xy[1] + resonator_config[f"meander_right_height_{section_no}"]
        right_side_meander_path_points.append([corner_x, corner_y])

        # Adding the horizontal section
        previous_corner_xy = right_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0] + resonator_config[f"meander_right_width_{section_no}"]
        corner_y = previous_corner_xy[1]
        right_side_meander_path_points.append([corner_x, corner_y])

    left_side_meander_path_points = [
        [-(meander_bot_width / 2), -(meander_lw / 2)],
    ]

    for section_no in range(1, number_of_left_meander_sections + 1):
        # Adding the vertical section
        previous_corner_xy = left_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0]
        corner_y = previous_corner_xy[1] + resonator_config[f"meander_left_height_{section_no}"]
        left_side_meander_path_points.append([corner_x, corner_y])

        # Adding the horizontal section
        previous_corner_xy = left_side_meander_path_points[-1]
        corner_x = previous_corner_xy[0] + resonator_config[f"meander_left_width_{section_no}"]
        corner_y = previous_corner_xy[1]
        left_side_meander_path_points.append([corner_x, corner_y])

    # Adding the final horizontal section to the left side. It has a last height but not a last width
    left_previous_corner_xy = left_side_meander_path_points[-1]
    corner_x = left_previous_corner_xy[0]
    corner_y = left_previous_corner_xy[1] + resonator_config[f"meander_left_height_10"]
    left_side_meander_path_points.append([corner_x, corner_y])

    meander_path_points = right_side_meander_path_points[::-1] + left_side_meander_path_points

    # Making the left and right side of the frame
    right_side_end_of_meander = meander_path_points[0]
    left_side_end_of_meander = meander_path_points[-1]

    coupler_arm_lw = resonator_config["coupler_arm_lw"]
    idc_arms_offset_from_meander = resonator_config["idc_arms_offset_from_meander"]
    top_idc_arm_to_coupler_arm = resonator_config["top_idc_arm_to_coupler_arm"]
    idc_arm_lw = resonator_config["idc_arm_lw"]
    idc_arm_spacing = resonator_config["idc_arm_spacing"]

    idc_arm_lengths, coupler_len = mux_func(f0)
    no_of_idc_arms = len(idc_arm_lengths)

    left_frame_height_overrun = resonator_config["left_frame_height_overrun"]
    left_frame_width = resonator_config["left_frame_width"]

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
    right_frame_width = resonator_config["right_frame_width"]

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
    coupler_fork_bot_arm_lw = resonator_config["coupler_fork_bot_arm_lw"]
    coupler_fork_top_arm_lw = resonator_config["coupler_fork_top_arm_lw"]
    coupler_fork_stub_lw = resonator_config["coupler_fork_stub_lw"]
    coupler_fork_x_offset_from_coupler_arm = resonator_config["coupler_fork_x_offset_from_coupler_arm"]
    coupler_fork_bot_arm_y_offset_from_coupler_arm = resonator_config["coupler_fork_bot_arm_y_offset_from_coupler_arm"]
    coupler_fork_top_arm_y_offset_from_coupler_arm = resonator_config["coupler_fork_top_arm_y_offset_from_coupler_arm"]
    coupler_fork_top_arm_offset_to_ground = resonator_config["coupler_fork_top_arm_offset_to_ground"]

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

    # Making the ground plane cutout
    grnd_cutout_gap_left = resonator_config["grnd_cutout_gap_left"]
    grnd_cutout_gap_right = resonator_config["grnd_cutout_gap_right"]
    grnd_cutout_gap_top = resonator_config["grnd_cutout_gap_top"]
    grnd_cutout_gap_bot = resonator_config["grnd_cutout_gap_bot"]

    ground_cutout_start_x = left_frame_points[0][0] - grnd_cutout_gap_left
    ground_cutout_end_x = right_frame_points[0][0] + grnd_cutout_gap_right

    ground_cutout_start_y = right_side_meander_path_points[0][1] - (meander_lw / 2) - grnd_cutout_gap_bot
    ground_cutout_end_y = coupler_fork_poly_points[0][1] + grnd_cutout_gap_top

    # Defined from the bot left going Clockwise
    ground_cutout_poly_points = [
        [ground_cutout_start_x, ground_cutout_start_y],
        [ground_cutout_end_x, ground_cutout_start_y],
        [ground_cutout_end_x, ground_cutout_end_y],
        [ground_cutout_start_x, ground_cutout_end_y],
    ]

    # Making the silicon nitride deposition cutout
    SiN_dep_cutout_gap_left = resonator_config["SiN_dep_cutout_gap_left"]
    SiN_dep_cutout_gap_right = resonator_config["SiN_dep_cutout_gap_right"]
    SiN_dep_cutout_gap_top = resonator_config["SiN_dep_cutout_gap_top"]
    SiN_dep_cutout_gap_bot = resonator_config["SiN_dep_cutout_gap_bot"]

    SiN_dep_cutout_start_x = left_frame_points[0][0] - SiN_dep_cutout_gap_left
    SiN_dep_cutout_end_x = right_frame_points[0][0] + SiN_dep_cutout_gap_right

    SiN_dep_cutout_start_y = right_side_meander_path_points[0][1] - (meander_lw / 2) - SiN_dep_cutout_gap_bot
    SiN_dep_cutout_end_y = coupler_fork_poly_points[0][1] + SiN_dep_cutout_gap_top

    # Defined from the bot left going Clockwise
    SiN_dep_cutout_poly_points = [
        [SiN_dep_cutout_start_x, SiN_dep_cutout_start_y],
        [SiN_dep_cutout_end_x, SiN_dep_cutout_start_y],
        [SiN_dep_cutout_end_x, SiN_dep_cutout_end_y],
        [SiN_dep_cutout_start_x, SiN_dep_cutout_end_y],
    ]

    # Making the silicon oxide cutout
    SiO_cutout_gap_left = resonator_config["SiO_cutout_gap_left"]
    SiO_cutout_gap_right = resonator_config["SiO_cutout_gap_right"]
    SiO_cutout_gap_top = resonator_config["SiO_cutout_gap_top"]
    SiO_cutout_gap_bot = resonator_config["SiO_cutout_gap_bot"]

    SiO_cutout_start_x = left_frame_points[0][0] - SiO_cutout_gap_left
    SiO_cutout_end_x = right_frame_points[0][0] + SiO_cutout_gap_right

    SiO_cutout_start_y = right_side_meander_path_points[0][1] - (meander_lw / 2) - SiO_cutout_gap_bot
    SiO_cutout_end_y = coupler_fork_poly_points[0][1] + SiO_cutout_gap_top

    # Defined from the bot left going Clockwise
    SiO_cutout_poly_points = [
        [SiO_cutout_start_x, SiO_cutout_start_y],
        [SiO_cutout_end_x, SiO_cutout_start_y],
        [SiO_cutout_end_x, SiO_cutout_end_y],
        [SiO_cutout_start_x, SiO_cutout_end_y],
    ]

    # Making the silicon nitride membrane cutout
    SiN_membrane_cutout_gap_left = resonator_config["SiN_membrane_cutout_gap_left"]
    SiN_membrane_cutout_gap_right = resonator_config["SiN_membrane_cutout_gap_right"]
    SiN_membrane_cutout_gap_top = resonator_config["SiN_membrane_cutout_gap_top"]
    SiN_membrane_cutout_gap_bot = resonator_config["SiN_membrane_cutout_gap_bot"]

    SiN_membrane_cutout_start_x = left_frame_points[0][0] - SiN_membrane_cutout_gap_left
    SiN_membrane_cutout_end_x = right_frame_points[0][0] + SiN_membrane_cutout_gap_right

    SiN_membrane_cutout_start_y = right_side_meander_path_points[0][1] - (meander_lw / 2) - SiN_membrane_cutout_gap_bot
    SiN_membrane_cutout_end_y = coupler_fork_poly_points[0][1] + SiN_membrane_cutout_gap_top

    # Defined from the bot left going Clockwise
    SiN_membrane_cutout_poly_points = [
        [SiN_membrane_cutout_start_x, SiN_membrane_cutout_start_y],
        [SiN_membrane_cutout_end_x, SiN_membrane_cutout_start_y],
        [SiN_membrane_cutout_end_x, SiN_membrane_cutout_end_y],
        [SiN_membrane_cutout_start_x, SiN_membrane_cutout_end_y],
    ]

    # Making the backside check
    backside_check_gap_left = resonator_config["backside_check_gap_left"]
    backside_check_gap_right = resonator_config["backside_check_gap_right"]
    backside_check_gap_top = resonator_config["backside_check_gap_top"]
    backside_check_gap_bot = resonator_config["backside_check_gap_bot"]

    backside_check_start_x = left_frame_points[0][0] - backside_check_gap_left
    backside_check_end_x = right_frame_points[0][0] + backside_check_gap_right

    backside_check_start_y = right_side_meander_path_points[0][1] - (meander_lw / 2) - backside_check_gap_bot
    backside_check_end_y = coupler_fork_poly_points[0][1] + backside_check_gap_top

    # Defined from the bot left going Clockwise
    backside_check_poly_points = [
        [backside_check_start_x, backside_check_start_y],
        [backside_check_end_x, backside_check_start_y],
        [backside_check_end_x, backside_check_end_y],
        [backside_check_start_x, backside_check_end_y],
    ]

    # Making the inductor cover
    inductor_cover_gap_left = resonator_config["inductor_cover_gap_left"]
    inductor_cover_gap_right = resonator_config["inductor_cover_gap_right"]
    inductor_cover_gap_top = resonator_config["inductor_cover_gap_top"]
    inductor_cover_gap_bot = resonator_config["inductor_cover_gap_bot"]

    inductor_cover_start_x = left_frame_points[0][0] - inductor_cover_gap_left
    inductor_cover_end_x = right_frame_points[0][0] + inductor_cover_gap_right

    inductor_cover_start_y = right_side_meander_path_points[0][1] - (meander_lw / 2) - inductor_cover_gap_bot
    inductor_cover_end_y = right_side_meander_path_points[-1][1] + (meander_lw / 2) + inductor_cover_gap_top

    # Defined from the bot left going Clockwise
    inductor_cover_poly_points = [
        [inductor_cover_start_x, inductor_cover_start_y],
        [inductor_cover_end_x, inductor_cover_start_y],
        [inductor_cover_end_x, inductor_cover_end_y],
        [inductor_cover_start_x, inductor_cover_end_y],
    ]

    # Aluminium patch and etch settings
    aluminium_patch_offset_bot = 15
    aluminium_patch_offset_top = 0
    aluminium_patch_offset_left = 20
    aluminium_patch_offset_right = 20

    aluminium_etch_offset_bot = 15
    aluminium_etch_offset_top = 20
    aluminium_etch_offset_left = 20
    aluminium_etch_offset_right = 20

    # Making the Aluminium patch under the meander
    aluminium_patch_meander_left_edge = left_frame_points[0][0] - aluminium_patch_offset_left
    aluminium_patch_meander_right_edge = right_side_end_of_meander[0] + aluminium_patch_offset_right
    aluminium_patch_meander_top_edge = left_side_end_of_meander[1] + aluminium_patch_offset_top
    aluminium_patch_meander_bot_edge = -meander_lw - aluminium_patch_offset_bot

    aluminium_etch_meander_left_edge = aluminium_patch_meander_left_edge - aluminium_etch_offset_left
    aluminium_etch_meander_right_edge = aluminium_patch_meander_right_edge + aluminium_etch_offset_right
    aluminium_etch_meander_top_edge = aluminium_patch_meander_top_edge + aluminium_etch_offset_top
    aluminium_etch_meander_bot_edge = aluminium_patch_meander_bot_edge - aluminium_etch_offset_bot

    # defined from bot left going anti-clockwise
    aluminium_patch_meander_points = [
        [aluminium_patch_meander_left_edge, aluminium_patch_meander_bot_edge],
        [aluminium_patch_meander_right_edge, aluminium_patch_meander_bot_edge],
        [aluminium_patch_meander_right_edge, aluminium_patch_meander_top_edge],
        [aluminium_patch_meander_left_edge, aluminium_patch_meander_top_edge],
    ]
    aluminium_etch_meander_points = [
        [aluminium_etch_meander_left_edge, aluminium_etch_meander_bot_edge],
        [aluminium_etch_meander_right_edge, aluminium_etch_meander_bot_edge],
        [aluminium_etch_meander_right_edge, aluminium_etch_meander_top_edge],
        [aluminium_etch_meander_left_edge, aluminium_etch_meander_top_edge],
    ]

    # Making the Aluminium patch under the idc section
    aluminium_patch_idc_left_edge = left_frame_points[0][0] - aluminium_patch_offset_left
    aluminium_patch_idc_right_edge = right_side_end_of_meander[0] + aluminium_patch_offset_right
    # aluminium_patch_idc_top_edge = right_frame_points[1][1] + aluminium_patch_offset_top
    aluminium_patch_idc_bot_edge = -(meander_lw / 2) - aluminium_patch_offset_bot

    aluminium_etch_idc_left_edge = aluminium_patch_idc_left_edge - aluminium_etch_offset_left
    aluminium_etch_idc_right_edge = aluminium_patch_idc_right_edge + aluminium_etch_offset_right
    # aluminium_etch_idc_top_edge = aluminium_patch_idc_top_edge + aluminium_etch_offset_top
    aluminium_etch_idc_bot_edge = aluminium_patch_idc_bot_edge - aluminium_etch_offset_bot

    # defined from bot left going anti-clockwise
    aluminium_patch_idc_points = [
        [aluminium_patch_idc_left_edge, aluminium_patch_idc_bot_edge],
        [aluminium_patch_idc_right_edge, aluminium_patch_idc_bot_edge],
        [aluminium_patch_idc_right_edge, right_frame_points[1][1] + coupler_fork_top_arm_y_offset_from_coupler_arm * 0.5],
        [
            right_frame_points[2][0] - (coupler_fork_x_offset_from_coupler_arm * 0.5),
            right_frame_points[2][1] + (coupler_fork_top_arm_y_offset_from_coupler_arm * 0.5),
        ],
        [
            right_frame_points[3][0] - (coupler_fork_x_offset_from_coupler_arm * 0.5),
            right_frame_points[3][1] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.5),
        ],
        [
            right_frame_points[4][0] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.5),
            right_frame_points[4][1] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.5),
        ],
        [
            right_frame_points[4][0] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.5),
            left_frame_points[1][1] + left_frame_width,
        ],
        [
            aluminium_patch_idc_left_edge,
            left_frame_points[1][1] + left_frame_width,
        ],
    ]
    aluminium_etch_idc_points = [
        [aluminium_etch_idc_left_edge, aluminium_etch_idc_bot_edge],
        [aluminium_etch_idc_right_edge, aluminium_etch_idc_bot_edge],
        [aluminium_etch_idc_right_edge, right_frame_points[1][1] + coupler_fork_top_arm_y_offset_from_coupler_arm * 0.75],
        [
            right_frame_points[2][0] - (coupler_fork_x_offset_from_coupler_arm * 0.75),
            right_frame_points[2][1] + (coupler_fork_top_arm_y_offset_from_coupler_arm * 0.75),
        ],
        [
            right_frame_points[3][0] - (coupler_fork_x_offset_from_coupler_arm * 0.75),
            right_frame_points[3][1] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.75),
        ],
        [
            right_frame_points[4][0] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.75),
            right_frame_points[4][1] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.75),
        ],
        [
            right_frame_points[4][0] - (coupler_fork_bot_arm_y_offset_from_coupler_arm * 0.75),
            left_frame_points[1][1] + (left_frame_width * 1.5),
        ],
        [
            aluminium_etch_idc_left_edge,
            left_frame_points[1][1] + (left_frame_width * 1.5),
        ],
    ]

    #################################################
    # Adding to the x and y transform to ensure the #
    # coupler arm sits at the cooardinates given.   #
    #################################################

    dx = -(coupler_fork_poly_points[0][0] + (coupler_fork_poly_points[-1][0] - coupler_fork_poly_points[0][0]) / 2)
    dy = -coupler_fork_poly_points[0][1]

    meander_path_points = mbu.move_points_list(meander_path_points, dx, dy)
    left_frame_points = mbu.move_points_list(left_frame_points, dx, dy)
    right_frame_points = mbu.move_points_list(right_frame_points, dx, dy)
    for key, val in idc_arm_poly_points_dict.items():
        idc_arm_poly_points_dict[key] = mbu.move_points_list(val, dx, dy)
    coupler_fork_poly_points = mbu.move_points_list(coupler_fork_poly_points, dx, dy)
    ground_cutout_poly_points = mbu.move_points_list(ground_cutout_poly_points, dx, dy)
    SiN_dep_cutout_poly_points = mbu.move_points_list(SiN_dep_cutout_poly_points, dx, dy)
    SiO_cutout_poly_points = mbu.move_points_list(SiO_cutout_poly_points, dx, dy)
    SiN_membrane_cutout_poly_points = mbu.move_points_list(SiN_membrane_cutout_poly_points, dx, dy)
    backside_check_poly_points = mbu.move_points_list(backside_check_poly_points, dx, dy)
    inductor_cover_poly_points = mbu.move_points_list(inductor_cover_poly_points, dx, dy)

    aluminium_patch_meander_points = mbu.move_points_list(aluminium_patch_meander_points, dx, dy)
    aluminium_etch_meander_points = mbu.move_points_list(aluminium_etch_meander_points, dx, dy)

    aluminium_patch_idc_points = mbu.move_points_list(aluminium_patch_idc_points, dx, dy)
    aluminium_etch_idc_points = mbu.move_points_list(aluminium_etch_idc_points, dx, dy)

    #################################
    # Adding everything to the mask #
    #################################

    # Adding the meander
    if mirror:
        new_meander_path_points = mbu.mirror_points_around_yaxis(meander_path_points)
        new_meander_path_points = mbu.rotate_and_move_points_list(new_meander_path_points, rot_angle, x, y)
    else:
        new_meander_path_points = mbu.rotate_and_move_points_list(meander_path_points, rot_angle, x, y)

    meander_path = gdspy.FlexPath(
        new_meander_path_points,
        meander_lw,
        corners="circular bend",
        bend_radius=meander_corner_bend_radius,
        layer=meander_material.number,
        datatype=meander_material.datatype,
    )

    mask_builder.make_flexpath_into_polygons_and_add_to_main(
        meander_path,
        layer=meander_material,
    )

    # Adding the left and right frame
    if mirror:
        new_left_frame_points = mbu.mirror_points_around_yaxis(left_frame_points)
        new_left_frame_points = mbu.rotate_and_move_points_list(new_left_frame_points, rot_angle, x, y)
    else:
        new_left_frame_points = mbu.rotate_and_move_points_list(left_frame_points, rot_angle, x, y)

    left_frame_polygon = gdspy.Polygon(
        new_left_frame_points,
        layer=IDC_and_frame_material.number,
        datatype=IDC_and_frame_material.datatype,
    )
    mask_builder.Main.add(left_frame_polygon)

    if mirror:
        new_right_frame_points = mbu.mirror_points_around_yaxis(right_frame_points)
        new_right_frame_points = mbu.rotate_and_move_points_list(new_right_frame_points, rot_angle, x, y)
    else:
        new_right_frame_points = mbu.rotate_and_move_points_list(right_frame_points, rot_angle, x, y)

    right_frame_polygon = gdspy.Polygon(
        new_right_frame_points,
        layer=IDC_and_frame_material.number,
        datatype=IDC_and_frame_material.datatype,
    )
    mask_builder.Main.add(right_frame_polygon)
    size_diff = (1.0, 1.0)
    # size_diff = (-1.0, -1.0)
    # size_diff = (1.0, 2.0)
    # size_diff = (1.0, 10.0)
    new_right_frame_points_expanded = mbu.size_polygon(new_right_frame_points, size_diff)

    new_right_frame_poly_expanded = gdspy.Polygon(
        new_right_frame_points_expanded,
        layer=mask_builder.layers.General_labeling.number,
        datatype=mask_builder.layers.General_labeling.datatype,
    )
    mask_builder.Main.add(new_right_frame_poly_expanded)

    # Adding the IDC arms
    for arm_poly_points in idc_arm_poly_points_dict.values():
        if mirror:
            new_arm_poly_points = mbu.mirror_points_around_yaxis(arm_poly_points)
            new_arm_poly_points = mbu.rotate_and_move_points_list(new_arm_poly_points, rot_angle, x, y)
        else:
            new_arm_poly_points = mbu.rotate_and_move_points_list(arm_poly_points, rot_angle, x, y)

        arm_polygon = gdspy.Polygon(
            new_arm_poly_points,
            layer=IDC_and_frame_material.number,
            datatype=IDC_and_frame_material.datatype,
        )
        mask_builder.Main.add(arm_polygon)

    # Adding the coupler fork
    if mirror:
        new_coupler_fork_poly_points = mbu.mirror_points_around_yaxis(coupler_fork_poly_points)
        new_coupler_fork_poly_points = mbu.rotate_and_move_points_list(new_coupler_fork_poly_points, rot_angle, x, y)
    else:
        new_coupler_fork_poly_points = mbu.rotate_and_move_points_list(coupler_fork_poly_points, rot_angle, x, y)

    coupler_fork_polygon = gdspy.Polygon(
        new_coupler_fork_poly_points,
        layer=coupler_fork_material.number,
        datatype=coupler_fork_material.datatype,
    )
    mask_builder.Main.add(coupler_fork_polygon)

    # Adding the ground plane cutout
    if add_grnd_cutout:
        if mirror:
            new_ground_cutout_poly_points = mbu.mirror_points_around_yaxis(ground_cutout_poly_points)
            new_ground_cutout_poly_points = mbu.rotate_and_move_points_list(new_ground_cutout_poly_points, rot_angle, x, y)
        else:
            new_ground_cutout_poly_points = mbu.rotate_and_move_points_list(ground_cutout_poly_points, rot_angle, x, y)

        ground_cutout_polygon = gdspy.Polygon(
            new_ground_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_cutouts.add(ground_cutout_polygon)

    # Adding the silicon nitride desposition cutout
    if add_SiN_dep_dielectric_cutout:
        if mirror:
            new_SiN_dep_cutout_poly_points = mbu.mirror_points_around_yaxis(SiN_dep_cutout_poly_points)
            new_SiN_dep_cutout_poly_points = mbu.rotate_and_move_points_list(new_SiN_dep_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiN_dep_cutout_poly_points = mbu.rotate_and_move_points_list(SiN_dep_cutout_poly_points, rot_angle, x, y)

        SiN_dep_cutout_polygon = gdspy.Polygon(
            new_SiN_dep_cutout_poly_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(SiN_dep_cutout_polygon)

    # Adding the silicon oxide cutout
    if add_SiO_cutout:
        if mirror:
            new_SiO_cutout_poly_points = mbu.mirror_points_around_yaxis(SiO_cutout_poly_points)
            new_SiO_cutout_poly_points = mbu.rotate_and_move_points_list(new_SiO_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiO_cutout_poly_points = mbu.rotate_and_move_points_list(SiO_cutout_poly_points, rot_angle, x, y)

        SiO_cutout_polygon = gdspy.Polygon(
            new_SiO_cutout_poly_points,
            layer=mask_builder.layers.SiO.number,
            datatype=mask_builder.layers.SiO.datatype,
        )
        mask_builder.silicon_oxide_cutouts.add(SiO_cutout_polygon)

    # Adding the silicon nitride membrane cutout
    if add_SiN_membrane_cutout:
        if mirror:
            new_SiN_membrane_cutout_poly_points = mbu.mirror_points_around_yaxis(SiN_membrane_cutout_poly_points)
            new_SiN_membrane_cutout_poly_points = mbu.rotate_and_move_points_list(new_SiN_membrane_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiN_membrane_cutout_poly_points = mbu.rotate_and_move_points_list(SiN_membrane_cutout_poly_points, rot_angle, x, y)

        SiN_membrane_cutout_polygon = gdspy.Polygon(
            new_SiN_membrane_cutout_poly_points,
            layer=mask_builder.layers.SiN_Membrane.number,
            datatype=mask_builder.layers.SiN_Membrane.datatype,
        )
        mask_builder.silicon_nitride_membrane_cutouts.add(SiN_membrane_cutout_polygon)

    # Adding the silicon nitride membrane cutout
    if add_backside_check:
        if mirror:
            new_backside_check_poly_points = mbu.mirror_points_around_yaxis(backside_check_poly_points)
            new_backside_check_poly_points = mbu.rotate_and_move_points_list(new_backside_check_poly_points, rot_angle, x, y)
        else:
            new_backside_check_poly_points = mbu.rotate_and_move_points_list(backside_check_poly_points, rot_angle, x, y)

        backside_check_polygon = gdspy.Polygon(
            new_backside_check_poly_points,
            layer=mask_builder.layers.Backside_Check.number,
            datatype=mask_builder.layers.Backside_Check.datatype,
        )
        mask_builder.Main.add(backside_check_polygon)

    # Adding the cover over the inductor
    if add_inductor_cover:
        if mirror:
            new_inductor_cover_poly_points = mbu.mirror_points_around_yaxis(inductor_cover_poly_points)
            new_inductor_cover_poly_points = mbu.rotate_and_move_points_list(new_inductor_cover_poly_points, rot_angle, x, y)
        else:
            new_inductor_cover_poly_points = mbu.rotate_and_move_points_list(inductor_cover_poly_points, rot_angle, x, y)

        inductor_cover_polygon = gdspy.Polygon(
            new_inductor_cover_poly_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.Main.add(inductor_cover_polygon)

    # if meander_material == "Al":
    if meander_material == mask_builder.layers.Aluminium:
        if mirror:
            new_aluminium_patch_meander_points = mbu.mirror_points_around_yaxis(aluminium_patch_meander_points)
            new_aluminium_patch_meander_points = mbu.rotate_and_move_points_list(new_aluminium_patch_meander_points, rot_angle, x, y)

            new_aluminium_etch_meander_points = mbu.mirror_points_around_yaxis(aluminium_etch_meander_points)
            new_aluminium_etch_meander_points = mbu.rotate_and_move_points_list(new_aluminium_etch_meander_points, rot_angle, x, y)
        else:
            new_aluminium_patch_meander_points = mbu.rotate_and_move_points_list(aluminium_patch_meander_points, rot_angle, x, y)
            new_aluminium_etch_meander_points = mbu.rotate_and_move_points_list(aluminium_etch_meander_points, rot_angle, x, y)

        aluminium_patch_meander = gdspy.Polygon(
            new_aluminium_patch_meander_points,
            layer=mask_builder.layers.Aluminium_Patch.number,
            datatype=mask_builder.layers.Aluminium_Patch.datatype,
        )
        mask_builder.Main.add(aluminium_patch_meander)

        aluminium_etch_meander = gdspy.Polygon(
            new_aluminium_etch_meander_points,
            layer=mask_builder.layers.Aluminium_Etch.number,
            datatype=mask_builder.layers.Aluminium_Etch.datatype,
        )
        mask_builder.aluminium_etch_positives.add(aluminium_etch_meander)

    if IDC_and_frame_material == mask_builder.layers.Aluminium:
        if mirror:
            new_aluminium_patch_idc_points = mbu.mirror_points_around_yaxis(aluminium_patch_idc_points)
            new_aluminium_patch_idc_points = mbu.rotate_and_move_points_list(new_aluminium_patch_idc_points, rot_angle, x, y)

            new_aluminium_etch_idc_points = mbu.mirror_points_around_yaxis(aluminium_etch_idc_points)
            new_aluminium_etch_idc_points = mbu.rotate_and_move_points_list(new_aluminium_etch_idc_points, rot_angle, x, y)
        else:
            new_aluminium_patch_idc_points = mbu.rotate_and_move_points_list(aluminium_patch_idc_points, rot_angle, x, y)
            new_aluminium_etch_idc_points = mbu.rotate_and_move_points_list(aluminium_etch_idc_points, rot_angle, x, y)

        aluminium_patch_idc = gdspy.Polygon(
            new_aluminium_patch_idc_points,
            layer=mask_builder.layers.Aluminium_Patch.number,
            datatype=mask_builder.layers.Aluminium_Patch.datatype,
        )
        mask_builder.Main.add(aluminium_patch_idc)

        aluminium_etch_idc = gdspy.Polygon(
            new_aluminium_etch_idc_points,
            layer=mask_builder.layers.Aluminium_Etch.number,
            datatype=mask_builder.layers.Aluminium_Etch.datatype,
        )
        mask_builder.aluminium_etch_positives.add(aluminium_etch_idc)

    return
