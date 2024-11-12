from collections.abc import Sequence

import gdspy
import numpy as np
from matplotlib.font_manager import FontProperties
from numpy import pi

from ... import amber_muxing, amber_resonators
from ... import mask_builder_utils as mbu
from ...amber_resonators import AmberResonatorType
from ...logging import TextColor, TextStyle, pretty_print
from ...souk_mask_builder import SoukMaskBuilder
from ...souk_mask_configs import config_keys
from .types import TestChipSettings


def _get_caller_info(return_full: bool = False) -> str:
    """Get the filename of the callee.

    KwArgs
    ------
    return_full: bool = False
        When True will return of filepath/filename, not just filename.

    Returns
    -------
    callee: str
        The file that called.
    """
    import inspect
    import os

    # first get the full filename (including path and file extension)
    caller_frame = inspect.stack()[2]  # This is 2 FrameInfos deep, 0, and 1 are this file.
    caller_filename_full = caller_frame.filename
    caller_filename_only = os.path.splitext(os.path.basename(caller_filename_full))[0]

    if return_full:
        return caller_filename_full
    else:
        return caller_filename_only


def fixed_amber_mux_func(kid_id):
    """."""
    max_idcl = 3960
    CCL = 990 + 1066  # (length of coupler in the fork) + (length to very right side of the right frame).

    kid_no_to_no_of_arms_and_top_arm_idcl = {
        0: (21, 3960),
        1: (21, 3180),
        2: (21, 2400),
        3: (21, 1620),
        4: (21, 840),
        5: (21, 60),
        6: (17, 3960),
        7: (17, 3180),
        8: (17, 2400),
        9: (17, 1620),
        10: (17, 840),
        11: (17, 60),
        12: (13, 3960),
        13: (13, 3180),
        14: (13, 2400),
        15: (13, 1620),
        16: (13, 840),
        17: (13, 60),
    }
    no_of_arms, top_arm_idcl = kid_no_to_no_of_arms_and_top_arm_idcl[kid_id]

    IDC_array = np.ones(no_of_arms) * max_idcl

    IDC_array[0] = top_arm_idcl

    return IDC_array, CCL


def add_test_chip(
    mask_builder: SoukMaskBuilder,
    test_chip_cetner_xy: Sequence[float | int],
    test_chip_settings: TestChipSettings,
) -> tuple[list[float], list[float]]:
    """Add souk test chip.

    Parameters
    ----------
    mask_builder: SoukMaskBuilder
        The mask_builder to add the test chip to.

    test_chip_settings: TestChipSettings
        The settings to be used in making this test chip.

    Returns
    -------
    feedline_ends: tuple[list[float], list[float]]
        This is ([x, y], [x, y]) where the first list is the feedline end at
        the top of the test chip, and the second list is the feedline end at
        the bottom of the test chip.
    """
    TEST_CHIP_HEIGHT: float = 27000.0
    TEST_CHIP_WIDTH: float = 27000.0

    FEEDLINE_CENTER_WIDTH: float = 16.0
    FEEDLINE_CUTOUT_WIDTH: float = 28.0

    X_OFFSET_BETWEEN_RESONATORS: float = 6000.0
    INITIAL_X_FOR_RESONATORS: float = 2900.0

    ###########################################################################
    # Setup.
    ###########################################################################

    mask_builder.add_test_chip_quad_tabbed_dicing_line(
        test_chip_cetner_xy,
        TEST_CHIP_WIDTH,
        TEST_CHIP_HEIGHT,
        tab_positions=["all"],
    )

    ground_plane_layer_str = test_chip_settings.get("ground_plane_layer", None)
    feedline_center_layer_str = test_chip_settings.get("feedline_center_layer", None)

    if ground_plane_layer_str is None:
        ground_plane_layer = mask_builder.layers.Nb_Groundplane
    else:
        ground_plane_layer = mask_builder.layers.get_layer_from_layer_name(ground_plane_layer_str)

    if feedline_center_layer_str is None:
        feedline_center_layer = None
    else:
        feedline_center_layer = mask_builder.layers.get_layer_from_layer_name(feedline_center_layer_str)

    ground_plane_rect = gdspy.Rectangle(
        (
            (test_chip_cetner_xy[0] - (TEST_CHIP_WIDTH / 2)),
            (test_chip_cetner_xy[1] - (TEST_CHIP_HEIGHT / 2)),
        ),
        (
            (test_chip_cetner_xy[0] + (TEST_CHIP_WIDTH / 2)),
            (test_chip_cetner_xy[1] + (TEST_CHIP_HEIGHT / 2)),
        ),
        layer=ground_plane_layer.number,
        datatype=ground_plane_layer.datatype,
        # layer=mask_builder.layers.Nb_Antenna.number,
        # datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    if ground_plane_layer == mask_builder.layers.Nb_Groundplane:
        mask_builder.ground_plane_positives.add(ground_plane_rect)
    else:
        mask_builder.Main.add(ground_plane_rect)

    ###########################################################################
    # Adding text to chip
    ###########################################################################
    top_left_text = test_chip_settings["top_left_text"]
    top_right_text = test_chip_settings["top_right_text"]
    bottom_left_text = test_chip_settings["bottom_left_text"]

    TEXT_CORNER_OFFSET = 1000
    TEXT_SIZE = 300
    TEXT_LAYER = mask_builder.layers.Aluminium

    text_sub: dict[str, str] = {
        "$filename": _get_caller_info(),
        "$full_filename": _get_caller_info(return_full=True),
    }

    if bottom_left_text != "":
        bottom_left_text = text_sub.get(bottom_left_text, bottom_left_text)
        bottom_left_chip_corner = [
            test_chip_cetner_xy[0] - (TEST_CHIP_WIDTH / 2),
            test_chip_cetner_xy[1] - (TEST_CHIP_HEIGHT / 2),
        ]
        bottom_left_text_x = bottom_left_chip_corner[0] + TEXT_CORNER_OFFSET
        bottom_left_text_y = bottom_left_chip_corner[1] + TEXT_CORNER_OFFSET

        mask_builder.add_fancy_text(
            bottom_left_text,
            bottom_left_text_x,
            bottom_left_text_y,
            TEXT_SIZE,
            TEXT_LAYER,
            bb_cutout_in_grnd=True,
            bb_cutout_in_sin_dep=True,
        )

    if top_left_text != "":
        top_left_text = text_sub.get(top_left_text, top_left_text)
        top_left_chip_corner = [
            test_chip_cetner_xy[0] - (TEST_CHIP_WIDTH / 2),
            test_chip_cetner_xy[1] + (TEST_CHIP_HEIGHT / 2),
        ]
        top_left_text_x = top_left_chip_corner[0] + TEXT_CORNER_OFFSET
        top_left_text_y = top_left_chip_corner[1] - TEXT_CORNER_OFFSET

        mask_builder.add_fancy_text(
            top_left_text,
            top_left_text_x,
            top_left_text_y,
            TEXT_SIZE,
            TEXT_LAYER,
            bb_cutout_in_grnd=True,
            bb_cutout_in_sin_dep=True,
            vertical_align="below",
        )

    # Adding the top_right_text if True
    if top_right_text != "":
        top_right_text = text_sub.get(top_right_text, top_right_text)
        top_right_chip_corner = [
            test_chip_cetner_xy[0] + (TEST_CHIP_WIDTH / 2),
            test_chip_cetner_xy[1] + (TEST_CHIP_HEIGHT / 2),
        ]
        top_right_text_x = top_right_chip_corner[0] - TEXT_CORNER_OFFSET
        top_right_text_y = top_right_chip_corner[1] - TEXT_CORNER_OFFSET

        mask_builder.add_fancy_text(
            top_right_text,
            top_right_text_x,
            top_right_text_y,
            TEXT_SIZE,
            TEXT_LAYER,
            bb_cutout_in_grnd=True,
            bb_cutout_in_sin_dep=True,
            vertical_align="below",
            horizontal_align="end",
        )

    ###########################################################################
    # Ports and feedline.
    ###########################################################################
    # [x, y] coordinates based on chip_center. Defined from bottom to top.
    feedline_inner_corners: list[list[float]] = [
        # Bottom Half
        [0.0, -10465.0],
        [12000.0, -10465.0],
        [12000.0, -7130.0],
        [-12000.0, -7130.0],
        [-12000.0, -3130.0],
        [12000.0, -3130.0],
        # Top Half
        [12000.0, 3130.0],
        [-12000.0, 3130.0],
        [-12000.0, 7130.0],
        [12000.0, 7130.0],
        [12000.0, 10465.0],
        [0.0, 10465.0],
    ]

    feedline_inner_corners = mbu.move_points_list(feedline_inner_corners, test_chip_cetner_xy[0], test_chip_cetner_xy[1])

    port_config: dict[str, float] = {
        config_keys.port.outer_feedline_width: 500.0,
        config_keys.port.outer_cutout_around_feedline_width: 875.0,
        config_keys.port.outer_dielectric_under_feedline_width: 0.0,
        config_keys.port.outer_back_length: 530.0,
        config_keys.port.taper_length: 1500.0,
        config_keys.port.dielectric_cutout_in_port_width: 0.0,
        config_keys.port.dielectric_cutout_in_port_length: 0.0,
    }

    feedline_config: dict[str, float] = {
        config_keys.cpw_feedline.feedline_width: FEEDLINE_CENTER_WIDTH,
        config_keys.cpw_feedline.cutout_around_feedline_width: FEEDLINE_CUTOUT_WIDTH,
        config_keys.cpw_feedline.dielectric_under_feedline_width: 0,
        config_keys.cpw_feedline.bend_radius: 60,
        config_keys.cpw_feedline.extra_straight_length: 0,
        config_keys.cpw_feedline.bridge_gap: 0,
        config_keys.cpw_feedline.bridge_width: 0,
    }

    top_port_end = [
        test_chip_cetner_xy[0],
        test_chip_cetner_xy[1] + (TEST_CHIP_HEIGHT / 2),
    ]

    top_port_connection = mask_builder.add_port_and_get_connection_point(
        top_port_end[0],
        top_port_end[1],
        -(pi / 2),
        cpw_feedline_config_override=feedline_config,
        port_config_override=port_config,
        add_extra_squares=False,
        dielectric_cutout_in_port=True,
    )

    bot_port_end = [
        test_chip_cetner_xy[0],
        test_chip_cetner_xy[1] - (TEST_CHIP_HEIGHT / 2),
    ]

    bot_port_connection = mask_builder.add_port_and_get_connection_point(
        bot_port_end[0],
        bot_port_end[1],
        (pi / 2),
        cpw_feedline_config_override=feedline_config,
        port_config_override=port_config,
        add_extra_squares=False,
        dielectric_cutout_in_port=True,
    )

    feedline_points: list[list[float]] = []

    feedline_points.append(
        [
            bot_port_connection[0],
            bot_port_connection[1],
        ]
    )

    for point in feedline_inner_corners:
        feedline_points.append(point)

    feedline_points.append(
        [
            top_port_connection[0],
            top_port_connection[1],
        ]
    )

    mask_builder.add_feedline_and_dielectric(
        feedline_points,
        cpw_feedline_config_override=feedline_config,
        center_material=feedline_center_layer,
        add_bridges=False,
    )

    ###########################################################################
    # Adding the amber resonators.
    ###########################################################################

    GROUND_PLANE_GAP = 20

    y_coords_for_kids = [
        # bot 3
        feedline_points[3][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[3][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[3][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        # next 3
        feedline_points[5][1] - FEEDLINE_CUTOUT_WIDTH - GROUND_PLANE_GAP,
        feedline_points[5][1] - FEEDLINE_CUTOUT_WIDTH - GROUND_PLANE_GAP,
        feedline_points[5][1] - FEEDLINE_CUTOUT_WIDTH - GROUND_PLANE_GAP,
        # next 3
        feedline_points[5][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[5][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[5][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        # next 3
        feedline_points[7][1] - FEEDLINE_CUTOUT_WIDTH - GROUND_PLANE_GAP,
        feedline_points[7][1] - FEEDLINE_CUTOUT_WIDTH - GROUND_PLANE_GAP,
        feedline_points[7][1] - FEEDLINE_CUTOUT_WIDTH - GROUND_PLANE_GAP,
        # next 3
        feedline_points[7][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[7][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[7][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        # top 3
        feedline_points[9][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[9][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
        feedline_points[9][1] + FEEDLINE_CUTOUT_WIDTH + GROUND_PLANE_GAP,
    ]

    x_coords_for_kids = [
        # bot 3
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS,
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (1 * X_OFFSET_BETWEEN_RESONATORS),
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (2 * X_OFFSET_BETWEEN_RESONATORS),
        # next 3
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS,
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (1 * X_OFFSET_BETWEEN_RESONATORS),
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (2 * X_OFFSET_BETWEEN_RESONATORS),
        # next 3
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS,
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (1 * X_OFFSET_BETWEEN_RESONATORS),
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (2 * X_OFFSET_BETWEEN_RESONATORS),
        # next 3
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS,
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (1 * X_OFFSET_BETWEEN_RESONATORS),
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (2 * X_OFFSET_BETWEEN_RESONATORS),
        # next 3
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS,
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (1 * X_OFFSET_BETWEEN_RESONATORS),
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (2 * X_OFFSET_BETWEEN_RESONATORS),
        # top 3
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS,
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (1 * X_OFFSET_BETWEEN_RESONATORS),
        test_chip_cetner_xy[0] + INITIAL_X_FOR_RESONATORS - (2 * X_OFFSET_BETWEEN_RESONATORS),
    ]
    rotation_for_kids = [
        # bot 3
        pi,
        pi,
        pi,
        # next 3
        0,
        0,
        0,
        # next 3
        pi,
        pi,
        pi,
        # next 3
        0,
        0,
        0,
        # next 3
        pi,
        pi,
        pi,
        # top 3
        pi,
        pi,
        pi,
    ]
    mirror_for_kids = [
        # bot 3
        True,
        True,
        True,
        # next 3
        False,
        False,
        False,
        # next 3
        True,
        True,
        True,
        # next 3
        False,
        False,
        False,
        # next 3
        True,
        True,
        True,
        # top 3
        True,
        True,
        True,
    ]

    for kid_number, (x, y, rot_angle, mirror) in enumerate(
        zip(
            x_coords_for_kids,
            y_coords_for_kids,
            rotation_for_kids,
            mirror_for_kids,
        )
    ):
        resonator_settings = test_chip_settings[f"resonator_{kid_number + 1}_settings"]

        resonator_type_str = resonator_settings.pop("resonator_type")
        resonator_type = AmberResonatorType(resonator_type_str)
        f0 = resonator_settings.pop("f0")

        if resonator_settings.get("mux_func_override", None) is not None:
            mux_func_override_str = resonator_settings.pop("mux_func_override")
            if mux_func_override_str == "fixed_amber_mux_func":
                mux_func_override = fixed_amber_mux_func
            else:
                mux_func_override = mux_func_override_str
        else:
            mux_func_override = None

        mask_builder.add_amber_resonator(
            resonator_type,
            x,
            y,
            rot_angle,
            f0,
            mirror=mirror,
            mux_func_override=mux_func_override,
            **resonator_settings,
        )

    pretty_print(
        f"Added Amber Test Chip to mask at:\n    x: {test_chip_cetner_xy[0]}\n    y: {test_chip_cetner_xy[1]}",
        color=TextColor.BRIGHT_BLUE,
    )

    return (top_port_end, bot_port_end)
