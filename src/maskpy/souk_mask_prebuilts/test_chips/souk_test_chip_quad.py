from collections.abc import Sequence

import numpy as np
from numpy import cos, pi, sin

from ... import souk_mask_components as smc
from ... import souk_muxing
from ...logging import TextColor, pretty_print
from ...souk_mask_builder import SoukMaskBuilder
from ...souk_resonators import SoukResonatorType
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
    NO_OF_SOUK_RESONATORS: int = 16
    TEST_CHIP_HEIGHT: float = 27000.0
    TEST_CHIP_WIDTH: float = 27000.0

    ###########################################################################
    # Making the boundary
    horn_positions: list[list[float]]

    text_keys = [
        "top_left_text",
        "top_right_text",
        "bottom_left_text",
    ]
    text_sub: dict[str, str] = {
        "$filename": _get_caller_info(),
        "$full_filename": _get_caller_info(return_full=True),
        "$chip_id": "Chip_" + str(test_chip_settings["chip_id"]),
    }
    for text_key in text_keys:
        if (new_val := text_sub.get(test_chip_settings[text_key], None)) is not None:
            test_chip_settings[text_key] = new_val

    (
        horn_positions,
        chip_outer_poly_points,
    ) = smc.add_test_chip_quadrent_boundary_and_get_horn_positions(
        mask_builder,
        test_chip_cetner_xy,
        test_chip_quad_config_override=None,
        # top_right_text=test_chip_settings["top_right_text"],
        # top_left_text=test_chip_settings["top_left_text"],
        # bottom_left_text=test_chip_settings["bot_left_text"],
        cardiff_logo=True,
        souk_logo=True,
        return_outer_poly_points=True,
        dice_line_tab_positions=["all"],
        **test_chip_settings,
        # time_stamp_position="TL",
        # add_groundplane_under_test_chip=True,
        # add_SiN_dep_under_test_chip=False,
    )

    ###########################################################################
    # Generating a mapping for KID IDs
    mapped_IDs = np.zeros(NO_OF_SOUK_RESONATORS, dtype=int)
    max_KID_ID = len(mapped_IDs)
    for i in range(int(max_KID_ID / 4)):
        mapped_IDs[i * 4 : (i * 4) + 4] = np.array(
            [
                0 + i,
                int(max_KID_ID / 4) + i,
                int(max_KID_ID / 2) + i,
                int(3 * max_KID_ID / 4) + i,
            ]
        )

    mapped_freqs = souk_muxing.get_evenly_spaced_freq_array(
        NO_OF_SOUK_RESONATORS,
        None,
        order=mapped_IDs,
    )

    ###########################################################################
    # building for all quadrents
    feedline_pass_through_points: list[list[float]] = []

    ###########################################################################
    # Making each quad
    for quad_id, (horn_center_x, horn_center_y) in enumerate(horn_positions):
        kid_numbers = mapped_IDs[quad_id * 4 : (quad_id + 1) * 4]
        f0s = mapped_freqs[quad_id * 4 : (quad_id + 1) * 4]

        quad_settings = test_chip_settings[f"quad_{quad_id}_settings"]

        under_text = quad_settings["text_under_quad"]
        mask_builder.add_text_under_horn_in_test_chip_quad(under_text, horn_center_x, horn_center_y)

        if isinstance(quad_settings["resonator_type"], str):
            resonator_types = [
                SoukResonatorType(quad_settings["resonator_type"]),
                SoukResonatorType(quad_settings["resonator_type"]),
                SoukResonatorType(quad_settings["resonator_type"]),
                SoukResonatorType(quad_settings["resonator_type"]),
            ]
        elif isinstance(quad_settings["resonator_type"], Sequence):
            resonator_types = [
                SoukResonatorType(quad_settings["resonator_type"][0]),
                SoukResonatorType(quad_settings["resonator_type"][1]),
                SoukResonatorType(quad_settings["resonator_type"][2]),
                SoukResonatorType(quad_settings["resonator_type"][3]),
            ]
        else:
            resonator_types = [
                SoukResonatorType(quad_settings["resonator_type"]),
                SoukResonatorType(quad_settings["resonator_type"]),
                SoukResonatorType(quad_settings["resonator_type"]),
                SoukResonatorType(quad_settings["resonator_type"]),
            ]

        relative_kid_positions = mask_builder.get_relative_kid_positions(
            resonator_types,
            **quad_settings,
        )

        mask_builder.add_4_resonators_around_horn(
            horn_center_x,
            horn_center_y,
            relative_kid_positions,
            kid_numbers,
            f0s,
            resonator_types,
            **quad_settings,
            # mux_func_overrides=quad_settings["mux_func_overrides"],
            # resonator_config_overrides=quad_settings["resonator_config_overrides"],
            # general_config_override=quad_settings["general_config_override"],
            # IDC_and_frame_materials=quad_settings["idc_and_frame_materials"],
            # meander_materials=quad_settings["meander_materials"],
            # trim_lengths=quad_settings["trim_lengths"],
            # add_grnd_cutout=quad_settings["add_grnd_cutout"],
            # add_SiN_dep_dielectric_around=quad_settings["add_SiN_dep_dielectric_around"],
            # add_SiN_dep_dielectric_cutout=quad_settings["add_SiN_dep_dielectric_cutout"],
            # add_SiO_cutout=quad_settings["add_SiO_cutout"],
            # add_SiN_membrane_cutout=quad_settings["add_SiN_membrane_cutout"],
            # add_backside_check=quad_settings["add_backside_check"],
            # add_grnd_cutout_over_inductor=quad_settings["add_grnd_cutout_over_inductor"],
            # add_SiN_dep_dielectric_cutout_over_inductor=quad_settings["add_SiN_dep_dielectric_cutout_over_inductor"],
            # add_Aluminium_Patch_and_Etch=quad_settings["add_Aluminium_Patch_and_Etch"],
        )

        for points in mask_builder.get_feedline_pass_through_points(
            horn_center_x,
            horn_center_y,
            relative_kid_positions,
            resonator_types,
            **test_chip_settings,
            **quad_settings,
        ):
            feedline_pass_through_points.append(points)

        if quad_settings["add_antenna"]:
            if quad_settings.get("antenna_rotation", None) is None:
                quad_settings["antenna_rotation"] = 0.0

            mask_builder.add_antenna(horn_center_x, horn_center_y, quad_settings["antenna_rotation"])
            relative_antena_conect_positions = mask_builder.get_relative_antenna_conect_positions(quad_settings["antenna_rotation"])

        if quad_settings["couple_KID_to_ANT"]:
            mask_builder.connect_ants_to_KIDs(
                horn_center_x,
                horn_center_y,
                relative_antena_conect_positions,
                relative_kid_positions,
                quad_settings["antenna_rotation"],
                **quad_settings,
            )

        if quad_settings["add_filter_bank"]:
            absolute_filter_bank_conect_points = mask_builder.add_filter_bank_and_get_conection_points(
                horn_center_x,
                horn_center_y,
                with_combiner=quad_settings["add_filter_bank"]["with_combiner"],
                with_crossover=quad_settings["add_filter_bank"]["with_crossover"],
                only_1_pol=quad_settings["add_filter_bank"]["only_1_pol"],
                **quad_settings,
            )

            # Connecting the filter bank to the KID Meanders.
            only_1_pol_no_comb = quad_settings["add_filter_bank"]["only_1_pol"] and (not quad_settings["add_filter_bank"]["with_combiner"])

            # 'temp'
            connect_filter_bank_to_KIDs_kwargs = {}
            if quad_settings["resonator_type"] == SoukResonatorType.CPW_COUPLED_V1:
                connect_filter_bank_to_KIDs_kwargs["KID_conection_linewidth"] = 5
                connect_filter_bank_to_KIDs_kwargs["kid_side_bend_radius"] = 20
                connect_filter_bank_to_KIDs_kwargs["filt_side_bend_radius"] = 20
                connect_filter_bank_to_KIDs_kwargs["KID_extend_out_distance"] = 20
                connect_filter_bank_to_KIDs_kwargs["filter_extend_out_distance"] = 20
            # 'temp'

            mask_builder.connect_filter_bank_to_KIDs(
                horn_center_x,
                horn_center_y,
                absolute_filter_bank_conect_points,
                relative_kid_positions,
                only_1_pol_no_comb=only_1_pol_no_comb,
                **connect_filter_bank_to_KIDs_kwargs,
            )

            # Connecting the end of the antennas to the filter bank.
            if quad_settings["add_filter_bank"]["only_1_pol"]:
                terminate_ants = ["L", "R"]
            else:
                terminate_ants = []

            mask_builder.connect_ants_to_filter_bank(
                horn_center_x,
                horn_center_y,
                relative_antena_conect_positions,
                quad_settings["antenna_rotation"],
                terminate_ants=terminate_ants,
                add_dielectric_under_conections=True,
                **quad_settings,
            )

        if quad_settings["add_top_choke_features"]:
            mask_builder.add_top_choke_features(
                horn_center_x,
                horn_center_y,
                **quad_settings,
            )

        if quad_settings["add_bot_choke_features"]:
            mask_builder.add_bottom_choke_features(
                horn_center_x,
                horn_center_y,
                relative_kid_positions,
                resonator_types,
                **quad_settings,
            )
        # End of quad loop

    ###########################################################################
    # Adding the feedline through the test quad chip and the ports at the top and bottom
    feedline_running_list = mask_builder.get_feedline_running_list(
        feedline_pass_through_points,
        init_direction="right",
        **test_chip_settings,
    )

    feedline_angle_top = -pi / 2
    feedline_angle_bot = pi / 2

    if test_chip_settings["add_ports"]:
        feedline_connect_point_top = mask_builder.add_port_and_get_connection_point(
            test_chip_cetner_xy[0],
            test_chip_cetner_xy[1] + (TEST_CHIP_HEIGHT / 2),
            feedline_angle_top,
        )
        feedline_connect_point_bot = mask_builder.add_port_and_get_connection_point(
            test_chip_cetner_xy[0],
            test_chip_cetner_xy[1] - (TEST_CHIP_HEIGHT / 2),
            feedline_angle_bot,
        )
    else:
        feedline_connect_point_top = [
            test_chip_cetner_xy[0],
            test_chip_cetner_xy[1] + (TEST_CHIP_HEIGHT / 2),
        ]
        feedline_connect_point_bot = [
            test_chip_cetner_xy[0],
            test_chip_cetner_xy[1] - (TEST_CHIP_HEIGHT / 2),
        ]

    # Snaking the feedline around to connect nicely
    # top
    delta_y_to_feedline_end_top = np.abs(feedline_running_list[-1][1] - feedline_connect_point_top[1])
    delta_y_to_feedline_end_bot = np.abs(feedline_running_list[0][1] - feedline_connect_point_bot[1])

    feedline_running_list.append(
        [
            feedline_running_list[-1][0],
            feedline_running_list[-1][1] + (delta_y_to_feedline_end_top / 2),
        ]
    )
    feedline_running_list.append(
        [
            feedline_connect_point_top[0] + (delta_y_to_feedline_end_top / 2) * cos(feedline_angle_top),
            feedline_connect_point_top[1] + (delta_y_to_feedline_end_top / 2) * sin(feedline_angle_top),
        ]
    )
    feedline_running_list.append(feedline_connect_point_top)

    # bot
    feedline_running_list.insert(0, feedline_connect_point_bot)
    feedline_running_list.insert(
        1,
        [
            feedline_connect_point_bot[0] + (delta_y_to_feedline_end_bot / 2) * cos(feedline_angle_bot),
            feedline_connect_point_bot[1] + (delta_y_to_feedline_end_bot / 2) * sin(feedline_angle_bot),
        ],
    )
    feedline_running_list.insert(
        2,
        [
            feedline_running_list[2][0],
            feedline_running_list[2][1] - delta_y_to_feedline_end_bot / 2,
        ],
    )

    mask_builder.add_feedline_and_dielectric(
        feedline_running_list,
        add_bridges=True,
        points_to_avoid_bridging=feedline_pass_through_points,
        **test_chip_settings,
    )

    pretty_print(
        f"Added SOUK Test Chip to mask at:\n    x: {test_chip_cetner_xy[0]}\n    y: {test_chip_cetner_xy[1]}",
        color=TextColor.BRIGHT_BLUE,
    )
    return (feedline_running_list[0], feedline_running_list[-1])
