from collections.abc import Sequence
from datetime import datetime
from typing import Literal

import gdspy
from numpy import pi

from .. import mask_builder_utils as mbu
from ..logging.pretty_print import TextColor, pretty_print, styled_type_error
from ..souk_mask_configs import SoukMaskConfig, get_mask_default_config
from . import add_fancy_text


def add_test_chip_quadrent_boundary_and_get_horn_positions(
    mask_builder,
    quadrent_center_xy: Sequence[float | int],
    test_chip_quad_config_override: dict[str, float | int] | None = None,
    bottom_left_text: str = "",
    top_right_text: str = "",
    top_left_text: str = "",
    time_stamp_position: Literal["BL", "BR", "ML", "MR", "TL", "TR"] | Sequence[Literal["BL", "BR", "ML", "MR", "TL", "TR"]] | None = (
        "TL",
        "TR",
        "BL",
        "BR",
    ),
    cardiff_logo: bool = True,
    souk_logo: bool = True,
    top_right_label_window: bool = True,
    return_outer_poly_points: bool = False,
    add_center_pin_cutout: bool = True,
    add_slotted_pin_cutout: bool = True,
    dice_line_tab_positions: list[str] = ["left"],
    add_groundplane_under_test_chip: bool = True,
    add_SiN_dep_under_test_chip: bool = False,
    SiN_dep_under_test_chip_edge_offset: float | int = 0,
    add_SiN_membrane_under_test_chip: bool = False,
    add_SiO_under_test_chip: bool = False,
    **kwargs,
) -> list[list[float]] | tuple[list[list[float]], list[list[float]]]:
    """Make the boundary for a test chip quad. Adds a centered pin hole and
    slotted pin hole.

    Parameters
    ----------
    quadrent_center_xy: Sequence[float | int]
        list containing the [x,y] coordinates for the center of the test
        chip quad.

    KwArgs
    ------
    test_chip_quad_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    bottom_left_text: str = \"\"
        Text to add the bottom left of the chip. Default is a blank str
        which will not add any text. When not empty this text will be
        placed offset from the bot left corner at a textsize and offset
        specified in the config.

    top_right_text: str = \"\"
        Text to add the top right of the chip. Default is a blank str
        which will not add any text. When not empty this text will be
        placed offset from the top right corner at a textsize and offset
        specified in the config.

    top_left_text: str = \"\"
        Text to add the top left of the chip. Default is a blank str which
        will not add any text. When not empty this text will be placed
        offset from the top left corner at a textsize and offset
        specified in the config.

    time_stamp_position: Literal["BL", "BR", "ML", "MR", "TL", "TR"] | None = "ML",
        This will by default add a time stamp to the MidLeft side of the test
        chip. This can be disabled by passing None as a value. The other
        positions availible are ["BL", "BR", "ML", "MR", "TL", "TR"] where the
        first letter refer to the verical and horizontal position respectively.

    cardiff_logo: bool = True
        Default True will add the cardiff logo to the nb groundplane in the
        bottom right corner of the test chip quad. The size of this is
        specified in the config.

    souk_logo: bool = True
        Default True will add the souk logo to the nb groundplane in the
        bottom right corner of the test chip quad. The size of this is
        specified in the config.

    top_right_label_window: bool = True
        Default True will add a window cutout in the nb groundplane in the
        top right corner of the test chip quad. The size of this window
        is specified in the config.

    return_outer_poly_points: bool = False
        Default False. If true this will return the outer polygon points
        around the test chip quadrent with extra exclusion.

    add_center_pin_cutout: bool = True
        Default True. If false the center pin hole for the test chip will
        not be added to the mask.

    add_slotted_pin_cutout: bool = True
        Default True. If false the slotted pin hole for the test chip will
        not be added to the mask.

    dice_line_tab_positions: list[str]
        The edge positions that the dicing line should be tabbed. Any edges
        not specified will be solid dicing lines. This should be a list of
        strings. The default is a list of just ["left"], but can take str
        values "top", "left", "bot", "right", or "all".

    add_groundplane_under_test_chip: bool = True
        Default True. If false no groundplane will be added below the test
        chip boundary.

    add_SiN_dep_under_test_chip: bool = False
        Default False. If True an SiN dep will be added below the test
        chip boundary.

    SiN_dep_under_test_chip_edge_offset: float | int = 0,
        Default 0. When non-zero the SiN dep will be added below the test
        chip boundary with an offset from the boundary by the amount specified.
        This will only draw if `add_SiN_dep_under_test_chip` is set to `True`.
        Value here can be negaitve for an oversize.

    add_SiN_membrane_under_test_chip: bool = False,
        Default False. If True an SiN membrane will be added below the test
        chip boundary.

    add_SiO_under_test_chip: bool = False,
        Default False. If True an SiO layer will be added below the test
        chip boundary.

    Returns
    -------
    horn_centers: list
        list of [x,y] lists that define the coordinates for the horn
        centers for the 4 horns on the test chip quad. The coords are for
        the top left, top right, bot left, bot right horns respectively.

    outer_poly_points: list
        list of [x,y] lists defining the coordinates of the outer polygon.
        **Note This is only returned if return_outer_poly_points KwArg is
        True**.
    """
    config = get_mask_default_config(SoukMaskConfig.TEST_CHIP_QUAD, config_override=test_chip_quad_config_override)

    test_chip_quad_width = config["test_chip_quad_width"]
    test_chip_quad_height = config["test_chip_quad_height"]

    horn_offset_from_center_x = config["horn_offset_from_center_x"]
    horn_offset_from_center_y = config["horn_offset_from_center_y"]

    top_right_label_window_offset_x = config["top_right_label_window_offset_x"]
    top_right_label_window_offset_y = config["top_right_label_window_offset_y"]
    top_right_label_window_width = config["top_right_label_window_width"]
    top_right_label_window_height = config["top_right_label_window_height"]

    top_right_chip_corner = [
        quadrent_center_xy[0] + (test_chip_quad_width / 2),
        quadrent_center_xy[1] + (test_chip_quad_height / 2),
    ]

    top_left_chip_corner = [
        quadrent_center_xy[0] - (test_chip_quad_width / 2),
        quadrent_center_xy[1] + (test_chip_quad_height / 2),
    ]

    bottom_left_chip_corner = [
        quadrent_center_xy[0] - (test_chip_quad_width / 2),
        quadrent_center_xy[1] - (test_chip_quad_height / 2),
    ]

    bottom_right_chip_corner = [
        quadrent_center_xy[0] + (test_chip_quad_width / 2),
        quadrent_center_xy[1] - (test_chip_quad_height / 2),
    ]

    middle_left_chip_side = [
        quadrent_center_xy[0] - (test_chip_quad_width / 2),
        quadrent_center_xy[1],
    ]
    middle_right_chip_side = [
        quadrent_center_xy[0] + (test_chip_quad_width / 2),
        quadrent_center_xy[1],
    ]

    if add_groundplane_under_test_chip:
        groundplane_rect = gdspy.Rectangle(
            top_right_chip_corner,
            bottom_left_chip_corner,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_positives.add(groundplane_rect)
        # mask_builder.Main.add(groundplane_rect)

    if add_SiN_dep_under_test_chip:
        if not isinstance(SiN_dep_under_test_chip_edge_offset, float | int):
            styled_type_error(
                SiN_dep_under_test_chip_edge_offset,
                "SiN_dep_under_test_chip_edge_offset",
                float | int,
            )
        SiN_dep_rect = gdspy.Rectangle(
            [
                top_right_chip_corner[0] - SiN_dep_under_test_chip_edge_offset,
                top_right_chip_corner[1] - SiN_dep_under_test_chip_edge_offset,
            ],
            [
                bottom_left_chip_corner[0] + SiN_dep_under_test_chip_edge_offset,
                bottom_left_chip_corner[1] + SiN_dep_under_test_chip_edge_offset,
            ],
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_positives.add(SiN_dep_rect)

    if add_SiN_membrane_under_test_chip:
        SiN_mem_rect = gdspy.Rectangle(
            top_right_chip_corner,
            bottom_left_chip_corner,
            layer=mask_builder.layers.SiN_Membrane.number,
            datatype=mask_builder.layers.SiN_Membrane.datatype,
        )
        mask_builder.silicon_nitride_membrane_positives.add(SiN_mem_rect)

    if add_SiO_under_test_chip:
        SiO_rect = gdspy.Rectangle(
            top_right_chip_corner,
            bottom_left_chip_corner,
            layer=mask_builder.layers.SiO.number,
            datatype=mask_builder.layers.SiO.datatype,
        )
        mask_builder.silicon_oxide_positives.add(SiO_rect)

    # Adding the bottom_left_text if True
    if bottom_left_text != "":
        bottom_left_text_offset_x = config["bottom_left_text_offset_x"]
        bottom_left_text_offset_y = config["bottom_left_text_offset_y"]
        bottom_left_text_size = config["bottom_left_text_size"]
        bottom_left_text_x = bottom_left_chip_corner[0] + bottom_left_text_offset_x
        bottom_left_text_y = bottom_left_chip_corner[1] + bottom_left_text_offset_y

        add_fancy_text(
            mask_builder,
            bottom_left_text,
            bottom_left_text_x,
            bottom_left_text_y,
            bottom_left_text_size,
            mask_builder.layers.SiN_dep,
            bb_cutout_in_grnd=True,
            bb_cutout_in_sin_dep=True,
            vertical_align="above",
        )

    # Adding the top_left_text if True
    if top_left_text != "":
        top_left_text_offset_x = config["top_left_text_offset_x"]
        top_left_text_offset_y = config["top_left_text_offset_y"]
        top_left_text_size = config["top_left_text_size"]
        top_left_text_x = top_left_chip_corner[0] + top_left_text_offset_x
        top_left_text_y = top_left_chip_corner[1] - top_left_text_offset_y

        add_fancy_text(
            mask_builder,
            top_left_text,
            top_left_text_x,
            top_left_text_y,
            top_left_text_size,
            [
                mask_builder.layers.Aluminium,
                mask_builder.layers.Aluminium_Direct,
            ],
            bb_cutout_in_grnd=True,
            bb_cutout_in_sin_dep=True,
            vertical_align="below",
        )

    # Adding the top_right_text if True
    if top_right_text != "":
        top_right_text_offset_x = config["top_right_text_offset_x"]
        top_right_text_offset_y = config["top_right_text_offset_y"]
        top_right_text_size = config["top_right_text_size"]
        top_right_text_x = top_right_chip_corner[0] - top_right_text_offset_x
        top_right_text_y = top_right_chip_corner[1] - top_right_text_offset_y

        add_fancy_text(
            mask_builder,
            top_right_text,
            top_right_text_x,
            top_right_text_y,
            top_right_text_size,
            mask_builder.layers.Nb_Antenna,
            bb_cutout_in_grnd=True,
            bb_cutout_in_sin_dep=True,
            vertical_align="below",
            horizontal_align="end",
        )

    # Adding the time_stamp_position if not None
    if time_stamp_position is not None:
        time_stamp_positions: list[str] = []
        valid_time_stamp_positions = ["BL", "BR", "ML", "MR", "TL", "TR"]
        DEFAULT_TIME_STAMP_POSITION = "ML"

        if isinstance(time_stamp_position, str):
            time_stamp_positions.append(time_stamp_position)
        elif isinstance(time_stamp_position, Sequence):
            for index, pos in enumerate(time_stamp_position):
                if pos not in valid_time_stamp_positions:
                    pretty_print(
                        f"time_stamp_position[{index}] `{time_stamp_position}` not recognised. Should be one of {valid_time_stamp_positions}. Defaulting to `ML`.",
                        color=TextColor.WARNING,
                    )
                    time_stamp_positions.append(DEFAULT_TIME_STAMP_POSITION)
                else:
                    time_stamp_positions.append(pos)

        else:
            styled_type_error(
                time_stamp_position,
                "time_stamp_position",
                Literal["BL", "BR", "ML", "MR", "TL", "TR"] | Sequence[Literal["BL", "BR", "ML", "MR", "TL", "TR"]] | None,
            )
        time_stamp_text = f"Generated: " + mbu.get_time(
            datetime.now(),
            time_format="%Y#-%m#-%d# %H#:%M",
        )

        for time_stamp_pos in time_stamp_positions:
            line_height = 1.2
            DEFAULT_TIME_STAMP_SIZE = 450
            match time_stamp_pos:
                case "BL":
                    bottom_left_text_offset_x = config["bottom_left_text_offset_x"]
                    bottom_left_text_offset_y = config["bottom_left_text_offset_y"]
                    time_stamp_size = config["bottom_left_text_size"]
                    time_stamp_x = bottom_left_chip_corner[0] + bottom_left_text_offset_x
                    if bottom_left_text != "":
                        time_stamp_y = bottom_left_chip_corner[1] + bottom_left_text_offset_y + (line_height * time_stamp_size)
                    else:
                        time_stamp_y = bottom_left_chip_corner[1] + bottom_left_text_offset_y
                    time_stamp_rot = 0
                    time_stamp_va = "above"
                    time_stamp_ha = "start"
                    time_stamp_layer = mask_builder.layers.SiN_dep
                case "BR":
                    bottom_right_text_offset_x = 1000
                    bottom_right_text_offset_y = 265
                    time_stamp_size = DEFAULT_TIME_STAMP_SIZE
                    time_stamp_x = bottom_right_chip_corner[0] - bottom_right_text_offset_x
                    time_stamp_y = bottom_right_chip_corner[1] + bottom_right_text_offset_y
                    time_stamp_rot = 0
                    time_stamp_va = "above"
                    time_stamp_ha = "end"
                    time_stamp_layer = mask_builder.layers.Nb_Groundplane
                case "ML":
                    middle_left_text_offset_x = 1000
                    middle_left_text_offset_y = 0
                    time_stamp_size = DEFAULT_TIME_STAMP_SIZE
                    time_stamp_x = middle_left_chip_side[0] + middle_left_text_offset_x
                    time_stamp_y = middle_left_chip_side[1] + middle_left_text_offset_y
                    time_stamp_rot = pi / 2
                    time_stamp_va = "below"
                    time_stamp_ha = "center"
                    time_stamp_layer = mask_builder.layers.Aluminium
                case "MR":
                    middle_right_text_offset_x = 1000
                    middle_right_text_offset_y = 0
                    time_stamp_size = DEFAULT_TIME_STAMP_SIZE
                    time_stamp_x = middle_right_chip_side[0] - middle_right_text_offset_x
                    time_stamp_y = middle_right_chip_side[1] + middle_right_text_offset_y
                    time_stamp_rot = -pi / 2
                    time_stamp_va = "below"
                    time_stamp_ha = "center"
                    time_stamp_layer = mask_builder.layers.Aluminium
                case "TL":
                    top_left_text_offset_x = config["top_left_text_offset_x"]
                    top_left_text_offset_y = config["top_left_text_offset_y"]
                    time_stamp_size = config["top_left_text_size"]
                    time_stamp_x = top_left_chip_corner[0] + top_left_text_offset_x
                    if top_left_text != "":
                        time_stamp_y = top_left_chip_corner[1] - top_left_text_offset_y - (line_height * time_stamp_size)
                    else:
                        time_stamp_y = top_left_chip_corner[1] - top_left_text_offset_y
                    time_stamp_rot = 0
                    time_stamp_va = "below"
                    time_stamp_ha = "start"
                    time_stamp_layer = [
                        mask_builder.layers.Aluminium,
                        mask_builder.layers.Aluminium_Direct,
                    ]

                case "TR":
                    top_right_text_offset_x = config["top_right_text_offset_x"]
                    top_right_text_offset_y = config["top_right_text_offset_y"]
                    time_stamp_size = config["top_right_text_size"]
                    time_stamp_x = top_right_chip_corner[0] - top_right_text_offset_x
                    if top_right_text != "":
                        time_stamp_y = top_right_chip_corner[1] - top_right_text_offset_y - (line_height * time_stamp_size)
                    else:
                        time_stamp_y = top_right_chip_corner[1] - top_right_text_offset_y
                    time_stamp_rot = 0
                    time_stamp_va = "below"
                    time_stamp_ha = "end"
                    time_stamp_layer = mask_builder.layers.Nb_Antenna
                case _:
                    raise ValueError("time_stamp_position not valid.")

            add_fancy_text(
                mask_builder,
                time_stamp_text,
                time_stamp_x,
                time_stamp_y,
                time_stamp_size,
                time_stamp_layer,
                rotation=time_stamp_rot,
                bb_cutout_in_grnd=True,
                bb_cutout_in_sin_dep=True,
                vertical_align=time_stamp_va,
                horizontal_align=time_stamp_ha,
                usetex=True,
            )

    # Adding the top_right_label_window if True
    if top_right_label_window:
        top_right_label_window_rect = gdspy.Rectangle(
            [top_right_chip_corner[0] - top_right_label_window_offset_x, top_right_chip_corner[1] - top_right_label_window_offset_y],
            [
                top_right_chip_corner[0] - top_right_label_window_offset_x - top_right_label_window_width,
                top_right_chip_corner[1] - top_right_label_window_offset_y - top_right_label_window_height,
            ],
        )
        mask_builder.ground_plane_cutouts.add(top_right_label_window_rect)

    # Adding the cardiff logo if True
    if cardiff_logo:
        cardiff_logo_size = config["cardiff_logo_size"]  # 2000
        cardiff_logo_offset_x = config["cardiff_logo_offset_x"]  # 1000
        cardiff_logo_offset_y = config["cardiff_logo_offset_y"]  # 1000

        logo_xy = [
            bottom_right_chip_corner[0] - cardiff_logo_offset_x - (cardiff_logo_size / 2),
            bottom_right_chip_corner[1] + cardiff_logo_offset_y + (cardiff_logo_size / 2),
        ]

        mask_builder.add_cardiff_logo(logo_xy, cardiff_logo_size)

    if souk_logo:
        souk_logo_size = config["souk_logo_size"]  # 2000
        souk_logo_offset_x = config["souk_logo_offset_x"]  # 3500
        souk_logo_offset_y = config["souk_logo_offset_y"]  # 3500

        logo_xy = [
            bottom_right_chip_corner[0] - souk_logo_offset_x - (souk_logo_size / 2),
            bottom_right_chip_corner[1] + souk_logo_offset_y + (souk_logo_size / 2),
        ]

        mask_builder.add_souk_logo(logo_xy, souk_logo_size, include_outer_ring_text=False, draw_all_in_one_layer=False)

    # Adding the dicing line in the groundplane and the tabbed dicing line
    mask_builder.add_test_chip_quad_tabbed_dicing_line(
        quadrent_center_xy,
        test_chip_quad_width,
        test_chip_quad_height,
        tab_positions=dice_line_tab_positions,
    )

    groundplane_edge_cuts_top_mid_gap = config["groundplane_edge_cuts_top_mid_gap"]  # 5500
    groundplane_edge_cuts_bot_mid_gap = config["groundplane_edge_cuts_bot_mid_gap"]  # 5500
    groundplane_edge_cuts_right_mid_gap = config["groundplane_edge_cuts_right_mid_gap"]  # 0
    groundplane_edge_cuts_left_mid_gap = config["groundplane_edge_cuts_left_mid_gap"]  # 0

    mask_builder.add_test_chip_quad_groundplane_edge_cuts(
        quadrent_center_xy,
        test_chip_quad_width,
        test_chip_quad_height,
        top_mid_gap=groundplane_edge_cuts_top_mid_gap,
        bot_mid_gap=groundplane_edge_cuts_bot_mid_gap,
        right_mid_gap=groundplane_edge_cuts_right_mid_gap,
        left_mid_gap=groundplane_edge_cuts_left_mid_gap,
    )

    # Adding the center pin and slotted pin
    center_pin_radius = config["center_pin_radius"]  # 2020/2

    center_pin_xy = quadrent_center_xy

    slotted_pin_radius = config["slotted_pin_radius"]  # 2020/2
    slotted_pin_length = config["slotted_pin_length"]  # 1000/2
    slotted_pin_offset_x_from_center_pin = config["slotted_pin_offset_x_from_center_pin"]  # 8000
    slotted_pin_offset_y_from_center_pin = config["slotted_pin_offset_y_from_center_pin"]  # 0

    slotted_pin_xy = [
        quadrent_center_xy[0] + slotted_pin_offset_x_from_center_pin,
        quadrent_center_xy[1] + slotted_pin_offset_y_from_center_pin,
    ]

    if add_center_pin_cutout:
        _ = mask_builder.add_center_pin_and_get_bounding_points(quadrent_center_xy, center_pin_radius)

    if add_slotted_pin_cutout:
        _ = mask_builder.add_slotted_pin_and_get_bounding_points(slotted_pin_xy, slotted_pin_length, slotted_pin_radius, center_pin_xy)

    horn_centers = [
        [quadrent_center_xy[0] - horn_offset_from_center_x, quadrent_center_xy[1] + horn_offset_from_center_y],
        [quadrent_center_xy[0] + horn_offset_from_center_x, quadrent_center_xy[1] + horn_offset_from_center_y],
        [quadrent_center_xy[0] - horn_offset_from_center_x, quadrent_center_xy[1] - horn_offset_from_center_y],
        [quadrent_center_xy[0] + horn_offset_from_center_x, quadrent_center_xy[1] - horn_offset_from_center_y],
    ]

    extra_exclusion_around_test_chip_quadrent = 2000  # config

    if not return_outer_poly_points:
        return horn_centers

    outer_poly_points = [
        [
            quadrent_center_xy[0] - (test_chip_quad_width / 2) - extra_exclusion_around_test_chip_quadrent,
            quadrent_center_xy[1] + (test_chip_quad_height / 2) + extra_exclusion_around_test_chip_quadrent,
        ],
        [
            quadrent_center_xy[0] + (test_chip_quad_width / 2) + extra_exclusion_around_test_chip_quadrent,
            quadrent_center_xy[1] + (test_chip_quad_height / 2) + extra_exclusion_around_test_chip_quadrent,
        ],
        [
            quadrent_center_xy[0] + (test_chip_quad_width / 2) + extra_exclusion_around_test_chip_quadrent,
            quadrent_center_xy[1] - (test_chip_quad_height / 2) - extra_exclusion_around_test_chip_quadrent,
        ],
        [
            quadrent_center_xy[0] - (test_chip_quad_width / 2) - extra_exclusion_around_test_chip_quadrent,
            quadrent_center_xy[1] - (test_chip_quad_height / 2) - extra_exclusion_around_test_chip_quadrent,
        ],
    ]

    return horn_centers, outer_poly_points
