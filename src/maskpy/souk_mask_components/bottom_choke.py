from collections.abc import Sequence

import gdspy
from numpy import cos, sin

from .. import mask_builder_utils as mbu
from ..souk_mask_configs import SoukMaskConfig, get_mask_default_config

# from ..souk_resonators import SoukResonatorType


def add_bottom_choke_backshort_hole(
    mask_builder,
    x: float | int,
    y: float | int,
    rel_kid_positions: Sequence[Sequence[float | int]],
    # resonator_types: Sequence[SoukResonatorType],
    resonator_types,
    bottom_choke_config_override: dict[str, float | int] | None = None,
    resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None],
) -> None:
    """Adds the bottom choke waveguide hole cutout for a horn block centered
    on the x,y given to the mask with dimensions from the default config or
    from the override if given.

    Parameters
    ----------
    x: float | int
        The x position to place the bottom choke features.

    y: float | int
        The y position to place the bottom choke features.

    rel_kid_positions: Sequence[Sequence[float | int]]
        Sequence of [x, y] lists that describe where the base of each KIDs meander
        is placed relative to the center of the antenna.
        Expected order is top_left, top_right, bot_left, bot_right

    resonator_types: Sequence[SoukResonatorType]
        This is the type of resonators drawn. The values accepted here are
        members of the SoukResonatorType enum. The order of the values
        passed in will be attributed to each KID and should be the same
        order as the rel_kid_positions, TL, TR, BL,
        BR.

    KwArgs
    ------
    bottom_choke_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None]
        This is a Sequence of 4 optional override dictionarys containing key
        value pairs for variable name and that variable's value respectively.
        Any keys required that do not exist in this dict will be got from the
        default config. If extra keys that are not expected are provided a
        warning will be printed but nothing is done with those.
    """
    bottom_choke_config = get_mask_default_config(SoukMaskConfig.BOTTOM_CHOKE, config_override=bottom_choke_config_override)

    backshort_circle_hole_radius = bottom_choke_config["backshort_circle_hole_radius"]

    backshort_circle_hole = gdspy.Round(
        [x, y],
        backshort_circle_hole_radius,
        layer=mask_builder.layers.Bottom_choke_backshort.number,
        datatype=mask_builder.layers.Bottom_choke_backshort.datatype,
    )

    mask_builder.Main.add(backshort_circle_hole)

    backshort_IDC_offset_top = bottom_choke_config["backshort_IDC_offset_top"]
    backshort_IDC_offset_bot = bottom_choke_config["backshort_IDC_offset_bot"]
    backshort_IDC_offset_left = bottom_choke_config["backshort_IDC_offset_left"]
    backshort_IDC_offset_right = bottom_choke_config["backshort_IDC_offset_right"]

    xy_signs = (
        (-1, +1),
        (+1, +1),
        (-1, -1),
        (+1, -1),
    )

    for rel_kid_position, resonator_type, resonator_config_override, (x_sign, y_sign) in zip(
        rel_kid_positions, resonator_types, resonator_config_overrides, xy_signs
    ):

        (IDC_cutout_width, IDC_cutout_height) = mask_builder.get_width_height_of_resonator_IDC_section(
            resonator_type,
            resonator_config_override=resonator_config_override,
        )
        tot_KID_height = mask_builder.get_total_height_of_resonator(
            resonator_type,
            resonator_config_override=resonator_config_override,
        )

        horizontal_offset_from_rel_kid_position_to_center_of_IDC_cutout = tot_KID_height - (IDC_cutout_height / 2)

        y_offset_between_rel_kid_pos_and_center_of_IDC_cutout = 0

        if resonator_type == "cpw_coupled_v1":
            y_offset_between_rel_kid_pos_and_center_of_IDC_cutout = -750  # TODO Temp adjustment for cpw coupled resonators.

        center_of_IDC_cutout_xy = [
            x + rel_kid_position[0] + x_sign * horizontal_offset_from_rel_kid_position_to_center_of_IDC_cutout,
            y + rel_kid_position[1] + y_sign * y_offset_between_rel_kid_pos_and_center_of_IDC_cutout,
        ]

        # Making the bottom choke IDC cutouts
        IDC_cutout = gdspy.Rectangle(
            [
                center_of_IDC_cutout_xy[0] - (IDC_cutout_height / 2) - backshort_IDC_offset_left,
                center_of_IDC_cutout_xy[1] - (IDC_cutout_width / 2) - backshort_IDC_offset_bot,
            ],
            [
                center_of_IDC_cutout_xy[0] + (IDC_cutout_height / 2) + backshort_IDC_offset_right,
                center_of_IDC_cutout_xy[1] + (IDC_cutout_width / 2) + backshort_IDC_offset_top,
            ],
            layer=mask_builder.layers.Bottom_choke_backshort.number,
            datatype=mask_builder.layers.Bottom_choke_backshort.datatype,
        )
        mask_builder.Main.add(IDC_cutout)

    return


def add_bottom_choke_wave_guide_hole(
    mask_builder,
    x: float | int,
    y: float | int,
    bottom_choke_config_override: dict[str, float | int] | None = None,
) -> None:
    """Adds the bottom choke waveguide hole cutout for a horn block centered
    on the x,y given to the mask with dimensions from the default config or
    from the override if given.

    Parameters
    ----------
    x: float | int
        The x position to place the bottom choke features.

    y: float | int
        The y position to place the bottom choke features.

    KwArgs
    ------
    bottom_choke_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.
    """

    bottom_choke_config = get_mask_default_config(SoukMaskConfig.BOTTOM_CHOKE, config_override=bottom_choke_config_override)

    # Making the bottom choke waveguide hole
    bottom_choke_wave_guide_hole = gdspy.Round(
        [x, y],
        bottom_choke_config["wave_guide_hole_radius"],
        layer=mask_builder.layers.Bottom_choke_waveguide_hole.number,
        datatype=mask_builder.layers.Bottom_choke_waveguide_hole.datatype,
    )
    mask_builder.Main.add(bottom_choke_wave_guide_hole)
    return


def add_bottom_choke_pads(
    mask_builder,
    x: float | int,
    y: float | int,
    bottom_choke_config_override: dict[str, float | int] | None = None,
) -> None:
    """Adds the bottom choke pads for a horn block centered on the x,y given
    to the mask with dimensions from the default config or from the override
    if given.

    Parameters
    ----------
    x: float | int
        The x position to place the bottom choke features.

    y: float | int
        The y position to place the bottom choke features.

    KwArgs
    ------
    bottom_choke_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.
    """

    bottom_choke_config = get_mask_default_config(SoukMaskConfig.BOTTOM_CHOKE, config_override=bottom_choke_config_override)

    pad_radius = bottom_choke_config["pad_radius"]
    pad_x_offset_from_center = bottom_choke_config["pad_x_offset_from_center"]
    pad_y_offset_from_center = bottom_choke_config["pad_y_offset_from_center"]

    xy_signs = (
        (-1, +1),
        (+1, +1),
        (-1, -1),
        (+1, -1),
    )

    for x_sign, y_sign in xy_signs:
        # Making the bottom choke pads
        pad_xy = [
            x + x_sign * pad_x_offset_from_center,
            y + y_sign * pad_y_offset_from_center,
        ]
        pad_round = gdspy.Round(
            pad_xy,
            pad_radius,
            layer=mask_builder.layers.Bottom_choke_pads.number,
            datatype=mask_builder.layers.Bottom_choke_pads.datatype,
        )
        mask_builder.Main.add(pad_round)

    return


def add_bottom_choke_IDC_holes(
    mask_builder,
    x: float | int,
    y: float | int,
    rel_kid_positions: Sequence[Sequence[float | int]],
    # resonator_types: Sequence[SoukResonatorType],
    resonator_types,
    bottom_choke_config_override: dict[str, float | int] | None = None,
    resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None],
) -> None:
    """Adds the bottom choke IDC hole cutouts for a horn block centered on the
    x,y given to the mask with dimensions from the default config or from the
    override if given.

    Parameters
    ----------
    x: float | int
        The x position to place the bottom choke features.

    y: float | int
        The y position to place the bottom choke features.

    rel_kid_positions: Sequence[Sequence[float | int]]
        Sequence of [x, y] lists that describe where the base of each KIDs meander
        is placed relative to the center of the antenna.
        Expected order is top_left, top_right, bot_left, bot_right

    resonator_types: Sequence[SoukResonatorType]
        This is the type of resonators drawn. The values accepted here are
        members of the SoukResonatorType enum. The order of the values
        passed in will be attributed to each KID and should be the same
        order as the rel_kid_positions, TL, TR, BL,
        BR.

    KwArgs
    ------
    bottom_choke_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None]
        This is a Sequence of 4 optional override dictionarys containing key
        value pairs for variable name and that variable's value respectively.
        Any keys required that do not exist in this dict will be got from the
        default config. If extra keys that are not expected are provided a
        warning will be printed but nothing is done with those.
    """

    bottom_choke_config = get_mask_default_config(SoukMaskConfig.BOTTOM_CHOKE, config_override=bottom_choke_config_override)

    IDC_cutout_offset_top = bottom_choke_config["IDC_cutout_offset_top"]
    IDC_cutout_offset_bot = bottom_choke_config["IDC_cutout_offset_bot"]
    IDC_cutout_offset_left = bottom_choke_config["IDC_cutout_offset_left"]
    IDC_cutout_offset_right = bottom_choke_config["IDC_cutout_offset_right"]

    xy_signs = (
        (-1, +1),
        (+1, +1),
        (-1, -1),
        (+1, -1),
    )

    for rel_kid_position, resonator_type, resonator_config_override, (x_sign, y_sign) in zip(
        rel_kid_positions, resonator_types, resonator_config_overrides, xy_signs
    ):

        (IDC_cutout_width, IDC_cutout_height) = mask_builder.get_width_height_of_resonator_IDC_section(
            resonator_type,
            resonator_config_override=resonator_config_override,
        )
        tot_KID_height = mask_builder.get_total_height_of_resonator(
            resonator_type,
            resonator_config_override=resonator_config_override,
        )

        horizontal_offset_from_rel_kid_position_to_center_of_IDC_cutout = tot_KID_height - (IDC_cutout_height / 2)

        y_offset_between_rel_kid_pos_and_center_of_IDC_cutout = 0

        if resonator_type == "cpw_coupled_v1":
            y_offset_between_rel_kid_pos_and_center_of_IDC_cutout = -750  # TODO Temp adjustment for cpw coupled resonators.

        center_of_IDC_cutout_xy = [
            x + rel_kid_position[0] + x_sign * horizontal_offset_from_rel_kid_position_to_center_of_IDC_cutout,
            y + rel_kid_position[1] + y_sign * y_offset_between_rel_kid_pos_and_center_of_IDC_cutout,
        ]

        # Making the bottom choke IDC cutouts
        IDC_cutout = gdspy.Rectangle(
            [
                center_of_IDC_cutout_xy[0] - (IDC_cutout_height / 2) - IDC_cutout_offset_left,
                center_of_IDC_cutout_xy[1] - (IDC_cutout_width / 2) - IDC_cutout_offset_bot,
            ],
            [
                center_of_IDC_cutout_xy[0] + (IDC_cutout_height / 2) + IDC_cutout_offset_right,
                center_of_IDC_cutout_xy[1] + (IDC_cutout_width / 2) + IDC_cutout_offset_top,
            ],
            layer=mask_builder.layers.Bottom_choke_IDC_hole.number,
            datatype=mask_builder.layers.Bottom_choke_IDC_hole.datatype,
        )
        mask_builder.Main.add(IDC_cutout)

    return
