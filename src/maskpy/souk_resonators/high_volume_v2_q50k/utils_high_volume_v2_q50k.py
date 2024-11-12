from ..resonator_types import SoukResonatorType
from ..utils.get_config import get_resonator_config


def _this_resonator_type() -> SoukResonatorType:
    return SoukResonatorType.HIGH_VOLUME_V2_Q50K


def get_total_height_of_resonator(
    resonator_config_override: dict[str, float | int] | None = None,
) -> float:
    """This will get the total height of the resonator from the base of the
    inductive meander to the end of the ground plane cutout at the top of
    the structure.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    total_resonator_height : float
        The total height of the resonator calculated from the config file.
    """

    config = get_resonator_config(_this_resonator_type(), resonator_config_override=resonator_config_override)

    total_resonator_height = (
        config["meander_lw"]
        + config["meander_left_height_1"]
        - config["meander_left_height_2"]
        + config["meander_left_height_3"]
        + config["meander_left_height_4"]
        + config["frame_left_height"]
        + config["coupler_frame_left_height"]
        + config["coupler_gap"]
        + config["coupler_lw"]
        + config["cutout_top_offset"]
        + (3 * config["step_down_distance_between_layers"])
    )
    return total_resonator_height


def get_horizontal_coupler_end_to_meander_base_distance(
    resonator_config_override: dict[str, float | int] | None = None,
) -> float:
    """This will calculate the the horizonatal distance between the end of the
    coupler arm (where it would connect to a feedline) and the center of the
    base of the resonator's inductive meander. This is calculated from the
    config.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    coupler_end_to_meander_base_distance : float | int
        The distance between the coupler end and the center of the base of the
        resonator's inductive meander.
    """
    config = get_resonator_config(_this_resonator_type(), resonator_config_override=resonator_config_override)

    coupler_end_to_meander_base_distance = (
        (config["meander_bot_width"] / 2)
        - (config["meander_lw"] / 2)
        + config["meander_left_width_1"]
        + config["meander_left_width_2"]
        - config["meander_left_width_3"]
        + config["frame_bot_left_offset"]
        - (config["frame_left_lw"] / 2)
        + (config["coupler_frame_left_lw"] / 2)
        + config["left_coupler_frame_to_feed_distance"]
    )
    return coupler_end_to_meander_base_distance


def get_vertical_coupler_center_to_meander_base_distance(
    resonator_config_override: dict[str, float | int] | None = None,
) -> float:
    """This will calculate the the vertical distance between the center of the
    coupler arm and the base of the resonator's inductive meander. This is
    calculated from the config.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    coupler_center_to_meander_base_distance : float | int
        The distance between the coupler center and the base of the resonator's
        inductive meander.
    """
    config = get_resonator_config(_this_resonator_type(), resonator_config_override=resonator_config_override)

    coupler_center_to_meander_base_distance = (
        config["meander_lw"]
        + config["meander_left_height_1"]
        - config["meander_left_height_2"]
        + config["meander_left_height_3"]
        + config["meander_left_height_4"]
        + config["frame_left_height"]
        + config["coupler_frame_left_height"]
        + config["coupler_gap"]
        + (config["coupler_lw"] / 2)
    )
    return coupler_center_to_meander_base_distance


def get_width_and_height_of_IDC_cutout_section(
    resonator_config_override: dict[str, float | int] | None = None,
) -> tuple[float | int, float | int]:
    """Get the total width and height of ground plane cutout around the IDC
    section of the resonator calculated from the config.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    width, height : tuple[float, float]
        This is a tuple containing, in order, the width and the height
        of the resonator's IDC section calculated from the config file.
    """
    config = get_resonator_config(_this_resonator_type(), resonator_config_override=resonator_config_override)

    left_side_width = (
        (config["meander_bot_width"] / 2)
        - (config["meander_lw"] / 2)
        + config["meander_left_width_1"]
        + config["meander_left_width_2"]
        - config["meander_left_width_3"]
        + config["frame_bot_left_offset"]
        + config["cutout_left_offset"]
        + (3 * config["step_down_distance_between_layers"])
    )

    right_side_width = (
        (config["meander_bot_width"] / 2)
        - (config["meander_lw"] / 2)
        + config["meander_right_width_1"]
        + config["meander_right_width_2"]
        - config["meander_right_width_3"]
        + config["frame_bot_right_offset"]
        + config["cutout_right_offset"]
        + (3 * config["step_down_distance_between_layers"])
    )

    width = left_side_width + right_side_width

    height = (
        (3 * config["step_down_distance_between_layers"])
        + config["cutout_bot_offset"]
        + config["frame_left_height"]
        + config["coupler_frame_left_height"]
        + config["coupler_gap"]
        + config["coupler_lw"]
        + config["cutout_top_offset"]
        + (3 * config["step_down_distance_between_layers"])
    )

    return width, height
