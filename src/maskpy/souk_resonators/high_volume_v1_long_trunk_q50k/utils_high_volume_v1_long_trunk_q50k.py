from .config_high_volume_v1_long_trunk_q50k import get_resonator_config


def _get_config_checking_override(config_override: dict[str, float | int] | None) -> dict[str, float | int]:

    if config_override is not None:
        if isinstance(config_override, dict):
            return config_override
        else:
            print(f"\033[93mWarning: config_override not of correct format. Defaulting to using base config\033[0m")
            return get_resonator_config()
    else:
        return get_resonator_config()


def get_total_height_of_resonator(config_override: dict[str, float | int] | None = None) -> float:
    """This will get the total height of the resonator from the base of the
    inductive meander to the end of the ground plane cutout at the top of
    the structure.

    Parameters
    ----------
    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    total_resonator_height : float
        The total height of the resonator calculated from the config file.
    """

    config = _get_config_checking_override(config_override)

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


def get_horizontal_coupler_end_to_meander_base_distance(config_override: dict[str, float | int] | None = None) -> float:
    """This will calculate the the horizonatal distance between the end of the
    coupler arm (where it would connect to a feedline) and the center of the
    base of the resonator's inductive meander. This is calculated from the
    config.

    Parameters
    ----------
    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    coupler_end_to_meander_base_distance : float | int
        The distance between the coupler end and the center of the base of the
        resonator's inductive meander.
    """
    config = _get_config_checking_override(config_override)

    coupler_end_to_meander_base_distance = (
        (config["meander_bot_width"] / 2)
        - (config["meander_lw"] / 2)
        + config["meander_left_width_1"]
        + config["meander_left_width_2"]
        - config["meander_left_width_3"]
        + config["frame_bot_left_width"]
        - (config["frame_left_lw"] / 2)
        + (config["coupler_frame_left_lw"] / 2)
        + config["left_coupler_frame_to_feed_distance"]
    )
    return coupler_end_to_meander_base_distance


def get_vertical_coupler_center_to_meander_base_distance(config_override: dict[str, float | int] | None = None) -> float:
    """This will calculate the the vertical distance between the center of the
    coupler arm and the base of the resonator's inductive meander. This is
    calculated from the config.

    Parameters
    ----------
    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    coupler_center_to_meander_base_distance : float | int
        The distance between the coupler center and the base of the resonator's
        inductive meander.
    """
    config = _get_config_checking_override(config_override)

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


def get_width_and_height_of_IDC_cutout_section(config_override: dict[str, float | int] | None = None) -> tuple[
    float | int,
    float | int,
]:
    """Get the total width and height of ground plane cutout around the IDC
    section of the resonator calculated from the config.

    Parameters
    ----------
    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    width, height : tuple[float, float]
        This is a tuple containing, in order, the width and the height
        of the resonator's IDC section calculated from the config file.
    """
    config = _get_config_checking_override(config_override)

    left_side_width = (
        (config["meander_bot_width"] / 2)
        - (config["meander_lw"] / 2)
        + config["meander_left_width_1"]
        + config["meander_left_width_2"]
        - config["meander_left_width_3"]
        + config["frame_bot_left_width"]
        + config["cutout_left_offset"]
        + (3 * config["step_down_distance_between_layers"])
    )

    right_side_width = (
        (config["meander_bot_width"] / 2)
        - (config["meander_lw"] / 2)
        + config["meander_right_width_1"]
        + config["meander_right_width_2"]
        - config["meander_right_width_3"]
        + config["frame_bot_right_width"]
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
