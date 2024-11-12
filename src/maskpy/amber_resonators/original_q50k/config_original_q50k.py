def get_default_resonator_config() -> dict:
    """Get the default config for the ORIGINAL_Q50K resonator.

    This is the config used in the muxing for the Q50k design. This means only
    this config is guaranteed to generate Q50k pixels reliably.

    Returns
    -------
    config : dict
        This is a config dictionary with key value pairs as variable name and
        variable values.
    """
    config = {
        "meander_lw": 2,
        "meander_corner_bend_radius": 8,
        #
        "meander_bot_width": 100,
        "meander_right_height_1": 150,
        "meander_right_width_1": 946,
        "meander_right_height_2": 34,
        "meander_right_width_2": -946,
        "meander_right_height_3": 34,
        "meander_right_width_3": 98,
        "meander_right_height_4": 104,
        "meander_right_width_4": -192,
        "meander_right_height_5": 24,
        "meander_right_width_5": 192,
        "meander_right_height_6": 60,
        "meander_right_width_6": -192,
        "meander_right_height_7": 24,
        "meander_right_width_7": 192,
        "meander_right_height_8": 60,
        "meander_right_width_8": -192,
        "meander_right_height_9": 24,
        "meander_right_width_9": 2044,
        #
        "meander_left_height_1": 150,
        "meander_left_width_1": -946,
        "meander_left_height_2": 34,
        "meander_left_width_2": 946,
        "meander_left_height_3": 34,
        "meander_left_width_3": -98,
        "meander_left_height_4": 62,
        "meander_left_width_4": 192,
        "meander_left_height_5": 24,
        "meander_left_width_5": -192,
        "meander_left_height_6": 60,
        "meander_left_width_6": 192,
        "meander_left_height_7": 24,
        "meander_left_width_7": -192,
        "meander_left_height_8": 60,
        "meander_left_width_8": 192,
        "meander_left_height_9": 24,
        "meander_left_width_9": -2040,
        "meander_left_height_10": 43,
        #
        "left_frame_width": 8,
        "left_frame_height_overrun": 8,
        "right_frame_width": 8,
        "top_idc_arm_to_coupler_arm": 62,
        #
        "coupler_arm_lw": 6,
        #
        "idc_arms_offset_from_meander": 38,
        "idc_arm_lw": 4,
        "idc_arm_spacing": 10,
        #
        "coupler_fork_bot_arm_lw": 6,
        "coupler_fork_top_arm_lw": 6,
        "coupler_fork_stub_lw": 6,
        "coupler_fork_x_offset_from_coupler_arm": 10,
        "coupler_fork_bot_arm_y_offset_from_coupler_arm": 10,
        "coupler_fork_top_arm_y_offset_from_coupler_arm": 10,
        "coupler_fork_top_arm_offset_to_ground": 70,
        #
        "grnd_cutout_gap_right": 60,
        "grnd_cutout_gap_left": 60,
        "grnd_cutout_gap_top": 0,
        "grnd_cutout_gap_bot": 45,
        #
        "SiN_dep_cutout_gap_right": 55,
        "SiN_dep_cutout_gap_left": 55,
        "SiN_dep_cutout_gap_top": -5,
        "SiN_dep_cutout_gap_bot": 40,
        #
        "SiO_cutout_gap_right": 45,
        "SiO_cutout_gap_left": 45,
        "SiO_cutout_gap_top": -15,
        "SiO_cutout_gap_bot": 30,
        #
        "SiN_membrane_cutout_gap_right": 50,
        "SiN_membrane_cutout_gap_left": 50,
        "SiN_membrane_cutout_gap_top": -10,
        "SiN_membrane_cutout_gap_bot": 35,
        #
        "backside_check_gap_right": 45,
        "backside_check_gap_left": 45,
        "backside_check_gap_top": -15,
        "backside_check_gap_bot": 30,
        #
        "inductor_cover_gap_right": 60,
        "inductor_cover_gap_left": 60,
        "inductor_cover_gap_top": 0,
        "inductor_cover_gap_bot": 45,
    }

    return config
