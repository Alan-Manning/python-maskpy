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
        "meander_corner_bend_radius": 6,
        "meander_bot_width": 18,
        "meander_right_height_1": 41,
        "meander_right_height_2": 24,
        "meander_right_height_3": 73,
        "meander_left_height_1": 24,
        "meander_left_height_2": 58,
        "meander_left_height_3": 56,
        "meander_left_width_1": 564,
        "meander_left_width_2": 564,
        "meander_right_width_1": 565,
        "meander_right_width_2": 565,
        "meander_step_back_from_frame": 0,  ## Added
        "ant_pad_box_width": 5,
        "ant_pad_box_height": 13,
        "frame_bot_lw": 8,
        "frame_bot_left_width": 996,
        "frame_bot_right_width": 996,
        "frame_left_lw": 8,
        "frame_left_height": 400,
        "frame_right_lw": 8,
        "frame_right_height": 400,
        "frame_meander_cover_box_width": 5,
        "frame_meander_cover_box_height": 48,
        "extra_frame_meander_cover_box_width_left": 0,
        "extra_frame_meander_cover_box_width_right": 0,
        "extra_frame_meander_cover_box_height_above": 0,
        "extra_frame_meander_cover_box_height_below": 0,
        "coupler_frame_left_lw": 10,
        "coupler_frame_left_height": 39,
        "coupler_frame_top_lw": 3,
        "IDC_bot_arm_gap": 30,
        "IDC_arm_gap": 8,
        "IDC_arm_lw": 3,
        "No_of_arms": 28,
        "trim_arm_offset_right_side": 380,
        "trim_arm_offset_left_side": 389,
        "trim_arm_lw": 3,
        "trim_arm_length_right_side": 1975,
        "trim_arm_length_left_side": 1975,
        "coupler_gap": 16,
        "coupler_lw": 3,
        "left_coupler_frame_to_feed_distance": 164,
        "text_size": 90,
        "text_x_offset": 800,
        "text_y_offset": 900,
        "text_underline_height": 3,
        "cutout_bot_offset": 15,
        "cutout_left_offset": 50,
        "cutout_right_offset": 50,
        "cutout_top_offset": 25,
        "grndpl_meander_cutout_width": 80,
        "grndpl_meander_cutout_height": 10,
        "grndpl_coupler_cutout_width": 10,
        "grndpln_gap_between_adjacent_resonators": 14,
        "grndpl_coupler_cutout_height": 30,
        "SiO_stepdown_cutout_width": 110,
        "SiO_stepdown_cutout_height": 0,
        "SiN_membrane_stepdown_cutout_width": 100,
        "SiN_membrane_stepdown_cutout_height": 0,
        "step_down_distance_between_layers": 5,
    }

    return config
