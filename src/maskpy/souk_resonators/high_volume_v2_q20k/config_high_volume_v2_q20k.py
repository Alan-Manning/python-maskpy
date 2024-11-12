def get_default_resonator_config() -> dict:
    """Get the default config for the HIGH_VOLUME_V2A_Q20K resonator.

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
        "meander_right_height_1": 49,
        "meander_right_height_2": 685,
        "meander_right_height_3": 704,
        "meander_right_height_4": 75,
        "meander_left_height_1": 30,
        "meander_left_height_2": 685,
        "meander_left_height_3": 742,
        "meander_left_height_4": 56,
        "meander_right_width_1": 966,
        "meander_right_width_2": 19,
        "meander_right_width_3": 985,
        "meander_left_width_1": 929,
        "meander_left_width_2": 57,
        "meander_left_width_3": 986,
        "frame_meander_lw": 2,
        "frame_meander_corner_bend_radius": 6,
        "left_frame_meander_width_1": 974,
        "left_frame_meander_width_2": 974,
        "left_frame_meander_height_1": 12,
        "left_frame_meander_height_2": 12,
        "right_frame_meander_width_1": 974,
        "right_frame_meander_width_2": 974,
        "right_frame_meander_height_1": 12,
        "right_frame_meander_height_2": 12,
        "meander_step_back_from_frame": 30,
        "ant_pad_box_width": 5,
        "ant_pad_box_height": 13,
        "frame_bot_lw": 8,
        "frame_bot_left_offset": 996,
        "frame_bot_right_offset": 996,
        "frame_left_lw": 8,
        "frame_left_height": 400,
        "frame_right_lw": 8,
        "frame_right_height": 400,
        "frame_meander_cover_box_width": 5,
        "frame_meander_cover_box_height": 48,
        "extra_frame_meander_cover_box_width": 0,
        "extra_frame_meander_cover_box_height": 0,
        "coupler_frame_left_lw": 10,
        "coupler_frame_left_height": 38,
        "coupler_frame_top_lw": 6,
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
        "coupler_lw": 6,
        "left_coupler_frame_to_feed_distance": 164,
        "text_size": 90,
        "text_x_offset": 800,
        "text_y_offset": 900,
        "cutout_bot_offset": 15,
        "cutout_left_offset": 50,
        "cutout_right_offset": 50,
        "cutout_top_offset": 23,
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
