def get_default_resonator_config() -> dict[str, float | int]:
    """Get the config for the CPW_COUPLED_V1 resonator.

    This is the config used in the muxing for the CPW_COUPLED_V1 design. This means only
    this config is guaranteed to generate pixels reliably.

    Returns
    -------
    config : dict[str, float | int]
        This is a config dictionary with the required key value pairs as
        variable name and variable values.
    """
    config: dict[str, float | int] = {
        "meander_cpw_center_width": 2,
        "meander_cpw_outer_width": 20,
        "meander_dielectric_under_cpw_width": 20,
        "meander_bend_radius": 11,
        "meander_cpw_bridge_width": 3,
        #
        "meander_right_width_1": +612,
        "meander_right_height_1": 58,
        "meander_right_width_2": -583,
        "meander_right_height_2": 58,
        "meander_right_width_3": +583,
        "meander_right_height_3": 58,
        "meander_right_width_4": -583,
        "meander_right_height_4": 43,
        #
        "meander_left_width_1": -612,
        "meander_left_height_1": 58,
        "meander_left_width_2": +583,
        "meander_left_height_2": 58,
        "meander_left_width_3": -583,
        "meander_left_height_3": 58,
        "meander_left_width_4": +583,
        "meander_left_height_4": 43,
        #
        "ant_pad_box_width": 6,
        "ant_pad_box_height": 13,
        #
        "ant_cpw_center_width": 5,
        "ant_cpw_outer_width": 11,
        "ant_cpw_to_transition_offset_y": 100,
        "ant_cpw_to_transition_offset_x": 20,
        "ant_cpw_to_transition_angle_diff": 0.0,
        #
        "ant_microstrip_lw": 5,
        "ant_microstrip_final_x_offset": 750,
        "ant_microstrip_final_y_offset": 130,
        #
        "frame_bot_lw": 8,
        # "frame_bot_left_width": 996,
        # "frame_bot_right_width": 996,
        "frame_bot_left_width": 975,
        "frame_bot_right_width": 975,
        "frame_left_lw": 8,
        "frame_left_height": 400,
        "frame_right_lw": 8,
        "frame_right_height": 400,
        "frame_meander_cover_box_width": 5,
        "frame_meander_cover_box_height": 28,
        #
        "coupler_frame_left_lw": 10,
        "coupler_frame_left_height": 39,
        "coupler_frame_top_lw": 3,
        "coupler_gap": 16,
        "coupler_lw": 3,
        "left_coupler_frame_to_feed_distance": 164,
        #
        "IDC_bot_arm_gap": 30,
        "IDC_arm_gap": 8,
        "IDC_arm_lw": 3,
        "No_of_arms": 28,
        #
        "trim_arm_offset_right_side": 380,
        "trim_arm_offset_left_side": 389,
        "trim_arm_lw": 3,
        "trim_arm_length_right_side": 1975,
        "trim_arm_length_left_side": 1975,
        "step_down_distance_between_layers": 5,
        #
        "cutout_bot_offset": 5,
        "cutout_left_offset": 50,
        "cutout_right_offset": 50,
        "cutout_top_offset": 25,
        #
        "grndpln_gap_between_adjacent_resonators": 14,
        #
        "SiO_stepdown_cutout_width": 110,
        "SiO_stepdown_cutout_height": 0,
        "SiN_membrane_stepdown_cutout_width": 100,
        "SiN_membrane_stepdown_cutout_height": 0,
        #
        "text_size": 90,
        "text_x_offset": 800,
        "text_y_offset": 200,
        "text_underline_height": 3,
    }

    return config
