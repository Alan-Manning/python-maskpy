from __future__ import annotations

from enum import StrEnum

# """
# To generate all the enums in this file.
# Uncomment this section and commment out line 5 importing from . module.
# """
#
# if __name__ == "__main__":
#     from maskpy.souk_mask_configs import SoukMaskConfig
#     from maskpy.souk_mask_configs import get_mask_default_config
#
#     def get_mask_default_config_keys(config_type: SoukMaskConfig) -> list[str]:
#         """."""
#         default_config = get_mask_default_config(config_type)
#         keys = list(default_config.keys())
#         return keys
#
#     def create_config_keys_for(config_type: SoukMaskConfig) -> None:
#         """Prints the config keys as a StrEnum for the config_type passed
#         in."""
#
#         config_type_name = config_type.name.lower()
#
#         default_config = get_mask_default_config(config_type)
#
#         print(f"class {config_type_name}(StrEnum):")
#         # docstring start
#         print(f'    """Config keys for {config_type_name}.')
#         print("")
#         print("    StrEnum Members")
#         print("    ---------------")
#         for key in default_config.keys():
#             print(f"    {key}")
#         print('    """')
#         print("")
#         # docstring end
#         for key, val in default_config.items():
#             print(f'    {key} = "{key}"')
#             print(f'    """Default: {type(val).__name__} = {val}"""')
#         print("")
#         print("")
#
#     for member in SoukMaskConfig:
#         print("################################################################################")
#         create_config_keys_for(member)
#
#     print("################################################################################")


################################################################################
################################################################################
class antenna(StrEnum):
    """Config keys for antenna.

    StrEnum Members
    ---------------
    distance_from_center
    base_width
    top_conect_width
    straight_height
    taper_height
    backside_check_circle_radius
    antenna_rotation
    sin_dep_cutout_circle_radius
    grnd_cutout_circle_radius
    """

    distance_from_center = "distance_from_center"
    """Default: int = 320"""
    base_width = "base_width"
    """Default: int = 390"""
    top_conect_width = "top_conect_width"
    """Default: int = 2"""
    straight_height = "straight_height"
    """Default: int = 600"""
    taper_height = "taper_height"
    """Default: int = 280"""
    backside_check_circle_radius = "backside_check_circle_radius"
    """Default: int = 1900"""
    antenna_rotation = "antenna_rotation"
    """Default: float = -22.5"""
    sin_dep_cutout_circle_radius = "sin_dep_cutout_circle_radius"
    """Default: int = 1900"""
    grnd_cutout_circle_radius = "grnd_cutout_circle_radius"
    """Default: int = 1205"""


################################################################################
class antenna_cpw_microstrip_trans(StrEnum):
    """Config keys for antenna_cpw_microstrip_trans.

    StrEnum Members
    ---------------
    CPW_tans_len1
    CPW_tans_len2
    CPW_tans_len3
    CPW_tans_len4
    CPW_tans_len5
    gap1
    gap2
    gap3
    gap4
    gap5
    microstrip_lw
    CPW_width
    extend_out_distance
    flex_path_bend_radius
    """

    CPW_tans_len1 = "CPW_tans_len1"
    """Default: int = 108"""
    CPW_tans_len2 = "CPW_tans_len2"
    """Default: int = 70"""
    CPW_tans_len3 = "CPW_tans_len3"
    """Default: int = 36"""
    CPW_tans_len4 = "CPW_tans_len4"
    """Default: int = 17"""
    CPW_tans_len5 = "CPW_tans_len5"
    """Default: int = 4"""
    gap1 = "gap1"
    """Default: int = 4"""
    gap2 = "gap2"
    """Default: int = 16"""
    gap3 = "gap3"
    """Default: int = 39"""
    gap4 = "gap4"
    """Default: int = 76"""
    gap5 = "gap5"
    """Default: int = 121"""
    microstrip_lw = "microstrip_lw"
    """Default: int = 5"""
    CPW_width = "CPW_width"
    """Default: int = 11"""
    extend_out_distance = "extend_out_distance"
    """Default: int = 100"""
    flex_path_bend_radius = "flex_path_bend_radius"
    """Default: int = 80"""


################################################################################
class bottom_choke(StrEnum):
    """Config keys for bottom_choke.

    StrEnum Members
    ---------------
    wave_guide_hole_radius
    pad_radius
    pad_x_offset_from_center
    pad_y_offset_from_center
    IDC_cutout_offset_top
    IDC_cutout_offset_bot
    IDC_cutout_offset_left
    IDC_cutout_offset_right
    backshort_circle_hole_radius
    backshort_IDC_offset_top
    backshort_IDC_offset_bot
    backshort_IDC_offset_left
    backshort_IDC_offset_right
    """

    wave_guide_hole_radius = "wave_guide_hole_radius"
    """Default: int = 1200"""
    pad_radius = "pad_radius"
    """Default: int = 150"""
    pad_x_offset_from_center = "pad_x_offset_from_center"
    """Default: int = 1150"""
    pad_y_offset_from_center = "pad_y_offset_from_center"
    """Default: int = 2030"""
    IDC_cutout_offset_top = "IDC_cutout_offset_top"
    """Default: int = 100"""
    IDC_cutout_offset_bot = "IDC_cutout_offset_bot"
    """Default: int = 100"""
    IDC_cutout_offset_left = "IDC_cutout_offset_left"
    """Default: int = 100"""
    IDC_cutout_offset_right = "IDC_cutout_offset_right"
    """Default: int = 100"""
    backshort_circle_hole_radius = "backshort_circle_hole_radius"
    """Default: int = 1200"""
    backshort_IDC_offset_top = "backshort_IDC_offset_top"
    """Default: int = 80"""
    backshort_IDC_offset_bot = "backshort_IDC_offset_bot"
    """Default: int = 80"""
    backshort_IDC_offset_left = "backshort_IDC_offset_left"
    """Default: int = 210"""
    backshort_IDC_offset_right = "backshort_IDC_offset_right"
    """Default: int = 210"""


################################################################################
class combiner_section_150ghz(StrEnum):
    """Config keys for combiner_section_150ghz.

    StrEnum Members
    ---------------
    first_linewidth
    distance_from_outer_ring
    first_section_height
    first_section_width
    first_section_back_height
    second_section_linewidth
    second_section_top_width1
    second_section_top_height1
    second_section_top_width2
    second_section_top_height2
    second_section_top_width3
    second_section_top_height3
    second_section_top_width4
    second_section_top_height4
    second_section_top_width5
    second_section_bot_width1
    second_section_bot_height1
    second_section_bot_width2
    second_section_bot_height2
    second_section_bot_width3
    second_section_bot_height3
    second_section_bot_width4
    second_section_bot_height4
    second_section_bot_width5
    start_of_second_section_vertical_linewidth
    end_of_second_section_vertical_linewidth
    third_section_linewidth
    third_section_width1
    third_section_height1
    third_section_width2
    third_section_height2
    third_section_width3
    meander_height
    meander_linewidth
    meander_initial_gap
    no_of_full_straights
    meander_final_gap
    meander_gap_spacing
    meander_to_frame_box_size
    meander_last_fork_top_box_size
    meander_last_fork_via_box_size_x
    meander_last_fork_via_box_size_y
    meander_last_fork_via_box_y_offset
    meander_conect_linewidth
    meander_conect_bot_width
    meander_conect_right_height
    meander_conect_final_width
    meander_last_fork_wdith
    meander_last_fork_height
    meander_last_fork_linewdith
    meander_last_fork_start_height
    conection_to_kid_path_linewidth
    conection_to_kid_path_start_piece_length
    """

    first_linewidth = "first_linewidth"
    """Default: int = 5"""
    distance_from_outer_ring = "distance_from_outer_ring"
    """Default: int = 50"""
    first_section_height = "first_section_height"
    """Default: int = 560"""
    first_section_width = "first_section_width"
    """Default: int = 155"""
    first_section_back_height = "first_section_back_height"
    """Default: float = 150.5"""
    second_section_linewidth = "second_section_linewidth"
    """Default: int = 4"""
    second_section_top_width1 = "second_section_top_width1"
    """Default: float = 24.5"""
    second_section_top_height1 = "second_section_top_height1"
    """Default: int = 29"""
    second_section_top_width2 = "second_section_top_width2"
    """Default: int = 29"""
    second_section_top_height2 = "second_section_top_height2"
    """Default: int = 29"""
    second_section_top_width3 = "second_section_top_width3"
    """Default: int = 27"""
    second_section_top_height3 = "second_section_top_height3"
    """Default: int = 0"""
    second_section_top_width4 = "second_section_top_width4"
    """Default: int = 0"""
    second_section_top_height4 = "second_section_top_height4"
    """Default: int = 0"""
    second_section_top_width5 = "second_section_top_width5"
    """Default: int = 0"""
    second_section_bot_width1 = "second_section_bot_width1"
    """Default: float = 24.5"""
    second_section_bot_height1 = "second_section_bot_height1"
    """Default: int = 29"""
    second_section_bot_width2 = "second_section_bot_width2"
    """Default: int = 29"""
    second_section_bot_height2 = "second_section_bot_height2"
    """Default: int = 29"""
    second_section_bot_width3 = "second_section_bot_width3"
    """Default: int = 27"""
    second_section_bot_height3 = "second_section_bot_height3"
    """Default: int = 0"""
    second_section_bot_width4 = "second_section_bot_width4"
    """Default: int = 0"""
    second_section_bot_height4 = "second_section_bot_height4"
    """Default: int = 0"""
    second_section_bot_width5 = "second_section_bot_width5"
    """Default: int = 0"""
    start_of_second_section_vertical_linewidth = "start_of_second_section_vertical_linewidth"
    """Default: float = 3.5"""
    end_of_second_section_vertical_linewidth = "end_of_second_section_vertical_linewidth"
    """Default: int = 6"""
    third_section_linewidth = "third_section_linewidth"
    """Default: int = 3"""
    third_section_width1 = "third_section_width1"
    """Default: float = 26.5"""
    third_section_height1 = "third_section_height1"
    """Default: int = 28"""
    third_section_width2 = "third_section_width2"
    """Default: int = 28"""
    third_section_height2 = "third_section_height2"
    """Default: int = 28"""
    third_section_width3 = "third_section_width3"
    """Default: int = 28"""
    meander_height = "meander_height"
    """Default: int = 525"""
    meander_linewidth = "meander_linewidth"
    """Default: int = 5"""
    meander_initial_gap = "meander_initial_gap"
    """Default: int = 10"""
    no_of_full_straights = "no_of_full_straights"
    """Default: int = 10"""
    meander_final_gap = "meander_final_gap"
    """Default: int = 20"""
    meander_gap_spacing = "meander_gap_spacing"
    """Default: int = 10"""
    meander_to_frame_box_size = "meander_to_frame_box_size"
    """Default: int = 10"""
    meander_last_fork_top_box_size = "meander_last_fork_top_box_size"
    """Default: int = 11"""
    meander_last_fork_via_box_size_x = "meander_last_fork_via_box_size_x"
    """Default: int = 13"""
    meander_last_fork_via_box_size_y = "meander_last_fork_via_box_size_y"
    """Default: int = 11"""
    meander_last_fork_via_box_y_offset = "meander_last_fork_via_box_y_offset"
    """Default: int = 1"""
    meander_conect_linewidth = "meander_conect_linewidth"
    """Default: int = 5"""
    meander_conect_bot_width = "meander_conect_bot_width"
    """Default: int = 30"""
    meander_conect_right_height = "meander_conect_right_height"
    """Default: int = 265"""
    meander_conect_final_width = "meander_conect_final_width"
    """Default: float = 7.5"""
    meander_last_fork_wdith = "meander_last_fork_wdith"
    """Default: int = 5"""
    meander_last_fork_height = "meander_last_fork_height"
    """Default: int = 168"""
    meander_last_fork_linewdith = "meander_last_fork_linewdith"
    """Default: float = 2.5"""
    meander_last_fork_start_height = "meander_last_fork_start_height"
    """Default: int = 95"""
    conection_to_kid_path_linewidth = "conection_to_kid_path_linewidth"
    """Default: int = 5"""
    conection_to_kid_path_start_piece_length = "conection_to_kid_path_start_piece_length"
    """Default: int = 10"""


################################################################################
class combiner_section_90ghz(StrEnum):
    """Config keys for combiner_section_90ghz.

    StrEnum Members
    ---------------
    first_linewidth
    distance_from_outer_ring
    first_section_height
    first_section_width
    first_section_back_height
    second_section_linewidth
    second_section_top_width1
    second_section_top_height1
    second_section_top_width2
    second_section_top_height2
    second_section_top_width3
    second_section_top_height3
    second_section_top_width4
    second_section_top_height4
    second_section_top_width5
    second_section_bot_width1
    second_section_bot_height1
    second_section_bot_width2
    second_section_bot_height2
    second_section_bot_width3
    second_section_bot_height3
    second_section_bot_width4
    second_section_bot_height4
    second_section_bot_width5
    start_of_second_section_vertical_linewidth
    end_of_second_section_vertical_linewidth
    third_section_linewidth
    third_section_width1
    third_section_height1
    third_section_width2
    third_section_height2
    third_section_width3
    meander_height
    meander_linewidth
    meander_initial_gap
    no_of_full_straights
    meander_final_gap
    meander_gap_spacing
    meander_to_frame_box_size
    meander_last_fork_top_box_size
    meander_last_fork_via_box_size_x
    meander_last_fork_via_box_size_y
    meander_last_fork_via_box_y_offset
    meander_conect_linewidth
    meander_conect_bot_width
    meander_conect_right_height
    meander_conect_final_width
    meander_last_fork_wdith
    meander_last_fork_height
    meander_last_fork_linewdith
    meander_last_fork_start_height
    conection_to_kid_path_linewidth
    conection_to_kid_path_start_piece_length
    """

    first_linewidth = "first_linewidth"
    """Default: int = 5"""
    distance_from_outer_ring = "distance_from_outer_ring"
    """Default: int = 50"""
    first_section_height = "first_section_height"
    """Default: int = 573"""
    first_section_width = "first_section_width"
    """Default: int = 155"""
    first_section_back_height = "first_section_back_height"
    """Default: int = 101"""
    second_section_linewidth = "second_section_linewidth"
    """Default: float = 4.5"""
    second_section_top_width1 = "second_section_top_width1"
    """Default: float = 22.25"""
    second_section_top_height1 = "second_section_top_height1"
    """Default: int = 61"""
    second_section_top_width2 = "second_section_top_width2"
    """Default: float = 44.5"""
    second_section_top_height2 = "second_section_top_height2"
    """Default: int = 61"""
    second_section_top_width3 = "second_section_top_width3"
    """Default: float = 22.25"""
    second_section_top_height3 = "second_section_top_height3"
    """Default: int = 0"""
    second_section_top_width4 = "second_section_top_width4"
    """Default: int = 0"""
    second_section_top_height4 = "second_section_top_height4"
    """Default: int = 0"""
    second_section_top_width5 = "second_section_top_width5"
    """Default: int = 0"""
    second_section_bot_width1 = "second_section_bot_width1"
    """Default: float = 22.25"""
    second_section_bot_height1 = "second_section_bot_height1"
    """Default: int = 61"""
    second_section_bot_width2 = "second_section_bot_width2"
    """Default: float = 44.5"""
    second_section_bot_height2 = "second_section_bot_height2"
    """Default: int = 61"""
    second_section_bot_width3 = "second_section_bot_width3"
    """Default: float = 22.25"""
    second_section_bot_height3 = "second_section_bot_height3"
    """Default: int = 0"""
    second_section_bot_width4 = "second_section_bot_width4"
    """Default: int = 0"""
    second_section_bot_height4 = "second_section_bot_height4"
    """Default: int = 0"""
    second_section_bot_width5 = "second_section_bot_width5"
    """Default: int = 0"""
    start_of_second_section_vertical_linewidth = "start_of_second_section_vertical_linewidth"
    """Default: float = 3.5"""
    end_of_second_section_vertical_linewidth = "end_of_second_section_vertical_linewidth"
    """Default: int = 5"""
    third_section_linewidth = "third_section_linewidth"
    """Default: int = 3"""
    third_section_width1 = "third_section_width1"
    """Default: float = 21.5"""
    third_section_height1 = "third_section_height1"
    """Default: int = 50"""
    third_section_width2 = "third_section_width2"
    """Default: int = 43"""
    third_section_height2 = "third_section_height2"
    """Default: int = 50"""
    third_section_width3 = "third_section_width3"
    """Default: int = 23"""
    meander_height = "meander_height"
    """Default: int = 525"""
    meander_linewidth = "meander_linewidth"
    """Default: int = 5"""
    meander_initial_gap = "meander_initial_gap"
    """Default: int = 20"""
    no_of_full_straights = "no_of_full_straights"
    """Default: int = 10"""
    meander_final_gap = "meander_final_gap"
    """Default: int = 15"""
    meander_gap_spacing = "meander_gap_spacing"
    """Default: int = 10"""
    meander_to_frame_box_size = "meander_to_frame_box_size"
    """Default: int = 10"""
    meander_last_fork_top_box_size = "meander_last_fork_top_box_size"
    """Default: int = 11"""
    meander_last_fork_via_box_size_x = "meander_last_fork_via_box_size_x"
    """Default: int = 13"""
    meander_last_fork_via_box_size_y = "meander_last_fork_via_box_size_y"
    """Default: int = 11"""
    meander_last_fork_via_box_y_offset = "meander_last_fork_via_box_y_offset"
    """Default: int = 1"""
    meander_conect_linewidth = "meander_conect_linewidth"
    """Default: int = 5"""
    meander_conect_bot_width = "meander_conect_bot_width"
    """Default: int = 25"""
    meander_conect_right_height = "meander_conect_right_height"
    """Default: int = 260"""
    meander_conect_final_width = "meander_conect_final_width"
    """Default: float = 7.5"""
    meander_last_fork_wdith = "meander_last_fork_wdith"
    """Default: int = 3"""
    meander_last_fork_height = "meander_last_fork_height"
    """Default: int = 250"""
    meander_last_fork_linewdith = "meander_last_fork_linewdith"
    """Default: float = 2.5"""
    meander_last_fork_start_height = "meander_last_fork_start_height"
    """Default: int = 14"""
    conection_to_kid_path_linewidth = "conection_to_kid_path_linewidth"
    """Default: int = 5"""
    conection_to_kid_path_start_piece_length = "conection_to_kid_path_start_piece_length"
    """Default: int = 10"""


################################################################################
class cpw_feedline(StrEnum):
    """Config keys for cpw_feedline.

    StrEnum Members
    ---------------
    feedline_width
    cutout_around_feedline_width
    dielectric_under_feedline_width
    bend_radius
    extra_straight_length
    bridge_gap
    bridge_width
    """

    feedline_width = "feedline_width"
    """Default: int = 36"""
    cutout_around_feedline_width = "cutout_around_feedline_width"
    """Default: int = 68"""
    dielectric_under_feedline_width = "dielectric_under_feedline_width"
    """Default: int = 188"""
    bend_radius = "bend_radius"
    """Default: int = 200"""
    extra_straight_length = "extra_straight_length"
    """Default: int = 200"""
    bridge_gap = "bridge_gap"
    """Default: int = 2750"""
    bridge_width = "bridge_width"
    """Default: int = 6"""


################################################################################
class filter_bank(StrEnum):
    """Config keys for filter_bank.

    StrEnum Members
    ---------------
    inner_ring_radius
    inner_ring_line_width
    outer_ring_radius
    outer_ring_line_width
    ring_overlap_distance_from_center
    inner_ring_overlap_gap
    outer_ring_overlap_gap
    outer_ring_conector_gap
    inner_straight_extend_distance
    inner_bend_radius
    outer_straight_extend_distance
    outer_bend_radius
    Hugh_filter_thin_length1
    Hugh_filter_thin_length2
    Hugh_filter_thick_length1
    Hugh_filter_thick_length2
    Hugh_filter_thin_width
    Hugh_filter_thick_width
    ring_overlap_0_rot
    ring_overlap_1_rot
    ring_overlap_2_rot
    ring_overlap_3_rot
    """

    inner_ring_radius = "inner_ring_radius"
    """Default: int = 1960"""
    inner_ring_line_width = "inner_ring_line_width"
    """Default: int = 5"""
    outer_ring_radius = "outer_ring_radius"
    """Default: float = 2082.5"""
    outer_ring_line_width = "outer_ring_line_width"
    """Default: int = 5"""
    ring_overlap_distance_from_center = "ring_overlap_distance_from_center"
    """Default: int = 2020"""
    inner_ring_overlap_gap = "inner_ring_overlap_gap"
    """Default: int = 325"""
    outer_ring_overlap_gap = "outer_ring_overlap_gap"
    """Default: int = 325"""
    outer_ring_conector_gap = "outer_ring_conector_gap"
    """Default: int = 20"""
    inner_straight_extend_distance = "inner_straight_extend_distance"
    """Default: int = 30"""
    inner_bend_radius = "inner_bend_radius"
    """Default: int = 25"""
    outer_straight_extend_distance = "outer_straight_extend_distance"
    """Default: int = 30"""
    outer_bend_radius = "outer_bend_radius"
    """Default: int = 25"""
    Hugh_filter_thin_length1 = "Hugh_filter_thin_length1"
    """Default: int = 40"""
    Hugh_filter_thin_length2 = "Hugh_filter_thin_length2"
    """Default: int = 85"""
    Hugh_filter_thick_length1 = "Hugh_filter_thick_length1"
    """Default: int = 38"""
    Hugh_filter_thick_length2 = "Hugh_filter_thick_length2"
    """Default: int = 46"""
    Hugh_filter_thin_width = "Hugh_filter_thin_width"
    """Default: int = 2"""
    Hugh_filter_thick_width = "Hugh_filter_thick_width"
    """Default: int = 20"""
    ring_overlap_0_rot = "ring_overlap_0_rot"
    """Default: int = -45"""
    ring_overlap_1_rot = "ring_overlap_1_rot"
    """Default: int = 45"""
    ring_overlap_2_rot = "ring_overlap_2_rot"
    """Default: int = 45"""
    ring_overlap_3_rot = "ring_overlap_3_rot"
    """Default: int = -45"""


################################################################################
class filter_bank_ring_overlap(StrEnum):
    """Config keys for filter_bank_ring_overlap.

    StrEnum Members
    ---------------
    outer_box_width
    outer_box_height
    outer_box_inner_cutout_width
    outer_box_inner_cutout_height
    linewidth
    half_conect_taper_extra_height
    half_conect_taper_stright_width
    half_conect_taper_diag_width
    half_conect_taper_inner_height
    half_conect_distacne_from_center
    half_conect_bridge_rect_width
    half_conect_bridge_rect_height
    half_conect_bridge_pad_width
    half_conect_bridge_pad_height
    half_conect_bridge_pad_offset_from_center
    full_conect_taper_extra_width
    full_conect_taper_stright_height
    full_conect_taper_diag_height
    full_conect_taper_start_from_center
    full_conect_center_width
    """

    outer_box_width = "outer_box_width"
    """Default: int = 100"""
    outer_box_height = "outer_box_height"
    """Default: int = 100"""
    outer_box_inner_cutout_width = "outer_box_inner_cutout_width"
    """Default: int = 27"""
    outer_box_inner_cutout_height = "outer_box_inner_cutout_height"
    """Default: int = 12"""
    linewidth = "linewidth"
    """Default: int = 5"""
    half_conect_taper_extra_height = "half_conect_taper_extra_height"
    """Default: int = 5"""
    half_conect_taper_stright_width = "half_conect_taper_stright_width"
    """Default: int = 5"""
    half_conect_taper_diag_width = "half_conect_taper_diag_width"
    """Default: float = 5.5"""
    half_conect_taper_inner_height = "half_conect_taper_inner_height"
    """Default: int = 4"""
    half_conect_distacne_from_center = "half_conect_distacne_from_center"
    """Default: int = 5"""
    half_conect_bridge_rect_width = "half_conect_bridge_rect_width"
    """Default: int = 21"""
    half_conect_bridge_rect_height = "half_conect_bridge_rect_height"
    """Default: int = 4"""
    half_conect_bridge_pad_width = "half_conect_bridge_pad_width"
    """Default: float = 5.5"""
    half_conect_bridge_pad_height = "half_conect_bridge_pad_height"
    """Default: int = 4"""
    half_conect_bridge_pad_offset_from_center = "half_conect_bridge_pad_offset_from_center"
    """Default: int = 5"""
    full_conect_taper_extra_width = "full_conect_taper_extra_width"
    """Default: int = 5"""
    full_conect_taper_stright_height = "full_conect_taper_stright_height"
    """Default: int = 5"""
    full_conect_taper_diag_height = "full_conect_taper_diag_height"
    """Default: float = 5.5"""
    full_conect_taper_start_from_center = "full_conect_taper_start_from_center"
    """Default: int = 4"""
    full_conect_center_width = "full_conect_center_width"
    """Default: int = 4"""


################################################################################
class general(StrEnum):
    """Config keys for general.

    StrEnum Members
    ---------------
    horizontal_pitch
    vertical_pitch
    chip_center_x
    chip_center_y
    """

    horizontal_pitch = "horizontal_pitch"
    """Default: int = 5306"""
    vertical_pitch = "vertical_pitch"
    """Default: int = 4590"""
    chip_center_x = "chip_center_x"
    """Default: int = 0"""
    chip_center_y = "chip_center_y"
    """Default: int = 0"""


################################################################################
class hi_pass_filters(StrEnum):
    """Config keys for hi_pass_filters.

    StrEnum Members
    ---------------
    Hi_pass_arm1_linewidth
    Hi_pass_arm1_height
    Hi_pass_arm1_length
    Hi_pass_arm1_offset
    Hi_pass_arm2_linewidth
    Hi_pass_arm2_height
    Hi_pass_arm2_length
    Hi_pass_arm2_offset
    Hi_pass_arm3_linewidth
    Hi_pass_arm3_height
    Hi_pass_arm3_length
    Hi_pass_arm3_offset
    Hi_pass_arm4_linewidth
    Hi_pass_arm4_height
    Hi_pass_arm4_length
    Hi_pass_arm4_offset
    Hi_pass_arm5_linewidth
    Hi_pass_arm5_height
    Hi_pass_arm5_length
    Hi_pass_arm5_offset
    via_overlap_length
    via_extend_length
    """

    Hi_pass_arm1_linewidth = "Hi_pass_arm1_linewidth"
    """Default: int = 22"""
    Hi_pass_arm1_height = "Hi_pass_arm1_height"
    """Default: int = 23"""
    Hi_pass_arm1_length = "Hi_pass_arm1_length"
    """Default: int = 97"""
    Hi_pass_arm1_offset = "Hi_pass_arm1_offset"
    """Default: int = 141"""
    Hi_pass_arm2_linewidth = "Hi_pass_arm2_linewidth"
    """Default: int = 42"""
    Hi_pass_arm2_height = "Hi_pass_arm2_height"
    """Default: int = 34"""
    Hi_pass_arm2_length = "Hi_pass_arm2_length"
    """Default: int = 77"""
    Hi_pass_arm2_offset = "Hi_pass_arm2_offset"
    """Default: int = 134"""
    Hi_pass_arm3_linewidth = "Hi_pass_arm3_linewidth"
    """Default: int = 42"""
    Hi_pass_arm3_height = "Hi_pass_arm3_height"
    """Default: int = 34"""
    Hi_pass_arm3_length = "Hi_pass_arm3_length"
    """Default: int = 77"""
    Hi_pass_arm3_offset = "Hi_pass_arm3_offset"
    """Default: int = 124"""
    Hi_pass_arm4_linewidth = "Hi_pass_arm4_linewidth"
    """Default: int = 42"""
    Hi_pass_arm4_height = "Hi_pass_arm4_height"
    """Default: int = 34"""
    Hi_pass_arm4_length = "Hi_pass_arm4_length"
    """Default: int = 77"""
    Hi_pass_arm4_offset = "Hi_pass_arm4_offset"
    """Default: int = 124"""
    Hi_pass_arm5_linewidth = "Hi_pass_arm5_linewidth"
    """Default: int = 22"""
    Hi_pass_arm5_height = "Hi_pass_arm5_height"
    """Default: int = 23"""
    Hi_pass_arm5_length = "Hi_pass_arm5_length"
    """Default: int = 97"""
    Hi_pass_arm5_offset = "Hi_pass_arm5_offset"
    """Default: int = 114"""
    via_overlap_length = "via_overlap_length"
    """Default: int = 5"""
    via_extend_length = "via_extend_length"
    """Default: int = 5"""


################################################################################
class lo_pass_filters(StrEnum):
    """Config keys for lo_pass_filters.

    StrEnum Members
    ---------------
    Lo_pass_arm1_linewidth
    Lo_pass_arm1_height
    Lo_pass_arm1_length
    Lo_pass_arm1_offset
    Lo_pass_arm2_linewidth
    Lo_pass_arm2_height
    Lo_pass_arm2_length
    Lo_pass_arm2_offset
    Lo_pass_arm3_linewidth
    Lo_pass_arm3_height
    Lo_pass_arm3_length
    Lo_pass_arm3_offset
    Lo_pass_arm4_linewidth
    Lo_pass_arm4_height
    Lo_pass_arm4_length
    Lo_pass_arm4_offset
    Lo_pass_arm5_linewidth
    Lo_pass_arm5_height
    Lo_pass_arm5_length
    Lo_pass_arm5_offset
    via_overlap_length
    via_extend_length
    """

    Lo_pass_arm1_linewidth = "Lo_pass_arm1_linewidth"
    """Default: int = 10"""
    Lo_pass_arm1_height = "Lo_pass_arm1_height"
    """Default: int = 20"""
    Lo_pass_arm1_length = "Lo_pass_arm1_length"
    """Default: int = 175"""
    Lo_pass_arm1_offset = "Lo_pass_arm1_offset"
    """Default: int = 139"""
    Lo_pass_arm2_linewidth = "Lo_pass_arm2_linewidth"
    """Default: int = 20"""
    Lo_pass_arm2_height = "Lo_pass_arm2_height"
    """Default: int = 41"""
    Lo_pass_arm2_length = "Lo_pass_arm2_length"
    """Default: int = 160"""
    Lo_pass_arm2_offset = "Lo_pass_arm2_offset"
    """Default: int = 255"""
    Lo_pass_arm3_linewidth = "Lo_pass_arm3_linewidth"
    """Default: int = 20"""
    Lo_pass_arm3_height = "Lo_pass_arm3_height"
    """Default: int = 41"""
    Lo_pass_arm3_length = "Lo_pass_arm3_length"
    """Default: int = 160"""
    Lo_pass_arm3_offset = "Lo_pass_arm3_offset"
    """Default: int = 215"""
    Lo_pass_arm4_linewidth = "Lo_pass_arm4_linewidth"
    """Default: int = 20"""
    Lo_pass_arm4_height = "Lo_pass_arm4_height"
    """Default: int = 41"""
    Lo_pass_arm4_length = "Lo_pass_arm4_length"
    """Default: int = 160"""
    Lo_pass_arm4_offset = "Lo_pass_arm4_offset"
    """Default: int = 215"""
    Lo_pass_arm5_linewidth = "Lo_pass_arm5_linewidth"
    """Default: int = 10"""
    Lo_pass_arm5_height = "Lo_pass_arm5_height"
    """Default: int = 41"""
    Lo_pass_arm5_length = "Lo_pass_arm5_length"
    """Default: int = 170"""
    Lo_pass_arm5_offset = "Lo_pass_arm5_offset"
    """Default: int = 210"""
    via_overlap_length = "via_overlap_length"
    """Default: int = 5"""
    via_extend_length = "via_extend_length"
    """Default: int = 5"""


################################################################################
class port(StrEnum):
    """Config keys for port.

    StrEnum Members
    ---------------
    outer_feedline_width
    outer_cutout_around_feedline_width
    outer_dielectric_under_feedline_width
    outer_back_length
    taper_length
    dielectric_cutout_in_port_width
    dielectric_cutout_in_port_length
    """

    outer_feedline_width = "outer_feedline_width"
    """Default: int = 500"""
    outer_cutout_around_feedline_width = "outer_cutout_around_feedline_width"
    """Default: int = 1620"""
    outer_dielectric_under_feedline_width = "outer_dielectric_under_feedline_width"
    """Default: int = 1740"""
    outer_back_length = "outer_back_length"
    """Default: int = 1000"""
    taper_length = "taper_length"
    """Default: int = 1500"""
    dielectric_cutout_in_port_width = "dielectric_cutout_in_port_width"
    """Default: int = 500"""
    dielectric_cutout_in_port_length = "dielectric_cutout_in_port_length"
    """Default: int = 1000"""


################################################################################
class sma_connector(StrEnum):
    """Config keys for sma_connector.

    StrEnum Members
    ---------------
    sma_square_offset_left
    sma_square_offset_right
    sma_square_offset_top
    sma_square_offset_bot
    sma_circle_radius
    sma_square_width
    sma_square_height
    taper_length
    central_linewidth
    gap
    dielectric_overlap_distance
    bend_offset_from_end
    extra_bend_center
    """

    sma_square_offset_left = "sma_square_offset_left"
    """Default: int = 0"""
    sma_square_offset_right = "sma_square_offset_right"
    """Default: int = 7000"""
    sma_square_offset_top = "sma_square_offset_top"
    """Default: int = 0"""
    sma_square_offset_bot = "sma_square_offset_bot"
    """Default: int = 0"""
    sma_circle_radius = "sma_circle_radius"
    """Default: int = 2160"""
    sma_square_width = "sma_square_width"
    """Default: int = 5000"""
    sma_square_height = "sma_square_height"
    """Default: int = 5000"""
    taper_length = "taper_length"
    """Default: int = 3000"""
    central_linewidth = "central_linewidth"
    """Default: int = 500"""
    gap = "gap"
    """Default: int = 500"""
    dielectric_overlap_distance = "dielectric_overlap_distance"
    """Default: int = 300"""
    bend_offset_from_end = "bend_offset_from_end"
    """Default: int = 200"""
    extra_bend_center = "extra_bend_center"
    """Default: int = 50"""


################################################################################
class test_chip_quad(StrEnum):
    """Config keys for test_chip_quad.

    StrEnum Members
    ---------------
    test_chip_quad_width
    test_chip_quad_height
    chip_to_holder_bend_rad
    chip_to_holder_feedline_width
    chip_to_holder_cutout_around_feedline_width
    chip_to_holder_dielectric_under_feedline_width
    chip_to_holder_port_extend
    chip_to_holder_holder_extend
    chip_to_holder_end_cutout_distance
    chip_to_holder_end_cutout_center_distance
    horn_offset_from_center_x
    horn_offset_from_center_y
    top_right_label_window_offset_x
    top_right_label_window_offset_y
    top_right_label_window_width
    top_right_label_window_height
    bottom_left_text_offset_x
    bottom_left_text_offset_y
    bottom_left_text_size
    top_right_text_offset_x
    top_right_text_offset_y
    top_right_text_size
    top_left_text_offset_x
    top_left_text_offset_y
    top_left_text_size
    cardiff_logo_size
    cardiff_logo_offset_x
    cardiff_logo_offset_y
    souk_logo_size
    souk_logo_offset_x
    souk_logo_offset_y
    groundplane_edge_cuts_top_mid_gap
    groundplane_edge_cuts_bot_mid_gap
    groundplane_edge_cuts_right_mid_gap
    groundplane_edge_cuts_left_mid_gap
    center_pin_radius
    slotted_pin_radius
    slotted_pin_length
    slotted_pin_offset_x_from_center_pin
    slotted_pin_offset_y_from_center_pin
    """

    test_chip_quad_width = "test_chip_quad_width"
    """Default: int = 27000"""
    test_chip_quad_height = "test_chip_quad_height"
    """Default: int = 27000"""
    chip_to_holder_bend_rad = "chip_to_holder_bend_rad"
    """Default: int = 2500"""
    chip_to_holder_feedline_width = "chip_to_holder_feedline_width"
    """Default: int = 500"""
    chip_to_holder_cutout_around_feedline_width = "chip_to_holder_cutout_around_feedline_width"
    """Default: int = 1620"""
    chip_to_holder_dielectric_under_feedline_width = "chip_to_holder_dielectric_under_feedline_width"
    """Default: int = 1740"""
    chip_to_holder_port_extend = "chip_to_holder_port_extend"
    """Default: int = 4000"""
    chip_to_holder_holder_extend = "chip_to_holder_holder_extend"
    """Default: int = 10000"""
    chip_to_holder_end_cutout_distance = "chip_to_holder_end_cutout_distance"
    """Default: int = 4000"""
    chip_to_holder_end_cutout_center_distance = "chip_to_holder_end_cutout_center_distance"
    """Default: int = 700"""
    horn_offset_from_center_x = "horn_offset_from_center_x"
    """Default: int = 6000"""
    horn_offset_from_center_y = "horn_offset_from_center_y"
    """Default: int = 6000"""
    top_right_label_window_offset_x = "top_right_label_window_offset_x"
    """Default: int = 2000"""
    top_right_label_window_offset_y = "top_right_label_window_offset_y"
    """Default: int = 1200"""
    top_right_label_window_width = "top_right_label_window_width"
    """Default: int = 4000"""
    top_right_label_window_height = "top_right_label_window_height"
    """Default: int = 1000"""
    bottom_left_text_offset_x = "bottom_left_text_offset_x"
    """Default: int = 1000"""
    bottom_left_text_offset_y = "bottom_left_text_offset_y"
    """Default: int = 1000"""
    bottom_left_text_size = "bottom_left_text_size"
    """Default: int = 300"""
    top_right_text_offset_x = "top_right_text_offset_x"
    """Default: int = 1000"""
    top_right_text_offset_y = "top_right_text_offset_y"
    """Default: int = 1000"""
    top_right_text_size = "top_right_text_size"
    """Default: int = 300"""
    top_left_text_offset_x = "top_left_text_offset_x"
    """Default: int = 1000"""
    top_left_text_offset_y = "top_left_text_offset_y"
    """Default: int = 1000"""
    top_left_text_size = "top_left_text_size"
    """Default: int = 300"""
    cardiff_logo_size = "cardiff_logo_size"
    """Default: int = 2000"""
    cardiff_logo_offset_x = "cardiff_logo_offset_x"
    """Default: int = 1000"""
    cardiff_logo_offset_y = "cardiff_logo_offset_y"
    """Default: int = 1000"""
    souk_logo_size = "souk_logo_size"
    """Default: int = 2600"""
    souk_logo_offset_x = "souk_logo_offset_x"
    """Default: int = 3500"""
    souk_logo_offset_y = "souk_logo_offset_y"
    """Default: int = 700"""
    groundplane_edge_cuts_top_mid_gap = "groundplane_edge_cuts_top_mid_gap"
    """Default: int = 5500"""
    groundplane_edge_cuts_bot_mid_gap = "groundplane_edge_cuts_bot_mid_gap"
    """Default: int = 5500"""
    groundplane_edge_cuts_right_mid_gap = "groundplane_edge_cuts_right_mid_gap"
    """Default: int = 0"""
    groundplane_edge_cuts_left_mid_gap = "groundplane_edge_cuts_left_mid_gap"
    """Default: int = 0"""
    center_pin_radius = "center_pin_radius"
    """Default: int = 1010"""
    slotted_pin_radius = "slotted_pin_radius"
    """Default: int = 1010"""
    slotted_pin_length = "slotted_pin_length"
    """Default: int = 500"""
    slotted_pin_offset_x_from_center_pin = "slotted_pin_offset_x_from_center_pin"
    """Default: int = 8000"""
    slotted_pin_offset_y_from_center_pin = "slotted_pin_offset_y_from_center_pin"
    """Default: int = 0"""


################################################################################
class top_choke(StrEnum):
    """Config keys for top_choke.

    StrEnum Members
    ---------------
    waveguide_hole_radius
    anulus_width
    """

    waveguide_hole_radius = "waveguide_hole_radius"
    """Default: int = 1200"""
    anulus_width = "anulus_width"
    """Default: int = 600"""


################################################################################
