from enum import StrEnum


class SoukResonatorType(StrEnum):
    """Type of a Resonator.

    Enum Members
    ------------
    ORIGINAL_Q10K
    ORIGINAL_Q20K
    ORIGINAL_Q50K

    ORIGINAL_LONG_TRUNK_Q10K
    ORIGINAL_LONG_TRUNK_Q20K
    ORIGINAL_LONG_TRUNK_Q50K

    HIGH_VOLUME_V1_Q20K
    HIGH_VOLUME_V1_Q50K

    HIGH_VOLUME_V1_LONG_TRUNK_Q20K
    HIGH_VOLUME_V1_LONG_TRUNK_Q50K

    HIGH_VOLUME_V2_Q20K
    HIGH_VOLUME_V2_Q50K

    HIGH_VOLUME_V2_LONG_TRUNK_Q20K
    HIGH_VOLUME_V2_LONG_TRUNK_Q50K
    """

    ORIGINAL_Q10K = "original_q10k"
    ORIGINAL_Q20K = "original_q20k"
    ORIGINAL_Q50K = "original_q50k"

    ORIGINAL_LONG_TRUNK_Q10K = "original_long_trunk_q10k"
    ORIGINAL_LONG_TRUNK_Q20K = "original_long_trunk_q20k"
    ORIGINAL_LONG_TRUNK_Q50K = "original_long_trunk_q50k"

    HIGH_VOLUME_V1_Q20K = "high_volume_v1_q20k"
    HIGH_VOLUME_V1_Q50K = "high_volume_v1_q50k"

    HIGH_VOLUME_V1_LONG_TRUNK_Q20K = "high_volume_v1_long_trunk_q20k"
    HIGH_VOLUME_V1_LONG_TRUNK_Q50K = "high_volume_v1_long_trunk_q50k"

    HIGH_VOLUME_V2_Q20K = "high_volume_v2_q20k"
    HIGH_VOLUME_V2_Q50K = "high_volume_v2_q50k"

    HIGH_VOLUME_V2_LONG_TRUNK_Q20K = "high_volume_v2_long_trunk_q20k"
    HIGH_VOLUME_V2_LONG_TRUNK_Q50K = "high_volume_v2_long_trunk_q50k"

    @classmethod
    def _missing_(cls, value):
        value = value.upper()
        for member in cls:
            if member == value:
                return member
        return None


"""
################################################################################
To Add a new SoukResonatorType:

    - The type needs to be registered here. e.g.
        >>> NEW_RES_V1_Q10k = "new_res_v1_q10k"
    This example name will be used as example in the below steps.

    - There needs to exist a muxing function for this new resonator type. This
    is by convenstion a file named in lowercase "mux_new_res_v1_q10k" in the
    'maskpy/souk_resonators/new_res_v1_q10k' directory. This could a function that
    points to an existing muxing function defined in another file. See existing
    mux files in exiting directorys for other resonator types for examples.

    - This muxing function needs to be imported and added to the
    get_mux_func_for_resonator_type() function in the souk_muxing.py file in
    the maskpy directory.

    - There needs to exist a base configuration dictionary returned by a
    get_resonator_config() function for this new resonator type.
    This is by convenstion a file named in lowercase "config_new_res_v1_q10k"
    in the 'maskpy/souk_resonators/new_res_v1_q10k' directory. See existing
    configs files in existing directorys for other resonator types for examples.

    - This new config getting function needs to be imported and added to the
    get_resonator_config() function in the resonator_base_utils.py file in the
    'maskpy/souk_resonators/utils' directory.

    - There needs to exist a 4 functins that get details about the new
    resonator made. These are:
        >>> get_total_height_of_resonator()
        >>> get_horizontal_coupler_end_to_meander_base_distance()
        >>> get_vertical_coupler_center_to_meander_base_distance()
        >>> get_width_and_height_of_IDC_cutout_section()
    These functions are by convenstion made in a file named in lowercase
    "utils_new_res_v1_q10k" in the 'maskpy/souk_resonators/new_res_v1_q10k'
    directory. See existing utils files in existing directorys for other
    resonator types for examples and detailed descriptions for each of
    the functions required in the respective docstrings.

    - These new functions need to be imported and added to their respective get
    functions in the file resonator_base_utils.py in the
    'maskpy/souk_resonators/utils' directory. See existing implimentations for
    details on how this should be done.

    - There should exist a function capable of drawing the resonator to the
    mask. This can be an existing function or a new function can be made. This
    function should be called in the add_4_resonators_around_horn() method in
    the SoukMaskBuilder class in the souk_mask_builder.py file in the maskpy
    directory. See existing calls to drawing methods, add_resonator_original(),
    add_resonator_high_volume_v1(), add_resonator_high_volume_v2() for examples
    for how to impliment a new function of hook into already existing functions.

    -- Congratulations.
"""
