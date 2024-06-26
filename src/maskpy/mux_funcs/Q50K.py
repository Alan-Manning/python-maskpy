import numpy as np


def get_resonator_config() -> dict:
    """Get the config for the Q50k resonator.

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
        "ant_pad_box_width": 5,
        "ant_pad_box_height": 10,
        "frame_bot_lw": 8,
        "frame_bot_left_width": 996,
        "frame_bot_right_width": 996,
        "frame_left_lw": 8,
        "frame_left_height": 400,
        "frame_right_lw": 8,
        "frame_right_height": 400,
        "frame_meander_cover_box_width": 5,
        "frame_meander_cover_box_height": 48,
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
        "cutout_bot_offset": 15,
        "cutout_left_offset": 50,
        "cutout_right_offset": 50,
        "cutout_top_offset": 25,
        "grndpl_meander_cutout_width": 80,
        "grndpl_meander_cutout_height": 10,
        "grndpl_coupler_cutout_width": 10,
        "grndpln_gap_between_adjacent_resonators": 44,
        "grndpl_coupler_cutout_height": 30,
        "SiO_stepdown_cutout_width": 110,
        "SiO_stepdown_cutout_height": 0,
        "SiN_membrane_stepdown_cutout_width": 100,
        "SiN_membrane_stepdown_cutout_height": 0,
    }

    return config


def get_IDCLs_and_CCL_from_f0(f0: float, rounding_precision: int | None = 2):
    """Get the IDC lengths and CC length for a resonator given a resonant
    frequency in Hz. These lengths are by default rounded to 2 decimal places.

    Parameters
    ----------
    f0 : float
        The desired frequency (**in Hz**) for the KID.

    KwArgs
    ------
    rounding_precision=2
        This should be an integer. This is the amount of decimal places to
        round the result of the IDC and CC lengths to. By default (recomended)
        this is 2 decimal places. For no rounding this value should be None.

    Returns
    -------
    IDC_array : list
        This is a list of length 28. Each element is a float representing the
        length of an IDC arm in um. The first element in this list is the
        length of arm 1, the second element is the length of the second arm
        and so on...

    CCL : float
        This is the length of the coupling capacitor in um.
    """
    arms_1_4_min_max_freqs = [1970120000.0, 2080770000.0]
    arms_5_8_min_max_freqs = [2080770000.0, 2239280000.0]
    arms_9_12_min_max_freqs = [2239280000.0, 2441040000.0]
    arms_13_16_min_max_freqs = [2441040000.0, 2705890000.0]
    arms_17_20_min_max_freqs = [2705890000.0, 3066980000.0]
    arms_21_24_min_max_freqs = [3066980000.0, 3571600000.0]
    arms_25_28_min_max_freqs = [3571600000.0, 4174530000.0]

    min_freq = arms_1_4_min_max_freqs[0]
    max_freq = arms_25_28_min_max_freqs[1]

    if f0 < min_freq or f0 > max_freq:
        print("error occured")
        raise (Exception(f"Desired frequency out of range. Should be within {min_freq} -> {max_freq}"))

    # Fits from data.
    arms_1_4_IDC_f0_polyfit = [4.40593480e-23, -2.67909667e-13, 5.34980654e-04, -3.49055119e05]
    arms_5_8_IDC_f0_polyfit = [1.87200482e-23, -1.21162839e-13, 2.55765811e-04, -1.74274315e05]
    arms_9_12_IDC_f0_polyfit = [3.12896341e-24, -2.11432074e-14, 4.31831099e-05, -2.38375201e04]
    arms_13_16_IDC_f0_polyfit = [9.30798289e-25, -6.46468470e-15, 1.14603758e-05, -1.01767002e03]
    arms_17_20_IDC_f0_polyfit = [1.14601699e-25, -4.38019706e-16, -2.76247107e-06, 1.03869882e04]
    arms_21_24_IDC_f0_polyfit = [-9.89518123e-26, 1.33651102e-15, -7.32943707e-06, 1.47377098e04]
    arms_25_28_IDC_f0_polyfit = [-2.78899336e-25, 3.28681434e-15, -1.43346369e-05, 2.39520423e04]

    arms_1_4_CC_f0_polyfit = [1.02325925e-22, -6.15094874e-13, 1.23020222e-03, -8.17665744e05]
    arms_5_8_CC_f0_polyfit = [-2.24898154e-23, 1.49655879e-13, -3.32533402e-04, 2.47412710e05]
    arms_9_12_CC_f0_polyfit = [-3.91010408e-24, 2.86521516e-14, -7.04177628e-05, 5.85968409e04]
    arms_13_16_CC_f0_polyfit = [-3.27444516e-24, 2.59330629e-14, -6.87761882e-05, 6.15470736e04]
    arms_17_20_CC_f0_polyfit = [-1.09047056e-24, 9.81438047e-15, -2.96571843e-05, 3.04441295e04]
    arms_21_24_CC_f0_polyfit = [-6.20264683e-25, 6.32823387e-15, -2.16963227e-05, 2.52534933e04]
    arms_25_28_CC_f0_polyfit = [-1.84259106e-25, 2.25955100e-15, -9.33504064e-06, 1.31424615e04]

    # Making fits into polynomial functions.
    arms_1_4_IDC_f0_polyfunc = np.poly1d(arms_1_4_IDC_f0_polyfit)
    arms_5_8_IDC_f0_polyfunc = np.poly1d(arms_5_8_IDC_f0_polyfit)
    arms_9_12_IDC_f0_polyfunc = np.poly1d(arms_9_12_IDC_f0_polyfit)
    arms_13_16_IDC_f0_polyfunc = np.poly1d(arms_13_16_IDC_f0_polyfit)
    arms_17_20_IDC_f0_polyfunc = np.poly1d(arms_17_20_IDC_f0_polyfit)
    arms_21_24_IDC_f0_polyfunc = np.poly1d(arms_21_24_IDC_f0_polyfit)
    arms_25_28_IDC_f0_polyfunc = np.poly1d(arms_25_28_IDC_f0_polyfit)

    arms_1_4_CC_f0_polyfunc = np.poly1d(arms_1_4_CC_f0_polyfit)
    arms_5_8_CC_f0_polyfunc = np.poly1d(arms_5_8_CC_f0_polyfit)
    arms_9_12_CC_f0_polyfunc = np.poly1d(arms_9_12_CC_f0_polyfit)
    arms_13_16_CC_f0_polyfunc = np.poly1d(arms_13_16_CC_f0_polyfit)
    arms_17_20_CC_f0_polyfunc = np.poly1d(arms_17_20_CC_f0_polyfit)
    arms_21_24_CC_f0_polyfunc = np.poly1d(arms_21_24_CC_f0_polyfit)
    arms_25_28_CC_f0_polyfunc = np.poly1d(arms_25_28_CC_f0_polyfit)

    # Getting the IDC and CC lengths from the relevant polynomial functions.
    IDC_array = np.zeros(28)

    if f0 >= arms_1_4_min_max_freqs[0] and f0 <= arms_1_4_min_max_freqs[1]:
        IDC_array[0:4] = arms_1_4_IDC_f0_polyfunc(f0)
        IDC_array[4:] = 1975
        CCL = arms_1_4_CC_f0_polyfunc(f0)

    elif f0 >= arms_5_8_min_max_freqs[0] and f0 <= arms_5_8_min_max_freqs[1]:
        IDC_array[0:4] = 1100
        IDC_array[4:8] = arms_5_8_IDC_f0_polyfunc(f0)
        IDC_array[8:] = 1975
        CCL = arms_5_8_CC_f0_polyfunc(f0)

    elif f0 >= arms_9_12_min_max_freqs[0] and f0 <= arms_9_12_min_max_freqs[1]:
        IDC_array[0:8] = 1100
        IDC_array[8:12] = arms_9_12_IDC_f0_polyfunc(f0)
        IDC_array[12:] = 1975
        CCL = arms_9_12_CC_f0_polyfunc(f0)

    elif f0 >= arms_13_16_min_max_freqs[0] and f0 <= arms_13_16_min_max_freqs[1]:
        IDC_array[0:12] = 1100
        IDC_array[12:16] = arms_13_16_IDC_f0_polyfunc(f0)
        IDC_array[16:] = 1975
        CCL = arms_13_16_CC_f0_polyfunc(f0)

    elif f0 >= arms_17_20_min_max_freqs[0] and f0 <= arms_17_20_min_max_freqs[1]:
        IDC_array[0:16] = 1100
        IDC_array[16:20] = arms_17_20_IDC_f0_polyfunc(f0)
        IDC_array[20:] = 1975
        CCL = arms_17_20_CC_f0_polyfunc(f0)

    elif f0 >= arms_21_24_min_max_freqs[0] and f0 <= arms_21_24_min_max_freqs[1]:
        IDC_array[0:20] = 1100
        IDC_array[20:24] = arms_21_24_IDC_f0_polyfunc(f0)
        IDC_array[24:] = 1975
        CCL = arms_21_24_CC_f0_polyfunc(f0)

    elif f0 >= arms_25_28_min_max_freqs[0] and f0 <= arms_25_28_min_max_freqs[1]:
        IDC_array[0:24] = 1100
        IDC_array[24:] = arms_25_28_IDC_f0_polyfunc(f0)
        CCL = arms_25_28_CC_f0_polyfunc(f0)

    else:
        print("error occured")
        return Exception("Error in obtaining IDCL and CCL.")

    if isinstance(rounding_precision, int):
        IDC_array = np.round(IDC_array, rounding_precision)
        CCL = round(CCL, rounding_precision)
        return IDC_array, CCL

    else:
        return IDC_array, CCL
