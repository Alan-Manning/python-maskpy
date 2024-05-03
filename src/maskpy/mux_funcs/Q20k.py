import numpy as np


def get_resonator_config() -> dict:
    """Get the config for the Q20k resonator.

    This is the config used in the muxing for the Q50k design. This means only
    this config is guaranteed to generate Q50k pixels reliably.

    Returns
    -------
    config : dict
        This is a config dictionary with key value pairs as variable name and
        variable values.
    """

    raise (NotImplementedError("Not done Yet"))
    # config = {
    #     "meander_lw": 2,
    #     "meander_corner_bend_radius": 6,
    #     "meander_bot_width": 18,
    #     "meander_right_height_1": 41,
    #     "meander_right_height_2": 24,
    #     "meander_right_height_3": 73,
    #     "meander_left_height_1": 24,
    #     "meander_left_height_2": 58,
    #     "meander_left_height_3": 56,
    #     "meander_left_width_1": 564,
    #     "meander_left_width_2": 564,
    #     "meander_right_width_1": 565,
    #     "meander_right_width_2": 565,
    #     "ant_pad_box_width": 5,
    #     "ant_pad_box_height": 10,
    #     "frame_bot_lw": 8,
    #     "frame_bot_left_width": 996,
    #     "frame_bot_right_width": 996,
    #     "frame_left_lw": 8,
    #     "frame_left_height": 400,
    #     "frame_right_lw": 8,
    #     "frame_right_height": 400,
    #     "frame_meander_cover_box_width": 5,
    #     "frame_meander_cover_box_height": 48,
    #     "coupler_frame_left_lw": 10,
    #     "coupler_frame_left_height": 39,
    #     "coupler_frame_top_lw": 3,
    #     "IDC_bot_arm_gap": 30,
    #     "IDC_arm_gap": 8,
    #     "IDC_arm_lw": 3,
    #     "No_of_arms": 28,
    #     "trim_arm_offset_right_side": 380,
    #     "trim_arm_offset_left_side": 389,
    #     "trim_arm_lw": 3,
    #     "trim_arm_length_right_side": 1975,
    #     "trim_arm_length_left_side": 1975,
    #     "coupler_gap": 16,
    #     "coupler_lw": 3,
    #     "left_coupler_frame_to_feed_distance": 164,
    #     "text_size": 90,
    #     "text_x_offset": 800,
    #     "text_y_offset": 900,
    #     "cutout_bot_offset": 15,
    #     "cutout_left_offset": 50,
    #     "cutout_right_offset": 50,
    #     "cutout_top_offset": 25,
    #     "grndpl_meander_cutout_width": 80,
    #     "grndpl_meander_cutout_height": 10,
    #     "grndpl_coupler_cutout_width": 10,
    #     "grndpln_gap_between_adjacent_resonators": 44,
    #     "grndpl_coupler_cutout_height": 30,
    #     "SiO_stepdown_cutout_width": 110,
    #     "SiO_stepdown_cutout_height": 0,
    #     "SiN_membrane_stepdown_cutout_width": 100,
    #     "SiN_membrane_stepdown_cutout_height": 0,
    # }
    # return config


def get_IDCL_and_CCL_from_f0(f0: float, rounding_precision: int | None = 2):
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
    min_freq = 1962200000.0
    max_freq = 4154170000.0

    if f0 < min_freq or f0 > max_freq:
        print("error occured")
        raise (Exception("Desired frequency out of range. Should be within {0} -> {1}".format(min_freq, max_freq)))

    # # Min Max freqs for each arm
    # for key in all_data.keys():
    #     print(f"{key}_min_max_freqs = [{ min(all_data[key]["f0s"]) }, { max(all_data[key]["f0s"]) }]")

    arms_1_4_min_max_freqs = [1968093809.30539, 2077841016.3035378]
    arms_5_8_min_max_freqs = [2077841016.3035378, 2236037336.485638]
    arms_9_12_min_max_freqs = [2236037336.485638, 2439068400.644538]
    arms_13_16_min_max_freqs = [2439068400.644538, 2704389674.018148]
    arms_17_20_min_max_freqs = [2704389674.018148, 3065341721.843039]
    arms_21_24_min_max_freqs = [3065341721.843039, 3568551825.367324]
    arms_25_28_min_max_freqs = [3568551825.367324, 4167558790.3447323]

    # # Fits from data.
    # for key in all_data.keys():
    #     polyfit = np.polyfit(all_data[key]["f0s"], all_data[key]["IDCLs"], deg=3)
    #     print(
    #         f"{key}_IDC_f0_polyfit = [{polyfit[0]}, {polyfit[1]}, {polyfit[2]}, {polyfit[3]}]"
    #     )

    arms_1_4_IDC_f0_polyfit = [
        5.2916752596151724e-23,
        -3.220289476824636e-13,
        0.0006451083239749925,
        -423708.60305550456,
    ]
    arms_5_8_IDC_f0_polyfit = [
        8.527444094481647e-24,
        -5.4573127445516176e-14,
        0.00011081942868031958,
        -69173.34484372247,
    ]
    arms_9_12_IDC_f0_polyfit = [
        3.48698148001402e-24,
        -2.3713586520950942e-14,
        4.935834482305299e-05,
        -28811.019128959422,
    ]
    arms_13_16_IDC_f0_polyfit = [
        9.657307171332723e-25,
        -6.7458134853809094e-15,
        1.2220931154257412e-05,
        -1713.9317068063656,
    ]
    arms_17_20_IDC_f0_polyfit = [
        1.490489999187265e-25,
        -7.349129329243865e-16,
        -1.909783546518452e-06,
        9567.117329169083,
    ]
    arms_21_24_IDC_f0_polyfit = [
        -8.974846669808551e-26,
        1.2423081228757256e-15,
        -7.012076328599261e-06,
        14381.709284953267,
    ]
    arms_25_28_IDC_f0_polyfit = [
        -2.8570074521003694e-25,
        3.3567644664989025e-15,
        -1.4579236234424436e-05,
        24238.548621322716,
    ]

    # for key in all_data.keys():
    #     polyfit = np.polyfit(all_data[key]["f0s"], all_data[key]["CCLs"], deg=3)
    #     print(
    #         f"{key}_CC_f0_polyfit = [{polyfit[0]}, {polyfit[1]}, {polyfit[2]}, {polyfit[3]}]"
    #     )

    arms_1_4_CC_f0_polyfit = [
        1.6292725296037252e-23,
        -8.381281854073881e-14,
        0.00013692658971607707,
        -67218.69104759378,
    ]
    arms_5_8_CC_f0_polyfit = [
        -7.8232162314757e-26,
        6.827895390098974e-15,
        -3.0276651896979475e-05,
        35729.389155015495,
    ]
    arms_9_12_CC_f0_polyfit = [
        -1.1346144061078278e-23,
        8.257609011508353e-14,
        -0.00020132221490716,
        165438.8200583322,
    ]
    arms_13_16_CC_f0_polyfit = [
        -3.252517023993467e-24,
        2.622938673119096e-14,
        -7.111319255092813e-05,
        65619.96856922799,
    ]
    arms_17_20_CC_f0_polyfit = [
        -2.1723388372763496e-24,
        1.922293669563417e-14,
        -5.711535372390779e-05,
        57642.999746188274,
    ]
    arms_21_24_CC_f0_polyfit = [
        -6.910787673760656e-25,
        7.103761446476217e-15,
        -2.4615627438850808e-05,
        29232.749179783306,
    ]
    arms_25_28_CC_f0_polyfit = [
        -1.3680505958180762e-25,
        1.7567007325093302e-15,
        -7.649066386263858e-06,
        11590.190841248128,
    ]

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

    if f0 >= arms_1_4_min_max_freqs[0] and f0 <= arms_1_4_min_max_freqs[1]:
        IDCL = arms_1_4_IDC_f0_polyfunc(f0)
        CCL = arms_1_4_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:4] = IDCL
        IDC_array[4:] = 1975

        return IDC_array, CCL

    elif f0 >= arms_5_8_min_max_freqs[0] and f0 <= arms_5_8_min_max_freqs[1]:
        IDCL = arms_5_8_IDC_f0_polyfunc(f0)
        CCL = arms_5_8_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:4] = 1100
        IDC_array[4:8] = IDCL
        IDC_array[8:] = 1975

        return IDC_array, CCL

    elif f0 >= arms_9_12_min_max_freqs[0] and f0 <= arms_9_12_min_max_freqs[1]:
        IDCL = arms_9_12_IDC_f0_polyfunc(f0)
        CCL = arms_9_12_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:8] = 1100
        IDC_array[8:12] = IDCL
        IDC_array[12:] = 1975

        return IDC_array, CCL

    elif f0 >= arms_13_16_min_max_freqs[0] and f0 <= arms_13_16_min_max_freqs[1]:
        IDCL = arms_13_16_IDC_f0_polyfunc(f0)
        CCL = arms_13_16_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:12] = 1100
        IDC_array[12:16] = IDCL
        IDC_array[16:] = 1975

        return IDC_array, CCL

    elif f0 >= arms_17_20_min_max_freqs[0] and f0 <= arms_17_20_min_max_freqs[1]:
        IDCL = arms_17_20_IDC_f0_polyfunc(f0)
        CCL = arms_17_20_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:16] = 1100
        IDC_array[16:20] = IDCL
        IDC_array[20:] = 1975

        return IDC_array, CCL

    elif f0 >= arms_21_24_min_max_freqs[0] and f0 <= arms_21_24_min_max_freqs[1]:
        IDCL = arms_21_24_IDC_f0_polyfunc(f0)
        CCL = arms_21_24_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:20] = 1100
        IDC_array[20:24] = IDCL
        IDC_array[24:] = 1975

        return IDC_array, CCL

    elif f0 >= arms_25_28_min_max_freqs[0] and f0 <= arms_25_28_min_max_freqs[1]:
        IDCL = arms_25_28_IDC_f0_polyfunc(f0)
        CCL = arms_25_28_CC_f0_polyfunc(f0)

        IDC_array = np.zeros(28)

        IDC_array[0:24] = 1100
        IDC_array[24:] = IDCL

        return IDC_array, CCL

    else:
        print("error occured")
        return Exception("error in obtaining IDCL and CCL")
