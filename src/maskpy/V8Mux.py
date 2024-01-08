import numpy as np
import pandas as pd


def V8_get_evenly_spaced_freq_array(
    number_of_resonators, min_freq=1970120000.0, max_freq=4174530000.0, order=None, csv_output=False, csv_filename=""
):
    """
    This returns a list of resonant frequencies for the resonators in an array.

    Parameters
    ----------
    number_of_resonators : int
        This is the number of resonators to generate

    KwArgs
    ------
    min_freq = 1970120000.0
        This is the minimum frequency for the array. Should be a float or int.
        The default value is the minumum for the V8 mux.

    max_freq = 4174530000.0
        This is the maximum frequency for the arrat. Should be a float or int.
        The default value is the maximum for the V8 mux.

    order = None
        This is the order for the resonant frequency array. When defined this
        should be a list or numpy array.

    csv_output = False
        When true will produce a CSV output of all the resonant frequencies,
        KID IDs, IDC lengths and CC length. This will be saved under the name
        specified in the csv_filename argument if specified or by default will
        be called "csv_mask_mux_output".

    csv_filename = ""
        The name of the csv_output when it is True. This has no effect when
        csv_output is not True. This filename should not include the ".csv"
        file extention, just include the name before this. If not specified
        the csv output will take the default name "csv_mask_mux_output".

    Returns
    -------

    f0s : np.array
        The list of resonant frequencies for the array.

    """

    f0s = np.linspace(min_freq, max_freq, number_of_resonators)

    if csv_output:
        if isinstance(csv_filename, str) and csv_filename != "":
            filename = csv_filename + ".csv"
        else:
            filename = "csv_mask_mux_output.csv"

        print(f"Making CSV Output with Resonator information as '{filename}'")

        column_names = [
            "Freq(Hz)",
            "KID_ID",
            "CCL",
            "IDC_ARM_1",
            "IDC_ARM_2",
            "IDC_ARM_3",
            "IDC_ARM_4",
            "IDC_ARM_5",
            "IDC_ARM_6",
            "IDC_ARM_7",
            "IDC_ARM_8",
            "IDC_ARM_9",
            "IDC_ARM_10",
            "IDC_ARM_11",
            "IDC_ARM_12",
            "IDC_ARM_13",
            "IDC_ARM_14",
            "IDC_ARM_15",
            "IDC_ARM_16",
            "IDC_ARM_17",
            "IDC_ARM_18",
            "IDC_ARM_19",
            "IDC_ARM_20",
            "IDC_ARM_21",
            "IDC_ARM_22",
            "IDC_ARM_23",
            "IDC_ARM_24",
            "IDC_ARM_25",
            "IDC_ARM_26",
            "IDC_ARM_27",
            "IDC_ARM_28",
        ]
        data_for_datafrane = np.zeros((number_of_resonators, len(column_names)))

        for i, f0 in enumerate(f0s):
            IDCLS, CCL = V8_get_IDCL_and_CCL_from_f0(f0)
            data_for_datafrane[i][0:3] = f0, i, CCL
            data_for_datafrane[i][3:] = IDCLS

        datafrane = pd.DataFrame(data=data_for_datafrane, columns=column_names)
        datafrane.to_csv(filename, index=False)

    if isinstance(order, (np.ndarray, list)):
        ordered_f0s = np.zeros(number_of_resonators)

        for i in range(number_of_resonators):
            ordered_f0s[i] = f0s[order[i]]

        return ordered_f0s

    return f0s


def V8_get_CC_from_IDCL(IDCL_array):
    """
    This obtains the coupling capacitor length from the lengths of the IDC arms
    using the V8 muxing function.

    Parameters
    ----------
    IDC_array : list
        This is a list of length 28. Each element is a float representing the
        length of an IDC arm in um. The first element in this list is the
        length of arm 1, the second element is the length of the second arm
        and so on...

    Returns
    -------
    CCL : float
        This is the length of the coupling capacitor in um.
    """

    IDCL_min = 1100

    Arms = 0
    arm_len = 0
    for i, length in enumerate(IDCL_array):
        if length == IDCL_min:
            Arms = "Arms_" + str(i + 2) + "_" + str(i + 5)
            arm_len = length

    # Fits from data.
    Arms_1_4_CCL_IDC_polyfit = [-1.81405896e-07, 9.44606414e-04, -1.35421121e00, 1.41799077e03]
    Arms_5_8_CCL_IDC_polyfit = [1.55490768e-07, -5.88921283e-04, 8.49902818e-01, 2.49672012e02]
    Arms_9_12_CCL_IDC_polyfit = [5.18302559e-08, -1.80758017e-04, 3.15759637e-01, 3.61562196e02]
    Arms_13_16_CCL_IDC_polyfit = [8.63837599e-08, -3.45966958e-04, 5.60328258e-01, 1.37678976e02]
    Arms_17_20_CCL_IDC_polyfit = [7.08346831e-08, -2.73663751e-04, 4.45935644e-01, 9.09340784e01]
    Arms_21_24_CCL_IDC_polyfit = [1.24392614e-07, -5.36443149e-04, 8.75024295e-01, -2.50854227e02]
    Arms_25_28_CCL_IDC_polyfit = [5.01025807e-08, -1.76870748e-04, 2.76162401e-01, -1.14356981e01]

    Arms_1_4_CCL_IDC_polyfunc = np.poly1d(Arms_1_4_CCL_IDC_polyfit)
    Arms_5_8_CCL_IDC_polyfunc = np.poly1d(Arms_5_8_CCL_IDC_polyfit)
    Arms_9_12_CCL_IDC_polyfunc = np.poly1d(Arms_9_12_CCL_IDC_polyfit)
    Arms_13_16_CCL_IDC_polyfunc = np.poly1d(Arms_13_16_CCL_IDC_polyfit)
    Arms_17_20_CCL_IDC_polyfunc = np.poly1d(Arms_17_20_CCL_IDC_polyfit)
    Arms_21_24_CCL_IDC_polyfunc = np.poly1d(Arms_21_24_CCL_IDC_polyfit)
    Arms_25_28_CCL_IDC_polyfunc = np.poly1d(Arms_25_28_CCL_IDC_polyfit)

    if Arms == "Arms_1_4":
        return Arms_1_4_CCL_IDC_polyfunc(arm_len)

    elif Arms == "Arms_5_8":
        return Arms_5_8_CCL_IDC_polyfunc(arm_len)

    elif Arms == "Arms_9_12":
        return Arms_9_12_CCL_IDC_polyfunc(arm_len)

    elif Arms == "Arms_13_16":
        return Arms_13_16_CCL_IDC_polyfunc(arm_len)

    elif Arms == "Arms_17_20":
        return Arms_17_20_CCL_IDC_polyfunc(arm_len)

    elif Arms == "Arms_21_24":
        return Arms_21_24_CCL_IDC_polyfunc(arm_len)

    elif Arms == "Arms_25_28":
        return Arms_25_28_CCL_IDC_polyfunc(arm_len)

    else:
        raise (Exception("Error obtaining CCL from IDCL array"))


def V8_get_IDCL_and_CCL_from_f0(f0, rounding_precision=2):
    """
    This obtains the IDC lengths and CC length of a resonator for a given
    resonant frequency. These lengths are by default rounded to 2 decimal
    places.

    Parameters
    ----------
    f0 : float, int
        The desired frequency (**in Hz**) for the KID.

    KwArgs
    ------
    rounding_precision=2
        This should be an integer. This is the amount of decimal places to
        round the result of the IDC and CC lengths to. By default (recomended)
        this is 2 decimal places. For no rounding this value should be very
        large eg 50.

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
    min_freq = 1970120000.0
    max_freq = 4174530000.0

    if f0 < min_freq or f0 > max_freq:
        print("error occured")
        raise (Exception(f"Desired frequency out of range. Should be within {min_freq} -> {max_freq}"))

    Arms_1_4_min_max_freqs = [1970120000.0, 2080770000.0]
    Arms_5_8_min_max_freqs = [2080770000.0, 2239280000.0]
    Arms_9_12_min_max_freqs = [2239280000.0, 2441040000.0]
    Arms_13_16_min_max_freqs = [2441040000.0, 2705890000.0]
    Arms_17_20_min_max_freqs = [2705890000.0, 3066980000.0]
    Arms_21_24_min_max_freqs = [3066980000.0, 3571600000.0]
    Arms_25_28_min_max_freqs = [3571600000.0, 4174530000.0]

    # Fits from data.
    Arms_1_4_IDC_f0_polyfit = [4.40593480e-23, -2.67909667e-13, 5.34980654e-04, -3.49055119e05]
    Arms_5_8_IDC_f0_polyfit = [1.87200482e-23, -1.21162839e-13, 2.55765811e-04, -1.74274315e05]
    Arms_9_12_IDC_f0_polyfit = [3.12896341e-24, -2.11432074e-14, 4.31831099e-05, -2.38375201e04]
    Arms_13_16_IDC_f0_polyfit = [9.30798289e-25, -6.46468470e-15, 1.14603758e-05, -1.01767002e03]
    Arms_17_20_IDC_f0_polyfit = [1.14601699e-25, -4.38019706e-16, -2.76247107e-06, 1.03869882e04]
    Arms_21_24_IDC_f0_polyfit = [-9.89518123e-26, 1.33651102e-15, -7.32943707e-06, 1.47377098e04]
    Arms_25_28_IDC_f0_polyfit = [-2.78899336e-25, 3.28681434e-15, -1.43346369e-05, 2.39520423e04]

    Arms_1_4_CC_f0_polyfit = [1.02325925e-22, -6.15094874e-13, 1.23020222e-03, -8.17665744e05]
    Arms_5_8_CC_f0_polyfit = [-2.24898154e-23, 1.49655879e-13, -3.32533402e-04, 2.47412710e05]
    Arms_9_12_CC_f0_polyfit = [-3.91010408e-24, 2.86521516e-14, -7.04177628e-05, 5.85968409e04]
    Arms_13_16_CC_f0_polyfit = [-3.27444516e-24, 2.59330629e-14, -6.87761882e-05, 6.15470736e04]
    Arms_17_20_CC_f0_polyfit = [-1.09047056e-24, 9.81438047e-15, -2.96571843e-05, 3.04441295e04]
    Arms_21_24_CC_f0_polyfit = [-6.20264683e-25, 6.32823387e-15, -2.16963227e-05, 2.52534933e04]
    Arms_25_28_CC_f0_polyfit = [-1.84259106e-25, 2.25955100e-15, -9.33504064e-06, 1.31424615e04]

    # Making fits into polynomial functions.
    Arms_1_4_IDC_f0_polyfunc = np.poly1d(Arms_1_4_IDC_f0_polyfit)
    Arms_5_8_IDC_f0_polyfunc = np.poly1d(Arms_5_8_IDC_f0_polyfit)
    Arms_9_12_IDC_f0_polyfunc = np.poly1d(Arms_9_12_IDC_f0_polyfit)
    Arms_13_16_IDC_f0_polyfunc = np.poly1d(Arms_13_16_IDC_f0_polyfit)
    Arms_17_20_IDC_f0_polyfunc = np.poly1d(Arms_17_20_IDC_f0_polyfit)
    Arms_21_24_IDC_f0_polyfunc = np.poly1d(Arms_21_24_IDC_f0_polyfit)
    Arms_25_28_IDC_f0_polyfunc = np.poly1d(Arms_25_28_IDC_f0_polyfit)

    Arms_1_4_CC_f0_polyfunc = np.poly1d(Arms_1_4_CC_f0_polyfit)
    Arms_5_8_CC_f0_polyfunc = np.poly1d(Arms_5_8_CC_f0_polyfit)
    Arms_9_12_CC_f0_polyfunc = np.poly1d(Arms_9_12_CC_f0_polyfit)
    Arms_13_16_CC_f0_polyfunc = np.poly1d(Arms_13_16_CC_f0_polyfit)
    Arms_17_20_CC_f0_polyfunc = np.poly1d(Arms_17_20_CC_f0_polyfit)
    Arms_21_24_CC_f0_polyfunc = np.poly1d(Arms_21_24_CC_f0_polyfit)
    Arms_25_28_CC_f0_polyfunc = np.poly1d(Arms_25_28_CC_f0_polyfit)

    # Getting the IDC and CC lengths from the relevant polynomial functions.
    IDC_array = np.zeros(28)

    if f0 >= Arms_1_4_min_max_freqs[0] and f0 <= Arms_1_4_min_max_freqs[1]:
        IDC_array[0:4] = Arms_1_4_IDC_f0_polyfunc(f0)
        IDC_array[4:] = 1975
        CCL = Arms_1_4_CC_f0_polyfunc(f0)

    elif f0 >= Arms_5_8_min_max_freqs[0] and f0 <= Arms_5_8_min_max_freqs[1]:
        IDC_array[0:4] = 1100
        IDC_array[4:8] = Arms_5_8_IDC_f0_polyfunc(f0)
        IDC_array[8:] = 1975
        CCL = Arms_5_8_CC_f0_polyfunc(f0)

    elif f0 >= Arms_9_12_min_max_freqs[0] and f0 <= Arms_9_12_min_max_freqs[1]:
        IDC_array[0:8] = 1100
        IDC_array[8:12] = Arms_9_12_IDC_f0_polyfunc(f0)
        IDC_array[12:] = 1975
        CCL = Arms_9_12_CC_f0_polyfunc(f0)

    elif f0 >= Arms_13_16_min_max_freqs[0] and f0 <= Arms_13_16_min_max_freqs[1]:
        IDC_array[0:12] = 1100
        IDC_array[12:16] = Arms_13_16_IDC_f0_polyfunc(f0)
        IDC_array[16:] = 1975
        CCL = Arms_13_16_CC_f0_polyfunc(f0)

    elif f0 >= Arms_17_20_min_max_freqs[0] and f0 <= Arms_17_20_min_max_freqs[1]:
        IDC_array[0:16] = 1100
        IDC_array[16:20] = Arms_17_20_IDC_f0_polyfunc(f0)
        IDC_array[20:] = 1975
        CCL = Arms_17_20_CC_f0_polyfunc(f0)

    elif f0 >= Arms_21_24_min_max_freqs[0] and f0 <= Arms_21_24_min_max_freqs[1]:
        IDC_array[0:20] = 1100
        IDC_array[20:24] = Arms_21_24_IDC_f0_polyfunc(f0)
        IDC_array[24:] = 1975
        CCL = Arms_21_24_CC_f0_polyfunc(f0)

    elif f0 >= Arms_25_28_min_max_freqs[0] and f0 <= Arms_25_28_min_max_freqs[1]:
        IDC_array[0:24] = 1100
        IDC_array[24:] = Arms_25_28_IDC_f0_polyfunc(f0)
        CCL = Arms_25_28_CC_f0_polyfunc(f0)

    else:
        print("error occured")
        return Exception("Error in obtaining IDCL and CCL.")

    if isinstance(rounding_precision, int):
        IDC_array = np.round(IDC_array, rounding_precision)
        CCL = round(CCL, rounding_precision)
        return IDC_array, CCL

    else:
        return IDC_array, CCL
