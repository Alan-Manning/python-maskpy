import numpy as np


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
    # Min Max freqs for each arm
    arms_1_4_min_max_freqs = [1962193055.745924, 2071546769.55142]
    arms_5_8_min_max_freqs = [2071546769.55142, 2229083971.505298]
    arms_9_12_min_max_freqs = [2229083971.505298, 2431651979.890808]
    arms_13_16_min_max_freqs = [2431651979.890808, 2696491178.751894]
    arms_17_20_min_max_freqs = [2696491178.751894, 3056186630.712749]
    arms_21_24_min_max_freqs = [3056186630.712749, 3557176493.872386]
    arms_25_28_min_max_freqs = [3557176493.872386, 4154184936.9751763]

    min_freq = arms_1_4_min_max_freqs[0]
    max_freq = arms_25_28_min_max_freqs[1]

    if f0 < min_freq or f0 > max_freq:
        print("error occured")
        raise (Exception(f"Desired frequency out of range. Should be within range {min_freq} -> {max_freq}"))

    # Fits from data.
    arms_1_4_IDC_f0_polyfit = [6.767696709327338e-23, -4.0998243440291944e-13, 0.000819677631910575, -539164.0136003464]
    arms_5_8_IDC_f0_polyfit = [6.730329291630004e-24, -4.268600455043413e-14, 8.462222866515004e-05, -49975.12714434396]
    arms_9_12_IDC_f0_polyfit = [3.485444963206079e-24, -2.3590757667167347e-14, 4.8811269516983334e-05, -28215.313711449875]
    arms_13_16_IDC_f0_polyfit = [9.845730863469365e-25, -6.877032261889292e-15, 1.2526532983062344e-05, -1977.6667747536535]
    arms_17_20_IDC_f0_polyfit = [9.443531070041871e-26, -2.654295099140603e-16, -3.2524344435929964e-06, 10823.993211480463]
    arms_21_24_IDC_f0_polyfit = [-1.1147316377050832e-25, 1.4570979929442625e-15, -7.719109582100416e-06, 15138.846451306348]
    arms_25_28_IDC_f0_polyfit = [-2.801791330399383e-25, 3.2864652121492133e-15, -1.4287781040336895e-05, 23825.298205884836]

    arms_1_4_CC_f0_polyfit = [-1.6060511535332507e-22, 9.798359652954358e-13, -0.0019941373079622334, 1355558.2923700206]
    arms_5_8_CC_f0_polyfit = [3.826678033412186e-23, -2.4359149810072817e-13, 0.0005147621732587452, -359542.0813002416]
    arms_9_12_CC_f0_polyfit = [-1.0365061470132613e-23, 7.48295766256644e-14, -0.00018113510144262927, 148136.87440822576]
    arms_13_16_CC_f0_polyfit = [-2.617278235412754e-24, 2.1390330519786914e-14, -5.8900447871742394e-05, 55487.868982187225]
    arms_17_20_CC_f0_polyfit = [-3.3011190814655043e-25, 3.4999958695322288e-15, -1.2482352491053908e-05, 15561.645803368623]
    arms_21_24_CC_f0_polyfit = [-2.3442630698128755e-25, 2.5830613479870064e-15, -9.739844107703402e-06, 13012.011320133017]
    arms_25_28_CC_f0_polyfit = [-2.138881161222542e-25, 2.619737848761559e-15, -1.0867933721816246e-05, 15636.17994038814]

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
        IDCL = arms_1_4_IDC_f0_polyfunc(f0)
        CCL = arms_1_4_CC_f0_polyfunc(f0)
        IDC_array[0:4] = IDCL
        IDC_array[4:] = 1975

    elif f0 >= arms_5_8_min_max_freqs[0] and f0 <= arms_5_8_min_max_freqs[1]:
        IDCL = arms_5_8_IDC_f0_polyfunc(f0)
        CCL = arms_5_8_CC_f0_polyfunc(f0)
        IDC_array[0:4] = 1100
        IDC_array[4:8] = IDCL
        IDC_array[8:] = 1975

    elif f0 >= arms_9_12_min_max_freqs[0] and f0 <= arms_9_12_min_max_freqs[1]:
        IDCL = arms_9_12_IDC_f0_polyfunc(f0)
        CCL = arms_9_12_CC_f0_polyfunc(f0)
        IDC_array[0:8] = 1100
        IDC_array[8:12] = IDCL
        IDC_array[12:] = 1975

    elif f0 >= arms_13_16_min_max_freqs[0] and f0 <= arms_13_16_min_max_freqs[1]:
        IDCL = arms_13_16_IDC_f0_polyfunc(f0)
        CCL = arms_13_16_CC_f0_polyfunc(f0)
        IDC_array[0:12] = 1100
        IDC_array[12:16] = IDCL
        IDC_array[16:] = 1975

    elif f0 >= arms_17_20_min_max_freqs[0] and f0 <= arms_17_20_min_max_freqs[1]:
        IDCL = arms_17_20_IDC_f0_polyfunc(f0)
        CCL = arms_17_20_CC_f0_polyfunc(f0)
        IDC_array[0:16] = 1100
        IDC_array[16:20] = IDCL
        IDC_array[20:] = 1975

    elif f0 >= arms_21_24_min_max_freqs[0] and f0 <= arms_21_24_min_max_freqs[1]:
        IDCL = arms_21_24_IDC_f0_polyfunc(f0)
        CCL = arms_21_24_CC_f0_polyfunc(f0)
        IDC_array[0:20] = 1100
        IDC_array[20:24] = IDCL
        IDC_array[24:] = 1975

    elif f0 >= arms_25_28_min_max_freqs[0] and f0 <= arms_25_28_min_max_freqs[1]:
        IDCL = arms_25_28_IDC_f0_polyfunc(f0)
        CCL = arms_25_28_CC_f0_polyfunc(f0)
        IDC_array[0:24] = 1100
        IDC_array[24:] = IDCL

    else:
        print("error occured")
        return Exception("error in obtaining IDCL and CCL")

    if isinstance(rounding_precision, int):
        IDC_array = np.round(IDC_array, rounding_precision)
        CCL = round(CCL, rounding_precision)
        return IDC_array, CCL

    else:
        return IDC_array, CCL
