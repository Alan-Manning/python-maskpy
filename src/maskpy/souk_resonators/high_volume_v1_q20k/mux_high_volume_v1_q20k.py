import numpy as np
from numpy import float64
from numpy._typing import NDArray


def get_IDCLs_and_CCL_from_f0(
    f0: float,
    rounding_precision: int | None = 2,
) -> tuple[NDArray[float64], float | int]:
    """Get the IDC lengths and CC length (both in microns) for a resonator
    given a resonant frequency in Hz. These lengths are by default rounded to
    2 decimal places.

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

    # for key in all_data.keys():
    #     print(f"arms_{key}_min_max_freqs = [{ min(all_data[key]["f0s"]) }, { max(all_data[key]["f0s"]) }]")

    # Min Max freqs for each arm
    arms_1_4_min_max_freqs = [1976460099.789582, 2086954601.9639728]
    arms_5_8_min_max_freqs = [2086954601.9639728, 2243852362.4006658]
    arms_9_12_min_max_freqs = [2243852362.4006658, 2444364619.4824467]
    arms_13_16_min_max_freqs = [2444364619.4824467, 2706314959.9834857]
    arms_17_20_min_max_freqs = [2706314959.9834857, 3060964703.562229]
    arms_21_24_min_max_freqs = [3060964703.562229, 3550726031.423563]
    arms_25_28_min_max_freqs = [3550726031.423563, 4129398782.8582683]

    min_freq = arms_1_4_min_max_freqs[0]
    max_freq = arms_25_28_min_max_freqs[1]

    if f0 < min_freq or f0 > max_freq:
        print("error occured")
        raise (Exception("Desired frequency out of range. Should be within {0} -> {1}".format(min_freq, max_freq)))

    # # Fits from data for IDC against f0.
    # for key in all_data.keys():
    #     polyfit = np.polyfit(all_data[key]["f0s"], all_data[key]["IDCLs"], deg=3)
    #     print(
    #         f"arms_{key}_IDC_f0_polyfit = [\n\t{polyfit[0]},\n\t{polyfit[1]},\n\t{polyfit[2]},\n\t{polyfit[3]}\n]"
    #     )

    # IDC against resonant frequency polyfits
    arms_1_4_IDC_f0_polyfit = [
        5.1746329481385594e-23,
        -3.154258628052129e-13,
        0.0006328301700505949,
        -416135.23674744234,
    ]
    arms_5_8_IDC_f0_polyfit = [
        9.45049223630743e-24,
        -6.0560659257401e-14,
        0.00012370364419712063,
        -78323.83135198422,
    ]
    arms_9_12_IDC_f0_polyfit = [
        3.3081758313347738e-24,
        -2.2450290673778175e-14,
        4.632197095301158e-05,
        -26304.016583370547,
    ]
    arms_13_16_IDC_f0_polyfit = [
        8.935722173751798e-25,
        -6.195928436331369e-15,
        1.0778566803574938e-05,
        -401.6795618478126,
    ]
    arms_17_20_IDC_f0_polyfit = [
        1.0212643731840563e-25,
        -3.3448850900694566e-16,
        -3.0887669467708895e-06,
        10760.182453763355,
    ]
    arms_21_24_IDC_f0_polyfit = [
        -1.0219700222179101e-25,
        1.3623609883752806e-15,
        -7.437175660056797e-06,
        14906.697303799336,
    ]
    arms_25_28_IDC_f0_polyfit = [
        -2.8473047005303024e-25,
        3.328914434742875e-15,
        -1.445836621777353e-05,
        24089.574628272563,
    ]

    # # Fits from data for IDC against f0.
    # for key in all_data.keys():
    #     polyfit = np.polyfit(all_data[key]["f0s"], all_data[key]["CCLs"], deg=3)
    #     print(
    #         f"arms_{key}_CC_f0_polyfit = [\n\t{polyfit[0]},\n\t{polyfit[1]},\n\t{polyfit[2]},\n\t{polyfit[3]}\n]"
    #     )
    # Coupler against resonant frequency polyfits
    arms_1_4_CC_f0_polyfit = [
        1.6487037725936303e-24,
        -7.316016579293542e-15,
        7.714968047933058e-06,
        1793.856757234163,
    ]
    arms_5_8_CC_f0_polyfit = [
        -2.5905088221176984e-24,
        1.7411022190506688e-14,
        -3.988762979929466e-05,
        31972.59549700681,
    ]
    arms_9_12_CC_f0_polyfit = [
        -6.950694114162767e-24,
        4.975624183187765e-14,
        -0.00011929808584302062,
        96565.17728704738,
    ]
    arms_13_16_CC_f0_polyfit = [
        -1.036673482048879e-24,
        8.600914469127921e-15,
        -2.4134204495861857e-05,
        23477.576426685464,
    ]
    arms_17_20_CC_f0_polyfit = [
        -4.047202255588042e-25,
        3.861223698216885e-15,
        -1.2511747109858653e-05,
        14211.10915839199,
    ]
    arms_21_24_CC_f0_polyfit = [
        -3.4381431098592545e-25,
        3.608323862967447e-15,
        -1.282382872726415e-05,
        15789.493403044815,
    ]
    arms_25_28_CC_f0_polyfit = [
        -2.501415819356193e-25,
        3.0066878531124996e-15,
        -1.217281449337219e-05,
        16869.90262788116,
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
        raise Exception("error in obtaining IDCL and CCL")

    if isinstance(rounding_precision, int):
        IDC_array = np.round(IDC_array, rounding_precision)
        CCL = round(CCL, rounding_precision)
        return IDC_array, CCL

    else:
        return IDC_array, CCL
