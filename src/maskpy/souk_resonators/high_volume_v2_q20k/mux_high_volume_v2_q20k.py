import numpy as np
from numpy import float64
from numpy._typing import NDArray

# def get_IDCLs_and_CCL_from_f0(
#     f0: float,
#     rounding_precision: int | None = 2,
# ) -> tuple[list[float | int], float | int]:
#     """Get the IDC lengths and CC length (both in microns) for a resonator
#     given a resonant frequency in Hz. These lengths are by default rounded to
#     2 decimal places.
#
#     Parameters
#     ----------
#     f0 : float
#         The desired frequency (**in Hz**) for the KID.
#
#     KwArgs
#     ------
#     rounding_precision=2
#         This should be an integer. This is the amount of decimal places to
#         round the result of the IDC and CC lengths to. By default (recomended)
#         this is 2 decimal places. For no rounding this value should be None.
#
#     Returns
#     -------
#     IDC_array : list
#         This is a list of length 28. Each element is a float representing the
#         length of an IDC arm in um. The first element in this list is the
#         length of arm 1, the second element is the length of the second arm
#         and so on...
#
#     CCL : float
#         This is the length of the coupling capacitor in um.
#     """
#     from ..original_q20k.mux_original_q20k import get_IDCLs_and_CCL_from_f0 as copy_of_mux_original_q20k
#
#     print(f"\033[93mWarning: using mux_original_q20k\033[0m")
#     return copy_of_mux_original_q20k(f0, rounding_precision=rounding_precision)
#     raise NotImplementedError("Not Done Yet.")


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
    #     print(
    #         f"arms_{key}_min_max_freqs = [{ min(all_data[key]["f0s"]) }, { max(all_data[key]["f0s"]) }]"
    #     )

    # Min Max freqs for each arm
    arms_1_4_min_max_freqs = [1416670868.347992, 1492629469.5994911]
    arms_5_8_min_max_freqs = [1492629469.5994911, 1600162956.332958]
    arms_9_12_min_max_freqs = [1600162956.332958, 1736949259.8367639]
    arms_13_16_min_max_freqs = [1736949259.8367639, 1915157896.415582]
    arms_17_20_min_max_freqs = [1915157896.415582, 2157641366.576265]
    arms_21_24_min_max_freqs = [2157641366.576265, 2506365009.793208]
    arms_25_28_min_max_freqs = [2506365009.793208, 3019632155.221455]

    min_freq = arms_1_4_min_max_freqs[0]
    max_freq = arms_25_28_min_max_freqs[1]

    if f0 < min_freq or f0 > max_freq:
        print("error occured")
        raise (Exception("Desired frequency out of range. Should be within {0} -> {1}".format(min_freq, max_freq)))

    # Fits from data for IDC against f0.
    # for key in all_data.keys():
    #     polyfit = np.polyfit(all_data[key]["f0s"], all_data[key]["IDCLs"], deg=3)
    #     print(
    #         f"arms_{key}_IDC_f0_polyfit = [\n\t{polyfit[0]},\n\t{polyfit[1]},\n\t{polyfit[2]},\n\t{polyfit[3]}\n]"
    #     )

    # IDC against resonant frequency polyfits
    arms_1_4_IDC_f0_polyfit = [
        1.2218649553017394e-22,
        -5.295982425128092e-13,
        0.0007534229590209816,
        -349896.68029857817,
    ]
    arms_5_8_IDC_f0_polyfit = [
        2.7038435123672166e-23,
        -1.2159576144077885e-13,
        0.00017388208496910745,
        -76573.88034322445,
    ]
    arms_9_12_IDC_f0_polyfit = [
        8.165797765648411e-24,
        -3.7786757425629434e-14,
        5.146137577715686e-05,
        -17074.637177324446,
    ]
    arms_13_16_IDC_f0_polyfit = [
        2.1206118590736494e-24,
        -9.197926947219147e-15,
        7.451929723699484e-06,
        5669.117048216304,
    ]
    arms_17_20_IDC_f0_polyfit = [
        -1.177398971739429e-25,
        2.46044703632626e-15,
        -1.2162518652390365e-05,
        17071.085505476054,
    ]
    arms_21_24_IDC_f0_polyfit = [
        -3.2161665835946403e-25,
        3.3724496557171584e-15,
        -1.2981370958985521e-05,
        17514.980031915504,
    ]
    arms_25_28_IDC_f0_polyfit = [
        -4.49038952427145e-25,
        4.299614224512425e-15,
        -1.5150703219304176e-05,
        20008.949454366193,
    ]

    # # Fits from data for IDC against f0.
    # for key in all_data.keys():
    #     polyfit = np.polyfit(all_data[key]["f0s"], all_data[key]["CCLs"], deg=3)
    #     print(
    #         f"arms_{key}_CC_f0_polyfit = [\n\t{polyfit[0]},\n\t{polyfit[1]},\n\t{polyfit[2]},\n\t{polyfit[3]}\n]"
    #     )

    # Coupler against resonant frequency polyfits
    arms_1_4_CC_f0_polyfit = [
        3.835913669090149e-22,
        -1.6674321173482434e-12,
        0.002413425959371915,
        -1161964.5799516723,
    ]
    arms_5_8_CC_f0_polyfit = [
        -7.130055583199499e-25,
        6.112418345208484e-15,
        -1.5186590563664362e-05,
        12494.230334047053,
    ]
    arms_9_12_CC_f0_polyfit = [
        6.4914182350310316e-24,
        -2.9985180461306897e-14,
        4.4797266040156876e-05,
        -20580.943428266182,
    ]
    arms_13_16_CC_f0_polyfit = [
        3.2493283201847e-24,
        -1.659926619194188e-14,
        2.735947809900348e-05,
        -13687.768399800281,
    ]
    arms_17_20_CC_f0_polyfit = [
        5.2977653967168984e-24,
        -3.1643323372542705e-14,
        6.236718538309709e-05,
        -39942.29825643871,
    ]
    arms_21_24_CC_f0_polyfit = [
        -1.2958139610422926e-24,
        9.396383377706894e-15,
        -2.302621997838463e-05,
        19479.347349623,
    ]
    arms_25_28_CC_f0_polyfit = [
        1.0014057069115719e-25,
        -6.659469006874031e-16,
        1.1405206329661366e-06,
        140.15134848460832,
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
