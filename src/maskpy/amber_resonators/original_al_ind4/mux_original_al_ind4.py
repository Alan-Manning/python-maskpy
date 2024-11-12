import numpy as np
from numpy import float64
from numpy.typing import NDArray


def get_freq_bands() -> tuple[
    tuple[float, float],
    tuple[float, float],
    tuple[float, float],
]:
    """Get the frequency low mid and high band bounds over which the resonator
    is muxed.

    Returns
    -------
    freq_bands: tuple[
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
    ]
        The frequency low, mid and high bounds for the frequency bands.
    """

    # Min Max freqs for each band
    low_band_min_max_freqs = (856392785.829775, 873086753.419883)
    mid_band_min_max_freqs = (952601900.03232, 976316686.3428609)
    high_band_min_max_freqs = (1088353533.869286, 1125251289.0911)

    return (low_band_min_max_freqs, mid_band_min_max_freqs, high_band_min_max_freqs)


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
    IDC_array : NDArray[float64]
        This is a list of length 13, 17 or 21. Each element is a float representing the
        length of an IDC arm in um. The first element in this list is the
        length of arm 1, the second element is the length of the second arm
        and so on...

    CCL : float
        This is the length of the coupling capacitor in um.
    """

    low_band_min_max_freqs, mid_band_min_max_freqs, high_band_min_max_freqs = get_freq_bands()

    max_idcl = 3960
    CCL = 990 + 1066  # (length of coupler in the fork) + (length to very right side of the right frame).

    # Check if in one of these bands else raise exception.
    if low_band_min_max_freqs[0] <= f0 <= low_band_min_max_freqs[1]:
        fit = np.poly1d(np.array([-1.15114570e-19, 2.95381222e-10, -2.52842282e-01, 7.22029639e07]))
        no_of_arms = 21
    elif mid_band_min_max_freqs[0] <= f0 <= mid_band_min_max_freqs[1]:
        fit = np.poly1d(np.array([-3.87758171e-20, 1.10598172e-10, -1.05288380e-01, 3.34588400e07]))
        no_of_arms = 17
    elif high_band_min_max_freqs[0] <= f0 <= high_band_min_max_freqs[1]:
        fit = np.poly1d(np.array([-9.78156192e-21, 3.18317344e-11, -3.46177814e-02, 1.25853135e07]))
        no_of_arms = 13
    else:
        raise ValueError(
            f"Desired frequency, {f0}, out of band ranges. Should be within:\n\tLow: {low_band_min_max_freqs[0]} -> {low_band_min_max_freqs[1]}\n\tMid: {mid_band_min_max_freqs[0]} -> {mid_band_min_max_freqs[1]}\n\tHigh: {low_band_min_max_freqs[0]} -> {low_band_min_max_freqs[1]}"
        )

    top_arm_idcl = fit(f0)
    IDC_array = np.ones(no_of_arms) * max_idcl
    IDC_array[-1] = top_arm_idcl

    if isinstance(rounding_precision, int):
        IDC_array = np.round(IDC_array, rounding_precision)
        CCL = round(CCL, rounding_precision)
        return IDC_array, CCL

    else:
        return IDC_array, CCL
