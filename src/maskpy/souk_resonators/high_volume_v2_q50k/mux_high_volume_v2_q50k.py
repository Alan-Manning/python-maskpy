import numpy as np


def get_IDCLs_and_CCL_from_f0(
    f0: float,
    rounding_precision: int | None = 2,
) -> tuple[list[float | int], float | int]:
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
    from ..original_q50k.mux_original_q50k import get_IDCLs_and_CCL_from_f0 as copy_of_mux_original_q50k

    print(f"\033[93mWarning: using mux_original_q50k\033[0m")
    return copy_of_mux_original_q50k(f0, rounding_precision=rounding_precision)
    raise NotImplementedError("Not Done Yet.")
