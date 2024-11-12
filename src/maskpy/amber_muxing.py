from typing import Callable

import numpy as np
import pandas as pd
from numpy import float64
from numpy.typing import NDArray

from .amber_resonators.original_al_ind2 import mux_original_al_ind2
from .amber_resonators.original_al_ind4 import mux_original_al_ind4
from .amber_resonators.original_nb_ind2 import mux_original_nb_ind2
from .amber_resonators.original_nb_ind4 import mux_original_nb_ind4

# from .amber_resonators.original_q50k import mux_original_q50k
from .amber_resonators.resonator_types import AmberResonatorType
from .logging import TextColor, styled_text


def get_mux_func_for_resonator_type(resonator_type: AmberResonatorType) -> Callable[[float], tuple[NDArray[float64], float | int]]:
    """Get the muxing function callable that returns the IDCLs and CCL given a
    resonant frequency in Hz."""
    if not isinstance(resonator_type, AmberResonatorType):
        raise TypeError(
            styled_text(
                f"resonator_type should be of type AmberResonatorType. Current type is: {type(resonator_type)}", color=TextColor.ERROR
            )
        )

    match resonator_type:
        case AmberResonatorType.ORIGINAL_AL_IND2:
            return mux_original_al_ind2.get_IDCLs_and_CCL_from_f0
        case AmberResonatorType.ORIGINAL_AL_IND4:
            return mux_original_al_ind4.get_IDCLs_and_CCL_from_f0
        case AmberResonatorType.ORIGINAL_NB_IND2:
            return mux_original_nb_ind2.get_IDCLs_and_CCL_from_f0
        case AmberResonatorType.ORIGINAL_NB_IND4:
            return mux_original_nb_ind4.get_IDCLs_and_CCL_from_f0
        case _:
            raise (ValueError(f"ResonatorType does not have an associated muxing function."))


def get_evenly_spaced_freq_array(
    number_of_resonators: int,
    resonator_type: AmberResonatorType,
    rounding_precision: int | None = 2,
    order: list | np.ndarray | None = None,
    csv_output: bool = False,
    csv_filename: str = "",
    return_dataframe: bool = False,
) -> NDArray[float64] | tuple[NDArray[float64], pd.DataFrame]:

    match resonator_type:
        case AmberResonatorType.ORIGINAL_AL_IND2:
            get_freq_band_func = mux_original_al_ind2.get_freq_bands
        case AmberResonatorType.ORIGINAL_AL_IND4:
            get_freq_band_func = mux_original_al_ind4.get_freq_bands
        case AmberResonatorType.ORIGINAL_NB_IND2:
            get_freq_band_func = mux_original_nb_ind2.get_freq_bands
        case AmberResonatorType.ORIGINAL_NB_IND4:
            get_freq_band_func = mux_original_nb_ind4.get_freq_bands
        case _:
            raise (ValueError(f"ResonatorType does not have an associated muxing function."))

    low_band_min_max_freqs, mid_band_min_max_freqs, high_band_min_max_freqs = get_freq_band_func()

    no_of_resonators_low_band = int(np.floor(number_of_resonators / 3))
    no_of_resonators_mid_band = int(np.floor(number_of_resonators / 3))
    no_of_resonators_high_band = number_of_resonators - no_of_resonators_mid_band - no_of_resonators_low_band

    low_band_f0s = np.linspace(low_band_min_max_freqs[0], low_band_min_max_freqs[1], no_of_resonators_low_band)
    mid_band_f0s = np.linspace(mid_band_min_max_freqs[0], mid_band_min_max_freqs[1], no_of_resonators_mid_band)
    high_band_f0s = np.linspace(high_band_min_max_freqs[0], high_band_min_max_freqs[1], no_of_resonators_high_band)

    all_f0s = np.concatenate((low_band_f0s, mid_band_f0s, high_band_f0s))

    return all_f0s


def get_resonator_IDCLs_and_CCL_from_f0(
    f0: float,
    resonator_type: AmberResonatorType,
    rounding_precision: int | None = 2,
) -> tuple[NDArray[float64], float | int]:
    """Get the IDC lengths and CC length (both in microns) for a resonator type
    given a resonant frequency in Hz. These lengths are by default rounded to 2
    decimal places.

    Parameters
    ----------
    f0 : float
        The desired frequency (**in Hz**) for the KID.

    resonator_type : AmberResonatorType
        The type of resonator. Only accepts members of Enum ResonatorType.

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
    mux_function = get_mux_func_for_resonator_type(resonator_type)
    IDC_array, CCL = mux_function(f0, rounding_precision=rounding_precision)

    return IDC_array, CCL
