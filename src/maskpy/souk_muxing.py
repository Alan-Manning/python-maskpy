from typing import Callable, Union

import numpy as np
import pandas as pd

from .souk_resonators import SoukResonatorType
from .souk_resonators.high_volume_v1_long_trunk_q20k import mux_high_volume_v1_long_trunk_q20k
from .souk_resonators.high_volume_v1_long_trunk_q50k import mux_high_volume_v1_long_trunk_q50k
from .souk_resonators.high_volume_v1_q20k import mux_high_volume_v1_q20k
from .souk_resonators.high_volume_v1_q50k import mux_high_volume_v1_q50k
from .souk_resonators.high_volume_v2_long_trunk_q20k import mux_high_volume_v2_long_trunk_q20k
from .souk_resonators.high_volume_v2_long_trunk_q50k import mux_high_volume_v2_long_trunk_q50k
from .souk_resonators.high_volume_v2_q20k import mux_high_volume_v2_q20k
from .souk_resonators.high_volume_v2_q50k import mux_high_volume_v2_q50k
from .souk_resonators.original_long_trunk_q10k import mux_original_long_trunk_q10k
from .souk_resonators.original_long_trunk_q20k import mux_original_long_trunk_q20k
from .souk_resonators.original_long_trunk_q50k import mux_original_long_trunk_q50k
from .souk_resonators.original_q10k import mux_original_q10k
from .souk_resonators.original_q20k import mux_original_q20k
from .souk_resonators.original_q50k import mux_original_q50k


def get_mux_func_for_resonator_type(resonator_type: SoukResonatorType) -> Callable:
    """Get the muxing function callable that returns the IDCLs and CCL given a
    resonant frequency in Hz."""
    if not isinstance(resonator_type, SoukResonatorType):
        raise TypeError(f"\033[31mresonator_type should be of type SoukResonatorType. Current type is: {type(resonator_type)}\033[0m")

    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return mux_original_q10k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.ORIGINAL_Q20K:
            return mux_original_q20k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.ORIGINAL_Q50K:
            return mux_original_q50k.get_IDCLs_and_CCL_from_f0

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return mux_original_long_trunk_q10k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return mux_original_long_trunk_q20k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return mux_original_long_trunk_q50k.get_IDCLs_and_CCL_from_f0

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return mux_high_volume_v1_q20k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return mux_high_volume_v1_q50k.get_IDCLs_and_CCL_from_f0

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return mux_high_volume_v1_long_trunk_q20k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return mux_high_volume_v1_long_trunk_q50k.get_IDCLs_and_CCL_from_f0

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return mux_high_volume_v2_q20k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return mux_high_volume_v2_q50k.get_IDCLs_and_CCL_from_f0

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return mux_high_volume_v2_long_trunk_q20k.get_IDCLs_and_CCL_from_f0
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return mux_high_volume_v2_long_trunk_q50k.get_IDCLs_and_CCL_from_f0
        case _:
            raise (ValueError(f"ResonatorType does not have an associated muxing function."))


def get_evenly_spaced_freq_array(
    number_of_resonators: int,
    resonator_type: SoukResonatorType,
    min_freq: float = 2000000000.0,
    max_freq: float = 4000000000.0,
    rounding_precision: int | None = 2,
    order: list | np.ndarray | None = None,
    csv_output: bool = False,
    csv_filename: str = "",
    return_dataframe: bool = False,
) -> np.ndarray | tuple[np.ndarray, pd.DataFrame]:
    """This returns a list of resonant frequencies for the resonators in an
    array.

    Parameters
    ----------
    number_of_resonators : int
        This is the number of resonators to generate

    resonator_type : ResonatorType
        The type of resonator. Only accepts members of Enum ResonatorType.

    KwArgs
    ------
    min_freq = 2000000000.0
        The minimum frequency for the array in Hz. Should be a float or int.

    max_freq = 4000000000.0
        The maximum frequency for the arrat in Hz. Should be a float or int.

    order = None
        This is the order for the resonant frequency array. When defined this
        should be a list or numpy array. This is an array of index values in
        order.
        e.g. if order = [0,9,1,8,2,7,4,6,4,5], the f0s returned will be in
        order lowest freq, highest freq, second lowest freq, second highest
        freq, third lowest freq, third highest freq...

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

    return_dataframe = False
        Whether or not to return the pandas DataFrame detailing the data for
        each resonator. See Returns for more.

    Returns
    -------

    f0s : np.array
        The list of resonant frequencies for the array.

    data: pd.DataFrame
        Only returned if return_dataframe parameter is True, Default is false.
        The DataFrame contains each of the resonators details. This data for
        each resonator is resonant frequency, KID_ID, coupler length, and all
        IDC lengths. This is stored in keys: "Freq(Hz)", "KID_ID", "CCL",
        "IDC_ARM_n" (where n = 1..28).
    """

    f0s = np.linspace(min_freq, max_freq, number_of_resonators)

    if csv_output or return_dataframe:
        if isinstance(csv_filename, str) and csv_filename != "":
            filename = csv_filename + ".csv"
        else:
            filename = "csv_mask_mux_output.csv"

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
        data_for_dataframe = np.zeros((number_of_resonators, len(column_names)))

        for i, f0 in enumerate(f0s):
            mux_function = get_mux_func_for_resonator_type(resonator_type)
            IDC_array, CCL = mux_function(f0, rounding_precision=rounding_precision)

            data_for_dataframe[i][0:3] = f0, i, CCL
            data_for_dataframe[i][3:] = IDC_array

        dataframe = pd.DataFrame(data=data_for_dataframe, columns=column_names)

        if csv_output:
            print(f"Making CSV Output with Resonator information as '{filename}'")
            dataframe.to_csv(filename, index=False)

    if isinstance(order, (np.ndarray, list)):
        ordered_f0s = np.zeros(number_of_resonators)

        for i in range(number_of_resonators):
            ordered_f0s[i] = f0s[order[i]]

        if return_dataframe:
            return ordered_f0s, dataframe
        return ordered_f0s

    if return_dataframe:
        return f0s, dataframe
    return f0s


def get_resonator_IDCLs_and_CCL_from_f0(
    f0: float,
    resonator_type: SoukResonatorType,
    rounding_precision: int | None = 2,
) -> tuple[list[float | int], float | int]:
    """Get the IDC lengths and CC length (both in microns) for a resonator type
    given a resonant frequency in Hz. These lengths are by default rounded to 2
    decimal places.

    Parameters
    ----------
    f0 : float
        The desired frequency (**in Hz**) for the KID.

    resonator_type : ResonatorType
        The type of resonator. Only accepts members of Enum ResonatorType.

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
    mux_function = get_mux_func_for_resonator_type(resonator_type)
    IDC_array, CCL = mux_function(f0, rounding_precision=rounding_precision)

    return IDC_array, CCL
