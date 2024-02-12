import numpy as np
import pandas as pd
from mux_funcs import Q10k, Q50k


def get_evenly_spaced_freq_array(
    number_of_resonators: int,
    Q_value: str,
    min_freq: float = 2000000000.0,
    max_freq: float = 4000000000.0,
    order: list | None = None,
    csv_output: bool = False,
    csv_filename: str = "",
) -> np.ndarray:
    """This returns a list of resonant frequencies for the resonators in an
    array.

    Parameters
    ----------
    number_of_resonators : int
        This is the number of resonators to generate

    Q_value : str
        The desired Q value for the resonators. Only accepts values "10k" or "50k".

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

    Returns
    -------

    f0s : np.array
        The list of resonant frequencies for the array.
    """

    acceptable_Q_value_inputs = ["10k", "50k"]

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
            match Q_value:
                case "10k":
                    IDC_array, CCL = Q10k.get_IDCLs_and_CCL_from_f0(f0, rounding_precision=rounding_precision)
                case "50k":
                    IDC_array, CCL = Q50k.get_IDCLs_and_CCL_from_f0(f0, rounding_precision=rounding_precision)
                case _:
                    raise (ValueError(f"Q_value argument is not acceptable. Valid inputs are: {acceptable_Q_value_inputs}"))
            data_for_datafrane[i][0:3] = f0, i, CCL
            data_for_datafrane[i][3:] = IDC_array

        datafrane = pd.DataFrame(data=data_for_datafrane, columns=column_names)
        datafrane.to_csv(filename, index=False)

    if isinstance(order, (np.ndarray, list)):
        ordered_f0s = np.zeros(number_of_resonators)

        for i in range(number_of_resonators):
            ordered_f0s[i] = f0s[order[i]]

        return ordered_f0s

    return f0s


def get_resonator_IDCLs_and_CCL_from_f0(f0: float, Q_value: str, rounding_precision: int = 2):
    """Get the IDC lengths and CC length for a resonator given a resonant
    frequency in Hz. These lengths are by default rounded to 2 decimal places.

    Parameters
    ----------
    f0 : float
        The desired frequency (**in Hz**) for the KID.

    Q_value : str
        The desired Q value for the resonator. Only accepts values "10k" or "50k".

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
    acceptable_Q_value_inputs = ["10k", "50k"]

    match Q_value:
        case "10k":
            IDC_array, CCL = Q10k.get_IDCLs_and_CCL_from_f0(f0, rounding_precision=rounding_precision)
        case "50k":
            IDC_array, CCL = Q50k.get_IDCLs_and_CCL_from_f0(f0, rounding_precision=rounding_precision)
        case _:
            raise (ValueError(f"Q_value argument is not acceptable. Valid inputs are: {acceptable_Q_value_inputs}"))

    return IDC_array, CCL
