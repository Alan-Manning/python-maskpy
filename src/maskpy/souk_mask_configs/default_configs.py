import os
from enum import StrEnum

import yaml

from ..logging import TextColor, pretty_print


class SoukMaskConfig(StrEnum):
    """Type of a Config.

    Enum Members
    ------------
    ANTENNA
    ANTENNA_CPW_MICROSTRIP_TRANS
    BOTTOM_CHOKE
    COMBINER_SECTION_150GHZ
    COMBINER_SECTION_90GHZ
    CPW_FEEDLINE
    FILTER_BANK
    FILTER_BANK_RING_OVERLAP
    GENERAL
    HI_PASS_FILTERS
    LO_PASS_FILTERS
    PORT
    SMA_CONNECTOR
    TEST_CHIP_QUAD
    TOP_CHOKE
    """

    ANTENNA = "antenna"
    ANTENNA_CPW_MICROSTRIP_TRANS = "antenna_cpw_microstrip_trans"
    BOTTOM_CHOKE = "bottom_choke"
    COMBINER_SECTION_150GHZ = "combiner_section_150GHz"
    COMBINER_SECTION_90GHZ = "combiner_section_90GHz"
    CPW_FEEDLINE = "cpw_feedline"
    FILTER_BANK = "filter_bank"
    FILTER_BANK_RING_OVERLAP = "filter_bank_ring_overlap"
    GENERAL = "general"
    HI_PASS_FILTERS = "Hi_pass_filters"
    LO_PASS_FILTERS = "Lo_pass_filters"
    PORT = "port"
    SMA_CONNECTOR = "sma_connector"
    TEST_CHIP_QUAD = "test_chip_quad"
    TOP_CHOKE = "top_choke"

    @classmethod
    def _missing_(cls, value):
        value = value.lower()
        for member in cls:
            if member == value:
                return member
        return None


def _get_default_config(config_type: SoukMaskConfig):
    """Get the default config for the config_type passed in from the yaml file.

    Parameters
    ----------
    config_type: Config
        The configuration for some type that is any member of the
        SoukMaskConfig Enum.

    Returns
    -------
    config: dict[str, float | int]
        This is a config dictionary with the required key value pairs which are
        variable name and variable value pairs.
    """
    file_path = os.path.dirname(os.path.realpath(__file__))
    config_file_path = os.path.join(file_path, "configs", f"{config_type.value}.yaml")

    with open(config_file_path, "r") as config_file:
        config = yaml.safe_load(config_file)

    return config


def get_mask_default_config(config_type: SoukMaskConfig, config_override: dict[str, float | int] | None = None) -> dict[str, float | int]:
    """Get the configuration variables and values with possible overrides for
    some values.

    Parameters
    ----------
    config_type: Config
        The configuration for some type that is any member of the
        SoukMaskConfig Enum.

    KwArgs
    ------
    config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    config: dict[str, float | int]
        This is a config dictionary with the required key value pairs which are
        variable name and variable value pairs.
    """

    if not isinstance(config_type, SoukMaskConfig):
        raise TypeError(f"config_type should be of type SoukMaskConfig, not '{type(config_type)}'.")

    if config_override is None:
        return _get_default_config(config_type)

    if not isinstance(config_override, dict):
        raise (TypeError(f"config_override should be of type dict, not of type {type(config_override)}"))

    default_config = _get_default_config(config_type)
    config_with_override: dict[str, float | int] = {}

    required_keys = default_config.keys()
    for key, value in config_override.items():
        if key not in required_keys:
            pretty_print(
                f"Warning: config_override contains superfluous key value pair '{key} : {value}'. This will be ignored.",
                color=TextColor.WARNING,
            )

    for key, default_value in default_config.items():
        config_with_override[key] = config_override.get(key, default_value)

    return config_with_override
