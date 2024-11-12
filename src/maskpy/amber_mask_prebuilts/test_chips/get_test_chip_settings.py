import os
import re
from enum import StrEnum
from typing import Any

import yaml

from ...logging import styled_type_error
from .types import ResonatorSettings, TestChipSettings


class TestChipSettingsPreset(StrEnum):
    """Type of a preset for TestChipSettings.

    Enum Members
    ------------
    NB_FEED_AL_RES_VARYING_LINEWIDTH
    GROUNDPLANE_FEED_VARYING_ORIGINAL_IND4_RES
    """

    NB_FEED_AL_RES_VARYING_LINEWIDTH = "nb_feed_al_res_varying_linewidth"
    GROUNDPLANE_FEED_VARYING_ORIGINAL_IND4_RES = "groundplane_feed_varying_original_ind4_res"

    @classmethod
    def _missing_(cls, value):
        value = value.lower()
        for member in cls:
            if member == value:
                return member
        return None


def _get_non_required_settings(
    settings: dict[str, Any],
    non_required_keys: list[str],
) -> dict[str, Any]:
    """."""
    non_required_settings: dict[str, Any] = {}
    for key in non_required_keys:
        try:
            non_required_settings[key] = settings.pop(key)
        except KeyError as _:
            pass
        except Exception as e:
            raise e
    return non_required_settings


def _get_resonator_settings(resonator_settings: dict[str, Any]) -> ResonatorSettings:
    """."""
    non_required_keys = [
        "mux_func_override",
        "resonator_config_override",
        "mirror",
        "IDC_and_frame_material",
        "meander_material",
        "coupler_fork_material",
        "add_grnd_cutout",
        "add_SiN_dep_dielectric_cutout",
        "add_SiO_cutout",
        "add_SiN_membrane_cutout",
        "add_backside_check",
        "add_inductor_cover",
    ]
    non_required_settings = _get_non_required_settings(resonator_settings, non_required_keys)

    settings = ResonatorSettings(
        resonator_type=resonator_settings["resonator_type"],
        f0=resonator_settings["f0"],
        **non_required_settings,
    )

    return settings


def _get_single_test_chip_settings(chip_settings: dict[str, Any]) -> TestChipSettings:
    """."""
    non_required_keys = [
        "ground_plane_layer",
        "feedline_center_layer",
    ]
    non_required_settings = _get_non_required_settings(chip_settings, non_required_keys)

    settings = TestChipSettings(
        top_left_text=chip_settings["top_left_text"],
        top_right_text=chip_settings["top_right_text"],
        bottom_left_text=chip_settings["bottom_left_text"],
        resonator_1_settings=_get_resonator_settings(chip_settings["resonator_1_settings"]),
        resonator_2_settings=_get_resonator_settings(chip_settings["resonator_2_settings"]),
        resonator_3_settings=_get_resonator_settings(chip_settings["resonator_3_settings"]),
        resonator_4_settings=_get_resonator_settings(chip_settings["resonator_4_settings"]),
        resonator_5_settings=_get_resonator_settings(chip_settings["resonator_5_settings"]),
        resonator_6_settings=_get_resonator_settings(chip_settings["resonator_6_settings"]),
        resonator_7_settings=_get_resonator_settings(chip_settings["resonator_7_settings"]),
        resonator_8_settings=_get_resonator_settings(chip_settings["resonator_8_settings"]),
        resonator_9_settings=_get_resonator_settings(chip_settings["resonator_9_settings"]),
        resonator_10_settings=_get_resonator_settings(chip_settings["resonator_10_settings"]),
        resonator_11_settings=_get_resonator_settings(chip_settings["resonator_11_settings"]),
        resonator_12_settings=_get_resonator_settings(chip_settings["resonator_12_settings"]),
        resonator_13_settings=_get_resonator_settings(chip_settings["resonator_13_settings"]),
        resonator_14_settings=_get_resonator_settings(chip_settings["resonator_14_settings"]),
        resonator_15_settings=_get_resonator_settings(chip_settings["resonator_15_settings"]),
        resonator_16_settings=_get_resonator_settings(chip_settings["resonator_16_settings"]),
        resonator_17_settings=_get_resonator_settings(chip_settings["resonator_17_settings"]),
        resonator_18_settings=_get_resonator_settings(chip_settings["resonator_18_settings"]),
        **non_required_settings,
    )

    return settings


def get_test_chip_settings_from_file(settings_file: str) -> TestChipSettings:
    """Get the test chip settings from a .yaml or a .yml file.
    Parameters
    ----------
    settings_file: str
        The yaml file including extention.

    Returns
    -------
    settings: TestChipSettings
        The settings for the chip in a TestChipSettings TypedDict or a
        list of TestChipSettings for each chip settings in the settings file
        if many are present.
    """

    yaml_exts = ["yaml", "yml"]
    if settings_file.split(".")[-1] not in yaml_exts:
        raise ValueError(f"settings_file should be a yaml file with extentions {yaml_exts}.")

    loader = yaml.SafeLoader
    loader.add_implicit_resolver(
        "tag:yaml.org,2002:float",
        re.compile(
            """^(?:
            [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$""",
            re.X,
        ),
        list("-+0123456789."),
    )

    try:
        with open(settings_file) as file:
            # settings: dict[str, Any] = yaml.safe_load(file)
            settings: dict[str, Any] = yaml.load(file, Loader=loader)
    except FileNotFoundError as fnfe:
        print(f"Unable to find file {str(settings_file)}.")
        raise fnfe
    except Exception as e:
        raise e

    test_chip_settings = _get_single_test_chip_settings(settings)

    return test_chip_settings


def _get_preset_file_path(preset_file_name: str) -> str:
    """Get the full file path for a preset file from the directory where all
    the amber test settings presets files are stored within maskpy.

    Returns
    -------
    full_file_path: str
        The full file path for the preset file.
    """
    preset_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_chip_settings_presets")

    full_file_path = os.path.join(preset_file_path, preset_file_name)

    if not os.path.isfile(full_file_path):
        print(f"preset file '{preset_file_name}' not found.")
        raise FileNotFoundError(f"File not found. Full file and path: {full_file_path}")

    return full_file_path


def get_test_chip_settings_preset(preset: TestChipSettingsPreset) -> TestChipSettings:
    """Get the test chip settings from a preset already put together.
    Parameters
    ----------
    preset:
        The type of preset.

    Returns
    -------
    settings: TestChipSettings
    """

    if not isinstance(preset, TestChipSettingsPreset):
        styled_type_error(preset, "preset", TestChipSettingsPreset)

    preset_file_name = f"{preset.value}.yaml"
    preset_full_file_path = _get_preset_file_path(preset_file_name)

    return get_test_chip_settings_from_file(preset_full_file_path)
