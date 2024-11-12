import os
from enum import StrEnum
from typing import Any

import yaml

from ...logging import TextColor, pretty_print, styled_text, styled_type_error
from .types import FilterBankStructureSettings, QuadSettings, TestChipSettings


def print_settings(settings: dict[str, Any]) -> None:
    print(f"\n\n\nsettings")
    for key, val in settings.items():
        print(f"{key} : {val}")
    print(f"\n\n")


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


def _get_filter_bank_settings(filter_bank_settings: dict[str, Any] | bool) -> FilterBankStructureSettings | bool:
    """."""
    # print_settings(filter_bank_settings)
    if not isinstance(filter_bank_settings, dict):
        return False

    settings = FilterBankStructureSettings(
        with_combiner=filter_bank_settings["with_combiner"],
        with_crossover=filter_bank_settings["with_crossover"],
        only_1_pol=filter_bank_settings["only_1_pol"],
    )
    return settings


def _get_quad_settings(quad_settings: dict[str, Any]) -> QuadSettings:
    """."""
    # print_settings(quad_settings)
    non_required_keys = [
        "antenna_rotation",
        "mux_func_overrides",
        "trim_lengths",
        "add_grnd_cutout",
        "add_SiN_dep_dielectric_around",
        "add_SiN_dep_dielectric_cutout",
        "add_SiO_cutout",
        "add_SiN_membrane_cutout",
        "add_backside_check",
        "add_grnd_cutout_over_inductor",
        "add_SiN_dep_dielectric_cutout_over_inductor",
        "add_Aluminium_Patch_and_Etch",
        "resonator_config_overrides",
        "general_config_override",
        "antenna_config_override",
        "antenna_cpw_microstrip_trans_config_override",
        "filter_bank_config_override",
        "filter_bank_ring_overlap_config_override",
        "Hi_pass_filters_config_override",
        "Lo_pass_filters_config_override",
        "combiner_section_90ghz_config_override",
        "combiner_section_150ghz_config_override",
        "top_choke_config_override",
        "bottom_choke_config_override",
    ]
    non_required_settings = _get_non_required_settings(quad_settings, non_required_keys)

    settings = QuadSettings(
        resonator_type=quad_settings["resonator_type"],
        add_antenna=quad_settings["add_antenna"],
        # antenna_rotation=quad_settings["antenna_rotation"], now no longer required
        text_under_quad=quad_settings["text_under_quad"],
        couple_KID_to_ANT=quad_settings["couple_KID_to_ANT"],
        add_filter_bank=_get_filter_bank_settings(quad_settings["add_filter_bank"]),
        add_top_choke_features=quad_settings["add_top_choke_features"],
        add_bot_choke_features=quad_settings["add_bot_choke_features"],
        meander_materials=quad_settings["meander_materials"],
        IDC_and_frame_materials=quad_settings["IDC_and_frame_materials"],
        **non_required_settings,
    )

    return settings


def _get_single_test_chip_settings(chip_settings: dict[str, Any]) -> TestChipSettings:
    """."""
    non_required_keys = [
        "time_stamp_position",
        "add_groundplane_under_test_chip",
        "add_SiN_dep_under_test_chip",
        "cpw_feedline_config_override",
        "port_config_override",
        "add_SiN_membrane_under_test_chip",
        "add_SiO_under_test_chip",
    ]
    non_required_settings = _get_non_required_settings(chip_settings, non_required_keys)

    settings = TestChipSettings(
        chip_id=chip_settings["chip_id"],
        top_right_text=chip_settings["top_right_text"],
        top_left_text=chip_settings["top_left_text"],
        bottom_left_text=chip_settings["bottom_left_text"],
        add_ports=chip_settings["add_ports"],
        quad_0_settings=_get_quad_settings(chip_settings["quad_0_settings"]),
        quad_1_settings=_get_quad_settings(chip_settings["quad_1_settings"]),
        quad_2_settings=_get_quad_settings(chip_settings["quad_2_settings"]),
        quad_3_settings=_get_quad_settings(chip_settings["quad_3_settings"]),
        **non_required_settings,
    )

    return settings


def get_test_chip_settings_from_file(settings_file: str) -> TestChipSettings | list[TestChipSettings]:
    """Get the test chip settings from a .yaml or a .yml file.

    Parameters
    ----------
    settings_file: str
        The yaml file including extention.

    Returns
    -------
    settings: TestChipSettings | list[TestChipSettings]
        The settings for the chip in a TestChipSettings TypedDict or a
        list of TestChipSettings for each chip settings in the settings file
        if many are present.
    """

    yaml_exts = ["yaml", "yml"]
    if settings_file.split(".")[-1] not in yaml_exts:
        raise ValueError(
            styled_text("settings_file should be a yaml file with extentions", color=TextColor.ERROR)
            + styled_text(f"`.{yaml_exts[0]}` or `.{yaml_exts[1]}`", color=TextColor.GREEN)
        )

    try:
        with open(settings_file) as file:
            settings: dict[str, Any] = yaml.safe_load(file)
    except FileNotFoundError as fnfe:
        pretty_print(f"Unable to find file:\n{str(settings_file)}.", color=TextColor.ERROR)
        raise fnfe
    except Exception as e:
        raise e

    if settings is None:
        raise ValueError(f"Unable to get settings from file: {settings_file}")

    # Check if there are multiple chips by checking for chip_id key which
    # should only be in a single test_chip_settings
    if "chip_id" in settings.keys():
        test_chip_settings = _get_single_test_chip_settings(settings)
        return test_chip_settings
    else:
        test_chips_settings: list[TestChipSettings] = []
        for _, single_chip_settings in settings.items():
            test_chips_settings.append(_get_single_test_chip_settings(single_chip_settings))
        return test_chips_settings


class TestChipSettingsPreset(StrEnum):
    """Type of a preset for TestChipSettings.

    Enum Members
    ------------
    ORIGINAL_Q50K_AL_OPTICAL_BUILDUP
    ORIGINAL_Q50K_AL_FULL_OPTICAL

    ORIGINAL_Q50K_AND_CPW_COUPLED_AL_BLIND_QUAD_SPLIT
    ORIGINAL_Q50K_AND_CPW_COUPLED_AL_BLIND_QUAD_SAME

    CPW_COUPLED_V1_AL_FULL_COUPLE_OPTICAL

    HIGH_VOLUME_V1_LONG_TRUNK_Q50K_AL_OPTICAL_BUILDUP
    HIGH_VOLUME_V1_LONG_TRUNK_Q50K_AL_FULL_OPTICAL

    HIGH_VOLUME_V2_LONG_TRUNK_Q50K_AL_OPTICAL_BUILDUP
    HIGH_VOLUME_V2_LONG_TRUNK_Q50K_AL_FULL_OPTICAL
    """

    ORIGINAL_Q50K_AL_OPTICAL_BUILDUP = "original_q50k_al_optical_buildup"
    ORIGINAL_Q50K_AL_FULL_OPTICAL = "original_q50k_al_full_optical"

    ORIGINAL_Q50K_AND_CPW_COUPLED_AL_BLIND_QUAD_SPLIT = "original_q50k_and_cpw_coupled_al_blind_quad_split"
    ORIGINAL_Q50K_AND_CPW_COUPLED_AL_BLIND_QUAD_SAME = "original_q50k_and_cpw_coupled_al_blind_quad_same"

    CPW_COUPLED_V1_AL_FULL_COUPLE_OPTICAL = "cpw_coupled_v1_al_full_couple_optical"

    HIGH_VOLUME_V1_LONG_TRUNK_Q50K_AL_OPTICAL_BUILDUP = "high_volume_v1_long_trunk_q50k_al_optical_buildup"
    HIGH_VOLUME_V1_LONG_TRUNK_Q50K_AL_FULL_OPTICAL = "high_volume_v1_long_trunk_q50k_al_full_optical"

    HIGH_VOLUME_V2_LONG_TRUNK_Q50K_AL_OPTICAL_BUILDUP = "high_volume_v2_long_trunk_q50k_al_optical_buildup"
    HIGH_VOLUME_V2_LONG_TRUNK_Q50K_AL_FULL_OPTICAL = "high_volume_v2_long_trunk_q50k_al_full_optical"

    @classmethod
    def _missing_(cls, value):
        value = value.lower()
        for member in cls:
            if member == value:
                return member
        return None


def get_test_chip_settings_preset(preset: TestChipSettingsPreset) -> TestChipSettings:
    """Get the test chip settings from a preset already put together.

    Parameters
    ----------
    preset: TestChipSettingsPreset
        A memeber of the enum TestChipSettingsPreset to specify the preset.

    Returns
    -------
    settings: TestChipSettings
    """

    if not isinstance(preset, TestChipSettingsPreset):
        styled_type_error(preset, "preset", TestChipSettingsPreset)

    full_preset_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_chip_settings_presets", f"{preset.value}.yaml")

    return get_test_chip_settings_from_file(full_preset_file_path)
