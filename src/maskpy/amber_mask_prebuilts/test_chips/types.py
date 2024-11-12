from collections.abc import Callable
from typing import NotRequired, Required, TypedDict

from ...amber_resonators import AmberResonatorType


class ResonatorSettings(TypedDict):
    resonator_type: Required[AmberResonatorType]
    f0: Required[float]
    #
    mux_func_override: NotRequired[Callable]
    resonator_config_override: NotRequired[dict[str, float | int]]
    mirror: NotRequired[bool]
    IDC_and_frame_material: NotRequired[str]
    meander_material: NotRequired[str]
    coupler_fork_material: NotRequired[str]
    add_grnd_cutout: NotRequired[bool]
    add_SiN_dep_dielectric_cutout: NotRequired[bool]
    add_SiO_cutout: NotRequired[bool]
    add_SiN_membrane_cutout: NotRequired[bool]
    add_backside_check: NotRequired[bool]
    add_inductor_cover: NotRequired[bool]


class TestChipSettings(TypedDict):
    top_left_text: Required[str]
    top_right_text: Required[str]
    bottom_left_text: Required[str]
    resonator_1_settings: Required[ResonatorSettings]
    resonator_2_settings: Required[ResonatorSettings]
    resonator_3_settings: Required[ResonatorSettings]
    resonator_4_settings: Required[ResonatorSettings]
    resonator_5_settings: Required[ResonatorSettings]
    resonator_6_settings: Required[ResonatorSettings]
    resonator_7_settings: Required[ResonatorSettings]
    resonator_8_settings: Required[ResonatorSettings]
    resonator_9_settings: Required[ResonatorSettings]
    resonator_10_settings: Required[ResonatorSettings]
    resonator_11_settings: Required[ResonatorSettings]
    resonator_12_settings: Required[ResonatorSettings]
    resonator_13_settings: Required[ResonatorSettings]
    resonator_14_settings: Required[ResonatorSettings]
    resonator_15_settings: Required[ResonatorSettings]
    resonator_16_settings: Required[ResonatorSettings]
    resonator_17_settings: Required[ResonatorSettings]
    resonator_18_settings: Required[ResonatorSettings]
    #
    ground_plane_layer: NotRequired[str]
    feedline_center_layer: NotRequired[str]
