from __future__ import annotations

from typing import Callable, Literal, NotRequired, Required, TypedDict

# from ...souk_resonators import SoukResonatorType
from maskpy.souk_resonators import SoukResonatorType


class FilterBankStructureSettings(TypedDict):
    """FilterBankStructureSettings.

    Required
    --------
    with_combiner: Required[bool]
    with_crossover: Required[bool]
    only_1_pol: Required[bool]
    """

    with_combiner: Required[bool]
    with_crossover: Required[bool]
    only_1_pol: Required[bool]


class QuadSettings(TypedDict):
    """QuadSettings.

    Required
    --------
    resonator_type: Required[SoukResonatorType]
    add_antenna: Required[bool]
    text_under_quad: Required[str]
    couple_KID_to_ANT: Required[bool]
    add_filter_bank: Required[FilterBankStructureSettings | Literal[False]]
    add_top_choke_features: Required[bool]
    add_bot_choke_features: Required[bool]
    meander_materials: Required[tuple[str, str, str, str]]
    IDC_and_frame_materials: Required[tuple[str, str, str, str]]

    Not Required
    ------------
    antenna_rotation: NotRequired[float]
    mux_func_overrides: NotRequired[tuple[Callable | None, Callable | None, Callable | None, Callable | None]]
    trim_lengths: NotRequired[tuple[float | None, float | None, float | None, float | None]]
    add_grnd_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiN_dep_dielectric_around: NotRequired[bool]
    add_SiN_dep_dielectric_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiO_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiN_membrane_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_backside_check: NotRequired[tuple[bool, bool, bool, bool]]
    add_grnd_cutout_over_inductor: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiN_dep_dielectric_cutout_over_inductor: NotRequired[tuple[bool, bool, bool, bool]]
    add_Aluminium_Patch_and_Etch: NotRequired[tuple[bool, bool, bool, bool]]

    # configs
    resonator_config_overrides: NotRequired[tuple[dict[str, float | int] | None, dict[str, float | int] | None, dict[str, float | int] | None, dict[str, float | int] | None]]
    general_config_override: NotRequired[dict[str, float | int] | None]
    antenna_config_override: NotRequired[dict[str, float | int] | None]
    antenna_cpw_microstrip_trans_config_override: NotRequired[dict[str, float | int] | None]
    filter_bank_config_override: NotRequired[dict[str, float | int] | None]
    filter_bank_ring_overlap_config_override: NotRequired[dict[str, float | int] | None]
    Hi_pass_filters_config_override: NotRequired[dict[str, float | int] | None]
    Lo_pass_filters_config_override: NotRequired[dict[str, float | int] | None]
    combiner_section_90ghz_config_override: NotRequired[dict[str, float | int] | None]
    combiner_section_150ghz_config_override: NotRequired[dict[str, float | int] | None]
    top_choke_config_override: NotRequired[dict[str, float | int] | None]
    bottom_choke_config_override: NotRequired[dict[str, float | int] | None]
    """

    # Required
    resonator_type: Required[SoukResonatorType]
    add_antenna: Required[bool]
    text_under_quad: Required[str]
    couple_KID_to_ANT: Required[bool]
    add_filter_bank: Required[FilterBankStructureSettings | Literal[False]]
    add_top_choke_features: Required[bool]
    add_bot_choke_features: Required[bool]
    meander_materials: Required[tuple[str, str, str, str]]
    IDC_and_frame_materials: Required[tuple[str, str, str, str]]
    # NotRequired
    antenna_rotation: NotRequired[float]
    mux_func_overrides: NotRequired[tuple[Callable | None, Callable | None, Callable | None, Callable | None]]
    trim_lengths: NotRequired[tuple[float | None, float | None, float | None, float | None]]
    add_grnd_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiN_dep_dielectric_around: NotRequired[bool]
    add_SiN_dep_dielectric_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiO_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiN_membrane_cutout: NotRequired[tuple[bool, bool, bool, bool]]
    add_backside_check: NotRequired[tuple[bool, bool, bool, bool]]
    add_grnd_cutout_over_inductor: NotRequired[tuple[bool, bool, bool, bool]]
    add_SiN_dep_dielectric_cutout_over_inductor: NotRequired[tuple[bool, bool, bool, bool]]
    add_Aluminium_Patch_and_Etch: NotRequired[tuple[bool, bool, bool, bool]]
    # configs
    resonator_config_overrides: NotRequired[
        tuple[dict[str, float | int] | None, dict[str, float | int] | None, dict[str, float | int] | None, dict[str, float | int] | None]
    ]
    general_config_override: NotRequired[dict[str, float | int] | None]
    antenna_config_override: NotRequired[dict[str, float | int] | None]
    antenna_cpw_microstrip_trans_config_override: NotRequired[dict[str, float | int] | None]
    filter_bank_config_override: NotRequired[dict[str, float | int] | None]
    filter_bank_ring_overlap_config_override: NotRequired[dict[str, float | int] | None]
    Hi_pass_filters_config_override: NotRequired[dict[str, float | int] | None]
    Lo_pass_filters_config_override: NotRequired[dict[str, float | int] | None]
    combiner_section_90ghz_config_override: NotRequired[dict[str, float | int] | None]
    combiner_section_150ghz_config_override: NotRequired[dict[str, float | int] | None]
    top_choke_config_override: NotRequired[dict[str, float | int] | None]
    bottom_choke_config_override: NotRequired[dict[str, float | int] | None]


class TestChipSettings(TypedDict):
    """TestChipSettings.

    Required
    --------
    chip_id: Required[int]
    top_right_text: Required[str]
    top_left_text: Required[str]
    bottom_left_text: Required[str]

    add_ports: Required[bool]
    quad_0_settings: Required[QuadSettings]
    quad_1_settings: Required[QuadSettings]
    quad_2_settings: Required[QuadSettings]
    quad_3_settings: Required[QuadSettings]

    Not Required
    ------------
    time_stamp_position: NotRequired[str]
    add_groundplane_under_test_chip: NotRequired[bool]
    add_SiN_dep_under_test_chip: NotRequired[bool]
    SiN_dep_under_test_chip_edge_offset: NotRequired[float | int]
    add_SiN_membrane_under_test_chip: NotRequired[bool]
    add_SiO_under_test_chip: NotRequired[bool]

    # configs
    cpw_feedline_config_override: NotRequired[dict[str, float | int] | None]
    port_config_override: NotRequired[dict[str, float | int] | None]
    """

    chip_id: Required[int]
    top_right_text: Required[str]
    top_left_text: Required[str]
    bottom_left_text: Required[str]
    time_stamp_position: NotRequired[str]
    add_groundplane_under_test_chip: NotRequired[bool]
    add_SiN_dep_under_test_chip: NotRequired[bool]
    SiN_dep_under_test_chip_edge_offset: NotRequired[float | int]
    add_SiN_membrane_under_test_chip: NotRequired[bool]
    add_SiO_under_test_chip: NotRequired[bool]
    add_ports: Required[bool]
    quad_0_settings: Required[QuadSettings]
    quad_1_settings: Required[QuadSettings]
    quad_2_settings: Required[QuadSettings]
    quad_3_settings: Required[QuadSettings]
    # configs
    cpw_feedline_config_override: NotRequired[dict[str, float | int] | None]
    port_config_override: NotRequired[dict[str, float | int] | None]
