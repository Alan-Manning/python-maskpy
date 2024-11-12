from collections.abc import Sequence
from typing import Any

import gdspy

from maskpy.souk_resonators.resonator_types import SoukResonatorType

from ...souk_mask_builder import SoukMaskBuilder

# type OptConfigOverride = dict[str, float | int] | None


def _add_general_labeling_square_around_horn_center(
    mask_builder: SoukMaskBuilder,
    horn_number: int,
    horn_center_x: float,
    horn_center_y: float,
) -> None:
    """Adding a small square showing pixel centers to general labeling."""
    mask_builder.Main.add(
        gdspy.Rectangle(
            [horn_center_x, horn_center_y],
            [horn_center_x + 10, horn_center_y + 10],
            layer=mask_builder.layers.General_labeling.number,
            datatype=mask_builder.layers.General_labeling.datatype,
        )
    )
    mask_builder.Main.add(
        gdspy.Rectangle(
            [horn_center_x, horn_center_y],
            [horn_center_x - 10, horn_center_y - 10],
            layer=mask_builder.layers.General_labeling.number,
            datatype=mask_builder.layers.General_labeling.datatype,
        )
    )

    mask_builder.add_fancy_text(
        str(horn_number),
        horn_center_x,
        horn_center_y,
        80,
        mask_builder.layers.General_labeling,
        horizontal_align="end",
        vertical_align="above",
    )


def add_full_horn_and_resonators(
    mask_builder: SoukMaskBuilder,
    horn_number: int,
    horn_center_x: float,
    horn_center_y: float,
    antenna_rotation: float,
    KID_numbers: Sequence[int],
    KID_freqs: Sequence[float],
    resonator_types: Sequence[SoukResonatorType],
    # Config Overrides.
    general_config_override: dict[str, float | int] | None = None,
    resonator_config_overrides: Sequence[dict[str, float | int] | None] = (None, None, None, None),
    antenna_config_override: dict[str, float | int] | None = None,
    cpw_feedline_config_override: dict[str, float | int] | None = None,
    antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None,
    top_choke_config_override: dict[str, float | int] | None = None,
    bottom_choke_config_override: dict[str, float | int] | None = None,
    # KwArgs
    add_antenna_kwargs: dict[str, Any] = {},
    add_4_resonators_around_horn_kwargs: dict[str, Any] = {},
    add_filter_bank_kwargs: dict[str, Any] = {},
) -> list[list[float]]:
    """Add a Full horn and resonators around the horn connecting then to the
    filter bank where needed. This finally returns the coordinates at the end
    of each the coupling capacitor where the feedline will connect with the
    resonaotrs. This is a list of lists containing the [x,y] for the TopLeft,
    TopRight, BotLeft, BotRight resonators respectively.

    This is a partial function which combines a slew of SoukMaskBuilder methods
    to generate the horn block. This aims to simplify mask making code
    abstracting some of the underlying setup away from the user. If more fine
    control over what is drawn is required then a lot of this code can be taken
    and used outside.

    The structure of how this code works is:
    >>> Adds the antenna.
    >>> Adds a small labeling square for the center of the horn.
    >>> Adds the 4 resonators around the horn block.
    >>> Gets the feedline_pass_through_points for those resonators.
    >>> Add the filter bank and combiners.
    >>> Connect antennas to the filter bank and the filter bank to the kids.
    >>> Add the top and bottom choke features.
    >>> return the feedline_pass_through_points.

    The SoukMaskBuilder methods called in this function are in order:
    >>> mask_builder.get_relative_antenna_conect_positions
    >>> mask_builder.get_relative_kid_positions
    >>> mask_builder.add_antenna
    >>> mask_builder.add_4_resonators_around_horn
    >>> mask_builder.get_feedline_pass_through_points
    >>> mask_builder.add_filter_bank_and_get_conection_points
    >>> mask_builder.connect_filter_bank_to_KIDs
    >>> mask_builder.connect_ants_to_filter_bank
    >>> mask_builder.add_top_choke_features
    >>> mask_builder.add_bottom_choke_features

    Parameters
    ----------
    mask_builder: SoukMaskBuilder
        The mask builder to add the full horn block and resonators to.

    horn_number: int,
        The number of this horn block, this is only really used for labeling.

    horn_center_x: float
        The center x coord for the full horn block and resonators assembly.

    horn_center_y: float
        The center y coord for the full horn block and resonators assembly.

    antenna_rotation: float,
        The angle (**in degrees**) the antenna structure should be rotated.

    KID_numbers: tuple[int, int, int, int] | list[int]
        This is a list or tuple of 4 integers which are the numbers each KID
        should have drawn next to it. This list of numbers for each KID should
        be in the order TopLeft, TopRigh, BotLeft, BotRight.

    KID_freqs: tuple[float, float, float, float] | list[float]
        The resonant frequencies of the resonators. Should be in the same unit
        that the mux_funcs function takes (default mux_funcs are in Hz).

    resonator_types: tuple[SoukResonatorType, SoukResonatorType, SoukResonatorType, SoukResonatorType] | list[SoukResonatorType],
        This is the type of resonators to be drawn. The values accepted
        here are members of the SoukResonatorType enum.
        The order of the values passed in will be attributed to each KID
        and should be in the order TopLeft, TopRigh, BotLeft, BotRight.

    Returns
    -------
    feedline_pass_through_points: list[list[float]]
        The points where the feedline should pass through to connect to the
        resonators coupling capacitors.
    """

    ###########################################################################
    # Getting the end points of the antenna pads where they will connect to the filter bank.
    relative_antena_conect_positions = mask_builder.get_relative_antenna_conect_positions(
        antenna_rotation,
        antenna_config_override=antenna_config_override,
    )

    ###########################################################################
    # Getting the relative KID postions for an antenna block.
    relative_kid_positions = mask_builder.get_relative_kid_positions(
        resonator_types,
        resonator_config_overrides=resonator_config_overrides,
        general_config_override=general_config_override,
    )

    ###########################################################################
    # adds the antennas

    antenna_config_override = add_antenna_kwargs.pop("antenna_config_override", None)
    add_grnd_cutout = add_antenna_kwargs.pop("add_grnd_cutout", True)
    add_SiN_dep_cutout = add_antenna_kwargs.pop("add_SiN_dep_cutout", True)
    add_backside_check = add_antenna_kwargs.pop("add_backside_check", True)

    mask_builder.add_antenna(
        horn_center_x,
        horn_center_y,
        antenna_rotation,
        add_grnd_cutout=add_grnd_cutout,
        add_SiN_dep_cutout=add_SiN_dep_cutout,
        add_backside_check=add_backside_check,
    )

    _add_general_labeling_square_around_horn_center(
        mask_builder,
        horn_number,
        horn_center_x,
        horn_center_y,
    )

    ###########################################################################
    # Adds the 4 KIDs to the mask and also adds the number for each KID

    mux_func_overrides = add_4_resonators_around_horn_kwargs.pop(
        "mux_func_overrides",
        [None, None, None, None],
    )
    IDC_and_frame_materials = add_4_resonators_around_horn_kwargs.pop(
        "IDC_and_frame_materials",
        ["IDC_Nb", "IDC_Nb", "IDC_Nb", "IDC_Nb"],
    )
    meander_materials = add_4_resonators_around_horn_kwargs.pop(
        "meander_materials",
        ["Al", "Al", "Al", "Al"],
    )
    trim_lengths = add_4_resonators_around_horn_kwargs.pop(
        "trim_lengths",
        [None, None, None, None],
    )
    add_grnd_cutout = add_4_resonators_around_horn_kwargs.pop(
        "add_grnd_cutout",
        [True, True, True, True],
    )
    add_SiN_dep_dielectric_around = add_4_resonators_around_horn_kwargs.pop(
        "add_SiN_dep_dielectric_around",
        True,
    )
    add_SiN_dep_dielectric_cutout = add_4_resonators_around_horn_kwargs.pop(
        "add_SiN_dep_dielectric_cutout",
        [True, True, True, True],
    )
    add_SiO_cutout = add_4_resonators_around_horn_kwargs.pop(
        "add_SiO_cutout",
        [True, True, True, True],
    )
    add_SiN_membrane_cutout = add_4_resonators_around_horn_kwargs.pop(
        "add_SiN_membrane_cutout",
        [True, True, True, True],
    )
    add_backside_check = add_4_resonators_around_horn_kwargs.pop(
        "add_backside_check",
        [False, False, False, False],
    )
    add_grnd_cutout_over_inductor = add_4_resonators_around_horn_kwargs.pop(
        "add_grnd_cutout_over_inductor",
        [False, False, False, False],
    )
    add_SiN_dep_dielectric_cutout_over_inductor = add_4_resonators_around_horn_kwargs.pop(
        "add_SiN_dep_dielectric_cutout_over_inductor",
        [False, False, False, False],
    )
    add_Aluminium_Patch_and_Etch = add_4_resonators_around_horn_kwargs.pop(
        "add_Aluminium_Patch_and_Etch",
        [False, False, False, False],
    )

    mask_builder.add_4_resonators_around_horn(
        horn_center_x,
        horn_center_y,
        relative_kid_positions,
        KID_numbers,
        KID_freqs,
        resonator_types,
        resonator_config_overrides=resonator_config_overrides,
        mux_func_overrides=mux_func_overrides,
        general_config_override=general_config_override,
        IDC_and_frame_materials=IDC_and_frame_materials,
        meander_materials=meander_materials,
        trim_lengths=trim_lengths,
        add_grnd_cutout=add_grnd_cutout,
        add_SiN_dep_dielectric_around=add_SiN_dep_dielectric_around,
        add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
        add_SiO_cutout=add_SiO_cutout,
        add_SiN_membrane_cutout=add_SiN_membrane_cutout,
        add_backside_check=add_backside_check,
        add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
        add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
        add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
    )

    ###########################################################################
    # Gets the coords of the end of the kids where they conect to the feedline and adds them to list
    feedline_pass_through_points = mask_builder.get_feedline_pass_through_points(
        horn_center_x,
        horn_center_y,
        relative_kid_positions,
        resonator_types,
        cpw_feedline_config_override=cpw_feedline_config_override,
        resonator_config_overrides=resonator_config_overrides,
    )

    # TODO return feedline_pass_through_points

    ###########################################################################
    # Adding the filter bank

    filter_bank_config_override = add_filter_bank_kwargs.pop("filter_bank_config_override", None)
    filter_bank_ring_overlap_config_override = add_filter_bank_kwargs.pop("filter_bank_ring_overlap_config_override", None)
    Hi_pass_filters_config_override = add_filter_bank_kwargs.pop("Hi_pass_filters_config_override", None)
    Lo_pass_filters_config_override = add_filter_bank_kwargs.pop("Lo_pass_filters_config_override", None)
    combiner_section_90ghz_config_override = add_filter_bank_kwargs.pop("combiner_section_90ghz_config_override", None)
    combiner_section_150ghz_config_override = add_filter_bank_kwargs.pop("combiner_section_150ghz_config_override", None)
    with_combiner = add_filter_bank_kwargs.pop("with_combiner", True)
    with_crossover = add_filter_bank_kwargs.pop("with_crossover", True)
    only_1_pol = add_filter_bank_kwargs.pop("only_1_pol", False)

    absolute_filter_bank_conect_points = mask_builder.add_filter_bank_and_get_conection_points(
        horn_center_x,
        horn_center_y,
        filter_bank_config_override=filter_bank_config_override,
        filter_bank_ring_overlap_config_override=filter_bank_ring_overlap_config_override,
        Hi_pass_filters_config_override=Hi_pass_filters_config_override,
        Lo_pass_filters_config_override=Lo_pass_filters_config_override,
        combiner_section_90ghz_config_override=combiner_section_90ghz_config_override,
        combiner_section_150ghz_config_override=combiner_section_150ghz_config_override,
        with_combiner=with_combiner,
        with_crossover=with_crossover,
        only_1_pol=only_1_pol,
    )  # add filter bank centered on x, y given and returns the points where to conect the KIDs to filter bank.

    only_1_pol_no_comb = only_1_pol and (not with_combiner)

    # Connecting the filter bank to the KIDs inductive meanders.
    mask_builder.connect_filter_bank_to_KIDs(
        horn_center_x,
        horn_center_y,
        absolute_filter_bank_conect_points,
        relative_kid_positions,
        only_1_pol_no_comb=only_1_pol_no_comb,
    )

    if only_1_pol_no_comb:
        terminate_ants = ["L", "R"]
    else:
        terminate_ants = []

    # Connecting the end of the antenna pads to the filter bank.
    mask_builder.connect_ants_to_filter_bank(
        horn_center_x,
        horn_center_y,
        relative_antena_conect_positions,
        antenna_rotation,
        filter_bank_config_override=filter_bank_config_override,
        antenna_config_override=antenna_config_override,
        antenna_cpw_microstrip_trans_config_override=antenna_cpw_microstrip_trans_config_override,
        terminate_ants=terminate_ants,
    )

    # Adding the top and bottom choke features.
    mask_builder.add_top_choke_features(
        horn_center_x,
        horn_center_y,
        top_choke_config_override=top_choke_config_override,
    )

    mask_builder.add_bottom_choke_features(
        horn_center_x,
        horn_center_y,
        relative_kid_positions,
        resonator_types,
        bottom_choke_config_override=bottom_choke_config_override,
        resonator_config_overrides=resonator_config_overrides,
    )

    return feedline_pass_through_points
