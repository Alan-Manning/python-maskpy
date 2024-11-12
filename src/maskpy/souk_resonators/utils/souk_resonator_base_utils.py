from .. import SoukResonatorType
from ..cpw_coupled_v1 import utils_cpw_coupled_v1
from ..high_volume_v1_long_trunk_q20k import utils_high_volume_v1_long_trunk_q20k
from ..high_volume_v1_long_trunk_q50k import utils_high_volume_v1_long_trunk_q50k
from ..high_volume_v1_q20k import utils_high_volume_v1_q20k
from ..high_volume_v1_q50k import utils_high_volume_v1_q50k
from ..high_volume_v2_long_trunk_q20k import utils_high_volume_v2_long_trunk_q20k
from ..high_volume_v2_long_trunk_q50k import utils_high_volume_v2_long_trunk_q50k
from ..high_volume_v2_q20k import utils_high_volume_v2_q20k
from ..high_volume_v2_q50k import utils_high_volume_v2_q50k
from ..original_long_trunk_q10k import utils_original_long_trunk_q10k
from ..original_long_trunk_q20k import utils_original_long_trunk_q20k
from ..original_long_trunk_q50k import utils_original_long_trunk_q50k
from ..original_q10k import utils_original_q10k
from ..original_q20k import utils_original_q20k
from ..original_q50k import utils_original_q50k


def get_total_height_of_resonator(
    resonator_type: SoukResonatorType,
    resonator_config_override: dict[str, float | int] | None = None,
) -> float:
    """This will get the total height of the resonator from the base of the
    inductive meander to the end of the ground plane cutout at the top of
    the structure.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    total_resonator_height : float
        The total height of the resonator calculated from the config file.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case SoukResonatorType.CPW_COUPLED_V1:
            return utils_cpw_coupled_v1.get_total_height_of_resonator(resonator_config_override=resonator_config_override)

        case _:
            raise (ValueError(f"resonator_type, '{resonator_type}', does not have an associated get_total_height_of_resonator function."))


def get_horizontal_coupler_end_to_meander_base_distance(
    resonator_type: SoukResonatorType,
    resonator_config_override: dict[str, float | int] | None = None,
) -> float:
    """This will calculate the the horizonatal distance between the end of the
    coupler arm (where it would connect to a feedline) and the center of the
    base of the resonator's inductive meander. This is calculated from the
    config.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    coupler_end_to_meander_base_distance : float | int
        The distance between the coupler end and the center of the base of the
        resonator's inductive meander.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.CPW_COUPLED_V1:
            return utils_cpw_coupled_v1.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case _:
            raise (
                ValueError(
                    f"resonator_type, '{resonator_type}', does not have an associated get_horizontal_coupler_end_to_meander_base_distance function."
                )
            )


def get_vertical_coupler_center_to_meander_base_distance(
    resonator_type: SoukResonatorType,
    resonator_config_override: dict[str, float | int] | None = None,
) -> float:
    """This will calculate the the vertical distance between the center of the
    coupler arm and the base of the resonator's inductive meander. This is
    calculated from the config.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    coupler_center_to_meander_base_distance : float | int
        The distance between the coupler center and the base of the resonator's
        inductive meander.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.CPW_COUPLED_V1:
            return utils_cpw_coupled_v1.get_vertical_coupler_center_to_meander_base_distance(
                resonator_config_override=resonator_config_override
            )

        case _:
            raise (
                ValueError(
                    f"resonator_type, '{resonator_type}', does not have an associated get_vertical_coupler_center_to_meander_base_distance function."
                )
            )


def get_width_and_height_of_IDC_cutout_section(
    resonator_type: SoukResonatorType,
    resonator_config_override: dict[str, float | int] | None = None,
) -> tuple[float | int, float | int]:
    """Get the total width and height of ground plane cutout around the IDC
    section of the resonator calculated from the config.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    KwArgs
    ------
    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    Returns
    -------
    width, height : tuple[float, float]
        This is a tuple containing, in order, the width and the height
        of the resonator's IDC section calculated from the config file.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(
                resonator_config_override=resonator_config_override
            )

        case SoukResonatorType.CPW_COUPLED_V1:
            return utils_cpw_coupled_v1.get_width_and_height_of_IDC_cutout_section(resonator_config_override=resonator_config_override)

        case _:
            raise (
                ValueError(
                    f"resonator_type, '{resonator_type}', does not have an associated get_width_and_height_of_IDC_cutout_section function."
                )
            )
