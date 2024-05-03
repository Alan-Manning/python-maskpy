from ..high_volume_v1_long_trunk_q20k import config_high_volume_v1_long_trunk_q20k, utils_high_volume_v1_long_trunk_q20k
from ..high_volume_v1_long_trunk_q50k import config_high_volume_v1_long_trunk_q50k, utils_high_volume_v1_long_trunk_q50k
from ..high_volume_v1_q20k import config_high_volume_v1_q20k, utils_high_volume_v1_q20k
from ..high_volume_v1_q50k import config_high_volume_v1_q50k, utils_high_volume_v1_q50k
from ..high_volume_v2_long_trunk_q20k import config_high_volume_v2_long_trunk_q20k, utils_high_volume_v2_long_trunk_q20k
from ..high_volume_v2_long_trunk_q50k import config_high_volume_v2_long_trunk_q50k, utils_high_volume_v2_long_trunk_q50k
from ..high_volume_v2_q20k import config_high_volume_v2_q20k, utils_high_volume_v2_q20k
from ..high_volume_v2_q50k import config_high_volume_v2_q50k, utils_high_volume_v2_q50k
from ..original_long_trunk_q10k import config_original_long_trunk_q10k, utils_original_long_trunk_q10k
from ..original_long_trunk_q20k import config_original_long_trunk_q20k, utils_original_long_trunk_q20k
from ..original_long_trunk_q50k import config_original_long_trunk_q50k, utils_original_long_trunk_q50k
from ..original_q10k import config_original_q10k, utils_original_q10k
from ..original_q20k import config_original_q20k, utils_original_q20k
from ..original_q50k import config_original_q50k, utils_original_q50k
from ..resonator_types import SoukResonatorType


def get_resonator_config(resonator_type: SoukResonatorType) -> dict:
    """Get the config for a resonator of a given ResonatorType.

    This is the config used in muxing these resonators. This means only the
    config returned here is guaranteed to generate pixels reliably.

    Parameters
    ----------
    resonator_type : ResonatorType
        The type of resonator. Only accepts members of Enum ResonatorType.

    Returns
    -------
    config : dict
        This is a config dictionary with key value pairs as variable name and
        variable values.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            config = config_original_q10k.get_resonator_config()
        case SoukResonatorType.ORIGINAL_Q20K:
            config = config_original_q20k.get_resonator_config()
        case SoukResonatorType.ORIGINAL_Q50K:
            config = config_original_q50k.get_resonator_config()

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            config = config_original_long_trunk_q10k.get_resonator_config()
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            config = config_original_long_trunk_q20k.get_resonator_config()
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            config = config_original_long_trunk_q50k.get_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            config = config_high_volume_v1_q20k.get_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            config = config_high_volume_v1_q50k.get_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            config = config_high_volume_v1_long_trunk_q20k.get_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            config = config_high_volume_v1_long_trunk_q50k.get_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            config = config_high_volume_v2_q20k.get_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            config = config_high_volume_v2_q50k.get_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            config = config_high_volume_v2_long_trunk_q20k.get_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            config = config_high_volume_v2_long_trunk_q50k.get_resonator_config()

        case _:
            raise (ValueError(f"ResonatorType does not have an associated config."))

    return config


def get_total_height_of_resonator(resonator_type: SoukResonatorType, config_override: dict[str, float | int] | None = None) -> float:
    """This will get the total height of the resonator from the base of the
    inductive meander to the end of the ground plane cutout at the top of
    the structure.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    total_resonator_height : float
        The total height of the resonator calculated from the config file.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_total_height_of_resonator(config_override=config_override)

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_total_height_of_resonator(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_total_height_of_resonator(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_total_height_of_resonator(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_total_height_of_resonator(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_total_height_of_resonator(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_total_height_of_resonator(config_override=config_override)

        case _:
            raise (ValueError(f"resonator_type, '{resonator_type}', does not have an associated get_total_height_of_resonator function."))

        # match resonator_type:
        #     case SoukResonatorType.ORIGINAL_Q10K:
        #         return utils_original_q10k.get_total_height_of_resonator(config_override=config_override)
        #     case SoukResonatorType.ORIGINAL_Q20K:
        #         return utils_original_q20k.get_total_height_of_resonator(config_override=config_override)
        #     case SoukResonatorType.ORIGINAL_Q50K:
        #         return utils_original_q50k.get_total_height_of_resonator(config_override=config_override)
        #
        #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
        #         return utils_original_long_trunk_q10k.get_total_height_of_resonator(config_override=config_override)
        #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
        #         return utils_original_long_trunk_q20k.get_total_height_of_resonator(config_override=config_override)
        #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
        #         return utils_original_long_trunk_q50k.get_total_height_of_resonator(config_override=config_override)
        #
        #     case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
        #         return utils_high_volume_v1_q20k.get_total_height_of_resonator(config_override=config_override)
        #     case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
        #         return utils_high_volume_v1_long_trunk_q20k.get_total_height_of_resonator(config_override=config_override)
        #
        #     case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
        #         return utils_high_volume_v2_q20k.get_total_height_of_resonator(config_override=config_override)
        #     case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
        #         return utils_high_volume_v2_long_trunk_q20k.get_total_height_of_resonator(config_override=config_override)
        #     case _:
        #         raise (ValueError(f"resonator_type does not have an associated get_total_height_of_resonator function."))


def get_horizontal_coupler_end_to_meander_base_distance(
    resonator_type: SoukResonatorType, config_override: dict[str, float | int] | None = None
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

    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    coupler_end_to_meander_base_distance : float | int
        The distance between the coupler end and the center of the base of the
        resonator's inductive meander.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)

        case _:
            raise (
                ValueError(
                    f"resonator_type, '{resonator_type}', does not have an associated get_horizontal_coupler_end_to_meander_base_distance function."
                )
            )
    # match resonator_type:
    #     case SoukResonatorType.ORIGINAL_Q10K:
    #         return utils_original_q10k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_Q20K:
    #         return utils_original_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_Q50K:
    #         return utils_original_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
    #         return utils_original_long_trunk_q10k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
    #         return utils_original_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
    #         return utils_original_long_trunk_q50k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #
    #     case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
    #         return utils_high_volume_v1_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
    #         return utils_high_volume_v1_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #
    #     case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
    #         return utils_high_volume_v2_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
    #         return utils_high_volume_v2_long_trunk_q20k.get_horizontal_coupler_end_to_meander_base_distance(config_override=config_override)


def get_vertical_coupler_center_to_meander_base_distance(
    resonator_type: SoukResonatorType, config_override: dict[str, float | int] | None = None
) -> float:
    """This will calculate the the vertical distance between the center of the
    coupler arm and the base of the resonator's inductive meander. This is
    calculated from the config.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    coupler_center_to_meander_base_distance : float | int
        The distance between the coupler center and the base of the resonator's
        inductive meander.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
                config_override=config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(
                config_override=config_override
            )

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
                config_override=config_override
            )
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(
                config_override=config_override
            )

        case _:
            raise (
                ValueError(
                    f"resonator_type, '{resonator_type}', does not have an associated get_vertical_coupler_center_to_meander_base_distance function."
                )
            )

    # match resonator_type:
    #     case SoukResonatorType.ORIGINAL_Q10K:
    #         return utils_original_q10k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_Q20K:
    #         return utils_original_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_Q50K:
    #         return utils_original_q50k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
    #         return utils_original_long_trunk_q10k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
    #         return utils_original_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
    #         return utils_original_long_trunk_q50k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #
    #     case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
    #         return utils_high_volume_v1_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
    #         return utils_high_volume_v1_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
    #             config_override=config_override
    #         )
    #
    #     case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
    #         return utils_high_volume_v2_q20k.get_vertical_coupler_center_to_meander_base_distance(config_override=config_override)
    #     case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
    #         return utils_high_volume_v2_long_trunk_q20k.get_vertical_coupler_center_to_meander_base_distance(
    #             config_override=config_override
    #         )
    #     case _:
    #         raise (ValueError(f"resonator_type does not have an associated get_vertical_coupler_center_to_meander_base_distance function."))


def get_width_and_height_of_IDC_cutout_section(
    resonator_type: SoukResonatorType, config_override: dict[str, float | int] | None = None
) -> tuple[
    float | int,
    float | int,
]:
    """Get the total width and height of ground plane cutout around the IDC
    section of the resonator calculated from the config.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        This is the type of resonators to add bottom choke features for.
        The value accepted here are members of the SoukResonatorType enum.

    config_override : dict | None
        This is an optional override to the base config for this resonator type.
        If nothing is provided the base config for this resonator will be used.
        When provided, this should be a dictionary conating all the keys
        required for this resonator type.

    Returns
    -------
    width, height : tuple[float, float]
        This is a tuple containing, in order, the width and the height
        of the resonator's IDC section calculated from the config file.
    """
    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            return utils_original_q10k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q20K:
            return utils_original_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.ORIGINAL_Q50K:
            return utils_original_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            return utils_original_long_trunk_q10k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            return utils_original_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            return utils_original_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            return utils_high_volume_v1_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            return utils_high_volume_v1_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            return utils_high_volume_v1_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            return utils_high_volume_v1_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            return utils_high_volume_v2_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            return utils_high_volume_v2_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            return utils_high_volume_v2_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            return utils_high_volume_v2_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)

        case _:
            raise (
                ValueError(
                    f"resonator_type, '{resonator_type}', does not have an associated get_width_and_height_of_IDC_cutout_section function."
                )
            )

    # match resonator_type:
    #     case SoukResonatorType.ORIGINAL_Q10K:
    #         return utils_original_q10k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_Q20K:
    #         return utils_original_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_Q50K:
    #         return utils_original_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
    #         return utils_original_long_trunk_q10k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
    #         return utils_original_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
    #         return utils_original_long_trunk_q50k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #
    #     case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
    #         return utils_high_volume_v1_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
    #         return utils_high_volume_v1_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #
    #     case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
    #         return utils_high_volume_v2_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
    #         return utils_high_volume_v2_long_trunk_q20k.get_width_and_height_of_IDC_cutout_section(config_override=config_override)
    #     case _:
    #         raise (ValueError(f"resonator_type does not have an associated get_width_and_height_of_IDC_cutout_section function."))
