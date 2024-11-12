from ...logging import TextColor, pretty_print
from .. import SoukResonatorType
from ..cpw_coupled_v1 import config_cpw_coupled_v1
from ..high_volume_v1_long_trunk_q20k import config_high_volume_v1_long_trunk_q20k
from ..high_volume_v1_long_trunk_q50k import config_high_volume_v1_long_trunk_q50k
from ..high_volume_v1_q20k import config_high_volume_v1_q20k
from ..high_volume_v1_q50k import config_high_volume_v1_q50k
from ..high_volume_v2_long_trunk_q20k import config_high_volume_v2_long_trunk_q20k
from ..high_volume_v2_long_trunk_q50k import config_high_volume_v2_long_trunk_q50k
from ..high_volume_v2_q20k import config_high_volume_v2_q20k
from ..high_volume_v2_q50k import config_high_volume_v2_q50k
from ..original_long_trunk_q10k import config_original_long_trunk_q10k
from ..original_long_trunk_q20k import config_original_long_trunk_q20k
from ..original_long_trunk_q50k import config_original_long_trunk_q50k
from ..original_q10k import config_original_q10k
from ..original_q20k import config_original_q20k
from ..original_q50k import config_original_q50k


def get_resonator_config(
    resonator_type: SoukResonatorType,
    resonator_config_override: dict[str, float | int] | None = None,
) -> dict[str, float | int]:
    """Get the config for a resonator of a given ResonatorType.

    This is the config used in muxing these resonators. This means only the
    config returned here is guaranteed to generate pixels reliably.

    Parameters
    ----------
    resonator_type : SoukResonatorType
        The type of resonator. Only accepts members of Enum ResonatorType.

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
    config : dict[str, float | int]
        This is a config dictionary with key value pairs as variable name and
        variable values.
    """

    match resonator_type:
        case SoukResonatorType.ORIGINAL_Q10K:
            default_config = config_original_q10k.get_default_resonator_config()
        case SoukResonatorType.ORIGINAL_Q20K:
            default_config = config_original_q20k.get_default_resonator_config()
        case SoukResonatorType.ORIGINAL_Q50K:
            default_config = config_original_q50k.get_default_resonator_config()

        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
            default_config = config_original_long_trunk_q10k.get_default_resonator_config()
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
            default_config = config_original_long_trunk_q20k.get_default_resonator_config()
        case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
            default_config = config_original_long_trunk_q50k.get_default_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
            default_config = config_high_volume_v1_q20k.get_default_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
            default_config = config_high_volume_v1_q50k.get_default_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
            default_config = config_high_volume_v1_long_trunk_q20k.get_default_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
            default_config = config_high_volume_v1_long_trunk_q50k.get_default_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
            default_config = config_high_volume_v2_q20k.get_default_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
            default_config = config_high_volume_v2_q50k.get_default_resonator_config()

        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
            default_config = config_high_volume_v2_long_trunk_q20k.get_default_resonator_config()
        case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
            default_config = config_high_volume_v2_long_trunk_q50k.get_default_resonator_config()

        case SoukResonatorType.CPW_COUPLED_V1:
            default_config = config_cpw_coupled_v1.get_default_resonator_config()

        case _:
            raise (ValueError(f"ResonatorType does not have an associated config."))

    if resonator_config_override is None:
        return default_config

    if not isinstance(resonator_config_override, dict):
        raise TypeError(f"resonator_config_override should be of type dict, not of type {type(resonator_config_override)}")

    config_with_override: dict[str, float | int] = {}

    required_keys = default_config.keys()
    for key, value in resonator_config_override.items():
        if key not in required_keys:
            pretty_print(
                f"Warning: config override contains superfluous key value pair '{key} : {value}'. This will be ignored.",
                color=TextColor.WARNING,
            )

    for key, default_value in default_config.items():
        config_with_override[key] = resonator_config_override.get(key, default_value)

    return config_with_override
