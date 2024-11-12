# from .config_original_q50k import get_resonator_config
from ..resonator_types import AmberResonatorType
from ..utils.get_config import get_resonator_config


def _this_resonator_type() -> AmberResonatorType:
    return AmberResonatorType.ORIGINAL_Q50K
