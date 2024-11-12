from .resonator_types import SoukResonatorType
from .utils.get_config import get_resonator_config
from .utils.souk_resonator_base_utils import (
    get_horizontal_coupler_end_to_meander_base_distance,
    get_total_height_of_resonator,
    get_vertical_coupler_center_to_meander_base_distance,
    get_width_and_height_of_IDC_cutout_section,
)

from .cpw_coupled_v1.cpw_coupled_v1 import CpwCoupledV1
from .original_q50k.original_q50k import OriginalQ50k
from .original_q20k.original_q20k import OriginalQ20k
from .original_q10k.original_q10k import OriginalQ10k

from .original_long_trunk_q50k.original_long_trunk_q50k import OriginalLongTrunkQ50k
from .original_long_trunk_q20k.original_long_trunk_q20k import OriginalLongTrunkQ20k
from .original_long_trunk_q10k.original_long_trunk_q10k import OriginalLongTrunkQ10k

from .high_volume_v1_q50k.high_volume_v1_q50k import HighVolumeV1Q50k
from .high_volume_v1_q20k.high_volume_v1_q20k import HighVolumeV1Q20k
from .high_volume_v1_long_trunk_q50k.high_volume_v1_long_trunk_q50k import HighVolumeV1LongTrunkQ50k
from .high_volume_v1_long_trunk_q20k.high_volume_v1_long_trunk_q20k import HighVolumeV1LongTrunkQ20k

from .high_volume_v2_q50k.high_volume_v2_q50k import HighVolumeV2Q50k
from .high_volume_v2_q20k.high_volume_v2_q20k import HighVolumeV2Q20k
from .high_volume_v2_long_trunk_q50k.high_volume_v2_long_trunk_q50k import HighVolumeV2LongTrunkQ50k
from .high_volume_v2_long_trunk_q20k.high_volume_v2_long_trunk_q20k import HighVolumeV2LongTrunkQ20k


type SoukResonator = CpwCoupledV1 | OriginalQ50k | OriginalQ20k | OriginalQ10k | OriginalLongTrunkQ50k | OriginalLongTrunkQ20k | OriginalLongTrunkQ10k | HighVolumeV1Q50k | HighVolumeV1Q20k | HighVolumeV1LongTrunkQ50k | HighVolumeV1LongTrunkQ20k | HighVolumeV2Q50k | HighVolumeV2Q20k | HighVolumeV2LongTrunkQ50k | HighVolumeV2LongTrunkQ20k
