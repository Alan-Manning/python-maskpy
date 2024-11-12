from typing import Callable

from ...layers import Layer
from ..ABCs import AmberOriginalTypeBase
from ..resonator_types import AmberResonatorType


class OriginalNbInd2(AmberOriginalTypeBase):
    """."""

    def __init__(
        self,
        mask_builder,
        resonator_type: AmberResonatorType,
        x: float,
        y: float,
        rot_angle: float,
        f0: float,
        mux_func_override: Callable | None = None,
        resonator_config_override: dict | None = None,
        mirror=False,
        IDC_and_frame_material: Layer | None = None,
        meander_material: Layer | None = None,
        coupler_fork_material: Layer | None = None,
        add_grnd_cutout=True,
        add_SiN_dep_dielectric_cutout=True,
        add_SiO_cutout=True,
        add_SiN_membrane_cutout=True,
        add_backside_check=True,
        add_inductor_cover=True,
    ):
        """."""
        super().__init__(
            mask_builder,
            resonator_type,
            x,
            y,
            rot_angle,
            f0,
            mux_func_override=mux_func_override,
            resonator_config_override=resonator_config_override,
            mirror=mirror,
            IDC_and_frame_material=IDC_and_frame_material,
            meander_material=meander_material,
            coupler_fork_material=coupler_fork_material,
            add_grnd_cutout=add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=add_SiO_cutout,
            add_SiN_membrane_cutout=add_SiN_membrane_cutout,
            add_backside_check=add_backside_check,
            add_inductor_cover=add_inductor_cover,
        )
        self.draw()
