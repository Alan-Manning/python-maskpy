from typing import Callable

from ..ABCs import SoukOriginalTypeBase
from ..resonator_types import SoukResonatorType


class OriginalQ20k(SoukOriginalTypeBase):
    """."""

    def __init__(
        self,
        mask_builder,
        resonator_type: SoukResonatorType,
        x: float | int,
        y: float | int,
        rot_angle: float | int,
        f0: float | int,
        mux_func_override: Callable | None = None,
        resonator_config_override: dict[str, float | int] | None = None,
        mirror: bool = False,
        IDC_and_frame_material: str = "IDC_Nb",
        meander_material: str = "Al",
        trim_length: float | int | None = None,
        add_grnd_cutout: bool = True,
        add_SiN_dep_dielectric_cutout: bool = True,
        add_SiO_cutout: bool = True,
        add_SiN_membrane_cutout: bool = True,
        add_backside_check: bool = False,
        add_grnd_cutout_over_inductor: bool = False,
        add_SiN_dep_dielectric_cutout_over_inductor: bool = False,
        add_Aluminium_Patch_and_Etch: bool = True,
        return_configurator_points: bool = False,
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
            trim_length=trim_length,
            add_grnd_cutout=add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=add_SiO_cutout,
            add_SiN_membrane_cutout=add_SiN_membrane_cutout,
            add_backside_check=add_backside_check,
            add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
            add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
            add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
            return_configurator_points=return_configurator_points,
        )
        self.draw()
