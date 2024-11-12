from abc import ABC
from typing import Callable, TypedDict

from ..logging import styled_type_error
from ..souk_muxing import get_mux_func_for_resonator_type
from .draw_funcs import draw_cpw_coupled_v1_type, draw_high_volume_v1_type, draw_high_volume_v2_type, draw_souk_original_type
from .resonator_types import SoukResonatorType
from .utils.get_config import get_resonator_config


class Details(TypedDict):
    KID_type: str
    KID_No: int | None
    x_coord: float | int
    y_coord: float | int
    rot: float | int
    f0: float | int
    mux_override: None | str
    mux_IDC: list[float | int]
    mux_CC: float | int
    trim: float | int | None
    config_override: None | dict[str, float | int]


class SoukOriginalTypeBase(ABC):
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

        # mask_builder
        self.mask_builder = mask_builder

        # resonator_type
        if not isinstance(resonator_type, SoukResonatorType):
            styled_type_error(resonator_type, "resonator_type", SoukResonatorType)

        accepted_resonator_types = [
            SoukResonatorType.ORIGINAL_Q10K,
            SoukResonatorType.ORIGINAL_Q20K,
            SoukResonatorType.ORIGINAL_Q50K,
            SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K,
            SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K,
            SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K,
        ]

        if resonator_type not in accepted_resonator_types:
            raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

        self.resonator_type = resonator_type

        self.x = x
        self.y = y
        self.rot_angle = rot_angle

        self.f0 = f0
        self.mux_func_override = mux_func_override

        self.resonator_config_override = resonator_config_override

        self.mirror = mirror
        self.IDC_and_frame_material = IDC_and_frame_material
        self.meander_material = meander_material
        self.trim_length = trim_length
        self.add_grnd_cutout = add_grnd_cutout
        self.add_SiN_dep_dielectric_cutout = add_SiN_dep_dielectric_cutout
        self.add_SiO_cutout = add_SiO_cutout
        self.add_SiN_membrane_cutout = add_SiN_membrane_cutout
        self.add_backside_check = add_backside_check
        self.add_grnd_cutout_over_inductor = add_grnd_cutout_over_inductor
        self.add_SiN_dep_dielectric_cutout_over_inductor = add_SiN_dep_dielectric_cutout_over_inductor
        self.add_Aluminium_Patch_and_Etch = add_Aluminium_Patch_and_Etch
        self.return_configurator_points = return_configurator_points

    def get_details(self) -> Details:
        """Get the details for the resonator."""

        if self.mux_func_override is None:
            IDC_and_CC_function = get_mux_func_for_resonator_type(self.resonator_type)
            mux_override = None
        else:
            IDC_and_CC_function = self.mux_func_override
            mux_override = getattr(IDC_and_CC_function, "__name__", str(IDC_and_CC_function))

        idcls, ccl = IDC_and_CC_function(self.f0)

        details: Details = {
            "KID_type": str(self.resonator_type),
            "KID_No": None,
            "x_coord": self.x,
            "y_coord": self.y,
            "rot": self.rot_angle,
            "f0": self.f0,
            "mux_override": mux_override,
            "mux_IDC": idcls,
            "mux_CC": ccl,
            "trim": self.trim_length,
            "config_override": self.resonator_config_override,
        }
        return details

    def get_config(self) -> dict[str, float | int]:
        """."""
        config = get_resonator_config(self.resonator_type, resonator_config_override=self.resonator_config_override)
        return config

    def draw(self):
        """."""
        draw_souk_original_type.draw(
            self.mask_builder,
            self.resonator_type,
            self.x,
            self.y,
            self.rot_angle,
            self.f0,
            mux_func_override=self.mux_func_override,
            resonator_config_override=self.resonator_config_override,
            mirror=self.mirror,
            IDC_and_frame_material=self.IDC_and_frame_material,
            meander_material=self.meander_material,
            trim_length=self.trim_length,
            add_grnd_cutout=self.add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=self.add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=self.add_SiO_cutout,
            add_SiN_membrane_cutout=self.add_SiN_membrane_cutout,
            add_backside_check=self.add_backside_check,
            add_grnd_cutout_over_inductor=self.add_grnd_cutout_over_inductor,
            add_SiN_dep_dielectric_cutout_over_inductor=self.add_SiN_dep_dielectric_cutout_over_inductor,
            add_Aluminium_Patch_and_Etch=self.add_Aluminium_Patch_and_Etch,
            return_configurator_points=self.return_configurator_points,
        )


class SoukCPWCoupledV1TypeBase(ABC):
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
        return_configurator_points: bool = False,
    ):

        # mask_builder
        self.mask_builder = mask_builder

        # resonator_type
        if not isinstance(resonator_type, SoukResonatorType):
            styled_type_error(resonator_type, "resonator_type", SoukResonatorType)

        accepted_resonator_types = [
            SoukResonatorType.CPW_COUPLED_V1,
        ]

        if resonator_type not in accepted_resonator_types:
            raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

        self.resonator_type = resonator_type

        self.x = x
        self.y = y
        self.rot_angle = rot_angle

        self.f0 = f0
        self.mux_func_override = mux_func_override

        self.resonator_config_override = resonator_config_override

        self.mirror = mirror
        self.IDC_and_frame_material = IDC_and_frame_material
        self.meander_material = meander_material
        self.trim_length = trim_length
        self.add_grnd_cutout = add_grnd_cutout
        self.add_SiN_dep_dielectric_cutout = add_SiN_dep_dielectric_cutout
        self.add_SiO_cutout = add_SiO_cutout
        self.add_SiN_membrane_cutout = add_SiN_membrane_cutout
        self.add_backside_check = add_backside_check
        self.add_grnd_cutout_over_inductor = add_grnd_cutout_over_inductor
        self.add_SiN_dep_dielectric_cutout_over_inductor = add_SiN_dep_dielectric_cutout_over_inductor
        self.return_configurator_points = return_configurator_points

    def get_details(self) -> Details:
        """Get the details for the resonator."""

        if self.mux_func_override is None:
            IDC_and_CC_function = get_mux_func_for_resonator_type(self.resonator_type)
            mux_override = None
        else:
            IDC_and_CC_function = self.mux_func_override
            mux_override = getattr(IDC_and_CC_function, "__name__", str(IDC_and_CC_function))

        idcls, ccl = IDC_and_CC_function(self.f0)

        details: Details = {
            "KID_type": str(self.resonator_type),
            "KID_No": None,
            "x_coord": self.x,
            "y_coord": self.y,
            "rot": self.rot_angle,
            "f0": self.f0,
            "mux_override": mux_override,
            "mux_IDC": idcls,
            "mux_CC": ccl,
            "trim": self.trim_length,
            "config_override": self.resonator_config_override,
        }
        return details

    def get_config(self) -> dict[str, float | int]:
        """."""
        config = get_resonator_config(self.resonator_type, resonator_config_override=self.resonator_config_override)
        return config

    def draw(self):
        """."""
        draw_cpw_coupled_v1_type.draw(
            self.mask_builder,
            self.resonator_type,
            self.x,
            self.y,
            self.rot_angle,
            self.f0,
            mux_func_override=self.mux_func_override,
            resonator_config_override=self.resonator_config_override,
            mirror=self.mirror,
            IDC_and_frame_material=self.IDC_and_frame_material,
            meander_material=self.meander_material,
            trim_length=self.trim_length,
            add_grnd_cutout=self.add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=self.add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=self.add_SiO_cutout,
            add_SiN_membrane_cutout=self.add_SiN_membrane_cutout,
            add_backside_check=self.add_backside_check,
            add_grnd_cutout_over_inductor=self.add_grnd_cutout_over_inductor,
            add_SiN_dep_dielectric_cutout_over_inductor=self.add_SiN_dep_dielectric_cutout_over_inductor,
            return_configurator_points=self.return_configurator_points,
        )


class SoukHighVolV1TypeBase(ABC):
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
        # mask_builder
        self.mask_builder = mask_builder

        # resonator_type
        if not isinstance(resonator_type, SoukResonatorType):
            styled_type_error(resonator_type, "resonator_type", SoukResonatorType)

        accepted_resonator_types = [
            SoukResonatorType.HIGH_VOLUME_V1_Q20K,
            SoukResonatorType.HIGH_VOLUME_V1_Q50K,
            SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K,
            SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K,
        ]

        if resonator_type not in accepted_resonator_types:
            raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

        self.resonator_type = resonator_type

        self.x = x
        self.y = y
        self.rot_angle = rot_angle

        self.f0 = f0
        self.mux_func_override = mux_func_override

        self.resonator_config_override = resonator_config_override

        self.mirror = mirror
        self.IDC_and_frame_material = IDC_and_frame_material
        self.meander_material = meander_material
        self.trim_length = trim_length
        self.add_grnd_cutout = add_grnd_cutout
        self.add_SiN_dep_dielectric_cutout = add_SiN_dep_dielectric_cutout
        self.add_SiO_cutout = add_SiO_cutout
        self.add_SiN_membrane_cutout = add_SiN_membrane_cutout
        self.add_backside_check = add_backside_check
        self.add_grnd_cutout_over_inductor = add_grnd_cutout_over_inductor
        self.add_SiN_dep_dielectric_cutout_over_inductor = add_SiN_dep_dielectric_cutout_over_inductor
        self.return_configurator_points = return_configurator_points

    def get_details(self) -> Details:
        """Get the details for the resonator."""

        if self.mux_func_override is None:
            IDC_and_CC_function = get_mux_func_for_resonator_type(self.resonator_type)
            mux_override = None
        else:
            IDC_and_CC_function = self.mux_func_override
            mux_override = getattr(IDC_and_CC_function, "__name__", str(IDC_and_CC_function))

        idcls, ccl = IDC_and_CC_function(self.f0)

        details: Details = {
            "KID_type": str(self.resonator_type),
            "KID_No": None,
            "x_coord": self.x,
            "y_coord": self.y,
            "rot": self.rot_angle,
            "f0": self.f0,
            "mux_override": mux_override,
            "mux_IDC": idcls,
            "mux_CC": ccl,
            "trim": self.trim_length,
            "config_override": self.resonator_config_override,
        }
        return details

    def get_config(self) -> dict[str, float | int]:
        """."""
        config = get_resonator_config(self.resonator_type, resonator_config_override=self.resonator_config_override)
        return config

    def draw(self):
        """."""
        draw_high_volume_v1_type.draw(
            self.mask_builder,
            self.resonator_type,
            self.x,
            self.y,
            self.rot_angle,
            self.f0,
            mux_func_override=self.mux_func_override,
            resonator_config_override=self.resonator_config_override,
            mirror=self.mirror,
            IDC_and_frame_material=self.IDC_and_frame_material,
            meander_material=self.meander_material,
            trim_length=self.trim_length,
            add_grnd_cutout=self.add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=self.add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=self.add_SiO_cutout,
            add_SiN_membrane_cutout=self.add_SiN_membrane_cutout,
            add_backside_check=self.add_backside_check,
            add_grnd_cutout_over_inductor=self.add_grnd_cutout_over_inductor,
            add_SiN_dep_dielectric_cutout_over_inductor=self.add_SiN_dep_dielectric_cutout_over_inductor,
            return_configurator_points=self.return_configurator_points,
        )


class SoukHighVolV2TypeBase(ABC):
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
        return_configurator_points: bool = False,
    ):
        # mask_builder
        self.mask_builder = mask_builder

        # resonator_type
        if not isinstance(resonator_type, SoukResonatorType):
            styled_type_error(resonator_type, SoukResonatorType)

        accepted_resonator_types = [
            SoukResonatorType.HIGH_VOLUME_V2_Q20K,
            SoukResonatorType.HIGH_VOLUME_V2_Q50K,
            SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K,
            SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K,
        ]

        if resonator_type not in accepted_resonator_types:
            raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

        self.resonator_type = resonator_type

        self.x = x
        self.y = y
        self.rot_angle = rot_angle

        self.f0 = f0
        self.mux_func_override = mux_func_override

        self.resonator_config_override = resonator_config_override

        self.mirror = mirror
        self.IDC_and_frame_material = IDC_and_frame_material
        self.meander_material = meander_material
        self.trim_length = trim_length
        self.add_grnd_cutout = add_grnd_cutout
        self.add_SiN_dep_dielectric_cutout = add_SiN_dep_dielectric_cutout
        self.add_SiO_cutout = add_SiO_cutout
        self.add_SiN_membrane_cutout = add_SiN_membrane_cutout
        self.add_backside_check = add_backside_check
        self.add_grnd_cutout_over_inductor = add_grnd_cutout_over_inductor
        self.add_SiN_dep_dielectric_cutout_over_inductor = add_SiN_dep_dielectric_cutout_over_inductor
        self.return_configurator_points = return_configurator_points

    def get_details(self) -> Details:
        """Get the details for the resonator."""

        if self.mux_func_override is None:
            IDC_and_CC_function = get_mux_func_for_resonator_type(self.resonator_type)
            mux_override = None
        else:
            IDC_and_CC_function = self.mux_func_override
            mux_override = getattr(IDC_and_CC_function, "__name__", str(IDC_and_CC_function))

        idcls, ccl = IDC_and_CC_function(self.f0)

        details: Details = {
            "KID_type": str(self.resonator_type),
            "KID_No": None,
            "x_coord": self.x,
            "y_coord": self.y,
            "rot": self.rot_angle,
            "f0": self.f0,
            "mux_override": mux_override,
            "mux_IDC": idcls,
            "mux_CC": ccl,
            "trim": self.trim_length,
            "config_override": self.resonator_config_override,
        }
        return details

    def get_config(self) -> dict[str, float | int]:
        """."""
        config = get_resonator_config(self.resonator_type, resonator_config_override=self.resonator_config_override)
        return config

    def draw(self):
        """."""
        draw_high_volume_v2_type.draw(
            self.mask_builder,
            self.resonator_type,
            self.x,
            self.y,
            self.rot_angle,
            self.f0,
            mux_func_override=self.mux_func_override,
            resonator_config_override=self.resonator_config_override,
            mirror=self.mirror,
            IDC_and_frame_material=self.IDC_and_frame_material,
            meander_material=self.meander_material,
            trim_length=self.trim_length,
            add_grnd_cutout=self.add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=self.add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=self.add_SiO_cutout,
            add_SiN_membrane_cutout=self.add_SiN_membrane_cutout,
            add_backside_check=self.add_backside_check,
            add_grnd_cutout_over_inductor=self.add_grnd_cutout_over_inductor,
            add_SiN_dep_dielectric_cutout_over_inductor=self.add_SiN_dep_dielectric_cutout_over_inductor,
            return_configurator_points=self.return_configurator_points,
        )
