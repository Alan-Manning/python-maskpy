from abc import ABC
from typing import Callable, TypedDict

from ..amber_muxing import get_mux_func_for_resonator_type
from ..layers import Layer
from .draw_funcs import draw_amber_origianl_type
from .resonator_types import AmberResonatorType
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


class AmberOriginalTypeBase(ABC):

    def __init__(
        self,
        mask_builder,
        resonator_type: AmberResonatorType,
        x: float | int,
        y: float | int,
        rot_angle: float | int,
        f0: float | int,
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

        # mask_builder
        self.mask_builder = mask_builder

        # resonator_type
        if not isinstance(resonator_type, AmberResonatorType):
            raise TypeError(f"resonator_type should be of type AmberResonatorType, not {type(resonator_type)}")

        accepted_resonator_types = [
            AmberResonatorType.ORIGINAL_NB_IND2,
            AmberResonatorType.ORIGINAL_NB_IND4,
            AmberResonatorType.ORIGINAL_AL_IND2,
            AmberResonatorType.ORIGINAL_AL_IND4,
        ]

        if resonator_type not in accepted_resonator_types:
            raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

        self.resonator_type = resonator_type

        self.x = x
        self.y = y
        self.rot_angle = rot_angle

        self.f0 = f0
        self.mux_func_override = mux_func_override

        # if mux_func_override is None:
        #     self.mux_func = get_mux_func_for_resonator_type(self.resonator_type)
        # else:
        #     self.mux_func = mux_func_override

        self.resonator_config_override = resonator_config_override

        self.mirror = mirror
        self.IDC_and_frame_material = IDC_and_frame_material
        self.meander_material = meander_material
        self.trim_length = None
        self.coupler_fork_material = coupler_fork_material
        self.add_grnd_cutout = add_grnd_cutout
        self.add_SiN_dep_dielectric_cutout = add_SiN_dep_dielectric_cutout
        self.add_SiO_cutout = add_SiO_cutout
        self.add_SiN_membrane_cutout = add_SiN_membrane_cutout
        self.add_backside_check = add_backside_check
        self.add_inductor_cover = add_inductor_cover

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
            "mux_IDC": list(idcls),
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
        draw_amber_origianl_type.draw(
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
            coupler_fork_material=self.coupler_fork_material,
            add_grnd_cutout=self.add_grnd_cutout,
            add_SiN_dep_dielectric_cutout=self.add_SiN_dep_dielectric_cutout,
            add_SiO_cutout=self.add_SiO_cutout,
            add_SiN_membrane_cutout=self.add_SiN_membrane_cutout,
            add_backside_check=self.add_backside_check,
            add_inductor_cover=self.add_inductor_cover,
        )
