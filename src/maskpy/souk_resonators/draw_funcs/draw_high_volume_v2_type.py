from typing import Callable

import gdspy

from ... import mask_builder_utils as mbu
from ...souk_muxing import get_mux_func_for_resonator_type
from ..resonator_types import SoukResonatorType
from ..utils.get_config import get_resonator_config


def draw(
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
    """Adds the KID geometry to the Main cell athe the x,y cooardinate
    given. The KID is placed where the base middle of the inductive meander
    is at this x,y. The KID geometry is defined by the dimensions within
    the Main_config_file_dict. By default it will, but optionally can
    choose not to, add all the neccessay cutouts for the structure.

    Parameters
    ----------
    mask_builder: SoukMaskBuilder
       The mask builder to add the resonator to.

    resonator_type : SoukResonatorType
        This is the type of resonator to be drawn. The values accepted
        here are a subset of members of the SoukResonatorType enum:
        - SoukResonatorType.HIGH_VOLUME_V2_Q20K
        - SoukResonatorType.HIGH_VOLUME_V2_Q50K
        - SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K
        - SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K

    x,y : float, int
        The x,y coordinates to place the KID. This is the very bottom center
        point of the inductive meander section.

    rot_angle : float, int
        The rotation angle (**in radians**) for the structure. Positive values
        are anti-clockwise, negative is clockwise. The default rotation is with
        the inductive meander at the bottom and coupler at the top with IDC
        arms running horizontally.

    f0 : float, int
        The resonant frequency of the resonator. Should be in the same unit

    KwArgs
    ------
    mux_func_override : Callable | None = None
        This is None or a callable function for getting the IDC and CC
        lengths from a given f0. When the None is provided, The resonator's
        default muxing function will be used. The function should take a
        frequency as an arguments and should return an array-like (28 long)
        and a single float **in this order**, the array-like should contain
        all the lengths for each of the 28 arms for the IDC and the float
        should be the CC length. A simple example funtion:

        >>> def example_IDC_CC_func(
        ...     f0: float | int
        ...     )->tuple[list[float | int], float | int]:
        ...     '''Example function for IDC and CC.
        ...
        ...     f0 : float | int
        ...         Resonant frequency in Hz.
        ...     Returns
        ...     -------
        ...     IDC_lengths: list[float | int]
        ...     CC_length: float | int
        ...     '''
        ...     if (f0 < 3.1e9):
        ...         IDC_lengths = np.ones(28)*1900.0
        ...         CC_length = 600.0
        ...     else:
        ...         IDC_lengths = np.ones(28)*1500.0
        ...         CC_length = 300.0
        ...
        ...     return IDC_lengths, CC_length

    resonator_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    mirror=False
        Whether the KID should be mirrored about the center vertical, **this
        mirroring is done before any roation is applied**. By default
        (with 0° ratation) the KID's coupler is attached on the left but when
        mirror=True the coupler on the right.

    IDC_and_frame_material = "IDC_Nb"
        The material to make the IDC and frame structure out of. By Default
        this is "IDC_Nb" which will make them out of the IDC_Nb material.
        This can take any of the values "IDC_Nb", "Nb", or "Al", which will
        make the frame and IDC out of IDC_Nb, Nb_Antenna, or Aluminium
        layers respectively.

    meander_material = "Al"
        The material to make the inductive meander structure out of. By
        Default this is "Al" which will make it out of the Aluminium
        material. This can take any of the values "Al", "IDC_Nb", or "Nb"
        which will make the inductive meander out of Aluminium, IDC_Nb, or
        Nb_Antenna layers respectively.

    trim_length=None
        When None nothing is done. When a float or int value is passed,
        there will be a trim layer added which will overlap the trim
        fingers to bring the trim fingers down to the length specified.
        For example, if you pass trim_length=1500. There will be a trim box
        added over both trim fingers that will make the length of the trim
        fingers equal this 1500um value.

    add_grnd_cutout=True
        Whether or not to add a cutout in the Nb_Groundplane layer in the
        neccary place for the KID structure.

    add_SiN_dep_dielectric_cutout=True
        Whether or not to add a cutout in the SiN depositon layer in the
        neccary place for the KID structure.

    add_SiO_cutout=True
        Whether or not to add a cutout in the Silicon Oxide layer in the
        neccary place for the KID structure.

    add_SiN_membrane_cutout=True
        Whether or not to add a cutout in the Silicon Nitride membrane layer
        in the neccary place for the KID structure.

    add_backside_check=False
        Whether or not to add a backside check cover in the neccary place
        for the KID structure.

    add_grnd_cutout_over_inductor=False
        Whether or not to add a groundplane cutout over the inductive
        meander. Defaulf is false, when True will create a cutout over the
        mander that is oversived by 30um in all directions relative to the
        center of the inductive meander.

    add_SiN_dep_dielectric_cutout_over_inductor=False
        Whether or not to add a SiN depositon cutout over the inductive
        meander. Defaulf is false, when True will create a cutout over the
        mander that is oversived by 20um in all directions relative to the
        center of the inductive meander.

    return_configurator_points=False
        return a the points for use in the configurator.
    """
    if not isinstance(resonator_type, SoukResonatorType):
        raise TypeError(f"resonator_type should be of type SoukResonatorType, not {type(resonator_type)}")

    accepted_resonator_types = [
        SoukResonatorType.HIGH_VOLUME_V2_Q20K,
        SoukResonatorType.HIGH_VOLUME_V2_Q50K,
        SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K,
        SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K,
    ]

    if resonator_type not in accepted_resonator_types:
        raise ValueError(f"resonator_type is not compatible, should be one of {accepted_resonator_types}.")

    resonator_config = get_resonator_config(resonator_type, resonator_config_override=resonator_config_override)

    if mux_func_override is None:
        IDC_and_CC_function = get_mux_func_for_resonator_type(resonator_type)
    else:
        IDC_and_CC_function = mux_func_override

    IDC_and_frame_material_lookup = {
        "IDC_Nb": mask_builder.layers.IDC_Nb,
        "Nb": mask_builder.layers.Nb_Antenna,
        "Al": mask_builder.layers.Aluminium,
    }
    material_idc_and_frame = IDC_and_frame_material_lookup[IDC_and_frame_material]

    meander_material_lookup = {
        "Al": mask_builder.layers.Aluminium,
        "IDC_Nb": mask_builder.layers.IDC_Nb,
        "Nb": mask_builder.layers.Nb_Antenna,
    }
    material_meander = meander_material_lookup[meander_material]

    # Making the meander section
    meander_lw = resonator_config["meander_lw"]
    meander_corner_bend_radius = resonator_config["meander_corner_bend_radius"]

    meander_bot_width = resonator_config["meander_bot_width"]

    meander_right_height_1 = resonator_config["meander_right_height_1"]
    meander_right_height_2 = resonator_config["meander_right_height_2"]
    meander_right_height_3 = resonator_config["meander_right_height_3"]
    meander_right_height_4 = resonator_config["meander_right_height_4"]

    meander_left_height_1 = resonator_config["meander_left_height_1"]
    meander_left_height_2 = resonator_config["meander_left_height_2"]
    meander_left_height_3 = resonator_config["meander_left_height_3"]
    meander_left_height_4 = resonator_config["meander_left_height_4"]

    meander_left_width_1 = resonator_config["meander_left_width_1"]
    meander_left_width_2 = resonator_config["meander_left_width_2"]
    meander_left_width_3 = resonator_config["meander_left_width_3"]

    meander_right_width_1 = resonator_config["meander_right_width_1"]
    meander_right_width_2 = resonator_config["meander_right_width_2"]
    meander_right_width_3 = resonator_config["meander_right_width_3"]

    meander_step_back_from_frame = resonator_config["meander_step_back_from_frame"]  # TODO Add this to the config.

    meander_path_points = [
        # TopRight to TopLeft clockwise
        [
            (meander_bot_width / 2) - meander_right_width_1 - meander_right_width_2 + meander_right_width_3,
            meander_lw
            + meander_right_height_1
            - meander_right_height_2
            + meander_right_height_3
            + meander_right_height_4
            - meander_step_back_from_frame,
        ],
        [
            (meander_bot_width / 2) - meander_right_width_1 - meander_right_width_2 + meander_right_width_3,
            (meander_lw / 2) + meander_right_height_1 - meander_right_height_2 + meander_right_height_3,
        ],
        [
            (meander_bot_width / 2) - meander_right_width_1 - meander_right_width_2,
            (meander_lw / 2) + meander_right_height_1 - meander_right_height_2 + meander_right_height_3,
        ],
        [
            (meander_bot_width / 2) - meander_right_width_1 - meander_right_width_2,
            (meander_lw / 2) + meander_right_height_1 - meander_right_height_2,
        ],
        [(meander_bot_width / 2) - meander_right_width_1, (meander_lw / 2) + meander_right_height_1 - meander_right_height_2],
        [(meander_bot_width / 2) - meander_right_width_1, (meander_lw / 2) + meander_right_height_1],
        [(meander_bot_width / 2), (meander_lw / 2) + meander_right_height_1],
        [(meander_bot_width / 2), (meander_lw / 2)],  # Mid right hand side
        [-(meander_bot_width / 2), (meander_lw / 2)],  # Mid left hand side
        [-(meander_bot_width / 2), (meander_lw / 2) + meander_left_height_1],
        [-(meander_bot_width / 2) - meander_left_width_1, (meander_lw / 2) + meander_left_height_1],
        [-(meander_bot_width / 2) - meander_left_width_1, (meander_lw / 2) + meander_left_height_1 - meander_left_height_2],
        [
            -(meander_bot_width / 2) - meander_left_width_1 - meander_left_width_2,
            (meander_lw / 2) + meander_left_height_1 - meander_left_height_2,
        ],
        [
            -(meander_bot_width / 2) - meander_left_width_1 - meander_left_width_2,
            (meander_lw / 2) + meander_left_height_1 - meander_left_height_2 + meander_left_height_3,
        ],
        [
            -(meander_bot_width / 2) - meander_left_width_1 - meander_left_width_2 + meander_left_width_3,
            (meander_lw / 2) + meander_left_height_1 - meander_left_height_2 + meander_left_height_3,
        ],
        [
            -(meander_bot_width / 2) - meander_left_width_1 - meander_left_width_2 + meander_left_width_3,
            meander_lw
            + meander_left_height_1
            - meander_left_height_2
            + meander_left_height_3
            + meander_left_height_4
            - meander_step_back_from_frame,
        ],
    ]

    # Making the meander ant pad
    ant_pad_box_width = resonator_config["ant_pad_box_width"]  # 5
    ant_pad_box_height = resonator_config["ant_pad_box_height"]  # 10

    ant_pad_box_points = [
        [-ant_pad_box_width / 2, 0],
        [-ant_pad_box_width / 2, -ant_pad_box_height],
        [ant_pad_box_width / 2, -ant_pad_box_height],
        [ant_pad_box_width / 2, 0],
    ]

    # Making the frame sections left and right
    frame_bot_lw = resonator_config["frame_bot_lw"]  # 8
    frame_bot_left_offset = resonator_config["frame_bot_left_offset"]
    frame_bot_right_offset = resonator_config["frame_bot_right_offset"]

    frame_left_lw = resonator_config["frame_left_lw"]  # 8
    frame_left_height = resonator_config["frame_left_height"]  # 400

    frame_right_lw = resonator_config["frame_right_lw"]  # 8
    frame_right_height = resonator_config["frame_right_height"]  # 400

    frame_start_y = meander_path_points[-1][1] + meander_step_back_from_frame
    frame_left_start_x = meander_path_points[-1][0] + (meander_lw / 2)
    frame_right_start_x = meander_path_points[0][0] - (meander_lw / 2)

    frame_meander_lw = resonator_config["frame_meander_lw"]
    frame_meander_corner_bend_radius = resonator_config["frame_meander_corner_bend_radius"]

    left_frame_meander_width_1 = resonator_config["left_frame_meander_width_1"]
    left_frame_meander_width_2 = resonator_config["left_frame_meander_width_2"]
    left_frame_meander_height_1 = resonator_config["left_frame_meander_height_1"]
    left_frame_meander_height_2 = resonator_config["left_frame_meander_height_2"]

    right_frame_meander_width_1 = resonator_config["right_frame_meander_width_1"]
    right_frame_meander_width_2 = resonator_config["right_frame_meander_width_2"]
    right_frame_meander_height_1 = resonator_config["right_frame_meander_height_1"]
    right_frame_meander_height_2 = resonator_config["right_frame_meander_height_2"]

    frame_left_points = [
        # [frame_left_start_x, frame_start_y + frame_bot_lw],
        [
            frame_left_start_x - frame_bot_left_offset + frame_left_lw,
            frame_start_y + left_frame_meander_height_1 + left_frame_meander_height_2,
        ],
        [frame_left_start_x - frame_bot_left_offset + frame_left_lw, frame_start_y + frame_left_height],
        [frame_left_start_x - frame_bot_left_offset, frame_start_y + frame_left_height],
        [frame_left_start_x - frame_bot_left_offset, frame_start_y + left_frame_meander_height_1 + left_frame_meander_height_2],
        # [frame_left_start_x, frame_start_y],
    ]

    frame_right_points = [
        # [frame_right_start_x, frame_start_y + frame_bot_lw],
        [
            frame_right_start_x + frame_bot_right_offset - frame_right_lw,
            frame_start_y + left_frame_meander_height_1 + left_frame_meander_height_2,
        ],
        [frame_right_start_x + frame_bot_right_offset - frame_right_lw, frame_start_y + frame_right_height],
        [frame_right_start_x + frame_bot_right_offset, frame_start_y + frame_right_height],
        [frame_right_start_x + frame_bot_right_offset, frame_start_y + left_frame_meander_height_1 + left_frame_meander_height_2],
        # [frame_right_start_x, frame_start_y],
    ]

    # Making the extra frame meander path sections left and right.
    bot_left_frame_path_points = [
        [frame_left_start_x, frame_start_y + (frame_meander_lw / 2)],
        [frame_left_start_x - left_frame_meander_width_1, frame_start_y + (frame_meander_lw / 2)],
        [frame_left_start_x - left_frame_meander_width_1, frame_start_y + left_frame_meander_height_1 + (frame_meander_lw / 2)],
        [
            frame_left_start_x - left_frame_meander_width_1 + left_frame_meander_width_2,
            frame_start_y + left_frame_meander_height_1 + (frame_meander_lw / 2),
        ],
        [
            frame_left_start_x - left_frame_meander_width_1 + left_frame_meander_width_2,
            frame_start_y + left_frame_meander_height_1 + left_frame_meander_height_2 + (frame_meander_lw / 2),
        ],
        [
            frame_left_start_x - left_frame_meander_width_1 + left_frame_meander_width_2 - frame_bot_left_offset + frame_left_lw,
            frame_start_y + left_frame_meander_height_1 + left_frame_meander_height_2 + (frame_meander_lw / 2),
        ],
    ]

    bot_right_frame_path_points = [
        [frame_right_start_x, frame_start_y + (frame_meander_lw / 2)],
        [frame_right_start_x + right_frame_meander_width_1, frame_start_y + (frame_meander_lw / 2)],
        [frame_right_start_x + right_frame_meander_width_1, frame_start_y + right_frame_meander_height_1 + (frame_meander_lw / 2)],
        [
            frame_right_start_x + right_frame_meander_width_1 - right_frame_meander_width_2,
            frame_start_y + right_frame_meander_height_1 + (frame_meander_lw / 2),
        ],
        [
            frame_right_start_x + right_frame_meander_width_1 - right_frame_meander_width_2,
            frame_start_y + right_frame_meander_height_1 + right_frame_meander_height_2 + (frame_meander_lw / 2),
        ],
        [
            frame_right_start_x + right_frame_meander_width_1 - right_frame_meander_width_2 + frame_bot_right_offset - frame_right_lw,
            frame_start_y + right_frame_meander_height_1 + right_frame_meander_height_2 + (frame_meander_lw / 2),
        ],
    ]

    # Making the pads the conect the meander to the frame
    frame_meander_cover_box_width = resonator_config["frame_meander_cover_box_width"]  # 5
    frame_meander_cover_box_height = resonator_config["frame_meander_cover_box_height"]  # 28

    extra_frame_meander_cover_box_width = resonator_config["extra_frame_meander_cover_box_width"]  # TODO add this to the config
    extra_frame_meander_cover_box_height = resonator_config["extra_frame_meander_cover_box_height"]  # TODO add this to the config

    meander_cover_box_right_points = [
        [meander_path_points[0][0] - (frame_meander_cover_box_width / 2), frame_start_y + frame_meander_lw],
        [
            meander_path_points[0][0] - (frame_meander_cover_box_width / 2),
            frame_start_y + frame_meander_lw - frame_meander_cover_box_height,
        ],
        [
            meander_path_points[0][0] + (frame_meander_cover_box_width / 2),
            frame_start_y + frame_meander_lw - frame_meander_cover_box_height,
        ],
        [meander_path_points[0][0] + (frame_meander_cover_box_width / 2), frame_start_y + frame_meander_lw],
    ]

    meander_cover_box_left_points = [
        [meander_path_points[-1][0] - (frame_meander_cover_box_width / 2), frame_start_y + frame_meander_lw],
        [
            meander_path_points[-1][0] - (frame_meander_cover_box_width / 2),
            frame_start_y + frame_meander_lw - frame_meander_cover_box_height,
        ],
        [
            meander_path_points[-1][0] + (frame_meander_cover_box_width / 2),
            frame_start_y + frame_meander_lw - frame_meander_cover_box_height,
        ],
        [meander_path_points[-1][0] + (frame_meander_cover_box_width / 2), frame_start_y + frame_meander_lw],
    ]

    meander_cover_box_right_end_center = [
        (meander_cover_box_right_points[1][0] + meander_cover_box_right_points[2][0]) / 2,
        (meander_cover_box_right_points[1][1] + meander_cover_box_right_points[2][1]) / 2,
    ]

    meander_cover_box_left_end_center = [
        (meander_cover_box_left_points[1][0] + meander_cover_box_left_points[2][0]) / 2,
        (meander_cover_box_left_points[1][1] + meander_cover_box_left_points[2][1]) / 2,
    ]

    extra_meander_cover_box_right_points = [
        [
            meander_cover_box_right_end_center[0] - (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_right_end_center[1] + (extra_frame_meander_cover_box_height / 2),
        ],
        [
            meander_cover_box_right_end_center[0] - (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_right_end_center[1] - (extra_frame_meander_cover_box_height / 2),
        ],
        [
            meander_cover_box_right_end_center[0] + (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_right_end_center[1] - (extra_frame_meander_cover_box_height / 2),
        ],
        [
            meander_cover_box_right_end_center[0] + (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_right_end_center[1] + (extra_frame_meander_cover_box_height / 2),
        ],
    ]

    extra_meander_cover_box_left_points = [
        [
            meander_cover_box_left_end_center[0] - (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_left_end_center[1] + (extra_frame_meander_cover_box_height / 2),
        ],
        [
            meander_cover_box_left_end_center[0] - (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_left_end_center[1] - (extra_frame_meander_cover_box_height / 2),
        ],
        [
            meander_cover_box_left_end_center[0] + (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_left_end_center[1] - (extra_frame_meander_cover_box_height / 2),
        ],
        [
            meander_cover_box_left_end_center[0] + (extra_frame_meander_cover_box_width / 2),
            meander_cover_box_left_end_center[1] + (extra_frame_meander_cover_box_height / 2),
        ],
    ]

    # Making the coupler attachement
    coupler_frame_left_lw = resonator_config["coupler_frame_left_lw"]  # 10
    coupler_frame_left_height = resonator_config["coupler_frame_left_height"]  # 39
    coupler_frame_top_lw = resonator_config["coupler_frame_top_lw"]  # 3

    coupler_frame_start_x = frame_left_start_x - frame_bot_left_offset + (frame_left_lw / 2)
    coupler_frame_start_y = frame_start_y + frame_left_height

    # Making the IDC and trim arms
    IDC_bot_arm_gap = resonator_config["IDC_bot_arm_gap"]  # 30
    IDC_arm_gap = resonator_config["IDC_arm_gap"]  # 8
    IDC_arm_lw = resonator_config["IDC_arm_lw"]  # 3
    No_of_arms = resonator_config["No_of_arms"]  # 28

    arm_start_x_left_side = frame_left_start_x - frame_bot_left_offset + frame_left_lw
    arm_start_x_right_side = frame_right_start_x + frame_bot_right_offset - frame_right_lw

    arm_start_y_right_side = frame_start_y + frame_bot_lw + IDC_bot_arm_gap + (IDC_arm_lw / 2)
    arm_start_y_left_side = arm_start_y_right_side + IDC_arm_gap + IDC_arm_lw

    trim_arm_offset_right_side = resonator_config["trim_arm_offset_right_side"]  # 380
    trim_arm_offset_left_side = resonator_config["trim_arm_offset_left_side"]  # 389
    trim_arm_lw = resonator_config["trim_arm_lw"]  # 3
    trim_arm_length_right_side = resonator_config["trim_arm_length_right_side"]  # 1975
    trim_arm_length_left_side = resonator_config["trim_arm_length_left_side"]  # 1975

    trim_arm_start_y_right_side = frame_start_y + frame_bot_lw + trim_arm_offset_right_side
    trim_arm_start_y_left_side = frame_start_y + frame_bot_lw + trim_arm_offset_left_side

    # Making the coupler ataching to feedline
    coupler_gap = resonator_config["coupler_gap"]  # 16
    coupler_lw = resonator_config["coupler_lw"]  # 3
    left_coupler_frame_to_feed_distance = resonator_config["left_coupler_frame_to_feed_distance"]  # 164

    # Making the ground plane cutout
    cutout_bot_offset = resonator_config["cutout_bot_offset"]  # 15
    cutout_left_offset = resonator_config["cutout_left_offset"]  # 50
    cutout_right_offset = resonator_config["cutout_right_offset"]  # 50
    cutout_top_offset = resonator_config["cutout_top_offset"]  # 25

    # grndpl_meander_cutout_width = config["grndpl_meander_cutout_width"]#80
    # grndpl_meander_cutout_height = config["grndpl_meander_cutout_height"]#10

    cutout_start_height = frame_start_y - cutout_bot_offset
    cutout_width = (frame_right_start_x + frame_bot_right_offset + cutout_right_offset) - (
        frame_left_start_x - frame_bot_left_offset - cutout_left_offset
    )
    cutout_height = (coupler_frame_start_y + coupler_frame_left_height + coupler_gap + coupler_lw + cutout_top_offset) - cutout_start_height

    step_down_distance_between_layers = resonator_config["step_down_distance_between_layers"]  # 5

    grnd_plane_cutout_width = cutout_width + (3 * (2 * step_down_distance_between_layers))
    grnd_plane_cutout_height = cutout_height + (3 * (2 * step_down_distance_between_layers))
    grnd_plane_cutout_start_height = cutout_start_height - (3 * step_down_distance_between_layers)

    grnd_plane_meander_cutout_poly_points = [
        [-grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height],
        [-grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height + grnd_plane_cutout_height],
        [grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height + grnd_plane_cutout_height],
        [grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height],
    ]

    grnd_plane_inductor_cutout_offset_right = 30
    grnd_plane_inductor_cutout_offset_left = 30
    grnd_plane_inductor_cutout_offset_top = 30
    grnd_plane_inductor_cutout_offset_bot = 30

    grnd_plane_inductor_cutout_poly_points = [
        # bot right, bot left, top left, top right
        [
            meander_path_points[5][0] + grnd_plane_inductor_cutout_offset_right,
            meander_path_points[5][1] - grnd_plane_inductor_cutout_offset_bot,
        ],
        [
            meander_path_points[8][0] - grnd_plane_inductor_cutout_offset_left,
            meander_path_points[5][1] - grnd_plane_inductor_cutout_offset_bot,
        ],
        [
            meander_path_points[8][0] - grnd_plane_inductor_cutout_offset_left,
            meander_path_points[0][1] + grnd_plane_inductor_cutout_offset_top,
        ],
        [
            meander_path_points[0][0] + grnd_plane_inductor_cutout_offset_right,
            meander_path_points[0][1] + grnd_plane_inductor_cutout_offset_top,
        ],
    ]

    SiN_dep_inductor_cutout_offset_right = 20
    SiN_dep_inductor_cutout_offset_left = 20
    SiN_dep_inductor_cutout_offset_top = 20
    SiN_dep_inductor_cutout_offset_bot = 20

    SiN_dep_inductor_cutout_poly_points = [
        # bot right, bot left, top left, top right
        [
            meander_path_points[5][0] + SiN_dep_inductor_cutout_offset_right,
            meander_path_points[5][1] - SiN_dep_inductor_cutout_offset_bot,
        ],
        [
            meander_path_points[8][0] - SiN_dep_inductor_cutout_offset_left,
            meander_path_points[5][1] - SiN_dep_inductor_cutout_offset_bot,
        ],
        [
            meander_path_points[8][0] - SiN_dep_inductor_cutout_offset_left,
            meander_path_points[0][1] + SiN_dep_inductor_cutout_offset_top,
        ],
        [
            meander_path_points[0][0] + SiN_dep_inductor_cutout_offset_right,
            meander_path_points[0][1] + SiN_dep_inductor_cutout_offset_top,
        ],
    ]

    SiN_dep_cutout_width = cutout_width + 20
    SiN_dep_cutout_height = cutout_height + 20
    SiN_dep_cutout_start_height = cutout_start_height - 10

    SiN_dep_cutout_poly_points = [
        [-SiN_dep_cutout_width / 2, SiN_dep_cutout_start_height],
        [-SiN_dep_cutout_width / 2, SiN_dep_cutout_start_height + SiN_dep_cutout_height],
        [SiN_dep_cutout_width / 2, SiN_dep_cutout_start_height + SiN_dep_cutout_height],
        [SiN_dep_cutout_width / 2, SiN_dep_cutout_start_height],
    ]

    # Making the SiO cutout rect
    SiO_stepdown_cutout_width = resonator_config["SiO_stepdown_cutout_width"]  # 110
    SiO_stepdown_cutout_height = resonator_config["SiO_stepdown_cutout_height"]  # 39
    SiO_cutout_width = cutout_width
    SiO_cutout_height = cutout_height
    SiO_cutout_start_height = cutout_start_height

    SiO_cutout_poly_points = [
        [-SiO_stepdown_cutout_width / 2, SiO_cutout_start_height + SiO_stepdown_cutout_height],
        [-SiO_stepdown_cutout_width / 2, SiO_cutout_start_height],
        [-SiO_cutout_width / 2, SiO_cutout_start_height],
        [-SiO_cutout_width / 2, SiO_cutout_start_height + SiO_cutout_height],
        [SiO_cutout_width / 2, SiO_cutout_start_height + SiO_cutout_height],
        [SiO_cutout_width / 2, SiO_cutout_start_height],
        [SiO_stepdown_cutout_width / 2, SiO_cutout_start_height],
        [SiO_stepdown_cutout_width / 2, SiO_cutout_start_height + SiO_stepdown_cutout_height],
    ]

    # Making the SiN membrane cutout
    SiN_membrane_stepdown_cutout_width = resonator_config["SiN_membrane_stepdown_cutout_width"]  # 100
    SiN_membrane_stepdown_cutout_height = resonator_config["SiN_membrane_stepdown_cutout_height"]  # 36
    SiN_membrane_cutout_width = cutout_width + 10
    SiN_membrane_cutout_height = cutout_height + 10
    SiN_membrane_cutout_start_height = cutout_start_height - 5

    SiN_membrane_cutout_poly_points = [
        [-SiN_membrane_stepdown_cutout_width / 2, SiN_membrane_cutout_start_height + SiN_membrane_stepdown_cutout_height],
        [-SiN_membrane_stepdown_cutout_width / 2, SiN_membrane_cutout_start_height],
        [-SiN_membrane_cutout_width / 2, SiN_membrane_cutout_start_height],
        [-SiN_membrane_cutout_width / 2, SiN_membrane_cutout_start_height + SiN_membrane_cutout_height],
        [SiN_membrane_cutout_width / 2, SiN_membrane_cutout_start_height + SiN_membrane_cutout_height],
        [SiN_membrane_cutout_width / 2, SiN_membrane_cutout_start_height],
        [SiN_membrane_stepdown_cutout_width / 2, SiN_membrane_cutout_start_height],
        [SiN_membrane_stepdown_cutout_width / 2, SiN_membrane_cutout_start_height + SiN_membrane_stepdown_cutout_height],
    ]

    # Making the backside check covers
    backside_check_cover_width = cutout_width
    backside_check_cover_height = cutout_height
    backside_check_cutout_start_height = cutout_start_height

    backside_check_cover_poly_points = [
        [-backside_check_cover_width / 2, backside_check_cutout_start_height],
        [-backside_check_cover_width / 2, backside_check_cutout_start_height + backside_check_cover_height],
        [backside_check_cover_width / 2, backside_check_cutout_start_height + backside_check_cover_height],
        [backside_check_cover_width / 2, backside_check_cutout_start_height],
    ]

    # Making the Aluminium Patch and Etch
    aluminium_patch_offset_bot = 0
    aluminium_patch_offset_top = 0
    aluminium_patch_offset_left = 0
    aluminium_patch_offset_right = 0
    aluminium_etch_offset_bot = 5
    aluminium_etch_offset_top = 5
    aluminium_etch_offset_left = 5
    aluminium_etch_offset_right = 5

    # defining the box that is the exact bounding box of the meander from the top left of the L shape then going counter colockwise
    box_around_meander = [
        [
            -(meander_bot_width / 2) - (meander_lw / 2) - meander_left_width_1 - meander_left_width_2,
            meander_lw
            + meander_left_height_1
            - meander_left_height_2
            + meander_left_height_3
            + meander_left_height_4
            - meander_step_back_from_frame,
        ],
        [
            -(meander_bot_width / 2) - (meander_lw / 2) - meander_left_width_1 - meander_left_width_2,
            meander_left_height_1 - meander_left_height_2,
        ],
        [
            -(meander_bot_width / 2) + (meander_lw / 2) - meander_left_width_1,
            meander_left_height_1 - meander_left_height_2,
        ],
        [
            -(meander_bot_width / 2) + (meander_lw / 2) - meander_left_width_1,
            0,
        ],
        [
            (meander_bot_width / 2)
            + (extra_frame_meander_cover_box_width / 2)
            - meander_right_width_1
            - meander_right_width_2
            + meander_right_width_3,
            0,
        ],  # bot right side of the meander base
        [
            (meander_bot_width / 2)
            + (extra_frame_meander_cover_box_width / 2)
            - meander_right_width_1
            - meander_right_width_2
            + meander_right_width_3,
            meander_lw
            + meander_right_height_1
            - meander_right_height_2
            + meander_right_height_3
            + meander_right_height_4
            - meander_step_back_from_frame,
        ],
    ]

    # defined from bot left going anti-clockwise
    aluminium_patch_under_meander_points = [
        [
            box_around_meander[0][0] - aluminium_patch_offset_left,
            box_around_meander[0][1] + aluminium_patch_offset_top,
        ],
        [
            box_around_meander[1][0] - aluminium_patch_offset_left,
            box_around_meander[1][1] - aluminium_patch_offset_bot,
        ],
        [
            box_around_meander[2][0] + aluminium_patch_offset_right,
            box_around_meander[2][1] - aluminium_patch_offset_bot,
        ],
        [
            box_around_meander[3][0] + aluminium_patch_offset_right,
            box_around_meander[3][1] - aluminium_patch_offset_bot,
        ],
        [
            box_around_meander[4][0] + aluminium_patch_offset_right,
            box_around_meander[4][1] - aluminium_patch_offset_bot,
        ],
        [
            box_around_meander[5][0] + aluminium_patch_offset_right,
            box_around_meander[5][1] + aluminium_patch_offset_top,
        ],
    ]

    # Making the aluminium etch under the meander
    # defined from bot left going anti-clockwise
    aluminium_etch_under_meander_points = [
        [
            box_around_meander[0][0] - aluminium_etch_offset_left,
            box_around_meander[0][1] + aluminium_etch_offset_top,
        ],
        [
            box_around_meander[1][0] - aluminium_etch_offset_left,
            box_around_meander[1][1] - aluminium_etch_offset_bot,
        ],
        [
            box_around_meander[2][0] + aluminium_etch_offset_right,
            box_around_meander[2][1] - aluminium_etch_offset_bot,
        ],
        [
            box_around_meander[3][0] + aluminium_etch_offset_right,
            box_around_meander[3][1] - aluminium_etch_offset_bot,
        ],
        [
            box_around_meander[4][0] + aluminium_etch_offset_right,
            box_around_meander[4][1] - aluminium_etch_offset_bot,
        ],
        [
            box_around_meander[5][0] + aluminium_etch_offset_right,
            box_around_meander[5][1] + aluminium_etch_offset_top,
        ],
    ]

    # Making the aluminium patch under the IDC
    aluminium_patch_IDC_left_edge = (
        frame_left_points[-2][0] - (frame_left_lw / 2) + (coupler_frame_left_lw / 2) - aluminium_patch_offset_left
    )
    aluminium_patch_IDC_right_edge = frame_right_points[-2][0] + aluminium_patch_offset_right
    aluminium_patch_IDC_top_edge = coupler_frame_start_y + coupler_frame_left_height
    aluminium_patch_IDC_bot_edge = frame_start_y + frame_bot_lw - frame_meander_cover_box_height

    # defined from bot left going anti-clockwise
    aluminium_patch_under_IDC_points = [
        [aluminium_patch_IDC_left_edge, aluminium_patch_IDC_bot_edge],
        [aluminium_patch_IDC_right_edge, aluminium_patch_IDC_bot_edge],
        [aluminium_patch_IDC_right_edge, aluminium_patch_IDC_top_edge],
        [aluminium_patch_IDC_left_edge, aluminium_patch_IDC_top_edge],
    ]

    aluminium_etch_under_IDC_points = [
        [
            aluminium_patch_under_IDC_points[0][0] - aluminium_etch_offset_left,
            aluminium_patch_under_IDC_points[0][1] - aluminium_etch_offset_bot,
        ],
        [
            aluminium_patch_under_IDC_points[1][0] + aluminium_etch_offset_right,
            aluminium_patch_under_IDC_points[1][1] - aluminium_etch_offset_bot,
        ],
        [
            aluminium_patch_under_IDC_points[2][0] + aluminium_etch_offset_right,
            aluminium_patch_under_IDC_points[2][1] + aluminium_etch_offset_top,
        ],
        [
            aluminium_patch_under_IDC_points[3][0] - aluminium_etch_offset_left,
            aluminium_patch_under_IDC_points[3][1] + aluminium_etch_offset_top,
        ],
    ]

    # Getting the IDC and CC lengths from the function
    IDCLs, CCL = IDC_and_CC_function(f0)

    # Adding the meander

    if mirror:
        new_meander_path_points = mbu.mirror_points_around_yaxis(meander_path_points)
        new_meander_path_points = mbu.rotate_and_move_points_list(new_meander_path_points, rot_angle, x, y)
    else:
        new_meander_path_points = mbu.rotate_and_move_points_list(meander_path_points, rot_angle, x, y)
    meander_path = gdspy.FlexPath(
        new_meander_path_points,
        meander_lw,
        corners="circular bend",
        bend_radius=meander_corner_bend_radius,
        layer=material_meander.number,
        datatype=material_meander.datatype,
    )
    meander_path_polygons = mbu.get_polys_from_flexpath(meander_path)
    for i in range(len(meander_path_polygons)):
        mask_builder.Main.add(
            gdspy.Polygon(
                meander_path_polygons[i],
                layer=material_meander.number,
                datatype=material_meander.datatype,
            )
        )

    # Adding the meander ant overlap box
    new_ant_pad_box_points = mbu.rotate_and_move_points_list(ant_pad_box_points, rot_angle, x, y)
    ant_pad_box = gdspy.Polygon(
        new_ant_pad_box_points,
        layer=material_meander.number,
        datatype=material_meander.datatype,
    )
    mask_builder.Main.add(ant_pad_box)

    # Adding the meander frame overlap boxes
    if mirror:
        new_meander_cover_box_right_points = mbu.mirror_points_around_yaxis(meander_cover_box_right_points)
        new_meander_cover_box_right_points = mbu.rotate_and_move_points_list(new_meander_cover_box_right_points, rot_angle, x, y)
    else:
        new_meander_cover_box_right_points = mbu.rotate_and_move_points_list(meander_cover_box_right_points, rot_angle, x, y)

    meander_cover_box_right = gdspy.Polygon(
        new_meander_cover_box_right_points,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    mask_builder.Main.add(meander_cover_box_right)

    if mirror:
        new_meander_cover_box_left_points = mbu.mirror_points_around_yaxis(meander_cover_box_left_points)
        new_meander_cover_box_left_points = mbu.rotate_and_move_points_list(new_meander_cover_box_left_points, rot_angle, x, y)
    else:
        new_meander_cover_box_left_points = mbu.rotate_and_move_points_list(meander_cover_box_left_points, rot_angle, x, y)

    meander_cover_box_left = gdspy.Polygon(
        new_meander_cover_box_left_points,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    mask_builder.Main.add(meander_cover_box_left)

    # Adding the extra frame meander left and right.
    if mirror:
        new_bot_left_frame_path_points = mbu.mirror_points_around_yaxis(bot_left_frame_path_points)
        new_bot_left_frame_path_points = mbu.rotate_and_move_points_list(new_bot_left_frame_path_points, rot_angle, x, y)
    else:
        new_bot_left_frame_path_points = mbu.rotate_and_move_points_list(bot_left_frame_path_points, rot_angle, x, y)

    bot_left_frame_path = gdspy.FlexPath(
        new_bot_left_frame_path_points,
        frame_meander_lw,
        corners="circular bend",
        bend_radius=frame_meander_corner_bend_radius,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )

    if mirror:
        new_bot_right_frame_path_points = mbu.mirror_points_around_yaxis(bot_right_frame_path_points)
        new_bot_right_frame_path_points = mbu.rotate_and_move_points_list(new_bot_right_frame_path_points, rot_angle, x, y)
    else:
        new_bot_right_frame_path_points = mbu.rotate_and_move_points_list(bot_right_frame_path_points, rot_angle, x, y)

    bot_right_frame_path = gdspy.FlexPath(
        new_bot_right_frame_path_points,
        frame_meander_lw,
        corners="circular bend",
        bend_radius=frame_meander_corner_bend_radius,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    mask_builder.make_flexpath_into_polygons_and_add_to_main(
        bot_left_frame_path,
        material_idc_and_frame,
    )
    mask_builder.make_flexpath_into_polygons_and_add_to_main(
        bot_right_frame_path,
        material_idc_and_frame,
    )

    # Adding the extra meander frame overlap boxes
    if mirror:
        new_extra_meander_cover_box_right_points = mbu.mirror_points_around_yaxis(extra_meander_cover_box_right_points)
        new_extra_meander_cover_box_right_points = mbu.rotate_and_move_points_list(
            new_extra_meander_cover_box_right_points, rot_angle, x, y
        )
    else:
        new_extra_meander_cover_box_right_points = mbu.rotate_and_move_points_list(extra_meander_cover_box_right_points, rot_angle, x, y)

    extra_meander_cover_box_right = gdspy.Polygon(
        new_extra_meander_cover_box_right_points,
        layer=material_meander.number,
        datatype=material_meander.datatype,
    )
    mask_builder.Main.add(extra_meander_cover_box_right)

    if mirror:
        new_extra_meander_cover_box_left_points = mbu.mirror_points_around_yaxis(extra_meander_cover_box_left_points)
        new_extra_meander_cover_box_left_points = mbu.rotate_and_move_points_list(new_extra_meander_cover_box_left_points, rot_angle, x, y)
    else:
        new_extra_meander_cover_box_left_points = mbu.rotate_and_move_points_list(extra_meander_cover_box_left_points, rot_angle, x, y)
    extra_meander_cover_box_left = gdspy.Polygon(
        new_extra_meander_cover_box_left_points,
        layer=material_meander.number,
        datatype=material_meander.datatype,
    )
    mask_builder.Main.add(extra_meander_cover_box_left)

    # Adding the frame left and frame right
    if mirror:
        new_frame_left_points = mbu.mirror_points_around_yaxis(frame_left_points)
        new_frame_left_points = mbu.rotate_and_move_points_list(new_frame_left_points, rot_angle, x, y)
    else:
        new_frame_left_points = mbu.rotate_and_move_points_list(frame_left_points, rot_angle, x, y)

    frame_left_poly = gdspy.Polygon(
        new_frame_left_points,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    mask_builder.Main.add(frame_left_poly)

    if mirror:
        new_frame_right_points = mbu.mirror_points_around_yaxis(frame_right_points)
        new_frame_right_points = mbu.rotate_and_move_points_list(new_frame_right_points, rot_angle, x, y)
    else:
        new_frame_right_points = mbu.rotate_and_move_points_list(frame_right_points, rot_angle, x, y)

    frame_right_poly = gdspy.Polygon(
        new_frame_right_points,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    mask_builder.Main.add(frame_right_poly)

    # Adding the coupler frame
    coupler_frame_points = [
        [coupler_frame_start_x - (coupler_frame_left_lw / 2), coupler_frame_start_y],
        [coupler_frame_start_x - (coupler_frame_left_lw / 2), coupler_frame_start_y + coupler_frame_left_height],
        [coupler_frame_start_x - (coupler_frame_left_lw / 2) + CCL, coupler_frame_start_y + coupler_frame_left_height],
        [
            coupler_frame_start_x - (coupler_frame_left_lw / 2) + CCL,
            coupler_frame_start_y + coupler_frame_left_height - coupler_frame_top_lw,
        ],
        [coupler_frame_start_x + (coupler_frame_left_lw / 2), coupler_frame_start_y + coupler_frame_left_height - coupler_frame_top_lw],
        [coupler_frame_start_x + (coupler_frame_left_lw / 2), coupler_frame_start_y],
    ]
    if mirror:
        new_coupler_frame_points = mbu.mirror_points_around_yaxis(coupler_frame_points)
        new_coupler_frame_points = mbu.rotate_and_move_points_list(new_coupler_frame_points, rot_angle, x, y)
    else:
        new_coupler_frame_points = mbu.rotate_and_move_points_list(coupler_frame_points, rot_angle, x, y)

    coupler_frame_poly = gdspy.Polygon(
        new_coupler_frame_points,
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    mask_builder.Main.add(coupler_frame_poly)

    # Adding the coupler arm
    coupler_arm_points = [
        [
            coupler_frame_start_x - (coupler_frame_left_lw / 2) - left_coupler_frame_to_feed_distance,
            coupler_frame_start_y + coupler_frame_left_height + coupler_gap,
        ],
        [coupler_frame_start_x - (coupler_frame_left_lw / 2) + CCL, coupler_frame_start_y + coupler_frame_left_height + coupler_gap],
        [
            coupler_frame_start_x - (coupler_frame_left_lw / 2) + CCL,
            coupler_frame_start_y + coupler_frame_left_height + coupler_gap + coupler_lw,
        ],
        [
            coupler_frame_start_x - (coupler_frame_left_lw / 2) - left_coupler_frame_to_feed_distance,
            coupler_frame_start_y + coupler_frame_left_height + coupler_gap + coupler_lw,
        ],
    ]
    if mirror:
        new_coupler_arm_points = mbu.mirror_points_around_yaxis(coupler_arm_points)
        new_coupler_arm_points = mbu.rotate_and_move_points_list(new_coupler_arm_points, rot_angle, x, y)
    else:
        new_coupler_arm_points = mbu.rotate_and_move_points_list(coupler_arm_points, rot_angle, x, y)

    coupler_arm_poly = gdspy.Polygon(
        new_coupler_arm_points,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(coupler_arm_poly)

    # Adding the IDC arms
    for i in range(0, No_of_arms, 2):
        right_arm = gdspy.Rectangle(
            [arm_start_x_right_side, arm_start_y_right_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            [arm_start_x_right_side - IDCLs[-(i + 1)], arm_start_y_right_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            layer=material_idc_and_frame.number,
            datatype=material_idc_and_frame.datatype,
        )
        right_arm.translate(x, y)
        if mirror:
            right_arm.mirror([x, y], [x, y + 10])
        right_arm.rotate(rot_angle, center=(x, y))
        mask_builder.Main.add(right_arm)

        left_arm = gdspy.Rectangle(
            [arm_start_x_left_side, arm_start_y_left_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            [arm_start_x_left_side + IDCLs[-(i + 2)], arm_start_y_left_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            layer=material_idc_and_frame.number,
            datatype=material_idc_and_frame.datatype,
        )
        left_arm.translate(x, y)
        if mirror:
            left_arm.mirror([x, y], [x, y + 10])
        left_arm.rotate(rot_angle, center=(x, y))
        mask_builder.Main.add(left_arm)

    # Adding the Trim arms.
    right_trim_arm = gdspy.Rectangle(
        [arm_start_x_right_side, trim_arm_start_y_right_side],
        [arm_start_x_right_side - trim_arm_length_right_side, trim_arm_start_y_right_side + trim_arm_lw],
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    right_trim_arm.translate(x, y)
    if mirror:
        right_trim_arm.mirror([x, y], [x, y + 10])
    right_trim_arm.rotate(rot_angle, center=(x, y))
    mask_builder.Main.add(right_trim_arm)

    left_trim_arm = gdspy.Rectangle(
        [arm_start_x_left_side, trim_arm_start_y_left_side],
        [arm_start_x_left_side + trim_arm_length_left_side, trim_arm_start_y_left_side + trim_arm_lw],
        layer=material_idc_and_frame.number,
        datatype=material_idc_and_frame.datatype,
    )
    left_trim_arm.translate(x, y)
    if mirror:
        left_trim_arm.mirror([x, y], [x, y + 10])
    left_trim_arm.rotate(rot_angle, center=(x, y))
    mask_builder.Main.add(left_trim_arm)

    # Adding the Trim boxes for the trim arms at Trim lengths specified if non-zero.
    if trim_length != None:
        if (
            trim_length < trim_arm_length_right_side and trim_length < trim_arm_length_left_side
        ):  # does not add a trim box if the trim length is longer than full length arm.
            inner_width = arm_start_x_right_side - arm_start_x_left_side

            trim_box_width = 3 * trim_arm_lw  # 6
            trim_box_length_overhang = (inner_width - trim_arm_length_right_side) / 2
            trim_box_for_right_trim_arm = gdspy.Rectangle(
                [
                    arm_start_x_right_side - trim_arm_length_right_side - trim_box_length_overhang,
                    trim_arm_start_y_right_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_left_side + (inner_width - trim_length),
                    trim_arm_start_y_right_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
                layer=mask_builder.layers.TrimLayer.number,
                datatype=mask_builder.layers.TrimLayer.datatype,
            )
            trim_box_for_right_trim_arm.translate(x, y)
            if mirror:
                trim_box_for_right_trim_arm.mirror([x, y], [x, y + 10])
            trim_box_for_right_trim_arm.rotate(rot_angle, center=(x, y))
            mask_builder.Main.add(trim_box_for_right_trim_arm)

            trim_box_for_left_trim_arm = gdspy.Rectangle(
                [
                    arm_start_x_left_side + trim_arm_length_left_side + trim_box_length_overhang,
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_right_side - (inner_width - trim_length),
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
                layer=mask_builder.layers.TrimLayer.number,
                datatype=mask_builder.layers.TrimLayer.datatype,
            )
            trim_box_for_left_trim_arm.translate(x, y)
            if mirror:
                trim_box_for_left_trim_arm.mirror([x, y], [x, y + 10])
            trim_box_for_left_trim_arm.rotate(rot_angle, center=(x, y))
            mask_builder.Main.add(trim_box_for_left_trim_arm)

    # Adding the cutout to the groundplane.
    if add_grnd_cutout:
        if mirror:
            new_grnd_plane_meander_cutout_poly_points = mbu.mirror_points_around_yaxis(grnd_plane_meander_cutout_poly_points)
            new_grnd_plane_meander_cutout_poly_points = mbu.rotate_and_move_points_list(
                new_grnd_plane_meander_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_grnd_plane_meander_cutout_poly_points = mbu.rotate_and_move_points_list(
                grnd_plane_meander_cutout_poly_points, rot_angle, x, y
            )

        grnd_plane_meander_cutout_poly = gdspy.Polygon(
            new_grnd_plane_meander_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_cutouts.add(grnd_plane_meander_cutout_poly)

    # Adding the cutout to the Silicon DiOxide membrane.
    if add_SiO_cutout:
        if mirror:
            new_SiO_cutout_poly_points = mbu.mirror_points_around_yaxis(SiO_cutout_poly_points)
            new_SiO_cutout_poly_points = mbu.rotate_and_move_points_list(new_SiO_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiO_cutout_poly_points = mbu.rotate_and_move_points_list(SiO_cutout_poly_points, rot_angle, x, y)

        SiO_cutout_poly = gdspy.Polygon(
            new_SiO_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.silicon_oxide_cutouts.add(SiO_cutout_poly)

    # Adding the cutout to the Silicon Nitride membrane.
    if add_SiN_membrane_cutout:
        if mirror:
            new_SiN_membrane_cutout_poly_points = mbu.mirror_points_around_yaxis(SiN_membrane_cutout_poly_points)
            new_SiN_membrane_cutout_poly_points = mbu.rotate_and_move_points_list(new_SiN_membrane_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiN_membrane_cutout_poly_points = mbu.rotate_and_move_points_list(SiN_membrane_cutout_poly_points, rot_angle, x, y)

        SiN_membrane_cutout_poly = gdspy.Polygon(
            new_SiN_membrane_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.silicon_nitride_membrane_cutouts.add(SiN_membrane_cutout_poly)

    # Adding the cutout to the SiN Dep layer.
    if add_SiN_dep_dielectric_cutout:
        if mirror:
            new_SiN_dep_cutout_poly_points = mbu.mirror_points_around_yaxis(SiN_dep_cutout_poly_points)
            new_SiN_dep_cutout_poly_points = mbu.rotate_and_move_points_list(new_SiN_dep_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiN_dep_cutout_poly_points = mbu.rotate_and_move_points_list(SiN_dep_cutout_poly_points, rot_angle, x, y)

        SiN_dep_cutout_poly = gdspy.Polygon(
            new_SiN_dep_cutout_poly_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(SiN_dep_cutout_poly)

    # Adding the backside check covers.
    if add_backside_check:
        if mirror:
            new_backside_check_cover_poly_points = mbu.mirror_points_around_yaxis(backside_check_cover_poly_points)
            new_backside_check_cover_poly_points = mbu.rotate_and_move_points_list(new_backside_check_cover_poly_points, rot_angle, x, y)
        else:
            new_backside_check_cover_poly_points = mbu.rotate_and_move_points_list(backside_check_cover_poly_points, rot_angle, x, y)

        backside_check_cover_poly = gdspy.Polygon(
            new_backside_check_cover_poly_points,
            layer=mask_builder.layers.Backside_Check.number,
            datatype=mask_builder.layers.Backside_Check.datatype,
        )
        mask_builder.Main.add(backside_check_cover_poly)

    # Adding the groundplane cutout over the inductive meander.
    if add_grnd_cutout_over_inductor:
        if mirror:
            new_grnd_plane_inductor_cutout_poly_points = mbu.mirror_points_around_yaxis(grnd_plane_inductor_cutout_poly_points)
            new_grnd_plane_inductor_cutout_poly_points = mbu.rotate_and_move_points_list(
                new_grnd_plane_inductor_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_grnd_plane_inductor_cutout_poly_points = mbu.rotate_and_move_points_list(
                grnd_plane_inductor_cutout_poly_points, rot_angle, x, y
            )

        grnd_plane_inductor_cutout_poly = gdspy.Polygon(
            new_grnd_plane_inductor_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_cutouts.add(grnd_plane_inductor_cutout_poly)

    # Adding the SiNdep cutout over the inductive meander.
    if add_SiN_dep_dielectric_cutout_over_inductor:
        if mirror:
            new_SiN_dep_inductor_cutout_poly_points = mbu.mirror_points_around_yaxis(SiN_dep_inductor_cutout_poly_points)
            new_SiN_dep_inductor_cutout_poly_points = mbu.rotate_and_move_points_list(
                new_SiN_dep_inductor_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_SiN_dep_inductor_cutout_poly_points = mbu.rotate_and_move_points_list(SiN_dep_inductor_cutout_poly_points, rot_angle, x, y)

        SiN_dep_inductor_cutout_poly = gdspy.Polygon(
            new_SiN_dep_inductor_cutout_poly_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(SiN_dep_inductor_cutout_poly)

    # # Adding the aluminium patch and etch under the meander
    if material_meander == mask_builder.layers.Aluminium:
        # Adding the aluminium patch
        # new_aluminium_patch_under_meander_points = handle_moving_points( aluminium_patch_under_meander_points, rot_angle, x, y, mirror)
        if mirror:
            new_aluminium_patch_under_meander_points = mbu.mirror_points_around_yaxis(aluminium_patch_under_meander_points)
            new_aluminium_patch_under_meander_points = mbu.rotate_and_move_points_list(
                new_aluminium_patch_under_meander_points, rot_angle, x, y
            )
        else:
            new_aluminium_patch_under_meander_points = mbu.rotate_and_move_points_list(
                aluminium_patch_under_meander_points, rot_angle, x, y
            )
        aluminium_patch_under_meander_polygon = gdspy.Polygon(
            new_aluminium_patch_under_meander_points,
            layer=mask_builder.layers.Aluminium_Patch.number,
            datatype=mask_builder.layers.Aluminium_Patch.datatype,
        )
        mask_builder.Main.add(aluminium_patch_under_meander_polygon)

        # Adding the aluminium etch
        # new_aluminium_etch_box_under_meander_points = handle_moving_points( aluminium_etch_under_meander_points, rot_angle, x, y, mirror)
        if mirror:
            new_aluminium_etch_box_under_meander_points = mbu.mirror_points_around_yaxis(aluminium_etch_under_meander_points)
            new_aluminium_etch_box_under_meander_points = mbu.rotate_and_move_points_list(
                new_aluminium_etch_box_under_meander_points, rot_angle, x, y
            )
        else:
            new_aluminium_etch_box_under_meander_points = mbu.rotate_and_move_points_list(
                aluminium_etch_under_meander_points, rot_angle, x, y
            )
        aluminium_etch_box_under_meander_polygon = gdspy.Polygon(
            new_aluminium_etch_box_under_meander_points,
            layer=mask_builder.layers.Aluminium_Etch.number,
            datatype=mask_builder.layers.Aluminium_Etch.datatype,
        )
        mask_builder.aluminium_etch_positives.add(aluminium_etch_box_under_meander_polygon)
        # mask_builder.Main.add(aluminium_etch_box_under_meander_polygon)

    # Adding the aluminium patch and etch under the IDC
    if material_idc_and_frame == mask_builder.layers.Aluminium:
        # Adding the aluminium patch
        # new_aluminium_patch_under_IDC_points = handle_moving_points(aluminium_patch_under_IDC_points, rot_angle, x, y, mirror)
        if mirror:
            new_aluminium_patch_under_IDC_points = mbu.mirror_points_around_yaxis(aluminium_patch_under_IDC_points)
            new_aluminium_patch_under_IDC_points = mbu.rotate_and_move_points_list(new_aluminium_patch_under_IDC_points, rot_angle, x, y)
        else:
            new_aluminium_patch_under_IDC_points = mbu.rotate_and_move_points_list(aluminium_patch_under_IDC_points, rot_angle, x, y)
        aluminium_patch_under_IDC_polygon = gdspy.Polygon(
            new_aluminium_patch_under_IDC_points,
            layer=mask_builder.layers.Aluminium_Patch.number,
            datatype=mask_builder.layers.Aluminium_Patch.datatype,
        )
        mask_builder.Main.add(aluminium_patch_under_IDC_polygon)

        # Adding the aluminium etch
        # new_aluminium_etch_under_IDC_points = handle_moving_points(aluminium_etch_under_IDC_points, rot_angle, x, y, mirror)
        if mirror:
            new_aluminium_etch_under_IDC_points = mbu.mirror_points_around_yaxis(aluminium_etch_under_IDC_points)
            new_aluminium_etch_under_IDC_points = mbu.rotate_and_move_points_list(new_aluminium_etch_under_IDC_points, rot_angle, x, y)
        else:
            new_aluminium_etch_under_IDC_points = mbu.rotate_and_move_points_list(aluminium_etch_under_IDC_points, rot_angle, x, y)
        aluminium_etch_under_IDC_polygon = gdspy.Polygon(
            new_aluminium_etch_under_IDC_points,
            layer=mask_builder.layers.Aluminium_Etch.number,
            datatype=mask_builder.layers.Aluminium_Etch.datatype,
        )
        mask_builder.aluminium_etch_positives.add(aluminium_etch_under_IDC_polygon)
        # mask_builder.Main.add(aluminium_etch_under_IDC_polygon)

    # Return the configurator_points if specified else just return.
    if not return_configurator_points:
        return

    configurator_points = {}

    ####################################################################### Meander
    configurator_points["meander_lw"] = {
        "text": "meander_lw",
        "start": [
            (new_meander_path_points[5][0] + new_meander_path_points[6][0]) / 2,
            new_meander_path_points[5][1] - (meander_lw / 2),
        ],
        "end": [(new_meander_path_points[5][0] + new_meander_path_points[6][0]) / 2, new_meander_path_points[5][1] + (meander_lw / 2)],
    }

    return configurator_points
