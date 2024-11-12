from collections.abc import Sequence
from typing import Callable, Literal

import gdspy
from numpy import cos, pi, sin

from ... import mask_builder_utils as mbu
from ...logging import TextColor, pretty_print
from ...mask_builder_utils import Layer
from ...souk_mask_components import add_cpw, add_microstrip_to_cpw_transition
from ...souk_muxing import get_mux_func_for_resonator_type
from ..resonator_types import SoukResonatorType
from ..utils.get_config import get_resonator_config


# from ...souk_mask_builder import SoukMaskBuilder  # isort:skip
def draw(
    # mask_builder: SoukMaskBuilder,
    mask_builder,
    resonator_type: SoukResonatorType,
    x: float,
    y: float,
    rot_angle: float,
    f0: float,
    mux_func_override: Callable | None = None,
    resonator_config_override: dict[str, float | int] | None = None,
    mirror: bool = False,
    IDC_and_frame_material: str = "IDC_Nb",
    meander_material: str = "Al",
    trim_length=None,
    add_grnd_cutout: bool = True,
    add_SiN_dep_dielectric_cutout: bool = True,
    add_SiO_cutout: bool = True,
    add_SiN_membrane_cutout: bool = True,
    add_backside_check: bool = False,
    add_grnd_cutout_over_inductor: bool = False,
    add_SiN_dep_dielectric_cutout_over_inductor: bool = False,
    return_configurator_points: bool = False,
) -> None | dict[str, dict]:
    """Adds the KID geometry to the Main cell at the x,y cooardinate
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
        - SoukResonatorType.CPW_COUPLED_V1

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
        that the IDC_and_CC_function function takes.

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

    config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

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

    mirror=False
        Whether the KID should be mirrored about the center vertical, **this
        mirroring is done before any roation is applied**. By default
        (with 0° ratation) the KID's coupler is attached on the left but when
        mirror=True the coupler on the right.

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

    add_Aluminium_Patch_and_Etch=True
        Whether of not to add an Aluminium patch and etch around aluminium
        elements.

    return_configurator_points=False
        return a the points for use in the configurator.
    """

    zero_marker = gdspy.Polygon(
        [
            [x - 5, y - 5],
            [x - 5, y],
            [x, y],
            [x, y + 5],
            [x + 5, y + 5],
            [x + 5, y],
            [x, y],
            [x, y - 5],
        ],
        layer=mask_builder.layers.General_labeling.number,
        datatype=mask_builder.layers.General_labeling.datatype,
    )
    mask_builder.Main.add(zero_marker)

    material_lookup = {
        "IDC_Nb": mask_builder.layers.IDC_Nb,
        "Nb": mask_builder.layers.Nb_Antenna,
        "Al": mask_builder.layers.Aluminium,
        "Nb_Groundplane": mask_builder.layers.Nb_Groundplane,
    }

    idc_and_frame_layer = material_lookup.get(IDC_and_frame_material, None)
    if idc_and_frame_layer is None:
        pretty_print(f"Could not find layer for `{IDC_and_frame_material}`. Defaulting to `IDC_Nb`.", TextColor.WARNING)
        idc_and_frame_layer = mask_builder.layers.IDC_Nb

    meander_layer = material_lookup.get(meander_material, None)
    if meander_layer is None:
        pretty_print(f"Could not find layer for `{IDC_and_frame_material}`. Defaulting to `Aluminium`.", TextColor.WARNING)
        meander_layer = mask_builder.layers.Aluminium

    resonator_config = get_resonator_config(resonator_type, resonator_config_override=resonator_config_override)

    if mux_func_override is None:
        IDC_and_CC_function = get_mux_func_for_resonator_type(resonator_type)
    else:
        IDC_and_CC_function = mux_func_override

    # Making the cpw meander section
    meander_cpw_center_width = resonator_config["meander_cpw_center_width"]
    meander_cpw_outer_width = resonator_config["meander_cpw_outer_width"]
    meander_dielectric_under_cpw_width = resonator_config["meander_dielectric_under_cpw_width"]
    meander_bend_radius = resonator_config["meander_bend_radius"]
    meander_cpw_bridge_width = resonator_config["meander_cpw_bridge_width"]

    meander_middle_base_x = 0
    meander_middle_base_y = meander_cpw_outer_width / 2

    meander_cpw_points_left = [[meander_middle_base_x, meander_middle_base_y]]
    meander_cpw_points_right = [[meander_middle_base_x, meander_middle_base_y]]

    number_of_right_meander_sections = 4
    number_of_left_meander_sections = 4

    for section_no in range(1, number_of_left_meander_sections + 1):
        # Adding the horizontal section
        previous_corner_xy = meander_cpw_points_left[-1]
        corner_x = previous_corner_xy[0] + resonator_config[f"meander_left_width_{section_no}"]
        corner_y = previous_corner_xy[1]
        meander_cpw_points_left.append([corner_x, corner_y])

        # Adding the vertical section
        previous_corner_xy = meander_cpw_points_left[-1]
        corner_x = previous_corner_xy[0]
        corner_y = previous_corner_xy[1] + resonator_config[f"meander_left_height_{section_no}"]
        meander_cpw_points_left.append([corner_x, corner_y])

    for section_no in range(1, number_of_right_meander_sections + 1):
        # Adding the horizontal section
        previous_corner_xy = meander_cpw_points_right[-1]
        corner_x = previous_corner_xy[0] + resonator_config[f"meander_right_width_{section_no}"]
        corner_y = previous_corner_xy[1]
        meander_cpw_points_right.append([corner_x, corner_y])

        # Adding the vertical section
        previous_corner_xy = meander_cpw_points_right[-1]
        corner_x = previous_corner_xy[0]
        corner_y = previous_corner_xy[1] + resonator_config[f"meander_right_height_{section_no}"]
        meander_cpw_points_right.append([corner_x, corner_y])

    # Making the meander ant pad
    ant_pad_box_width = resonator_config["ant_pad_box_width"]
    ant_pad_box_height = resonator_config["ant_pad_box_height"]

    ant_pad_box_points = [
        [meander_middle_base_x - (ant_pad_box_width / 2), meander_middle_base_y],
        [-(ant_pad_box_width / 2), meander_middle_base_y - ant_pad_box_height],
        [(ant_pad_box_width / 2), meander_middle_base_y - ant_pad_box_height],
        [meander_middle_base_x + (ant_pad_box_width / 2), meander_middle_base_y],
    ]

    ant_cpw_center_width = resonator_config["ant_cpw_center_width"]
    ant_cpw_outer_width = resonator_config["ant_cpw_outer_width"]
    ant_cpw_to_transition_offset_y = resonator_config["ant_cpw_to_transition_offset_y"]
    ant_cpw_to_transition_offset_x = resonator_config["ant_cpw_to_transition_offset_x"]
    ant_cpw_to_transition_angle_diff = resonator_config["ant_cpw_to_transition_angle_diff"]

    ant_microstrip_lw = resonator_config["ant_microstrip_lw"]
    ant_microstrip_final_x_offset = resonator_config["ant_microstrip_final_x_offset"]
    ant_microstrip_final_y_offset = resonator_config["ant_microstrip_final_y_offset"]

    ant_cpw_to_transition_path_points = [
        [meander_middle_base_x, meander_middle_base_y],
        [meander_middle_base_x, meander_middle_base_y - ant_cpw_to_transition_offset_y],
        [
            meander_middle_base_x - ant_cpw_to_transition_offset_x * cos(ant_cpw_to_transition_angle_diff),
            meander_middle_base_y - ant_cpw_to_transition_offset_y - ant_cpw_to_transition_offset_x * sin(ant_cpw_to_transition_angle_diff),
        ],
    ]

    # Making the frame sections left and right
    frame_bot_lw = resonator_config["frame_bot_lw"]

    frame_bot_left_width = resonator_config["frame_bot_left_width"]
    frame_bot_right_width = resonator_config["frame_bot_right_width"]

    frame_left_lw = resonator_config["frame_left_lw"]
    frame_left_height = resonator_config["frame_left_height"]

    frame_right_lw = resonator_config["frame_right_lw"]
    frame_right_height = resonator_config["frame_right_height"]

    frame_start_y = meander_cpw_points_left[-1][1]
    frame_left_start_x = meander_cpw_points_left[-1][0]
    frame_right_start_x = meander_cpw_points_right[-1][0]

    frame_left_points = [
        [frame_left_start_x, frame_start_y + frame_bot_lw],
        [frame_left_start_x - frame_bot_left_width + frame_left_lw, frame_start_y + frame_bot_lw],
        [frame_left_start_x - frame_bot_left_width + frame_left_lw, frame_start_y + frame_left_height],
        [frame_left_start_x - frame_bot_left_width, frame_start_y + frame_left_height],
        [frame_left_start_x - frame_bot_left_width, frame_start_y],
        [frame_left_start_x, frame_start_y],
    ]

    frame_right_points = [
        [frame_right_start_x, frame_start_y + frame_bot_lw],
        [frame_right_start_x + frame_bot_right_width - frame_right_lw, frame_start_y + frame_bot_lw],
        [frame_right_start_x + frame_bot_right_width - frame_right_lw, frame_start_y + frame_right_height],
        [frame_right_start_x + frame_bot_right_width, frame_start_y + frame_right_height],
        [frame_right_start_x + frame_bot_right_width, frame_start_y],
        [frame_right_start_x, frame_start_y],
    ]

    # Making the pads the conect the meander to the frame
    frame_meander_cover_box_width = resonator_config["frame_meander_cover_box_width"]
    frame_meander_cover_box_height = resonator_config["frame_meander_cover_box_height"]

    meander_cover_box_left_points = [
        [meander_cpw_points_left[-1][0] - (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
        [
            meander_cpw_points_left[-1][0] - (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [
            meander_cpw_points_left[-1][0] + (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [meander_cpw_points_left[-1][0] + (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
    ]

    meander_cover_box_right_points = [
        [meander_cpw_points_right[-1][0] - (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
        [
            meander_cpw_points_right[-1][0] - (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [
            meander_cpw_points_right[-1][0] + (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [meander_cpw_points_right[-1][0] + (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
    ]

    meander_cover_box_right_end_center = [
        (meander_cover_box_right_points[1][0] + meander_cover_box_right_points[2][0]) / 2,
        (meander_cover_box_right_points[1][1] + meander_cover_box_right_points[2][1]) / 2,
    ]

    meander_cover_box_left_end_center = [
        (meander_cover_box_left_points[1][0] + meander_cover_box_left_points[2][0]) / 2,
        (meander_cover_box_left_points[1][1] + meander_cover_box_left_points[2][1]) / 2,
    ]

    # Making the coupler attachement
    coupler_frame_left_lw = resonator_config["coupler_frame_left_lw"]  # 10
    coupler_frame_left_height = resonator_config["coupler_frame_left_height"]  # 39
    coupler_frame_top_lw = resonator_config["coupler_frame_top_lw"]  # 3

    coupler_frame_start_x = frame_left_start_x - frame_bot_left_width + (frame_left_lw / 2)
    coupler_frame_start_y = frame_start_y + frame_left_height

    # Making the IDC and trim arms
    IDC_bot_arm_gap = resonator_config["IDC_bot_arm_gap"]  # 30
    IDC_arm_gap = resonator_config["IDC_arm_gap"]  # 8
    IDC_arm_lw = resonator_config["IDC_arm_lw"]  # 3
    No_of_arms = resonator_config["No_of_arms"]  # 28

    arm_start_x_left_side = frame_left_start_x - frame_bot_left_width + frame_left_lw
    arm_start_x_right_side = frame_right_start_x + frame_bot_right_width - frame_right_lw

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

    cutout_start_height = (
        (meander_cpw_outer_width / 2)
        + resonator_config["meander_right_height_1"]
        + resonator_config["meander_right_height_2"]
        + resonator_config["meander_right_height_3"]
        + resonator_config["meander_right_height_4"]
        - resonator_config["cutout_bot_offset"]
    )
    cutout_width = (frame_right_start_x + frame_bot_right_width + cutout_right_offset) - (
        frame_left_start_x - frame_bot_left_width - cutout_left_offset
    )
    cutout_height = (coupler_frame_start_y + coupler_frame_left_height + coupler_gap + coupler_lw + cutout_top_offset) - cutout_start_height

    step_down_distance_between_layers = resonator_config["step_down_distance_between_layers"]  # 5

    grnd_plane_cutout_width = cutout_width + (3 * (2 * step_down_distance_between_layers))
    grnd_plane_cutout_height = cutout_height + (3 * (2 * step_down_distance_between_layers))
    grnd_plane_cutout_start_height = cutout_start_height - (3 * step_down_distance_between_layers)
    # grnd_plane_cutout_start_height = cutout_start_height  # - (3 * step_down_distance_between_layers)

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
            meander_cpw_points_right[1][0] + grnd_plane_inductor_cutout_offset_right,
            meander_cpw_points_right[1][1] - grnd_plane_inductor_cutout_offset_bot,
        ],
        [
            meander_cpw_points_left[1][0] - grnd_plane_inductor_cutout_offset_left,
            meander_cpw_points_left[1][1] - grnd_plane_inductor_cutout_offset_bot,
        ],
        [
            meander_cpw_points_left[6][0] - grnd_plane_inductor_cutout_offset_left,
            meander_cpw_points_left[6][1] + grnd_plane_inductor_cutout_offset_top,
        ],
        [
            meander_cpw_points_right[6][0] + grnd_plane_inductor_cutout_offset_right,
            meander_cpw_points_right[6][1] + grnd_plane_inductor_cutout_offset_top,
        ],
    ]

    SiN_dep_inductor_cutout_offset_right = 20
    SiN_dep_inductor_cutout_offset_left = 20
    SiN_dep_inductor_cutout_offset_top = 20
    SiN_dep_inductor_cutout_offset_bot = 20

    SiN_dep_inductor_cutout_poly_points = [
        # bot right, bot left, top left, top right
        [
            meander_cpw_points_right[1][0] + SiN_dep_inductor_cutout_offset_right,
            meander_cpw_points_right[1][1] - SiN_dep_inductor_cutout_offset_bot,
        ],
        [
            meander_cpw_points_left[1][0] - SiN_dep_inductor_cutout_offset_left,
            meander_cpw_points_left[1][1] - SiN_dep_inductor_cutout_offset_bot,
        ],
        [
            meander_cpw_points_left[6][0] - SiN_dep_inductor_cutout_offset_left,
            meander_cpw_points_left[6][1] + SiN_dep_inductor_cutout_offset_top,
        ],
        [
            meander_cpw_points_right[6][0] + SiN_dep_inductor_cutout_offset_right,
            meander_cpw_points_right[6][1] + SiN_dep_inductor_cutout_offset_top,
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

    aluminium_patch_meander_width_left = min([meander_cpw_points_left[i][0] for i in range(len(meander_cpw_points_left))]) - (
        meander_cpw_center_width / 2
    )
    aluminium_patch_meander_width_right = max([meander_cpw_points_right[i][0] for i in range(len(meander_cpw_points_right))]) + (
        meander_cpw_center_width / 2
    )
    aluminium_patch_meander_height_above = max(meander_cpw_points_left[-1][1], meander_cpw_points_right[-1][1])
    aluminium_patch_meander_height_below = -ant_pad_box_height - aluminium_patch_offset_bot

    # defined from bot left going anti-clockwise
    aluminium_patch_under_meander_points = [
        [aluminium_patch_meander_width_left, aluminium_patch_meander_height_below],
        [aluminium_patch_meander_width_right, aluminium_patch_meander_height_below],
        [aluminium_patch_meander_width_right, aluminium_patch_meander_height_above],
        [aluminium_patch_meander_width_left, aluminium_patch_meander_height_above],
    ]

    aluminium_etch_under_meander_points = [
        [
            aluminium_patch_under_meander_points[0][0] - aluminium_etch_offset_left,
            aluminium_patch_under_meander_points[0][1] - aluminium_etch_offset_bot,
        ],
        [
            aluminium_patch_under_meander_points[1][0] + aluminium_etch_offset_right,
            aluminium_patch_under_meander_points[1][1] - aluminium_etch_offset_bot,
        ],
        [
            aluminium_patch_under_meander_points[2][0] + aluminium_etch_offset_right,
            aluminium_patch_under_meander_points[2][1] + aluminium_etch_offset_top,
        ],
        [
            aluminium_patch_under_meander_points[3][0] - aluminium_etch_offset_left,
            aluminium_patch_under_meander_points[3][1] + aluminium_etch_offset_top,
        ],
    ]

    # Making the aluminium patch under the IDC
    aluminium_patch_IDC_left_edge = (
        frame_left_points[-2][0] + (frame_left_lw / 2) - (coupler_frame_left_lw / 2) - aluminium_patch_offset_left
    )
    aluminium_patch_IDC_right_edge = frame_right_points[-2][0] + aluminium_patch_offset_right
    aluminium_patch_IDC_top_edge = coupler_frame_start_y + coupler_frame_left_height
    aluminium_patch_IDC_bot_edge = frame_start_y + frame_bot_lw - frame_meander_cover_box_height

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

    ###########################################################################
    # Adding eveything to mask
    ###########################################################################
    # Adjusting the x, y to place the resonator at the end of the microstip connection point.
    temp_sign = -1 if mirror else +1

    # dx = temp_sign * (ant_microstrip_final_x_offset * cos(rot_angle))
    # dy = ant_microstrip_final_y_offset * (cos(rot_angle))

    dx = temp_sign * ant_microstrip_final_x_offset
    dy = ant_microstrip_final_y_offset

    # x = x + dx
    # y = y + dy

    # x = x + ant_microstrip_final_x_offset
    # y = y + ant_microstrip_final_y_offset

    def handle_moving_points(
        points: Sequence[Sequence[float | int]],
        offset_x: float | int,
        offset_y: float | int,
        angle: float | int,
        final_x_pos: float | int,
        final_y_pos: float | int,
        should_mirror: bool,
    ):
        """This will handle moving the points to thier final location and
        rotation to be placed on the mask.

         - This will first mirror the points across the y axis if needed.
         - Then it will move the points by an offset to accomidate for the cpw
         to meander transition.
         - Then it will rotate and then move the points to thier final
         lcoation on the mask.
        """
        if should_mirror:
            new_points = mbu.mirror_points_around_yaxis(points)
            new_points = mbu.move_points_list(new_points, offset_x, offset_y)
        else:
            new_points = mbu.move_points_list(points, offset_x, offset_y)

        new_points = mbu.rotate_and_move_points_list(
            new_points,
            angle,
            final_x_pos,
            final_y_pos,
        )

        return new_points

    marker_points = [
        [0 - 5, 0 - 5],
        [0 - 5, 0],
        [0, 0],
        [0, 0 + 5],
        [0 + 5, 0 + 5],
        [0 + 5, 0],
        [0, 0],
        [0, 0 - 5],
    ]
    base_marker_points = handle_moving_points(marker_points, 0, 0, rot_angle, x, y, mirror)
    zero_marker = gdspy.Polygon(
        base_marker_points,
        layer=mask_builder.layers.General_labeling.number,
        datatype=mask_builder.layers.General_labeling.datatype,
    )
    mask_builder.Main.add(zero_marker)

    shifted_marker_points = handle_moving_points(marker_points, dx, dy, rot_angle, x, y, mirror)
    zero_marker = gdspy.Polygon(
        shifted_marker_points,
        layer=mask_builder.layers.General_labeling.number,
        datatype=mask_builder.layers.General_labeling.datatype,
    )
    mask_builder.Main.add(zero_marker)

    # print(f"    x,     y = {x}, {y}")
    # print(f"new_x, new_y = {shifted_marker_points[2][0]}, {shifted_marker_points[2][1]}")

    # Adding the cpw meander
    new_meander_cpw_points_left = handle_moving_points(meander_cpw_points_left, dx, dy, rot_angle, x, y, mirror)
    new_meander_cpw_points_right = handle_moving_points(meander_cpw_points_right, dx, dy, rot_angle, x, y, mirror)

    def _get_bridge_lengths(side: Literal["left", "right"]) -> list[float]:
        """get the mid way between the height of the leftmost/rightmost
        upright paths of the meander"""

        lengths_to_add_bridges: list[float] = []

        corner_correction = -meander_bend_radius * (2 - (pi / 2))

        # 0th bridge at start of path
        lengths_to_add_bridges.append(ant_cpw_outer_width + meander_cpw_bridge_width)
        lengths_to_add_bridges.append(abs(resonator_config[f"meander_{side}_width_1"]) * (1 / 3))
        lengths_to_add_bridges.append(abs(resonator_config[f"meander_{side}_width_1"]) * (2 / 3))

        # Adding a bridge in the middle of each of the vertical
        for no_of_vertical_sections in range(1, 3 + 1):
            len_into_path = 0.0
            no_of_corners = (2 * no_of_vertical_sections) - 1

            for no in range(1, no_of_vertical_sections + 1):
                len_into_path += abs(resonator_config[f"meander_{side}_width_{no}"]) + abs(resonator_config[f"meander_{side}_height_{no}"])
            len_into_path -= abs(resonator_config[f"meander_{side}_height_{no_of_vertical_sections}"] / 2)
            len_into_path += no_of_corners * corner_correction

            # getting the bridge position in the middle of the vertical
            lengths_to_add_bridges.append(len_into_path)

            # Going back and filling in two more bridges inbetween every vertical
            len_into_path += abs(resonator_config[f"meander_{side}_height_{no_of_vertical_sections}"] / 2)
            len_into_path += corner_correction

            len_into_path += abs(resonator_config[f"meander_{side}_width_{no_of_vertical_sections}"]) * (1 / 3)
            lengths_to_add_bridges.append(len_into_path)
            len_into_path += abs(resonator_config[f"meander_{side}_width_{no_of_vertical_sections}"]) * (1 / 3)
            lengths_to_add_bridges.append(len_into_path)

        lengths_to_add_bridges.sort()
        return lengths_to_add_bridges

    bridge_lengths_left = _get_bridge_lengths("left")
    bridge_lengths_right = _get_bridge_lengths("right")

    add_cpw(
        mask_builder,
        new_meander_cpw_points_left,
        meander_cpw_center_width,
        meander_cpw_outer_width,
        center_material=meander_layer,
        dielectric_in_cutout_or_positive="cutout",
        center_corners_type="circular bend",
        center_bend_radius=meander_bend_radius,
        cutout_corners_type="circular bend",
        cutout_bend_radius=meander_bend_radius,
        dielectric_width=meander_dielectric_under_cpw_width,
        dielectric_corners_type="circular bend",
        dielectric_bend_radius=meander_bend_radius,
        cutout_bridge_gap=bridge_lengths_left,
        cutout_bridge_width=meander_cpw_bridge_width,
        dielectric_bridge_gap=bridge_lengths_left,
        dielectric_bridge_width=meander_cpw_bridge_width + 4,
    )

    add_cpw(
        mask_builder,
        new_meander_cpw_points_right,
        meander_cpw_center_width,
        meander_cpw_outer_width,
        center_material=meander_layer,
        dielectric_in_cutout_or_positive="cutout",
        center_corners_type="circular bend",
        center_bend_radius=meander_bend_radius,
        cutout_corners_type="circular bend",
        cutout_bend_radius=meander_bend_radius,
        dielectric_width=meander_dielectric_under_cpw_width,
        dielectric_corners_type="circular bend",
        dielectric_bend_radius=meander_bend_radius,
        cutout_bridge_gap=bridge_lengths_right,
        cutout_bridge_width=meander_cpw_bridge_width,
        dielectric_bridge_gap=bridge_lengths_right,
        dielectric_bridge_width=meander_cpw_bridge_width + 4,
    )

    # Adding the meander ant overlap box
    new_ant_pad_box_points = handle_moving_points(ant_pad_box_points, dx, dy, rot_angle, x, y, mirror)

    ant_pad_box = gdspy.Polygon(
        new_ant_pad_box_points,
        layer=meander_layer.number,
        datatype=meander_layer.datatype,
    )
    mask_builder.Main.add(ant_pad_box)

    # Adding the ant to cpw transition
    new_ant_cpw_to_transition_path_points = handle_moving_points(ant_cpw_to_transition_path_points, dx, dy, rot_angle, x, y, mirror)
    ant_cpw_to_transition_path = gdspy.FlexPath(
        new_ant_cpw_to_transition_path_points,
        ant_cpw_center_width,
        corners="circular bend",
        bend_radius=meander_bend_radius,
        gdsii_path=True,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.make_flexpath_into_polygons_and_add_to_main(ant_cpw_to_transition_path, mask_builder.layers.Nb_Antenna)

    cutout_ant_cpw_to_transition_path = gdspy.FlexPath(
        new_ant_cpw_to_transition_path_points,
        ant_cpw_outer_width,
        corners="circular bend",
        bend_radius=meander_bend_radius,
        gdsii_path=True,
        layer=mask_builder.layers.Nb_Groundplane.number,
        datatype=mask_builder.layers.Nb_Groundplane.datatype,
    )
    for poly_points in mbu.get_polys_from_flexpath(cutout_ant_cpw_to_transition_path):
        mask_builder.ground_plane_cutouts.add(
            gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
        )

    if mirror:
        microstip_transition_angle = rot_angle + ant_cpw_to_transition_angle_diff
    else:
        microstip_transition_angle = rot_angle + ant_cpw_to_transition_angle_diff - pi

    add_microstrip_to_cpw_transition(
        mask_builder,
        ant_microstrip_lw,
        [
            new_ant_cpw_to_transition_path_points[-1][0],
            new_ant_cpw_to_transition_path_points[-1][1],
            microstip_transition_angle,
        ],
        ant_cpw_center_width,
        ant_cpw_outer_width,
        cpw_end_x_y_rot=[x, y, rot_angle - (pi / 2)],
    )

    # Adding the meander frame overlap boxes
    new_meander_cover_box_left_points = handle_moving_points(meander_cover_box_left_points, dx, dy, rot_angle, x, y, mirror)
    new_meander_cover_box_right_points = handle_moving_points(meander_cover_box_right_points, dx, dy, rot_angle, x, y, mirror)

    meander_cover_box_left = gdspy.Polygon(
        new_meander_cover_box_left_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
    )
    mask_builder.Main.add(meander_cover_box_left)
    meander_cover_box_right = gdspy.Polygon(
        new_meander_cover_box_right_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
    )
    mask_builder.Main.add(meander_cover_box_right)

    # Adding the frame left and frame right
    new_frame_left_points = handle_moving_points(frame_left_points, dx, dy, rot_angle, x, y, mirror)
    new_frame_right_points = handle_moving_points(frame_right_points, dx, dy, rot_angle, x, y, mirror)

    frame_left_poly = gdspy.Polygon(
        new_frame_left_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
    )
    mask_builder.Main.add(frame_left_poly)

    frame_right_poly = gdspy.Polygon(
        new_frame_right_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
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

    new_coupler_frame_points = handle_moving_points(coupler_frame_points, dx, dy, rot_angle, x, y, mirror)

    coupler_frame_poly = gdspy.Polygon(
        new_coupler_frame_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
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
    new_coupler_arm_points = handle_moving_points(coupler_arm_points, dx, dy, rot_angle, x, y, mirror)

    coupler_arm_poly = gdspy.Polygon(
        new_coupler_arm_points,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(coupler_arm_poly)

    # Adding the IDC arms
    for i in range(0, No_of_arms, 2):
        right_arm_points = [
            [
                arm_start_x_right_side,
                arm_start_y_right_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
            [
                arm_start_x_right_side - IDCLs[-(i + 1)],
                arm_start_y_right_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
            [
                arm_start_x_right_side - IDCLs[-(i + 1)],
                arm_start_y_right_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
            [
                arm_start_x_right_side,
                arm_start_y_right_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
        ]
        new_right_arm_points = handle_moving_points(right_arm_points, dx, dy, rot_angle, x, y, mirror)
        right_arm = gdspy.Polygon(
            new_right_arm_points,
            layer=idc_and_frame_layer.number,
            datatype=idc_and_frame_layer.datatype,
        )
        mask_builder.Main.add(right_arm)

        left_arm_points = [
            [
                arm_start_x_left_side,
                arm_start_y_left_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
            [
                arm_start_x_left_side + IDCLs[-(i + 2)],
                arm_start_y_left_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
            [
                arm_start_x_left_side + IDCLs[-(i + 2)],
                arm_start_y_left_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
            [
                arm_start_x_left_side,
                arm_start_y_left_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw)),
            ],
        ]
        new_left_arm_points = handle_moving_points(left_arm_points, dx, dy, rot_angle, x, y, mirror)
        left_arm = gdspy.Polygon(
            new_left_arm_points,
            layer=idc_and_frame_layer.number,
            datatype=idc_and_frame_layer.datatype,
        )
        mask_builder.Main.add(left_arm)

    # Adding the Trim arms.
    right_arm_trim_points = [
        [
            arm_start_x_right_side,
            trim_arm_start_y_right_side,
        ],
        [
            arm_start_x_right_side - trim_arm_length_right_side,
            trim_arm_start_y_right_side,
        ],
        [
            arm_start_x_right_side - trim_arm_length_right_side,
            trim_arm_start_y_right_side + trim_arm_lw,
        ],
        [
            arm_start_x_right_side,
            trim_arm_start_y_right_side + trim_arm_lw,
        ],
    ]
    new_right_arm_trim_points = handle_moving_points(right_arm_trim_points, dx, dy, rot_angle, x, y, mirror)
    right_trim_arm = gdspy.Polygon(
        new_right_arm_trim_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
    )
    mask_builder.Main.add(right_trim_arm)

    left_arm_trim_points = [
        [
            arm_start_x_left_side,
            trim_arm_start_y_left_side,
        ],
        [
            arm_start_x_left_side + trim_arm_length_left_side,
            trim_arm_start_y_left_side,
        ],
        [
            arm_start_x_left_side + trim_arm_length_left_side,
            trim_arm_start_y_left_side + trim_arm_lw,
        ],
        [
            arm_start_x_left_side,
            trim_arm_start_y_left_side + trim_arm_lw,
        ],
    ]

    new_left_arm_trim_points = handle_moving_points(left_arm_trim_points, dx, dy, rot_angle, x, y, mirror)
    left_trim_arm = gdspy.Polygon(
        new_left_arm_trim_points,
        layer=idc_and_frame_layer.number,
        datatype=idc_and_frame_layer.datatype,
    )
    mask_builder.Main.add(left_trim_arm)

    # Adding the Trim boxes for the trim arms at Trim lengths specified if non-zero.
    if trim_length != None:
        if (
            trim_length < trim_arm_length_right_side and trim_length < trim_arm_length_left_side
        ):  # does not add a trim box if the trim length is longer than full length arm.
            inner_width = arm_start_x_right_side - arm_start_x_left_side

            trim_box_width = 3 * trim_arm_lw  # 6
            trim_box_length_overhang = (inner_width - trim_arm_length_right_side) / 2
            trim_box_for_right_trim_arm_points = [
                [
                    arm_start_x_right_side - trim_arm_length_right_side - trim_box_length_overhang,
                    trim_arm_start_y_right_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_left_side + (inner_width - trim_length),
                    trim_arm_start_y_right_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_left_side + (inner_width - trim_length),
                    trim_arm_start_y_right_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
                [
                    arm_start_x_right_side - trim_arm_length_right_side - trim_box_length_overhang,
                    trim_arm_start_y_right_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
            ]
            new_trim_box_for_right_trim_arm_points = handle_moving_points(
                trim_box_for_right_trim_arm_points, dx, dy, rot_angle, x, y, mirror
            )
            trim_box_for_right_trim_arm = gdspy.Polygon(
                new_trim_box_for_right_trim_arm_points,
                layer=mask_builder.layers.TrimLayer.number,
                datatype=mask_builder.layers.TrimLayer.datatype,
            )
            mask_builder.Main.add(trim_box_for_right_trim_arm)

            trim_box_for_left_trim_arm_points = [
                [
                    arm_start_x_left_side + trim_arm_length_left_side + trim_box_length_overhang,
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_right_side - (inner_width - trim_length),
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_right_side - (inner_width - trim_length),
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
                [
                    arm_start_x_left_side + trim_arm_length_left_side + trim_box_length_overhang,
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
            ]
            new_trim_box_for_left_trim_arm_points = handle_moving_points(trim_box_for_left_trim_arm_points, dx, dy, rot_angle, x, y, mirror)
            trim_box_for_left_trim_arm = gdspy.Polygon(
                new_trim_box_for_left_trim_arm_points,
                layer=mask_builder.layers.TrimLayer.number,
                datatype=mask_builder.layers.TrimLayer.datatype,
            )
            mask_builder.Main.add(trim_box_for_left_trim_arm)

    # Adding the cutout to the groundplane.
    if add_grnd_cutout:
        new_grnd_plane_meander_cutout_poly_points = handle_moving_points(
            grnd_plane_meander_cutout_poly_points, dx, dy, rot_angle, x, y, mirror
        )
        grnd_plane_meander_cutout_poly = gdspy.Polygon(
            new_grnd_plane_meander_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_cutouts.add(grnd_plane_meander_cutout_poly)

    # Adding the cutout to the Silicon DiOxide membrane.
    if add_SiO_cutout:
        new_SiO_cutout_poly_points = handle_moving_points(SiO_cutout_poly_points, dx, dy, rot_angle, x, y, mirror)
        SiO_cutout_poly = gdspy.Polygon(
            new_SiO_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.silicon_oxide_cutouts.add(SiO_cutout_poly)

    # Adding the cutout to the Silicon Nitride membrane.
    if add_SiN_membrane_cutout:
        new_SiN_membrane_cutout_poly_points = handle_moving_points(SiN_membrane_cutout_poly_points, dx, dy, rot_angle, x, y, mirror)
        SiN_membrane_cutout_poly = gdspy.Polygon(
            new_SiN_membrane_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.silicon_nitride_membrane_cutouts.add(SiN_membrane_cutout_poly)

    # Adding the cutout to the SiN Dep layer.
    if add_SiN_dep_dielectric_cutout:
        new_SiN_dep_cutout_poly_points = handle_moving_points(SiN_dep_cutout_poly_points, dx, dy, rot_angle, x, y, mirror)
        SiN_dep_cutout_poly = gdspy.Polygon(
            new_SiN_dep_cutout_poly_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(SiN_dep_cutout_poly)

    # Adding the backside check covers.
    if add_backside_check:
        new_backside_check_cover_poly_points = handle_moving_points(backside_check_cover_poly_points, dx, dy, rot_angle, x, y, mirror)
        backside_check_cover_poly = gdspy.Polygon(
            new_backside_check_cover_poly_points,
            layer=mask_builder.layers.Backside_Check.number,
            datatype=mask_builder.layers.Backside_Check.datatype,
        )
        mask_builder.Main.add(backside_check_cover_poly)

    # Adding the groundplane cutout over the inductive meander.
    if add_grnd_cutout_over_inductor:
        new_grnd_plane_inductor_cutout_poly_points = handle_moving_points(
            grnd_plane_inductor_cutout_poly_points, dx, dy, rot_angle, x, y, mirror
        )
        grnd_plane_inductor_cutout_poly = gdspy.Polygon(
            new_grnd_plane_inductor_cutout_poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_cutouts.add(grnd_plane_inductor_cutout_poly)

    # Adding the SiNdep cutout over the inductive meander.
    if add_SiN_dep_dielectric_cutout_over_inductor:
        new_SiN_dep_inductor_cutout_poly_points = handle_moving_points(SiN_dep_inductor_cutout_poly_points, dx, dy, rot_angle, x, y, mirror)
        SiN_dep_inductor_cutout_poly = gdspy.Polygon(
            new_SiN_dep_inductor_cutout_poly_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(SiN_dep_inductor_cutout_poly)

    # Adding the meander Aluminium Patch and Etch
    if meander_layer == mask_builder.layers.Aluminium:
        # Adding the aluminium patch
        new_aluminium_patch_under_meander_points = handle_moving_points(
            aluminium_patch_under_meander_points,
            dx,
            dy,
            rot_angle,
            x,
            y,
            mirror,
        )
        aluminium_patch_under_meander_polygon = gdspy.Polygon(
            new_aluminium_patch_under_meander_points,
            layer=mask_builder.layers.Aluminium_Patch.number,
            datatype=mask_builder.layers.Aluminium_Patch.datatype,
        )
        mask_builder.Main.add(aluminium_patch_under_meander_polygon)

        # Adding the aluminium etch
        new_aluminium_etch_box_under_meander_points = handle_moving_points(
            aluminium_etch_under_meander_points,
            dx,
            dy,
            rot_angle,
            x,
            y,
            mirror,
        )
        aluminium_etch_box_under_meander_polygon = gdspy.Polygon(
            new_aluminium_etch_box_under_meander_points,
            layer=mask_builder.layers.Aluminium_Etch.number,
            datatype=mask_builder.layers.Aluminium_Etch.datatype,
        )
        mask_builder.aluminium_etch_positives.add(aluminium_etch_box_under_meander_polygon)
        # mask_builder.Main.add(aluminium_etch_box_under_meander_polygon)

    # Adding the IDC Aluminium Patch and Etch
    if idc_and_frame_layer == mask_builder.layers.Aluminium:
        # Adding the aluminium patch
        new_aluminium_patch_under_IDC_points = handle_moving_points(
            aluminium_patch_under_IDC_points,
            dx,
            dy,
            rot_angle,
            x,
            y,
            mirror,
        )
        aluminium_patch_under_IDC_polygon = gdspy.Polygon(
            new_aluminium_patch_under_IDC_points,
            layer=mask_builder.layers.Aluminium_Patch.number,
            datatype=mask_builder.layers.Aluminium_Patch.datatype,
        )
        mask_builder.Main.add(aluminium_patch_under_IDC_polygon)

        # Adding the aluminium etch
        new_aluminium_etch_box_under_IDC_points = handle_moving_points(
            aluminium_etch_under_IDC_points,
            dx,
            dy,
            rot_angle,
            x,
            y,
            mirror,
        )
        aluminium_etch_box_under_IDC_polygon = gdspy.Polygon(
            new_aluminium_etch_box_under_IDC_points,
            layer=mask_builder.layers.Aluminium_Etch.number,
            datatype=mask_builder.layers.Aluminium_Etch.datatype,
        )
        mask_builder.aluminium_etch_positives.add(aluminium_etch_box_under_IDC_polygon)
        # mask_builder.Main.add(aluminium_etch_box_under_IDC_polygon)

    if not return_configurator_points:
        return None

    # TODO configurator_points
    configurator_points = {}
    ####################################################################### Meander
    configurator_points["TODO"] = {
        "text": "TODO",
        "start": [0, 0],
        "end": [0, 0],
    }

    return configurator_points
