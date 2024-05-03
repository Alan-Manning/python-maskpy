from typing import Callable

from pandas import test

from ...souk_mask_builder import SoukMaskBuilder
from ..configs.config_original_q50k import get_resonator_config
from ..mux_funcs.mux_original_q50k import get_IDCLs_and_CCL_from_f0


def add_resonator(
    souk_mask_builder: SoukMaskBuilder,
    x: float,
    y: float,
    rot_angle: float,
    f0: float,
    IDC_and_CC_function: Callable = get_IDCLs_and_CCL_from_f0,
    Main_config_file_dict: dict = get_resonator_config(),
    mirror=False,
    IDC_and_frame_material="IDC_Nb",
    meander_material="Al",
    trim_length=None,
    add_grnd_cutout=True,
    add_SiN_dep_dielectric_cutout=True,
    add_SiO_cutout=True,
    add_SiN_membrane_cutout=True,
    add_backside_check=True,
    add_grnd_cutout_over_inductor=False,
    add_SiN_dep_dielectric_cutout_over_inductor=False,
    return_configurator_points=False,
):
    """Adds the KID geometry to the Main cell athe the x,y cooardinate
    given. The KID is placed where the base middle of the inductive meander
    is at this x,y. The KID geometry is defined by the dimensions within
    the Main_config_file_dict. By default it will, but optionally can
    choose not to, add all the neccessay cutouts for the structure.

    Parameters
    ----------
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

    IDC_and_CC_function : function
        function for getting the IDC and CC lengths from a given f0. The
        function should take a frequency and QR value string as arguments
        and should return an array-like (28 long) and a single float
        **in this order**, the array-like should contain all the lengths
        for each of the 28 arms for the IDC and the float should be the CC
        length. A simple example funtion:

        >>> def example_IDC_CC_func(f0, Q_value):
        ...     '''
        ...     Example function for IDC and CC
        ...     f0 : float
        ...         Resonant frequency in Hz.
        ...     QR : str
        ...         The QR value
        ...     '''
        ...     if (f0 < 3.1e9) and (Q_value=="50k"):
        ...         IDC_lengths = np.ones(28)*1900.0
        ...         CC_length = 600.0
        ...     else:
        ...         IDC_lengths = np.ones(28)*1500.0
        ...         CC_length = 300.0
        ...
        ...     # return IDC array, then CC seperately
        ...     return IDC_lengths, CC_length


    Main_config_file_dict : dict
        dictionary containing individual dictionarys of config settings.
        Requires "resonator".

    KwArgs
    ------
    Q_value = "50k"
        This is the QR value to use for the resonator in the
        IDC_and_CC_function. By default this is "50k" string.

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

    add_backside_check=True
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

    IDC_and_frame_material_lookup = {
        "IDC_Nb": souk_mask_builder.IDC_Nb,
        "Nb": souk_mask_builder.Nb_Antenna,
        "Al": souk_mask_builder.Aluminium,
    }
    material_idc_and_frame = IDC_and_frame_material_lookup[IDC_and_frame_material]

    meander_material_lookup = {"Al": souk_mask_builder.Aluminium, "IDC_Nb": souk_mask_builder.IDC_Nb, "Nb": souk_mask_builder.Nb_Antenna}
    material_meander = meander_material_lookup[meander_material]

    config = Main_config_file_dict["resonator"]

    # config = resonator_config_dict

    # Making the meander section
    meander_lw = config["meander_lw"]  # 2
    meander_corner_bend_radius = config["meander_corner_bend_radius"]  # 6

    meander_bot_width = config["meander_bot_width"]  # 18

    meander_right_height_1 = config["meander_right_height_1"]  # 41
    meander_right_height_2 = config["meander_right_height_2"]  # 24
    meander_right_height_3 = config["meander_right_height_3"]  # 53

    meander_left_height_1 = config["meander_left_height_1"]  # 24
    meander_left_height_2 = config["meander_left_height_2"]  # 58
    meander_left_height_3 = config["meander_left_height_3"]  # 36

    meander_left_width_1 = config["meander_left_width_1"]  # 564
    meander_left_width_2 = config["meander_left_width_2"]  # 564

    meander_right_width_1 = config["meander_right_width_1"]  # 565
    meander_right_width_2 = config["meander_right_width_2"]  # 565

    meander_step_back_from_frame = config["meander_step_back_from_frame"]  # TODO Add this to the config.

    meander_path_points = [
        [
            (meander_bot_width / 2) - meander_right_width_1 + meander_right_width_2,
            meander_lw
            + meander_right_height_1
            + meander_right_height_2
            + meander_right_height_3
            - meander_step_back_from_frame,  # TODO THIS MIGHT BE AN ERROR WHERE THE MEANDER_LW SHOULD BE (MEANDER_LW/2)
        ],
        [
            (meander_bot_width / 2) - meander_right_width_1 + meander_right_width_2,
            (meander_lw / 2) + meander_right_height_1 + meander_right_height_2,
        ],
        [(meander_bot_width / 2) - meander_right_width_1, (meander_lw / 2) + meander_right_height_1 + meander_right_height_2],
        [(meander_bot_width / 2) - meander_right_width_1, (meander_lw / 2) + meander_right_height_1],
        [(meander_bot_width / 2), (meander_lw / 2) + meander_right_height_1],
        [(meander_bot_width / 2), (meander_lw / 2)],
        [-(meander_bot_width / 2), (meander_lw / 2)],
        [-(meander_bot_width / 2), (meander_lw / 2) + meander_left_height_1],
        [-(meander_bot_width / 2) - meander_left_width_1, (meander_lw / 2) + meander_left_height_1],
        [-(meander_bot_width / 2) - meander_left_width_1, (meander_lw / 2) + meander_left_height_1 + meander_left_height_2],
        [
            -(meander_bot_width / 2) - meander_left_width_1 + meander_left_width_2,
            (meander_lw / 2) + meander_left_height_1 + meander_left_height_2,
        ],
        [
            -(meander_bot_width / 2) - meander_left_width_1 + meander_left_width_2,
            meander_lw
            + meander_left_height_1
            + meander_left_height_2
            + meander_left_height_3
            - meander_step_back_from_frame,  # TODO THIS MIGHT BE AN ERROR WHERE THE MEANDER_LW SHOULD BE (MEANDER_LW/2)
        ],
    ]
    # Making the meander ant pad
    ant_pad_box_width = config["ant_pad_box_width"]  # 5
    ant_pad_box_height = config["ant_pad_box_height"]  # 10

    ant_pad_box_points = [
        [-ant_pad_box_width / 2, 0],
        [-ant_pad_box_width / 2, -ant_pad_box_height],
        [ant_pad_box_width / 2, -ant_pad_box_height],
        [ant_pad_box_width / 2, 0],
    ]

    # Making the frame sections left and right
    frame_bot_lw = config["frame_bot_lw"]  # 8
    frame_bot_left_width = config["frame_bot_left_width"]  # 996
    frame_bot_right_width = config["frame_bot_right_width"]  # 996

    frame_left_lw = config["frame_left_lw"]  # 8
    frame_left_height = config["frame_left_height"]  # 400

    frame_right_lw = config["frame_right_lw"]  # 8
    frame_right_height = config["frame_right_height"]  # 400

    frame_start_y = meander_right_height_1 + meander_right_height_2 + meander_right_height_3 + meander_lw
    frame_left_start_x = meander_path_points[-1][0] + (meander_lw / 2)
    frame_right_start_x = meander_path_points[0][0] - (meander_lw / 2)

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
    frame_meander_cover_box_width = config["frame_meander_cover_box_width"]  # 5
    frame_meander_cover_box_height = config["frame_meander_cover_box_height"]  # 28

    extra_frame_meander_cover_box_width = config["extra_frame_meander_cover_box_width"]  # TODO add this to the config
    extra_frame_meander_cover_box_height = config["extra_frame_meander_cover_box_height"]  # TODO add this to the config

    meander_cover_box_right_points = [
        [meander_path_points[0][0] - (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
        [
            meander_path_points[0][0] - (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [
            meander_path_points[0][0] + (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [meander_path_points[0][0] + (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
    ]

    meander_cover_box_left_points = [
        [meander_path_points[-1][0] - (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
        [
            meander_path_points[-1][0] - (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [
            meander_path_points[-1][0] + (frame_meander_cover_box_width / 2),
            frame_start_y + frame_bot_lw - frame_meander_cover_box_height,
        ],
        [meander_path_points[-1][0] + (frame_meander_cover_box_width / 2), frame_start_y + frame_bot_lw],
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
    coupler_frame_left_lw = config["coupler_frame_left_lw"]  # 10
    coupler_frame_left_height = config["coupler_frame_left_height"]  # 39
    coupler_frame_top_lw = config["coupler_frame_top_lw"]  # 3

    coupler_frame_start_x = frame_left_start_x - frame_bot_left_width + (frame_left_lw / 2)
    coupler_frame_start_y = frame_start_y + frame_left_height

    # Making the IDC and trim arms
    IDC_bot_arm_gap = config["IDC_bot_arm_gap"]  # 30
    IDC_arm_gap = config["IDC_arm_gap"]  # 8
    IDC_arm_lw = config["IDC_arm_lw"]  # 3
    No_of_arms = config["No_of_arms"]  # 28

    arm_start_x_left_side = frame_left_start_x - frame_bot_left_width + frame_left_lw
    arm_start_x_right_side = frame_right_start_x + frame_bot_right_width - frame_right_lw

    arm_start_y_right_side = frame_start_y + frame_bot_lw + IDC_bot_arm_gap + (IDC_arm_lw / 2)
    arm_start_y_left_side = arm_start_y_right_side + IDC_arm_gap + IDC_arm_lw

    trim_arm_offset_right_side = config["trim_arm_offset_right_side"]  # 380
    trim_arm_offset_left_side = config["trim_arm_offset_left_side"]  # 389
    trim_arm_lw = config["trim_arm_lw"]  # 3
    trim_arm_length_right_side = config["trim_arm_length_right_side"]  # 1975
    trim_arm_length_left_side = config["trim_arm_length_left_side"]  # 1975

    trim_arm_start_y_right_side = frame_start_y + frame_bot_lw + trim_arm_offset_right_side
    trim_arm_start_y_left_side = frame_start_y + frame_bot_lw + trim_arm_offset_left_side

    # Making the coupler ataching to feedline
    coupler_gap = config["coupler_gap"]  # 16
    coupler_lw = config["coupler_lw"]  # 3
    left_coupler_frame_to_feed_distance = config["left_coupler_frame_to_feed_distance"]  # 164

    # Making the ground plane cutout
    cutout_bot_offset = config["cutout_bot_offset"]  # 15
    cutout_left_offset = config["cutout_left_offset"]  # 50
    cutout_right_offset = config["cutout_right_offset"]  # 50
    cutout_top_offset = config["cutout_top_offset"]  # 25

    # grndpl_meander_cutout_width = config["grndpl_meander_cutout_width"]#80
    # grndpl_meander_cutout_height = config["grndpl_meander_cutout_height"]#10

    cutout_start_height = meander_lw + meander_right_height_1 + meander_right_height_2 + meander_right_height_3 - cutout_bot_offset
    cutout_width = (frame_right_start_x + frame_bot_right_width + cutout_right_offset) - (
        frame_left_start_x - frame_bot_left_width - cutout_left_offset
    )
    cutout_height = (coupler_frame_start_y + coupler_frame_left_height + coupler_gap + coupler_lw + cutout_top_offset) - cutout_start_height

    grnd_plane_cutout_width = cutout_width + 30
    grnd_plane_cutout_height = cutout_height + 30
    grnd_plane_cutout_start_height = cutout_start_height - 15

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
    SiO_stepdown_cutout_width = config["SiO_stepdown_cutout_width"]  # 110
    SiO_stepdown_cutout_height = config["SiO_stepdown_cutout_height"]  # 39
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
    SiN_membrane_stepdown_cutout_width = config["SiN_membrane_stepdown_cutout_width"]  # 100
    SiN_membrane_stepdown_cutout_height = config["SiN_membrane_stepdown_cutout_height"]  # 36
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

    # Getting the IDC and CC lengths from the function
    IDCLs, CCL = IDC_and_CC_function(f0)

    # Adding the meander

    if mirror:
        new_meander_path_points = souk_mask_builder.mirror_points_around_yaxis(meander_path_points)
        new_meander_path_points = souk_mask_builder.rotate_and_move_points_list(new_meander_path_points, rot_angle, x, y)
    else:
        new_meander_path_points = souk_mask_builder.rotate_and_move_points_list(meander_path_points, rot_angle, x, y)
    meander_path = gdspy.FlexPath(
        new_meander_path_points, meander_lw, corners="circular bend", bend_radius=meander_corner_bend_radius, **material_meander
    )
    meander_path_polygons = souk_mask_builder.get_polys_from_flexpath(meander_path)
    for i in range(len(meander_path_polygons)):
        souk_mask_builder.Main.add(gdspy.Polygon(meander_path_polygons[i], **material_meander))

    # Adding the meander ant overlap box
    new_ant_pad_box_points = souk_mask_builder.rotate_and_move_points_list(ant_pad_box_points, rot_angle, x, y)
    ant_pad_box = gdspy.Polygon(new_ant_pad_box_points, **material_meander)
    souk_mask_builder.Main.add(ant_pad_box)

    # Adding the meander frame overlap boxes
    if mirror:
        new_meander_cover_box_right_points = souk_mask_builder.mirror_points_around_yaxis(meander_cover_box_right_points)
        new_meander_cover_box_right_points = souk_mask_builder.rotate_and_move_points_list(
            new_meander_cover_box_right_points, rot_angle, x, y
        )
    else:
        new_meander_cover_box_right_points = souk_mask_builder.rotate_and_move_points_list(meander_cover_box_right_points, rot_angle, x, y)

    meander_cover_box_right = gdspy.Polygon(new_meander_cover_box_right_points, **material_idc_and_frame)
    souk_mask_builder.Main.add(meander_cover_box_right)

    if mirror:
        new_meander_cover_box_left_points = souk_mask_builder.mirror_points_around_yaxis(meander_cover_box_left_points)
        new_meander_cover_box_left_points = souk_mask_builder.rotate_and_move_points_list(
            new_meander_cover_box_left_points, rot_angle, x, y
        )
    else:
        new_meander_cover_box_left_points = souk_mask_builder.rotate_and_move_points_list(meander_cover_box_left_points, rot_angle, x, y)

    meander_cover_box_left = gdspy.Polygon(new_meander_cover_box_left_points, **material_idc_and_frame)
    souk_mask_builder.Main.add(meander_cover_box_left)

    # Adding the extra meander frame overlap boxes
    if mirror:
        new_extra_meander_cover_box_right_points = souk_mask_builder.mirror_points_around_yaxis(extra_meander_cover_box_right_points)
        new_extra_meander_cover_box_right_points = souk_mask_builder.rotate_and_move_points_list(
            new_extra_meander_cover_box_right_points, rot_angle, x, y
        )
    else:
        new_extra_meander_cover_box_right_points = souk_mask_builder.rotate_and_move_points_list(
            extra_meander_cover_box_right_points, rot_angle, x, y
        )

    extra_meander_cover_box_right = gdspy.Polygon(new_extra_meander_cover_box_right_points, **material_meander)
    souk_mask_builder.Main.add(extra_meander_cover_box_right)

    if mirror:
        new_extra_meander_cover_box_left_points = souk_mask_builder.mirror_points_around_yaxis(extra_meander_cover_box_left_points)
        new_extra_meander_cover_box_left_points = souk_mask_builder.rotate_and_move_points_list(
            new_extra_meander_cover_box_left_points, rot_angle, x, y
        )
    else:
        new_extra_meander_cover_box_left_points = souk_mask_builder.rotate_and_move_points_list(
            extra_meander_cover_box_left_points, rot_angle, x, y
        )
    extra_meander_cover_box_left = gdspy.Polygon(new_extra_meander_cover_box_left_points, **material_meander)
    souk_mask_builder.Main.add(extra_meander_cover_box_left)

    # Adding the frame left and frame right
    if mirror:
        new_frame_left_points = souk_mask_builder.mirror_points_around_yaxis(frame_left_points)
        new_frame_left_points = souk_mask_builder.rotate_and_move_points_list(new_frame_left_points, rot_angle, x, y)
    else:
        new_frame_left_points = souk_mask_builder.rotate_and_move_points_list(frame_left_points, rot_angle, x, y)

    frame_left_poly = gdspy.Polygon(new_frame_left_points, **material_idc_and_frame)
    souk_mask_builder.Main.add(frame_left_poly)

    if mirror:
        new_frame_right_points = souk_mask_builder.mirror_points_around_yaxis(frame_right_points)
        new_frame_right_points = souk_mask_builder.rotate_and_move_points_list(new_frame_right_points, rot_angle, x, y)
    else:
        new_frame_right_points = souk_mask_builder.rotate_and_move_points_list(frame_right_points, rot_angle, x, y)

    frame_right_poly = gdspy.Polygon(new_frame_right_points, **material_idc_and_frame)
    souk_mask_builder.Main.add(frame_right_poly)

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
        new_coupler_frame_points = souk_mask_builder.mirror_points_around_yaxis(coupler_frame_points)
        new_coupler_frame_points = souk_mask_builder.rotate_and_move_points_list(new_coupler_frame_points, rot_angle, x, y)
    else:
        new_coupler_frame_points = souk_mask_builder.rotate_and_move_points_list(coupler_frame_points, rot_angle, x, y)

    coupler_frame_poly = gdspy.Polygon(new_coupler_frame_points, **material_idc_and_frame)
    souk_mask_builder.Main.add(coupler_frame_poly)

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
        new_coupler_arm_points = souk_mask_builder.mirror_points_around_yaxis(coupler_arm_points)
        new_coupler_arm_points = souk_mask_builder.rotate_and_move_points_list(new_coupler_arm_points, rot_angle, x, y)
    else:
        new_coupler_arm_points = souk_mask_builder.rotate_and_move_points_list(coupler_arm_points, rot_angle, x, y)

    coupler_arm_poly = gdspy.Polygon(new_coupler_arm_points, **souk_mask_builder.Nb_Antenna)
    souk_mask_builder.Main.add(coupler_arm_poly)

    # Adding the IDC arms
    for i in range(0, No_of_arms, 2):
        right_arm = gdspy.Rectangle(
            [arm_start_x_right_side, arm_start_y_right_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            [arm_start_x_right_side - IDCLs[-(i + 1)], arm_start_y_right_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            **material_idc_and_frame,
        )
        right_arm.translate(x, y)
        if mirror:
            right_arm.mirror([x, y], [x, y + 10])
        right_arm.rotate(rot_angle, center=(x, y))
        souk_mask_builder.Main.add(right_arm)

        left_arm = gdspy.Rectangle(
            [arm_start_x_left_side, arm_start_y_left_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            [arm_start_x_left_side + IDCLs[-(i + 2)], arm_start_y_left_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
            **material_idc_and_frame,
        )
        left_arm.translate(x, y)
        if mirror:
            left_arm.mirror([x, y], [x, y + 10])
        left_arm.rotate(rot_angle, center=(x, y))
        souk_mask_builder.Main.add(left_arm)

    # Adding the Trim arms.
    right_trim_arm = gdspy.Rectangle(
        [arm_start_x_right_side, trim_arm_start_y_right_side],
        [arm_start_x_right_side - trim_arm_length_right_side, trim_arm_start_y_right_side + trim_arm_lw],
        **material_idc_and_frame,
    )
    right_trim_arm.translate(x, y)
    if mirror:
        right_trim_arm.mirror([x, y], [x, y + 10])
    right_trim_arm.rotate(rot_angle, center=(x, y))
    souk_mask_builder.Main.add(right_trim_arm)

    left_trim_arm = gdspy.Rectangle(
        [arm_start_x_left_side, trim_arm_start_y_left_side],
        [arm_start_x_left_side + trim_arm_length_left_side, trim_arm_start_y_left_side + trim_arm_lw],
        **material_idc_and_frame,
    )
    left_trim_arm.translate(x, y)
    if mirror:
        left_trim_arm.mirror([x, y], [x, y + 10])
    left_trim_arm.rotate(rot_angle, center=(x, y))
    souk_mask_builder.Main.add(left_trim_arm)

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
                **souk_mask_builder.TrimLayer,
            )
            trim_box_for_right_trim_arm.translate(x, y)
            if mirror:
                trim_box_for_right_trim_arm.mirror([x, y], [x, y + 10])
            trim_box_for_right_trim_arm.rotate(rot_angle, center=(x, y))
            souk_mask_builder.Main.add(trim_box_for_right_trim_arm)

            trim_box_for_left_trim_arm = gdspy.Rectangle(
                [
                    arm_start_x_left_side + trim_arm_length_left_side + trim_box_length_overhang,
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                ],
                [
                    arm_start_x_right_side - (inner_width - trim_length),
                    trim_arm_start_y_left_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                ],
                **souk_mask_builder.TrimLayer,
            )
            trim_box_for_left_trim_arm.translate(x, y)
            if mirror:
                trim_box_for_left_trim_arm.mirror([x, y], [x, y + 10])
            trim_box_for_left_trim_arm.rotate(rot_angle, center=(x, y))
            souk_mask_builder.Main.add(trim_box_for_left_trim_arm)

    # Adding the cutout to the groundplane.
    if add_grnd_cutout:
        if mirror:
            new_grnd_plane_meander_cutout_poly_points = souk_mask_builder.mirror_points_around_yaxis(grnd_plane_meander_cutout_poly_points)
            new_grnd_plane_meander_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                new_grnd_plane_meander_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_grnd_plane_meander_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                grnd_plane_meander_cutout_poly_points, rot_angle, x, y
            )

        grnd_plane_meander_cutout_poly = gdspy.Polygon(new_grnd_plane_meander_cutout_poly_points, **souk_mask_builder.Nb_Groundplane)
        souk_mask_builder.ground_plane_cutouts.add(grnd_plane_meander_cutout_poly)

    # Adding the cutout to the Silicon DiOxide membrane.
    if add_SiO_cutout:
        if mirror:
            new_SiO_cutout_poly_points = souk_mask_builder.mirror_points_around_yaxis(SiO_cutout_poly_points)
            new_SiO_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(new_SiO_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiO_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(SiO_cutout_poly_points, rot_angle, x, y)

        SiO_cutout_poly = gdspy.Polygon(new_SiO_cutout_poly_points, **souk_mask_builder.Nb_Groundplane)
        souk_mask_builder.silicon_oxide_cutouts.add(SiO_cutout_poly)

    # Adding the cutout to the Silicon Nitride membrane.
    if add_SiN_membrane_cutout:
        if mirror:
            new_SiN_membrane_cutout_poly_points = souk_mask_builder.mirror_points_around_yaxis(SiN_membrane_cutout_poly_points)
            new_SiN_membrane_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                new_SiN_membrane_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_SiN_membrane_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                SiN_membrane_cutout_poly_points, rot_angle, x, y
            )

        SiN_membrane_cutout_poly = gdspy.Polygon(new_SiN_membrane_cutout_poly_points, **souk_mask_builder.Nb_Groundplane)
        souk_mask_builder.silicon_nitride_membrane_cutouts.add(SiN_membrane_cutout_poly)

    # Adding the cutout to the SiN Dep layer.
    if add_SiN_dep_dielectric_cutout:
        if mirror:
            new_SiN_dep_cutout_poly_points = souk_mask_builder.mirror_points_around_yaxis(SiN_dep_cutout_poly_points)
            new_SiN_dep_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(new_SiN_dep_cutout_poly_points, rot_angle, x, y)
        else:
            new_SiN_dep_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(SiN_dep_cutout_poly_points, rot_angle, x, y)

        SiN_dep_cutout_poly = gdspy.Polygon(new_SiN_dep_cutout_poly_points, **souk_mask_builder.SiN_dep)
        souk_mask_builder.silicon_nitride_cutouts.add(SiN_dep_cutout_poly)

    # Adding the backside check covers.
    if add_backside_check:
        if mirror:
            new_backside_check_cover_poly_points = souk_mask_builder.mirror_points_around_yaxis(backside_check_cover_poly_points)
            new_backside_check_cover_poly_points = souk_mask_builder.rotate_and_move_points_list(
                new_backside_check_cover_poly_points, rot_angle, x, y
            )
        else:
            new_backside_check_cover_poly_points = souk_mask_builder.rotate_and_move_points_list(
                backside_check_cover_poly_points, rot_angle, x, y
            )

        backside_check_cover_poly = gdspy.Polygon(new_backside_check_cover_poly_points, **souk_mask_builder.Backside_Check)
        souk_mask_builder.Main.add(backside_check_cover_poly)

    # Adding the groundplane cutout over the inductive meander.
    if add_grnd_cutout_over_inductor:
        if mirror:
            new_grnd_plane_inductor_cutout_poly_points = souk_mask_builder.mirror_points_around_yaxis(
                grnd_plane_inductor_cutout_poly_points
            )
            new_grnd_plane_inductor_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                new_grnd_plane_inductor_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_grnd_plane_inductor_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                grnd_plane_inductor_cutout_poly_points, rot_angle, x, y
            )

        grnd_plane_inductor_cutout_poly = gdspy.Polygon(new_grnd_plane_inductor_cutout_poly_points, **souk_mask_builder.Nb_Groundplane)
        souk_mask_builder.ground_plane_cutouts.add(grnd_plane_inductor_cutout_poly)

    # Adding the SiNdep cutout over the inductive meander.
    if add_SiN_dep_dielectric_cutout_over_inductor:
        if mirror:
            new_SiN_dep_inductor_cutout_poly_points = souk_mask_builder.mirror_points_around_yaxis(SiN_dep_inductor_cutout_poly_points)
            new_SiN_dep_inductor_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                new_SiN_dep_inductor_cutout_poly_points, rot_angle, x, y
            )
        else:
            new_SiN_dep_inductor_cutout_poly_points = souk_mask_builder.rotate_and_move_points_list(
                SiN_dep_inductor_cutout_poly_points, rot_angle, x, y
            )

        SiN_dep_inductor_cutout_poly = gdspy.Polygon(new_SiN_dep_inductor_cutout_poly_points, **souk_mask_builder.SiN_dep)
        souk_mask_builder.silicon_nitride_cutouts.add(SiN_dep_inductor_cutout_poly)

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

    configurator_points["ant_pad_box_width"] = {
        "text": "ant_pad_box_width",
        "start": [new_ant_pad_box_points[0][0], (new_ant_pad_box_points[0][1] + new_ant_pad_box_points[1][1]) / 3],
        "end": [new_ant_pad_box_points[3][0], (new_ant_pad_box_points[0][1] + new_ant_pad_box_points[1][1]) / 3],
    }

    configurator_points["ant_pad_box_height"] = {
        "text": "ant_pad_box_height",
        "start": [
            (new_ant_pad_box_points[0][0] + new_ant_pad_box_points[3][0]) / 2,
            new_ant_pad_box_points[0][1],
        ],
        "end": [
            (new_ant_pad_box_points[1][0] + new_ant_pad_box_points[2][0]) / 2,
            new_ant_pad_box_points[1][1],
        ],
    }

    configurator_points["meander_corner_bend_radius"] = {
        "text": "meander_corner_bend_radius",
        "start": [
            new_meander_path_points[6][0] + meander_corner_bend_radius,
            new_meander_path_points[6][1] + meander_corner_bend_radius,
        ],
        "end": [
            new_meander_path_points[6][0] + meander_corner_bend_radius - (np.sqrt(2) / 2) * meander_corner_bend_radius,
            new_meander_path_points[6][1] + meander_corner_bend_radius - (np.sqrt(2) / 2) * meander_corner_bend_radius,
        ],
    }

    configurator_points["meander_bot_width"] = {
        "text": "meander_bot_width",
        "start": [
            new_meander_path_points[6][0],
            new_meander_path_points[6][1] + meander_left_height_1 / 2,
        ],
        "end": [
            new_meander_path_points[5][0],
            new_meander_path_points[5][1] + meander_left_height_1 / 2,
        ],
    }

    configurator_points["meander_left_height_1"] = {
        "text": "meander_left_height_1",
        "start": [
            new_meander_path_points[6][0] - meander_corner_bend_radius,
            new_meander_path_points[6][1],
        ],
        "end": [
            new_meander_path_points[7][0] - meander_corner_bend_radius,
            new_meander_path_points[7][1],
        ],
    }

    configurator_points["meander_right_height_1"] = {
        "text": "meander_right_height_1",
        "start": [
            new_meander_path_points[5][0],
            new_meander_path_points[5][1],
        ],
        "end": [
            new_meander_path_points[4][0],
            new_meander_path_points[4][1],
        ],
    }

    configurator_points["meander_left_width_1"] = {
        "text": "meander_left_width_1",
        "start": [
            new_meander_path_points[7][0],
            new_meander_path_points[7][1],
        ],
        "end": [
            new_meander_path_points[8][0],
            new_meander_path_points[8][1],
        ],
    }

    configurator_points["meander_right_width_1"] = {
        "text": "meander_right_width_1",
        "start": [
            new_meander_path_points[4][0],
            new_meander_path_points[4][1],
        ],
        "end": [
            new_meander_path_points[3][0],
            new_meander_path_points[3][1],
        ],
    }

    configurator_points["meander_left_height_2"] = {
        "text": "meander_left_height_2",
        "start": [
            new_meander_path_points[8][0],
            new_meander_path_points[8][1],
        ],
        "end": [
            new_meander_path_points[9][0],
            new_meander_path_points[9][1],
        ],
    }

    configurator_points["meander_right_height_2"] = {
        "text": "meander_right_height_2",
        "start": [
            new_meander_path_points[3][0],
            new_meander_path_points[3][1],
        ],
        "end": [
            new_meander_path_points[2][0],
            new_meander_path_points[2][1],
        ],
    }

    configurator_points["meander_left_width_2"] = {
        "text": "meander_left_width_2",
        "start": [
            new_meander_path_points[9][0],
            new_meander_path_points[9][1],
        ],
        "end": [
            new_meander_path_points[10][0],
            new_meander_path_points[10][1],
        ],
    }

    configurator_points["meander_right_width_2"] = {
        "text": "meander_right_width_2",
        "start": [
            new_meander_path_points[2][0],
            new_meander_path_points[2][1],
        ],
        "end": [
            new_meander_path_points[1][0],
            new_meander_path_points[1][1],
        ],
    }

    configurator_points["meander_left_height_3"] = {
        "text": "meander_left_height_3",
        "start": [
            new_meander_path_points[10][0],
            new_meander_path_points[10][1],
        ],
        "end": [
            new_meander_path_points[11][0],
            new_meander_path_points[11][1],
        ],
    }

    configurator_points["meander_right_height_3"] = {
        "text": "meander_right_height_3",
        "start": [
            new_meander_path_points[1][0],
            new_meander_path_points[1][1],
        ],
        "end": [
            new_meander_path_points[0][0],
            new_meander_path_points[0][1],
        ],
    }

    ####################################################################### IDC Frame
    # left side
    configurator_points["frame_meander_cover_box_width_left"] = {
        "text": "frame_meander_cover_box_width",
        "start": [
            new_meander_cover_box_left_points[1][0],
            new_meander_cover_box_left_points[1][1]
            + (new_meander_cover_box_left_points[0][1] - new_meander_cover_box_left_points[1][1]) / 1.5,
        ],
        "end": [
            new_meander_cover_box_left_points[2][0],
            new_meander_cover_box_left_points[2][1]
            + (new_meander_cover_box_left_points[3][1] - new_meander_cover_box_left_points[2][1]) / 1.5,
        ],
    }
    # right side
    configurator_points["frame_meander_cover_box_width_right"] = {
        "text": "frame_meander_cover_box_width",
        "start": [
            new_meander_cover_box_right_points[1][0],
            new_meander_cover_box_right_points[1][1]
            + (new_meander_cover_box_right_points[0][1] - new_meander_cover_box_right_points[1][1]) / 1.5,
        ],
        "end": [
            new_meander_cover_box_right_points[2][0],
            new_meander_cover_box_right_points[2][1]
            + (new_meander_cover_box_right_points[3][1] - new_meander_cover_box_right_points[2][1]) / 1.5,
        ],
    }

    # left side
    configurator_points["frame_meander_cover_box_height_left"] = {
        "text": "frame_meander_cover_box_height",
        "start": [new_meander_cover_box_left_points[3][0], new_meander_cover_box_left_points[3][1]],
        "end": [new_meander_cover_box_left_points[2][0], new_meander_cover_box_left_points[2][1]],
    }
    # right side
    configurator_points["frame_meander_cover_box_height_right"] = {
        "text": "frame_meander_cover_box_height",
        "start": [new_meander_cover_box_right_points[3][0], new_meander_cover_box_right_points[3][1]],
        "end": [new_meander_cover_box_right_points[2][0], new_meander_cover_box_right_points[2][1]],
    }

    # left side
    configurator_points["frame_bot_lw_left"] = {
        "text": "frame_bot_lw",
        "start": [
            new_frame_left_points[0][0] + (new_frame_left_points[1][0] - new_frame_left_points[0][0]) / 2,
            new_frame_left_points[5][1],
        ],
        "end": [
            new_frame_left_points[0][0] + (new_frame_left_points[1][0] - new_frame_left_points[0][0]) / 2,
            new_frame_left_points[0][1],
        ],
    }
    # right side
    configurator_points["frame_bot_lw_right"] = {
        "text": "frame_bot_lw",
        "start": [
            new_frame_right_points[0][0] + (new_frame_right_points[1][0] - new_frame_right_points[0][0]) / 2,
            new_frame_right_points[5][1],
        ],
        "end": [
            new_frame_right_points[0][0] + (new_frame_right_points[1][0] - new_frame_right_points[0][0]) / 2,
            new_frame_right_points[0][1],
        ],
    }

    configurator_points["frame_bot_left_width"] = {
        "text": "frame_bot_left_width",
        "start": [new_frame_left_points[5][0], new_frame_left_points[5][1]],
        "end": [new_frame_left_points[4][0], new_frame_left_points[4][1]],
    }

    configurator_points["frame_bot_right_width"] = {
        "text": "frame_bot_right_width",
        "start": [new_frame_right_points[5][0], new_frame_right_points[5][1]],
        "end": [new_frame_right_points[4][0], new_frame_right_points[4][1]],
    }

    configurator_points["frame_left_lw"] = {
        "text": "frame_left_lw",
        "start": [new_frame_left_points[1][0], new_frame_left_points[1][1]],
        "end": [new_frame_left_points[4][0], new_frame_left_points[1][1]],
    }

    configurator_points["frame_right_lw"] = {
        "text": "frame_right_lw",
        "start": [new_frame_right_points[1][0], new_frame_right_points[1][1]],
        "end": [new_frame_right_points[4][0], new_frame_right_points[1][1]],
    }

    configurator_points["frame_left_height"] = {
        "text": "frame_left_height",
        "start": [new_frame_left_points[3][0], new_frame_left_points[3][1]],
        "end": [new_frame_left_points[4][0], new_frame_left_points[4][1]],
    }

    configurator_points["frame_right_height"] = {
        "text": "frame_right_height",
        "start": [new_frame_right_points[3][0], new_frame_right_points[3][1]],
        "end": [new_frame_right_points[4][0], new_frame_right_points[4][1]],
    }

    ####################################################################### Coupler Frame
    configurator_points["coupler_frame_left_lw"] = {
        "text": "coupler_frame_left_lw",
        "start": [
            new_coupler_frame_points[0][0],
            new_coupler_frame_points[0][1] + (new_coupler_frame_points[1][1] - new_coupler_frame_points[0][1]) / 1.5,
        ],
        "end": [
            new_coupler_frame_points[5][0],
            new_coupler_frame_points[5][1] + (new_coupler_frame_points[1][1] - new_coupler_frame_points[0][1]) / 1.5,
        ],
    }

    configurator_points["coupler_frame_left_height"] = {
        "text": "coupler_frame_left_height",
        "start": [new_coupler_frame_points[0][0], new_coupler_frame_points[0][1]],
        "end": [new_coupler_frame_points[1][0], new_coupler_frame_points[1][1]],
    }

    configurator_points["coupler_frame_top_lw"] = {
        "text": "coupler_frame_top_lw",
        "start": [
            new_coupler_frame_points[1][0] + (new_coupler_frame_points[2][0] - new_coupler_frame_points[1][0]) / 2.5,
            new_coupler_frame_points[1][1],
        ],
        "end": [
            new_coupler_frame_points[1][0] + (new_coupler_frame_points[2][0] - new_coupler_frame_points[1][0]) / 2.5,
            new_coupler_frame_points[4][1],
        ],
    }

    ####################################################################### IDC Arms
    configurator_points["IDC_bot_arm_gap"] = {
        "text": "IDC_bot_arm_gap",
        "start": [
            new_frame_left_points[0][0] + (new_frame_left_points[1][0] - new_frame_left_points[0][0]) / 2.5,
            new_frame_left_points[0][1],
        ],
        "end": [
            new_frame_left_points[0][0] + (new_frame_left_points[1][0] - new_frame_left_points[0][0]) / 2.5,
            new_frame_left_points[0][1] + IDC_bot_arm_gap,
        ],
    }

    configurator_points["IDC_arm_gap"] = {
        "text": "IDC_arm_gap",
        "start": [
            (arm_start_x_left_side + arm_start_x_right_side) / 2,
            arm_start_y_right_side + IDC_arm_lw / 2,
        ],
        "end": [
            (arm_start_x_left_side + arm_start_x_right_side) / 2,
            arm_start_y_right_side + IDC_arm_gap + IDC_arm_lw / 2,
        ],
    }

    configurator_points["IDC_arm_lw"] = {
        "text": "IDC_arm_lw",
        "start": [
            (arm_start_x_left_side + arm_start_x_right_side) / 2,
            arm_start_y_right_side - IDC_arm_lw / 2,
        ],
        "end": [
            (arm_start_x_left_side + arm_start_x_right_side) / 2,
            arm_start_y_right_side + IDC_arm_lw / 2,
        ],
    }

    configurator_points["No_of_arms"] = {
        "text": "No_of_arms",
        "start": [
            ((arm_start_x_left_side + arm_start_x_right_side) / 2) + (arm_start_x_right_side / 10),
            arm_start_y_right_side - IDC_arm_lw / 2,
        ],
        "end": [
            ((arm_start_x_left_side + arm_start_x_right_side) / 2) + (arm_start_x_right_side / 10),
            arm_start_y_right_side + No_of_arms * (IDC_arm_lw + IDC_arm_gap) - IDC_arm_gap,
        ],
    }

    ####################################################################### Trim Arms
    configurator_points["trim_arm_offset_left_side"] = {
        "text": "trim_arm_offset_left_side",
        "start": [
            new_frame_left_points[1][0] + (frame_bot_left_width / 2),
            new_frame_left_points[1][1],
        ],
        "end": [
            new_frame_left_points[1][0] + (frame_bot_left_width / 2),
            trim_arm_start_y_left_side,
        ],
    }

    configurator_points["trim_arm_offset_right_side"] = {
        "text": "trim_arm_offset_right_side",
        "start": [
            new_frame_right_points[1][0] - (frame_bot_right_width / 2),
            new_frame_right_points[1][1],
        ],
        "end": [
            new_frame_right_points[1][0] - (frame_bot_right_width / 2),
            trim_arm_start_y_right_side,
        ],
    }

    configurator_points["trim_arm_lw"] = {
        "text": "trim_arm_lw",
        "start": [
            (arm_start_x_left_side + arm_start_x_right_side) / 2,
            trim_arm_start_y_right_side,
        ],
        "end": [
            (arm_start_x_left_side + arm_start_x_right_side) / 2,
            trim_arm_start_y_right_side + trim_arm_lw,
        ],
    }

    configurator_points["trim_arm_length_left_side"] = {
        "text": "trim_arm_length_left_side",
        "start": [
            arm_start_x_left_side,
            trim_arm_start_y_left_side,
        ],
        "end": [
            arm_start_x_left_side + trim_arm_length_left_side,
            trim_arm_start_y_left_side,
        ],
    }

    configurator_points["trim_arm_length_right_side"] = {
        "text": "trim_arm_length_right_side",
        "start": [
            arm_start_x_right_side,
            trim_arm_start_y_right_side,
        ],
        "end": [
            arm_start_x_right_side - trim_arm_length_right_side,
            trim_arm_start_y_right_side,
        ],
    }

    ####################################################################### Coupler Arm
    configurator_points["coupler_gap"] = {
        "text": "coupler_gap",
        "start": [new_coupler_frame_points[4][0], new_coupler_frame_points[1][1]],
        "end": [new_coupler_frame_points[4][0], (new_coupler_frame_points[1][1] + coupler_gap)],
    }

    configurator_points["coupler_lw"] = {
        "text": "coupler_lw",
        "start": [new_coupler_frame_points[4][0], (new_coupler_frame_points[1][1] + coupler_gap)],
        "end": [new_coupler_frame_points[4][0], (new_coupler_frame_points[1][1] + coupler_gap + coupler_lw)],
    }

    configurator_points["left_coupler_frame_to_feed_distance"] = {
        "text": "left_coupler_frame_to_feed_distance",
        "start": [new_coupler_frame_points[1][0], new_coupler_frame_points[1][1]],
        "end": [(new_coupler_frame_points[1][0] - left_coupler_frame_to_feed_distance), new_coupler_frame_points[1][1]],
    }

    ####################################################################### KID Number Text
    # TODO
    # configurator_points["text_size"] = { }
    #
    # configurator_points["text_x_offset"] = { }
    #
    # configurator_points["text_y_offset"] = { }

    ####################################################################### Cutouts
    configurator_points["cutout_bot_offset"] = {
        "text": "cutout_bot_offset",
        "start": [
            new_frame_left_points[0][0] + (new_frame_left_points[1][0] - new_frame_left_points[0][0]) / 2.5,
            cutout_start_height + cutout_bot_offset,
        ],
        "end": [new_frame_left_points[0][0] + (new_frame_left_points[1][0] - new_frame_left_points[0][0]) / 2.5, cutout_start_height],
    }

    configurator_points["cutout_left_offset"] = {
        "text": "cutout_left_offset",
        "start": [new_frame_left_points[4][0], new_frame_left_points[4][1]],
        "end": [new_frame_left_points[4][0] - cutout_left_offset, new_frame_left_points[4][1]],
    }
    configurator_points["cutout_right_offset"] = {
        "text": "cutout_right_offset",
        "start": [new_frame_right_points[4][0], new_frame_right_points[4][1]],
        "end": [new_frame_right_points[4][0] + cutout_right_offset, new_frame_right_points[4][1]],
    }
    configurator_points["cutout_top_offset"] = {
        "text": "cutout_top_offset",
        "start": [new_coupler_frame_points[4][0], (new_coupler_frame_points[1][1] + coupler_gap + coupler_lw)],
        "end": [new_coupler_frame_points[4][0], (new_coupler_frame_points[1][1] + coupler_gap + coupler_lw + cutout_top_offset)],
    }

    configurator_points["SiO_stepdown_cutout_width"] = {
        "text": "SiO_stepdown_cutout_width",
        "start": [new_SiO_cutout_poly_points[0][0], new_SiO_cutout_poly_points[0][1]],
        "end": [new_SiO_cutout_poly_points[7][0], new_SiO_cutout_poly_points[7][1]],
    }

    configurator_points["SiO_stepdown_cutout_height"] = {
        "text": "SiO_stepdown_cutout_height",
        "start": [new_SiO_cutout_poly_points[0][0], new_SiO_cutout_poly_points[0][1]],
        "end": [new_SiO_cutout_poly_points[1][0], new_SiO_cutout_poly_points[1][1]],
    }

    configurator_points["SiN_membrane_stepdown_cutout_width"] = {
        "text": "SiN_membrane_stepdown_cutout_width",
        "start": [new_SiN_membrane_cutout_poly_points[0][0], new_SiN_membrane_cutout_poly_points[0][1]],
        "end": [new_SiN_membrane_cutout_poly_points[7][0], new_SiN_membrane_cutout_poly_points[7][1]],
    }

    configurator_points["SiN_membrane_stepdown_cutout_height"] = {
        "text": "SiN_membrane_stepdown_cutout_height",
        "start": [new_SiN_membrane_cutout_poly_points[0][0], new_SiN_membrane_cutout_poly_points[0][1]],
        "end": [new_SiN_membrane_cutout_poly_points[1][0], new_SiN_membrane_cutout_poly_points[1][1]],
    }

    return configurator_points
