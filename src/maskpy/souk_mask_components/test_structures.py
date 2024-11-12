from collections.abc import Sequence

import gdspy
import numpy as np
from numpy import cos, pi, sin

from .. import mask_builder_utils as mbu
from ..layers import Layer
from ..logging.pretty_print import styled_type_error


def add_test_H_pad_connected_box_section(
    mask_builder,
    x: float | int,
    y: float | int,
    pad_materials: Layer | Sequence[Layer],
    line_materials: Layer | Sequence[Layer],
    rot: float | int = 0,
    linewidths: float | int = 2.0,
    TL_text: str = "",
    add_text_for_linewidth: bool = True,
    add_dielectric_square: bool = True,
    add_groundplane_square: bool = True,
    add_dielectric_cutout_in_ports: bool | Sequence[bool] = True,
    add_dielectric_cutout_over_H: bool | Sequence[bool] = (False, False, False),
    add_dielectric_cutout_over_line: bool | Sequence[bool] = (False, False, False),
) -> None:
    """Adds a box of size 5500 centered on the x,y given that contains a series
    of 6 pads where 2 of those connect to a line connecting them.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test structure.

    pad_materials: Layer | Sequence[Layer]
        This is the layers that each of the pad pairs should be drawn in.
        This is an instance of Layer of a Sequence of many Layers. If a single
        Layer is given then all the ports will be drawn in this layer.
        If a Sequence of Layers is given then it should be of length 6 and each
        pad pair will take the respective layer from left to right. For more info
        see maskpy.layers.Layer. Usually this is within the
        SoukMaskBuilder.layers.xxx. e.g. `mask_builder.layers.Aluminium`

    line_materials: Layer | Sequence[Layer]
        This is the layers that each of the lines should be drawn in.
        This is an instance of Layer of a Sequence of many Layers. If a single
        Layer is given then all the lines will be drawn in this layer.
        If a Sequence of Layers is given then it should be of length 3 and each
        of the lines will take the respective layer from left to right. For more
        info see maskpy.layers.Layer. Usually this is within the
        SoukMaskBuilder.layers.xxx. e.g. `mask_builder.layers.Aluminium`

    KwArgs
    ------
    rot: float | int = 0
        The rotation for the entire structure given in radians.

    linewidths: float | int = 2.0,
        This is the linewidths to use for each of the lines. Default is 2.0
        which will draw the lines all of this linewidth. This wil be the same
        for the U shape connections leading up to the connecting line.

    TL_text: str = "",
        The optional string to write in the top left corner of the structure.

    add_text_for_linewidth: bool = True,
        This is whether to add text denoting each linewidth next to the line.
        Default True will add the linewidth of the current line for each line.
        This is disabled with False. When True this is drawn in the same layer
        as the linewidths and pads.

    add_dielectric_square: bool = True,
        Default True will add a dielectric square covering the entire
        structure.

    add_groundplane_square: bool = True,
        Default True will add a groundplane square covering the entire
        structure.

    add_dielectric_cutout_in_ports: bool | Sequence[bool] = True,
        This is whether to add a dielectric cutout in the ports for the
        linewidths pads. Default True will add a dielectric cutout in the pads
        for all the linewidths. If False then no cutouts will be added. If a
        Sequence of bool is given then it should be of length 6 and each
        line will take the respective cutout from left to right.

    add_dielectric_cutout_over_H: bool | Sequence[bool] = (False, False, False),
        This is whether to add a dielectric cutout over the entire H section
        with the 4 pads total and the line connecting them. Default all False
        will not add this dielectric cutout. If a Sequence of bool is given
        then it should be of length 3 and each H will take the respective
        cutout.

    add_dielectric_cutout_over_line: bool | Sequence[bool] = (False, False, False),
        This is whether to add a dielectric cutout over the line section.
        Default all False will not add this dielectric cutout. If a Sequence of
        bool is given then it should be of length 3 and each line will take the
        respective cutout.
    """
    BOX_SIZE = 5500
    NUMBER_OF_PADS = 6
    PAD_OFFSET_FROM_EDGE = 250

    # Adding the outer square
    ###########################################################################
    outer_lw = 100

    outer_square_points = [
        [x - BOX_SIZE / 2 - outer_lw / 2, y + BOX_SIZE / 2 + outer_lw],
        [x - BOX_SIZE / 2 - outer_lw / 2, y - BOX_SIZE / 2 - outer_lw / 2],
        [x + BOX_SIZE / 2 + outer_lw / 2, y - BOX_SIZE / 2 - outer_lw / 2],
        [x + BOX_SIZE / 2 + outer_lw / 2, y + BOX_SIZE / 2 + outer_lw / 2],
        [x - BOX_SIZE / 2 - outer_lw, y + BOX_SIZE / 2 + outer_lw / 2],
    ]
    outer_square = gdspy.FlexPath(
        outer_square_points,
        outer_lw,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(outer_square)

    # Adding the main square
    ###########################################################################
    # defined bot left going clockwise
    main_square_points = [
        [x - (BOX_SIZE / 2), y - (BOX_SIZE / 2)],
        [x - (BOX_SIZE / 2), y + (BOX_SIZE / 2)],
        [x + (BOX_SIZE / 2), y + (BOX_SIZE / 2)],
        [x + (BOX_SIZE / 2), y - (BOX_SIZE / 2)],
    ]

    if add_dielectric_square:
        SiN_square = gdspy.Polygon(
            main_square_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        SiN_square.rotate(rot, [x, y])
        mask_builder.silicon_nitride_positives.add(SiN_square)

    if add_groundplane_square:
        groundplane_square = gdspy.Polygon(
            main_square_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        groundplane_square.rotate(rot, [x, y])
        mask_builder.ground_plane_positives.add(groundplane_square)

    # settings up all the parameters for the pads
    ###########################################################################

    # Gettings the linewidths
    if not isinstance(linewidths, float | int):
        styled_type_error(linewidths, "linewidths", float | int)

    # Gettings the layers
    pad_layers: list[Layer] = []
    if isinstance(pad_materials, Layer):
        for _ in range(NUMBER_OF_PADS):
            pad_layers.append(pad_materials)

    elif isinstance(pad_materials, Sequence):
        if len(pad_materials) != (NUMBER_OF_PADS):
            raise ValueError(f"The length of `pad_materials` Sequence is {len(pad_materials)} but should be should be {NUMBER_OF_PADS}.")
        if all(isinstance(mat, Layer) for mat in pad_materials):
            for material in pad_materials:
                pad_layers.append(material)
        else:
            styled_type_error(pad_materials, "pad_materials", Layer | Sequence[Layer])
    else:
        styled_type_error(pad_materials, "pad_materials", Layer | Sequence[Layer])

    line_layers: list[Layer] = []
    if isinstance(line_materials, Layer):
        for _ in range(int(NUMBER_OF_PADS / 2)):
            line_layers.append(line_materials)

    elif isinstance(line_materials, Sequence):
        if len(line_materials) != (int(NUMBER_OF_PADS / 2)):
            raise ValueError(
                f"The length of `line_materials` Sequence is {len(line_materials)} but should be should be {NUMBER_OF_PADS/2}."
            )
        if all(isinstance(mat, Layer) for mat in line_materials):
            for material in line_materials:
                line_layers.append(material)
        else:
            styled_type_error(line_materials, "line_materials", Layer | Sequence[Layer])
    else:
        styled_type_error(line_materials, "line_materials", Layer | Sequence[Layer])

    # Getting the cutouts over the H sections
    dielectric_cutout_over_Hs: list[bool] = []
    if isinstance(add_dielectric_cutout_over_H, bool):
        for _ in range(int(NUMBER_OF_PADS / 2)):
            dielectric_cutout_over_Hs.append(add_dielectric_cutout_over_H)
    elif isinstance(add_dielectric_cutout_over_H, Sequence):
        if len(add_dielectric_cutout_over_H) != (int(NUMBER_OF_PADS / 2)):
            raise ValueError(
                f"The length of `add_dielectric_cutout_over_H` Sequence is {len(add_dielectric_cutout_over_H)} but should be should be {NUMBER_OF_PADS/2}."
            )
        if all(isinstance(cut, bool) for cut in add_dielectric_cutout_over_H):
            for cut in add_dielectric_cutout_over_H:
                dielectric_cutout_over_Hs.append(cut)
        else:
            styled_type_error(add_dielectric_cutout_over_H, "add_dielectric_cutout_over_H", bool | Sequence[bool])

    else:
        styled_type_error(add_dielectric_cutout_over_H, "add_dielectric_cutout_over_H", bool | Sequence[bool])

    # Getting the cutouts over the line sections
    dielectric_cutout_over_lines: list[bool] = []
    if isinstance(add_dielectric_cutout_over_line, bool):
        for _ in range(int(NUMBER_OF_PADS / 2)):
            dielectric_cutout_over_lines.append(add_dielectric_cutout_over_line)
    elif isinstance(add_dielectric_cutout_over_H, Sequence):
        if len(add_dielectric_cutout_over_line) != (int(NUMBER_OF_PADS / 2)):
            raise ValueError(
                f"The length of `add_dielectric_cutout_over_line` Sequence is {len(add_dielectric_cutout_over_line)} but should be should be {NUMBER_OF_PADS/2}."
            )
        if all(isinstance(cut, bool) for cut in add_dielectric_cutout_over_line):
            for cut in add_dielectric_cutout_over_line:
                dielectric_cutout_over_lines.append(cut)
        else:
            styled_type_error(add_dielectric_cutout_over_line, "add_dielectric_cutout_over_line", bool | Sequence[bool])

    else:
        styled_type_error(add_dielectric_cutout_over_line, "add_dielectric_cutout_over_line", bool | Sequence[bool])

    # Gettings the dielectric cutouts
    dielectric_cutout_in_ports: list[bool] = []
    if isinstance(add_dielectric_cutout_in_ports, bool):
        for _ in range(NUMBER_OF_PADS):
            dielectric_cutout_in_ports.append(add_dielectric_cutout_in_ports)

    elif isinstance(add_dielectric_cutout_in_ports, Sequence):
        if len(add_dielectric_cutout_in_ports) != NUMBER_OF_PADS:
            raise ValueError(
                f"The length of `add_dielectric_cutout_in_ports` Sequence is {len(add_dielectric_cutout_in_ports)} but should be should be {NUMBER_OF_PADS}."
            )
        if all(isinstance(cutout, bool) for cutout in add_dielectric_cutout_in_ports):
            for cutout in add_dielectric_cutout_in_ports:
                dielectric_cutout_in_ports.append(cutout)
        else:
            styled_type_error(add_dielectric_cutout_in_ports, "add_dielectric_cutout_in_ports", bool | Sequence[bool])
    else:
        styled_type_error(add_dielectric_cutout_in_ports, "add_dielectric_cutout_in_ports", bool | Sequence[bool])

    x_offset_between_pad_cetners = (BOX_SIZE - (2 * PAD_OFFSET_FROM_EDGE)) / (NUMBER_OF_PADS + 0)

    x_offsets = [
        -2.5 * x_offset_between_pad_cetners,
        -1.5 * x_offset_between_pad_cetners,
        -0.5 * x_offset_between_pad_cetners,
        +0.5 * x_offset_between_pad_cetners,
        +1.5 * x_offset_between_pad_cetners,
        +2.5 * x_offset_between_pad_cetners,
    ]

    # Adding the top right text
    ###########################################################################
    if TL_text != "":
        text_x = x - (BOX_SIZE / 2) * cos(rot)
        text_y = y + (BOX_SIZE / 2) * cos(rot)
        text_size = 0.9 * PAD_OFFSET_FROM_EDGE
        text_layer = pad_layers[0]
        mask_builder.add_fancy_text(
            TL_text,
            text_x,
            text_y,
            text_size,
            text_layer,
            rotation=rot,
            vertical_align="below",
        )

    # Adding the pads
    ###########################################################################
    PAD_BACK_HEIGHT = 1000
    PAD_BACK_WIDTH = 500
    PAD_TAPER_HEIGHT = 750

    pad_y_top = y + (BOX_SIZE / 2) - PAD_OFFSET_FROM_EDGE
    pad_y_bot = y - (BOX_SIZE / 2) + PAD_OFFSET_FROM_EDGE

    pad_top_connects: list[tuple[float, float]] = []
    pad_bot_connects: list[tuple[float, float]] = []

    for i, (pad_offset, layer, dielectric_cutout) in enumerate(zip(x_offsets, pad_layers, dielectric_cutout_in_ports, strict=True)):
        pad_x = x + pad_offset

        # defined from top right going clockwise
        top_pad_points = [
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_top],
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_top - PAD_BACK_HEIGHT],
            [pad_x + (linewidths / 2), pad_y_top - PAD_BACK_HEIGHT - PAD_TAPER_HEIGHT],
            [pad_x - (linewidths / 2), pad_y_top - PAD_BACK_HEIGHT - PAD_TAPER_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_top - PAD_BACK_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_top],
        ]
        pad_top_connects.append((pad_x, pad_y_top - PAD_BACK_HEIGHT - PAD_TAPER_HEIGHT))

        pad_top = gdspy.Polygon(
            top_pad_points,
            layer=layer.number,
            datatype=layer.datatype,
        )
        pad_top.rotate(rot, [x, y])
        mask_builder.Main.add(pad_top)

        # defined from bot right going anti-clockwise
        bot_pad_points = [
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_bot],
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_bot + PAD_BACK_HEIGHT],
            [pad_x + (linewidths / 2), pad_y_bot + PAD_BACK_HEIGHT + PAD_TAPER_HEIGHT],
            [pad_x - (linewidths / 2), pad_y_bot + PAD_BACK_HEIGHT + PAD_TAPER_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_bot + PAD_BACK_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_bot],
        ]
        pad_bot_connects.append((pad_x, pad_y_bot + PAD_BACK_HEIGHT + PAD_TAPER_HEIGHT))
        pad_bot = gdspy.Polygon(
            bot_pad_points,
            layer=layer.number,
            datatype=layer.datatype,
        )
        pad_bot.rotate(rot, [x, y])
        mask_builder.Main.add(pad_bot)

        if add_text_for_linewidth:
            if i % 2 == 0:  # Only draw every other pad.
                lw_text = f"{linewidths}um"
                lw_text_x = bot_pad_points[2][0] + (PAD_BACK_WIDTH / 3)
                lw_text_y = bot_pad_points[2][1] + (PAD_BACK_WIDTH / 3)
                lw_text_size = PAD_BACK_WIDTH / 10
                lw_text_layer = layer
                mask_builder.add_fancy_text(
                    lw_text,
                    lw_text_x,
                    lw_text_y,
                    lw_text_size,
                    lw_text_layer,
                    rotation=rot - (pi / 2),
                    vertical_align="above",
                    horizontal_align="end",
                )

        if dielectric_cutout:
            SiN_cutout_offset = 25
            # defined from top right going clockwise
            SiN_top_pad_cutout_points = [
                [
                    top_pad_points[0][0] - SiN_cutout_offset,
                    top_pad_points[0][1] - SiN_cutout_offset,
                ],
                [
                    top_pad_points[1][0] - SiN_cutout_offset,
                    top_pad_points[1][1] + SiN_cutout_offset,
                ],
                [
                    top_pad_points[4][0] + SiN_cutout_offset,
                    top_pad_points[4][1] + SiN_cutout_offset,
                ],
                [
                    top_pad_points[5][0] + SiN_cutout_offset,
                    top_pad_points[5][1] - SiN_cutout_offset,
                ],
            ]
            SiN_top_cutout = gdspy.Polygon(
                SiN_top_pad_cutout_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            SiN_top_cutout.rotate(rot, [x, y])
            mask_builder.silicon_nitride_cutouts.add(SiN_top_cutout)
            SiN_bot_pad_cutout_points = [
                [
                    bot_pad_points[0][0] - SiN_cutout_offset,
                    bot_pad_points[0][1] + SiN_cutout_offset,
                ],
                [
                    bot_pad_points[1][0] - SiN_cutout_offset,
                    bot_pad_points[1][1] - SiN_cutout_offset,
                ],
                [
                    bot_pad_points[4][0] + SiN_cutout_offset,
                    bot_pad_points[4][1] - SiN_cutout_offset,
                ],
                [
                    bot_pad_points[5][0] + SiN_cutout_offset,
                    bot_pad_points[5][1] + SiN_cutout_offset,
                ],
            ]
            SiN_bot_cutout = gdspy.Polygon(
                SiN_bot_pad_cutout_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            SiN_bot_cutout.rotate(rot, [x, y])
            mask_builder.silicon_nitride_cutouts.add(SiN_bot_cutout)

    U_LOOP_HEIGHT = 100
    U_LOOP_BEND_RAD = 70
    JUT_OUT_HEIGHT = 100

    for i in range(int(NUMBER_OF_PADS / 2)):
        # Adding the pads U connect path
        #######################################################################
        top_left_pad = pad_top_connects[2 * i]
        top_right_pad = pad_top_connects[(2 * i) + 1]
        bot_left_pad = pad_bot_connects[2 * i]
        bot_right_pad = pad_bot_connects[(2 * i) + 1]

        U_path_points_for_top = [
            [top_left_pad[0], top_left_pad[1]],
            [top_left_pad[0], top_left_pad[1] - U_LOOP_HEIGHT],
            [top_right_pad[0], top_right_pad[1] - U_LOOP_HEIGHT],
            [top_right_pad[0], top_right_pad[1]],
        ]
        U_path_points_for_bot = [
            [bot_left_pad[0], bot_left_pad[1]],
            [bot_left_pad[0], bot_left_pad[1] + U_LOOP_HEIGHT],
            [bot_right_pad[0], bot_right_pad[1] + U_LOOP_HEIGHT],
            [bot_right_pad[0], bot_right_pad[1]],
        ]
        U_path_top = gdspy.FlexPath(
            U_path_points_for_top,
            linewidths,
            corners="circular bend",
            bend_radius=U_LOOP_BEND_RAD,
            layer=pad_layers[i].number,
            datatype=pad_layers[i].datatype,
        )
        U_path_top.rotate(rot, (x, y))
        mask_builder.make_flexpath_into_polygons_and_add_to_main(U_path_top, pad_layers[i])

        U_path_bot = gdspy.FlexPath(
            U_path_points_for_bot,
            linewidths,
            corners="circular bend",
            bend_radius=U_LOOP_BEND_RAD,
            layer=pad_layers[i].number,
            datatype=pad_layers[i].datatype,
        )
        U_path_bot.rotate(rot, (x, y))
        mask_builder.make_flexpath_into_polygons_and_add_to_main(U_path_bot, pad_layers[i])

        # Adding the U connect path jut out straight
        #######################################################################
        jut_out_path_top_points = [
            [
                top_right_pad[0] + (top_left_pad[0] - top_right_pad[0]) / 2,
                top_left_pad[1] - U_LOOP_HEIGHT - (linewidths / 2),
            ],
            [
                top_right_pad[0] + (top_left_pad[0] - top_right_pad[0]) / 2,
                top_left_pad[1] - U_LOOP_HEIGHT - JUT_OUT_HEIGHT - (linewidths / 2),
            ],
        ]
        jut_out_path_bot_points = [
            [
                bot_right_pad[0] + (bot_left_pad[0] - bot_right_pad[0]) / 2,
                bot_left_pad[1] + U_LOOP_HEIGHT + (linewidths / 2),
            ],
            [
                bot_right_pad[0] + (bot_left_pad[0] - bot_right_pad[0]) / 2,
                bot_left_pad[1] + U_LOOP_HEIGHT + JUT_OUT_HEIGHT + (linewidths / 2),
            ],
        ]
        jut_out_path_top = gdspy.FlexPath(
            jut_out_path_top_points,
            linewidths,
            layer=pad_layers[i].number,
            datatype=pad_layers[i].datatype,
        )
        jut_out_path_top.rotate(rot, (x, y))
        mask_builder.make_flexpath_into_polygons_and_add_to_main(jut_out_path_top, pad_layers[i])

        jut_out_path_bot = gdspy.FlexPath(
            jut_out_path_bot_points,
            linewidths,
            layer=pad_layers[i].number,
            datatype=pad_layers[i].datatype,
        )
        jut_out_path_bot.rotate(rot, (x, y))
        mask_builder.make_flexpath_into_polygons_and_add_to_main(jut_out_path_bot, pad_layers[i])

        # Adding the line connecting the jut outs on U paths
        #######################################################################
        pad_size_x = 2 * linewidths
        pad_size_y_non_overlap = 2 * linewidths
        pad_size_y_overlap = 6 * linewidths

        top_con = [
            top_right_pad[0] + (top_left_pad[0] - top_right_pad[0]) / 2,
            top_left_pad[1] - U_LOOP_HEIGHT - JUT_OUT_HEIGHT - (linewidths / 2),
        ]
        bot_con = [
            bot_right_pad[0] + (bot_left_pad[0] - bot_right_pad[0]) / 2,
            bot_left_pad[1] + U_LOOP_HEIGHT + JUT_OUT_HEIGHT + (linewidths / 2),
        ]

        # defined from top right going clockwise
        line_connecting_points = [
            [top_con[0] + (pad_size_x / 2), top_con[1] + (pad_size_y_overlap / 2)],
            [top_con[0] + (pad_size_x / 2), top_con[1] - (pad_size_y_non_overlap / 2)],
            [top_con[0] + (linewidths / 2), top_con[1] - (pad_size_y_non_overlap / 2)],
            [bot_con[0] + (linewidths / 2), bot_con[1] + (pad_size_y_non_overlap / 2)],
            [bot_con[0] + (pad_size_x / 2), bot_con[1] + (pad_size_y_non_overlap / 2)],
            [bot_con[0] + (pad_size_x / 2), bot_con[1] - (pad_size_y_overlap / 2)],
            [bot_con[0] - (pad_size_x / 2), bot_con[1] - (pad_size_y_overlap / 2)],
            [bot_con[0] - (pad_size_x / 2), bot_con[1] + (pad_size_y_non_overlap / 2)],
            [bot_con[0] - (linewidths / 2), bot_con[1] + (pad_size_y_non_overlap / 2)],
            [top_con[0] - (linewidths / 2), top_con[1] - (pad_size_y_non_overlap / 2)],
            [top_con[0] - (pad_size_x / 2), top_con[1] - (pad_size_y_non_overlap / 2)],
            [top_con[0] - (pad_size_x / 2), top_con[1] + (pad_size_y_overlap / 2)],
        ]
        line_connecting = gdspy.Polygon(
            line_connecting_points,
            layer=line_layers[i].number,
            datatype=line_layers[i].datatype,
        )
        line_connecting.rotate(rot, (x, y))
        mask_builder.Main.add(line_connecting)

        # Adding the dielectric cutout over line section
        #######################################################################
        if dielectric_cutout_over_lines[i]:
            LINE_CUTOUT_OFFSET_Y = 110
            line_cutout_offset_x = 7 * linewidths
            # defined from top left going clockwise
            line_cutout_points = [
                [top_con[0] - line_cutout_offset_x, top_con[1] + LINE_CUTOUT_OFFSET_Y],
                [top_con[0] + line_cutout_offset_x, top_con[1] + LINE_CUTOUT_OFFSET_Y],
                [bot_con[0] + line_cutout_offset_x, bot_con[1] - LINE_CUTOUT_OFFSET_Y],
                [bot_con[0] - line_cutout_offset_x, bot_con[1] - LINE_CUTOUT_OFFSET_Y],
            ]
            line_cutout = gdspy.Polygon(
                line_cutout_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            line_cutout.rotate(rot, (x, y))
            mask_builder.silicon_nitride_cutouts.add(line_cutout)

        # Adding the dielectric cutout over whole H section
        #######################################################################
        if dielectric_cutout_over_Hs[i]:
            # defined from top left going clockwise
            H_CUTOUT_OFFSET_FROM_PORT = 20
            x_offset = (linewidths / 2) + (PAD_BACK_WIDTH / 2) + H_CUTOUT_OFFSET_FROM_PORT
            y_offset = PAD_TAPER_HEIGHT + PAD_BACK_HEIGHT + H_CUTOUT_OFFSET_FROM_PORT
            dielectric_cutout_over_H_points = [
                [
                    top_left_pad[0] - x_offset,
                    top_left_pad[1] + y_offset,
                ],
                [
                    top_right_pad[0] + x_offset,
                    top_right_pad[1] + y_offset,
                ],
                [
                    bot_right_pad[0] + x_offset,
                    bot_right_pad[1] - y_offset,
                ],
                [
                    bot_left_pad[0] - x_offset,
                    bot_left_pad[1] - y_offset,
                ],
            ]
            dielectric_cutout_over_H = gdspy.Polygon(
                dielectric_cutout_over_H_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            dielectric_cutout_over_H.rotate(rot, (x, y))
            mask_builder.silicon_nitride_cutouts.add(dielectric_cutout_over_H)
    return


def add_test_linewidths_pad_connected_box_section(
    mask_builder,
    x: float | int,
    y: float | int,
    materials: Layer | Sequence[Layer],
    rot: float | int = 0,
    linewidths: Sequence[float | int] | None = None,
    TL_text: str = "",
    add_text_for_linewidth: bool = True,
    add_dielectric_square: bool = True,
    add_groundplane_square: bool = True,
    add_dielectric_cutout_in_ports: bool | Sequence[bool] = True,
) -> None:
    """Adds a box of size 5500 centered on the x,y given that contains a series
    of 6 lines and pads connecting them. The test linewidths by default are (1.0,
    2.0, 3.0, 4.0, 5.0, 10.0).

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test structure.

    materials: Layer | Sequence[Layer]
        This is the layers that each of the lines and pads should be drawn in.
        This is an instance of Layer of a Sequence of many Layers. If a single
        Layer is given then all the lines and ports will be drawn in this layer.
        If a Sequence of Layers is given then it should be of length 6 and each
        line will take the respective layer from left to right. For more info
        see maskpy.layers.Layer. Usually this is within the
        SoukMaskBuilder.layers.xxx. e.g. `mask_builder.layers.Aluminium`

    KwArgs
    ------
    rot: float | int = 0
        The rotation for the entire structure given in radians.

    linewidths: Sequence[float | int] | None = None,
        This is the linewidths to use for each of the lines from left to right.
        Default None will assign the default (1.0, 2.0, 3.0, 4.0, 5.0, 10.0).
        When defined this should be a Sequence of float or int of length 6.

    TL_text: str = "",
        The optional string to write in the top left corner of the structure.

    add_text_for_linewidth: bool = True,
        This is whether to add text denoting each linewidth next to the line.
        Default True will add the linewidth of the current line for each line.
        This is disabled with False. When True this is drawn in the same layer
        as the linewidths and pads.

    add_dielectric_square: bool = True,
        Default True will add a dielectric square covering the entire
        structure.

    add_groundplane_square: bool = True,
        Default True will add a groundplane square covering the entire
        structure.

    add_dielectric_cutout_in_ports: bool | Sequence[bool] = True,
        This is whether to add a dielectric cutout in the ports for the
        linewidths pads. Default True will add a dielectric cutout in the pads
        for all the linewidths. If False then no cutouts will be added. If a
        Sequence of bool is given then it should be of length 6 and each
        line will take the respective cutout from left to right.
    """
    BOX_SIZE = 5500
    NUMBER_OF_PADS = 6
    PAD_OFFSET_FROM_EDGE = 250
    DEFAULT_LINEWIDTHS = (1.0, 2.0, 3.0, 4.0, 5.0, 10.0)

    # Adding the outer square
    ###########################################################################
    outer_lw = 100

    outer_square_points = [
        [x - BOX_SIZE / 2 - outer_lw / 2, y + BOX_SIZE / 2 + outer_lw],
        [x - BOX_SIZE / 2 - outer_lw / 2, y - BOX_SIZE / 2 - outer_lw / 2],
        [x + BOX_SIZE / 2 + outer_lw / 2, y - BOX_SIZE / 2 - outer_lw / 2],
        [x + BOX_SIZE / 2 + outer_lw / 2, y + BOX_SIZE / 2 + outer_lw / 2],
        [x - BOX_SIZE / 2 - outer_lw, y + BOX_SIZE / 2 + outer_lw / 2],
    ]
    outer_square = gdspy.FlexPath(
        outer_square_points,
        outer_lw,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(outer_square)

    # Adding the main square
    ###########################################################################
    # defined bot left going clockwise
    main_square_points = [
        [x - (BOX_SIZE / 2), y - (BOX_SIZE / 2)],
        [x - (BOX_SIZE / 2), y + (BOX_SIZE / 2)],
        [x + (BOX_SIZE / 2), y + (BOX_SIZE / 2)],
        [x + (BOX_SIZE / 2), y - (BOX_SIZE / 2)],
    ]

    if add_dielectric_square:
        SiN_square = gdspy.Polygon(
            main_square_points,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        SiN_square.rotate(rot, [x, y])
        mask_builder.silicon_nitride_positives.add(SiN_square)

    if add_groundplane_square:
        groundplane_square = gdspy.Polygon(
            main_square_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        groundplane_square.rotate(rot, [x, y])
        mask_builder.ground_plane_positives.add(groundplane_square)

    # settings up all the parameters for the pads
    ###########################################################################

    # Gettings the linewidths
    if linewidths is None:
        linewidths = DEFAULT_LINEWIDTHS
    elif isinstance(linewidths, Sequence):
        if len(linewidths) != NUMBER_OF_PADS:
            raise ValueError(f"The length of `linewidths` Sequence is {len(linewidths)} but should be should be {NUMBER_OF_PADS}.")
        if not all(isinstance(lw, float | int) for lw in linewidths):
            styled_type_error(linewidths, "linewidths", Sequence[float | int])
    else:
        styled_type_error(linewidths, "linewidths", Sequence[float | int])

    # Gettings the layers
    layers: list[Layer] = []
    if isinstance(materials, Layer):
        for _ in range(NUMBER_OF_PADS):
            layers.append(materials)

    elif isinstance(materials, Sequence):
        if len(materials) != NUMBER_OF_PADS:
            raise ValueError(f"The length of `materials` Sequence is {len(materials)} but should be should be {NUMBER_OF_PADS}.")
        if all(isinstance(mat, Layer) for mat in materials):
            for material in materials:
                layers.append(material)
        else:
            styled_type_error(materials, "materials", Layer | Sequence[Layer])
    else:
        styled_type_error(materials, "materials", Layer | Sequence[Layer])

    # Gettings the dielectric cutouts
    dielectric_cutout_in_ports: list[bool] = []
    if isinstance(add_dielectric_cutout_in_ports, bool):
        for _ in range(NUMBER_OF_PADS):
            dielectric_cutout_in_ports.append(add_dielectric_cutout_in_ports)

    elif isinstance(add_dielectric_cutout_in_ports, Sequence):
        if len(add_dielectric_cutout_in_ports) != NUMBER_OF_PADS:
            raise ValueError(
                f"The length of `add_dielectric_cutout_in_ports` Sequence is {len(add_dielectric_cutout_in_ports)} but should be should be {NUMBER_OF_PADS}."
            )
        if all(isinstance(cutout, bool) for cutout in add_dielectric_cutout_in_ports):
            for cutout in add_dielectric_cutout_in_ports:
                dielectric_cutout_in_ports.append(cutout)
        else:
            styled_type_error(add_dielectric_cutout_in_ports, "add_dielectric_cutout_in_ports", bool | Sequence[bool])
    else:
        styled_type_error(add_dielectric_cutout_in_ports, "add_dielectric_cutout_in_ports", bool | Sequence[bool])

    x_offset_between_pad_cetners = (BOX_SIZE - (2 * PAD_OFFSET_FROM_EDGE)) / (NUMBER_OF_PADS + 0)

    x_offsets = [
        -2.5 * x_offset_between_pad_cetners,
        -1.5 * x_offset_between_pad_cetners,
        -0.5 * x_offset_between_pad_cetners,
        +0.5 * x_offset_between_pad_cetners,
        +1.5 * x_offset_between_pad_cetners,
        +2.5 * x_offset_between_pad_cetners,
    ]

    # Adding the top right text
    ###########################################################################
    if TL_text != "":
        text_x = x - (BOX_SIZE / 2) * cos(rot)
        text_y = y + (BOX_SIZE / 2) * cos(rot)
        text_size = 0.9 * PAD_OFFSET_FROM_EDGE
        text_layer = layers[0]
        mask_builder.add_fancy_text(
            TL_text,
            text_x,
            text_y,
            text_size,
            text_layer,
            rotation=rot,
            vertical_align="below",
        )

    # Adding the pads
    ###########################################################################
    PAD_BACK_HEIGHT = 1000
    PAD_BACK_WIDTH = 500
    PAD_TAPER_HEIGHT = 750

    pad_y_top = y + (BOX_SIZE / 2) - PAD_OFFSET_FROM_EDGE
    pad_y_bot = y - (BOX_SIZE / 2) + PAD_OFFSET_FROM_EDGE

    for pad_offset, layer, lw, dielectric_cutout in zip(x_offsets, layers, linewidths, dielectric_cutout_in_ports, strict=True):
        pad_x = x + pad_offset

        # defined from top right going clockwise
        top_pad_points = [
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_top],
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_top - PAD_BACK_HEIGHT],
            [pad_x + (lw / 2), pad_y_top - PAD_BACK_HEIGHT - PAD_TAPER_HEIGHT],
            [pad_x - (lw / 2), pad_y_top - PAD_BACK_HEIGHT - PAD_TAPER_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_top - PAD_BACK_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_top],
        ]
        pad_top = gdspy.Polygon(
            top_pad_points,
            layer=layer.number,
            datatype=layer.datatype,
        )
        pad_top.rotate(rot, [x, y])
        mask_builder.Main.add(pad_top)

        # defined from bot right going anti-clockwise
        bot_pad_points = [
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_bot],
            [pad_x + (PAD_BACK_WIDTH / 2), pad_y_bot + PAD_BACK_HEIGHT],
            [pad_x + (lw / 2), pad_y_bot + PAD_BACK_HEIGHT + PAD_TAPER_HEIGHT],
            [pad_x - (lw / 2), pad_y_bot + PAD_BACK_HEIGHT + PAD_TAPER_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_bot + PAD_BACK_HEIGHT],
            [pad_x - (PAD_BACK_WIDTH / 2), pad_y_bot],
        ]
        pad_bot = gdspy.Polygon(
            bot_pad_points,
            layer=layer.number,
            datatype=layer.datatype,
        )
        pad_bot.rotate(rot, [x, y])
        mask_builder.Main.add(pad_bot)

        # defined from the top right going anti-clockwise
        line_connecting_pads_points = [
            top_pad_points[2],
            top_pad_points[3],
            bot_pad_points[3],
            bot_pad_points[2],
        ]
        line_connecting_pads = gdspy.Polygon(
            line_connecting_pads_points,
            layer=layer.number,
            datatype=layer.datatype,
        )
        line_connecting_pads.rotate(rot, [x, y])
        mask_builder.Main.add(line_connecting_pads)

        if add_text_for_linewidth:
            lw_text = f"{lw}um"
            lw_text_x = bot_pad_points[2][0] + (PAD_BACK_WIDTH / 3)
            lw_text_y = bot_pad_points[2][1] + (PAD_BACK_WIDTH / 3)
            lw_text_size = PAD_BACK_WIDTH / 10
            lw_text_layer = layer
            mask_builder.add_fancy_text(
                lw_text,
                lw_text_x,
                lw_text_y,
                lw_text_size,
                lw_text_layer,
                rotation=rot - (pi / 2),
                vertical_align="above",
                horizontal_align="end",
            )

        if dielectric_cutout:
            SiN_cutout_offset = 25
            # defined from top right going clockwise
            SiN_top_pad_cutout_points = [
                [
                    top_pad_points[0][0] - SiN_cutout_offset,
                    top_pad_points[0][1] - SiN_cutout_offset,
                ],
                [
                    top_pad_points[1][0] - SiN_cutout_offset,
                    top_pad_points[1][1] + SiN_cutout_offset,
                ],
                [
                    top_pad_points[4][0] + SiN_cutout_offset,
                    top_pad_points[4][1] + SiN_cutout_offset,
                ],
                [
                    top_pad_points[5][0] + SiN_cutout_offset,
                    top_pad_points[5][1] - SiN_cutout_offset,
                ],
            ]
            SiN_top_cutout = gdspy.Polygon(
                SiN_top_pad_cutout_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            SiN_top_cutout.rotate(rot, [x, y])
            mask_builder.silicon_nitride_cutouts.add(SiN_top_cutout)
            SiN_bot_pad_cutout_points = [
                [
                    bot_pad_points[0][0] - SiN_cutout_offset,
                    bot_pad_points[0][1] + SiN_cutout_offset,
                ],
                [
                    bot_pad_points[1][0] - SiN_cutout_offset,
                    bot_pad_points[1][1] - SiN_cutout_offset,
                ],
                [
                    bot_pad_points[4][0] + SiN_cutout_offset,
                    bot_pad_points[4][1] - SiN_cutout_offset,
                ],
                [
                    bot_pad_points[5][0] + SiN_cutout_offset,
                    bot_pad_points[5][1] + SiN_cutout_offset,
                ],
            ]
            SiN_bot_cutout = gdspy.Polygon(
                SiN_bot_pad_cutout_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            SiN_bot_cutout.rotate(rot, [x, y])
            mask_builder.silicon_nitride_cutouts.add(SiN_bot_cutout)

    return


def add_test_linewidth_structure_box_section(
    mask_builder,
    x: float | int,
    y: float | int,
    linewidths: Sequence[float | int] = (0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10),
) -> None:
    """
    Adds a box of size 5500 centered on the x,y given that contains a
    series of lines to test the linewidth. Test linewidths are [0.25, 0.5,
    0.75, 1, 1.5, 2, 3, 5, 10]. Materials are ["Aluminium", "Nb_Antenna",
    "Nb_grnd", "SiN_Dep","IDC_Nb"].

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test linewidth structure.

    KwArgs
    ------
    linewidths: Sequence[float | int]
        This is a Sequence of float or ints that define the linewidths in the
        test structure. By default these linewiths are:
        >>> [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10].
    """

    box_w_h = 5500
    main_square = gdspy.Rectangle(
        [x - box_w_h / 2, y - box_w_h / 2],
        [x + box_w_h / 2, y + box_w_h / 2],
        layer=mask_builder.layers.Nb_Groundplane.number,
        datatype=mask_builder.layers.Nb_Groundplane.datatype,
    )

    mask_builder.ground_plane_cutouts.add(main_square)

    outer_square_side_length = box_w_h
    outer_lw = 100

    outer_square_points = [
        [x - outer_square_side_length / 2 - outer_lw / 2, y + outer_square_side_length / 2 + outer_lw],
        [x - outer_square_side_length / 2 - outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y + outer_square_side_length / 2 + outer_lw / 2],
        [x - outer_square_side_length / 2 - outer_lw, y + outer_square_side_length / 2 + outer_lw / 2],
    ]
    outer_square = gdspy.FlexPath(
        outer_square_points,
        outer_lw,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(outer_square)

    GAP = 50
    BIG_GAP = 200

    # linewidths = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10]

    HEIGHT = 4000
    LW_TEXT_SIZE = 36
    LABEL_TEXT_HEIGHT = 5 * LW_TEXT_SIZE

    materials = [
        mask_builder.layers.Aluminium,
        mask_builder.layers.Nb_Antenna,
        mask_builder.layers.Nb_Groundplane,
        mask_builder.layers.SiN_dep,
        mask_builder.layers.IDC_Nb,
    ]

    y_pos = y
    init_x = (
        x
        - (
            len(materials) * (np.sum(linewidths) + (len(linewidths) - 1) * GAP)
            + (len(materials) - 1) * (BIG_GAP)
            + len(materials) * (9.5 / 9) * LABEL_TEXT_HEIGHT
        )
        / 2
        + ((9.5 / 9) * LABEL_TEXT_HEIGHT + 0.125)
    )

    for j, material in enumerate(materials):
        x_pos = init_x + j * (BIG_GAP + (len(linewidths) - 1) * 50 + np.sum(linewidths) + (9.5 / 9) * LABEL_TEXT_HEIGHT)

        material_name = material.name

        for i, lw in enumerate(linewidths):
            # add text label for first linewidth
            if i == 0:
                mask_builder.add_fancy_text(
                    material_name,
                    x_pos - GAP,
                    y_pos - HEIGHT / 2,
                    LABEL_TEXT_HEIGHT,
                    material,
                    rotation=pi / 2,
                    vertical_align="above",
                )

            mask_builder.Main.add(
                gdspy.Rectangle(
                    [x_pos - lw / 2, y_pos - HEIGHT / 2],
                    [x_pos + lw / 2, y_pos + HEIGHT / 2],
                    layer=material.number,
                    datatype=material.datatype,
                )
            )

            mask_builder.add_fancy_text(
                f"{lw}um",
                x_pos,
                y_pos - HEIGHT / 2 - (LW_TEXT_SIZE / 2),
                LW_TEXT_SIZE,
                material,
                rotation=-pi / 2,
                vertical_align="center",
            )

            if i < len(linewidths) - 1:
                x_pos += GAP + lw / 2 + linewidths[i + 1] / 2

    return


def add_test_crossover_structure_box_section(
    mask_builder,
    x: float | int,
    y: float | int,
    filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None,
    rot: float | int = 0.0,
) -> None:
    """Adds a box of size 5500 centered on the x,y given that contains a series
    of 3 crossover sections that connect to pads.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test crossover structure.

    KwArgs
    ------
    filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    rot: float, int
        The angle (**in Radians**) the structure should be rotated.
        Default is 0 which has the pads to connect to the crossovers on
        the left and right. This will only take the values in the range
        0 to pi/2. Other values will cause issues with the connections
        to the crossovers.
    """

    box_w_h = 5500

    main_square = gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2])
    main_square.rotate(rot, [x, y])

    working_groundplane_cutout = main_square
    # working_sin_dep_cutout = main_square
    mask_builder.silicon_nitride_positives.add(
        gdspy.Rectangle(
            [x - box_w_h / 2, y - box_w_h / 2],
            [x + box_w_h / 2, y + box_w_h / 2],
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
    )

    mask_builder.ground_plane_positives.add(
        gdspy.Rectangle(
            [x - box_w_h / 2, y - box_w_h / 2],
            [x + box_w_h / 2, y + box_w_h / 2],
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
    )

    outer_square_side_length = box_w_h
    outer_lw = 100

    outer_square_points = [
        [x - outer_square_side_length / 2 - outer_lw / 2, y + outer_square_side_length / 2 + outer_lw],
        [x - outer_square_side_length / 2 - outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y + outer_square_side_length / 2 + outer_lw / 2],
        [x - outer_square_side_length / 2 - outer_lw, y + outer_square_side_length / 2 + outer_lw / 2],
    ]
    outer_square = gdspy.FlexPath(
        outer_square_points,
        outer_lw,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    outer_square.rotate(rot, [x, y])
    mask_builder.Main.add(outer_square)

    pad_back_length = 1000
    pad_back_width = 500
    pad_taper_length = 750
    conect_linewidth = 5

    pad_cutout_length = 950
    pad_cutout_width = 450

    pad_left_poly_points = [
        [0, pad_back_width / 2],
        [pad_back_length, pad_back_width / 2],
        [pad_back_length + pad_taper_length, conect_linewidth / 2],
        [pad_back_length + pad_taper_length, -conect_linewidth / 2],
        [pad_back_length, -pad_back_width / 2],
        [0, -pad_back_width / 2],
    ]

    pad_left_connect = [
        (pad_left_poly_points[2][0] + pad_left_poly_points[3][0]) / 2,
        (pad_left_poly_points[2][1] + pad_left_poly_points[3][1]) / 2,
    ]
    pad_left_center = [(pad_back_length + 0) / 2, 0]

    pad_left_cutout_poly_points = [
        [pad_left_center[0] - (pad_cutout_length / 2), pad_left_center[1] - (pad_cutout_width / 2)],
        [pad_left_center[0] - (pad_cutout_length / 2), pad_left_center[1] + (pad_cutout_width / 2)],
        [pad_left_center[0] + (pad_cutout_length / 2), pad_left_center[1] + (pad_cutout_width / 2)],
        [pad_left_center[0] + (pad_cutout_length / 2), pad_left_center[1] - (pad_cutout_width / 2)],
    ]

    pad_right_poly_points = [
        [0, pad_back_width / 2],
        [-pad_back_length, pad_back_width / 2],
        [-pad_back_length - pad_taper_length, conect_linewidth / 2],
        [-pad_back_length - pad_taper_length, -conect_linewidth / 2],
        [-pad_back_length, -pad_back_width / 2],
        [0, -pad_back_width / 2],
    ]

    pad_right_connect = [
        (pad_right_poly_points[2][0] + pad_right_poly_points[3][0]) / 2,
        (pad_right_poly_points[2][1] + pad_right_poly_points[3][1]) / 2,
    ]
    pad_right_center = [(-pad_back_length - 0) / 2, 0]

    pad_right_cutout_poly_points = [
        [pad_right_center[0] - (pad_cutout_length / 2), pad_right_center[1] - (pad_cutout_width / 2)],
        [pad_right_center[0] - (pad_cutout_length / 2), pad_right_center[1] + (pad_cutout_width / 2)],
        [pad_right_center[0] + (pad_cutout_length / 2), pad_right_center[1] + (pad_cutout_width / 2)],
        [pad_right_center[0] + (pad_cutout_length / 2), pad_right_center[1] - (pad_cutout_width / 2)],
    ]

    top_offset = 1950
    mid_offset = 0
    bot_offset = -1950

    set_y_offset = [top_offset, mid_offset, bot_offset]

    ground_under_crossover_w_h = 1000
    cross_angle = pi / 4

    vert_spacing = 400
    left_x_offset = -2500
    right_x_offset = 2500

    TL_offset = [x + left_x_offset, y + (vert_spacing / 2 + pad_back_width / 2)]
    TR_offset = [x + right_x_offset, y + (vert_spacing / 2 + pad_back_width / 2)]
    BL_offset = [x + left_x_offset, y - (vert_spacing / 2 + pad_back_width / 2)]
    BR_offset = [x + right_x_offset, y - (vert_spacing / 2 + pad_back_width / 2)]

    pads_offset = [TL_offset, TR_offset, BL_offset, BR_offset]

    pad_no = 0

    for set_no in range(3):
        cross_xy = [x, y + set_y_offset[set_no]]
        cross_xy = mbu.rotate(x, y, cross_xy[0], cross_xy[1], rot)
        cross_connections = mask_builder.add_filter_bank_ring_overlap_and_get_conections(
            cross_xy[0],
            cross_xy[1],
            cross_xy[0] + 100,
            cross_xy[1],
            0,
            cross_angle + rot,
            filter_bank_ring_overlap_config_override=filter_bank_ring_overlap_config_override,
        )
        ground_under_cross = gdspy.Rectangle(
            [cross_xy[0] - (ground_under_crossover_w_h / 2), cross_xy[1] - (ground_under_crossover_w_h / 2)],
            [cross_xy[0] + (ground_under_crossover_w_h / 2), cross_xy[1] + (ground_under_crossover_w_h / 2)],
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        ground_under_cross.rotate(rot, cross_xy)

        working_groundplane_cutout = gdspy.boolean(working_groundplane_cutout, ground_under_cross, "not")

        if rot < pi / 4:
            cross_connects = [
                cross_connections["inner_conect_1"],
                cross_connections["outer_conect_1"],
                cross_connections["inner_conect_0"],
                cross_connections["outer_conect_0"],
            ]
        elif rot < pi / 2:
            cross_connects = [
                cross_connections["outer_conect_0"],
                cross_connections["outer_conect_1"],
                cross_connections["inner_conect_0"],
                cross_connections["inner_conect_1"],
            ]
        elif rot < 3 * pi / 4:
            cross_connects = [
                cross_connections["inner_conect_0"],
                cross_connections["inner_conect_1"],
                cross_connections["outer_conect_0"],
                cross_connections["outer_conect_1"],
            ]
        else:
            raise RuntimeError

        cross_connect_rots = [rot + 3 * pi / 4, rot + pi / 4, rot - 3 * pi / 4, rot - pi / 4]

        for i in range(4):
            pad_no += 1

            pad_points = pad_right_poly_points if pad_no % 2 == 0 else pad_left_poly_points
            pad_cutout_points = pad_right_cutout_poly_points if pad_no % 2 == 0 else pad_left_cutout_poly_points
            pad_connect = pad_right_connect if pad_no % 2 == 0 else pad_left_connect
            pad_connect_rot = rot + pi if pad_no % 2 == 0 else rot

            pad_points = mbu.move_points_list(pad_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
            pad_cutout_points = mbu.move_points_list(pad_cutout_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
            pad_connect = mbu.move_single_point(pad_connect, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])

            pad_points = mbu.rotate_points_list(pad_points, rot, x, y)
            pad_cutout_points = mbu.rotate_points_list(pad_cutout_points, rot, x, y)
            pad_connect = mbu.rotate(x, y, pad_connect[0], pad_connect[1], rot)

            pad_poly = gdspy.Polygon(
                pad_points,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            pad_cutout_poly = gdspy.Polygon(
                pad_cutout_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )

            mask_builder.Main.add(pad_poly)
            mask_builder.silicon_nitride_cutouts.add(pad_cutout_poly)

            pad1_path = gdspy.FlexPath(
                [
                    [pad_connect[0], pad_connect[1]],
                    [pad_connect[0] + 100 * cos(pad_connect_rot), pad_connect[1] + 100 * sin(pad_connect_rot)],
                    [cross_connects[i][0] + 100 * cos(cross_connect_rots[i]), cross_connects[i][1] + 100 * sin(cross_connect_rots[i])],
                    [cross_connects[i][0], cross_connects[i][1]],
                ],
                5,
                corners="circular bend",
                bend_radius=50,
                gdsii_path=True,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.Main.add(pad1_path)

    final_groundplane_cutout = gdspy.boolean(working_groundplane_cutout, main_square, "and")
    mask_builder.ground_plane_cutouts.add(final_groundplane_cutout)

    return


def add_test_straight_line_structure_box_section(
    mask_builder,
    x: float | int,
    y: float | int,
    rot: float | int = 0.0,
):
    """Adds a box of size 5500 centered on the x,y given that contains a series
    of 6 stright line sections that connect to pads.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test crossover structure.

    KwArgs
    ------
    rot: float, int
        The angle (**in Radians**) the structure should be rotated.
        Default is 0 which has the pads on the left and right. This will
        only take the values in the range 0 to pi/2.
    """

    box_w_h = 5500

    main_square = gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2])
    main_square.rotate(rot, [x, y])

    working_groundplane_cutout = main_square
    # working_sin_dep_cutout = main_square
    mask_builder.silicon_nitride_positives.add(
        gdspy.Rectangle(
            [x - box_w_h / 2, y - box_w_h / 2],
            [x + box_w_h / 2, y + box_w_h / 2],
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
    )

    outer_square_side_length = box_w_h
    outer_lw = 100

    outer_square_points = [
        [x - outer_square_side_length / 2 - outer_lw / 2, y + outer_square_side_length / 2 + outer_lw],
        [x - outer_square_side_length / 2 - outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y + outer_square_side_length / 2 + outer_lw / 2],
        [x - outer_square_side_length / 2 - outer_lw, y + outer_square_side_length / 2 + outer_lw / 2],
    ]
    outer_square = gdspy.FlexPath(
        outer_square_points,
        outer_lw,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    outer_square.rotate(rot, [x, y])
    mask_builder.Main.add(outer_square)

    pad_back_length = 1000
    pad_back_width = 500
    pad_taper_length = 750
    conect_linewidth = 5

    pad_cutout_length = 950
    pad_cutout_width = 450

    pad_left_poly_points = [
        [0, pad_back_width / 2],
        [pad_back_length, pad_back_width / 2],
        [pad_back_length + pad_taper_length, conect_linewidth / 2],
        [pad_back_length + pad_taper_length, -conect_linewidth / 2],
        [pad_back_length, -pad_back_width / 2],
        [0, -pad_back_width / 2],
    ]

    pad_left_connect = [
        (pad_left_poly_points[2][0] + pad_left_poly_points[3][0]) / 2,
        (pad_left_poly_points[2][1] + pad_left_poly_points[3][1]) / 2,
    ]
    pad_left_center = [(pad_back_length + 0) / 2, 0]

    pad_left_cutout_poly_points = [
        [pad_left_center[0] - (pad_cutout_length / 2), pad_left_center[1] - (pad_cutout_width / 2)],
        [pad_left_center[0] - (pad_cutout_length / 2), pad_left_center[1] + (pad_cutout_width / 2)],
        [pad_left_center[0] + (pad_cutout_length / 2), pad_left_center[1] + (pad_cutout_width / 2)],
        [pad_left_center[0] + (pad_cutout_length / 2), pad_left_center[1] - (pad_cutout_width / 2)],
    ]

    pad_right_poly_points = [
        [0, pad_back_width / 2],
        [-pad_back_length, pad_back_width / 2],
        [-pad_back_length - pad_taper_length, conect_linewidth / 2],
        [-pad_back_length - pad_taper_length, -conect_linewidth / 2],
        [-pad_back_length, -pad_back_width / 2],
        [0, -pad_back_width / 2],
    ]

    pad_right_connect = [
        (pad_right_poly_points[2][0] + pad_right_poly_points[3][0]) / 2,
        (pad_right_poly_points[2][1] + pad_right_poly_points[3][1]) / 2,
    ]
    pad_right_center = [(-pad_back_length - 0) / 2, 0]

    pad_right_cutout_poly_points = [
        [pad_right_center[0] - (pad_cutout_length / 2), pad_right_center[1] - (pad_cutout_width / 2)],
        [pad_right_center[0] - (pad_cutout_length / 2), pad_right_center[1] + (pad_cutout_width / 2)],
        [pad_right_center[0] + (pad_cutout_length / 2), pad_right_center[1] + (pad_cutout_width / 2)],
        [pad_right_center[0] + (pad_cutout_length / 2), pad_right_center[1] - (pad_cutout_width / 2)],
    ]

    top_offset = 1950
    mid_offset = 0
    bot_offset = -1950

    set_y_offset = [top_offset, mid_offset, bot_offset]

    vert_spacing = 400
    left_x_offset = -2500
    right_x_offset = 2500

    TL_offset = [x + left_x_offset, y + (vert_spacing / 2 + pad_back_width / 2)]
    TR_offset = [x + right_x_offset, y + (vert_spacing / 2 + pad_back_width / 2)]
    BL_offset = [x + left_x_offset, y - (vert_spacing / 2 + pad_back_width / 2)]
    BR_offset = [x + right_x_offset, y - (vert_spacing / 2 + pad_back_width / 2)]

    pads_offset = [TL_offset, TR_offset, BL_offset, BR_offset]

    pad_no = 0

    straight_line_length = 200
    straight_line_width = 5

    al_pad_length = 20
    al_pad_width = 10
    al_straight_line_length = (right_x_offset - left_x_offset) - 2 * (pad_back_length + pad_taper_length) - 2 * (straight_line_length)
    al_straight_line_width = 2
    for set_no in range(3):
        for i in range(4):
            pad_no += 1

            pad_points = pad_right_poly_points if pad_no % 2 == 0 else pad_left_poly_points
            pad_cutout_points = pad_right_cutout_poly_points if pad_no % 2 == 0 else pad_left_cutout_poly_points
            pad_connect = pad_right_connect if pad_no % 2 == 0 else pad_left_connect
            pad_connect_rot = rot + pi if pad_no % 2 == 0 else rot

            pad_points = mbu.move_points_list(pad_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
            pad_cutout_points = mbu.move_points_list(pad_cutout_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
            pad_connect = mbu.move_single_point(pad_connect, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])

            pad_points = mbu.rotate_points_list(pad_points, rot, x, y)
            pad_cutout_points = mbu.rotate_points_list(pad_cutout_points, rot, x, y)
            pad_connect = mbu.rotate(x, y, pad_connect[0], pad_connect[1], rot)

            straight_line_path_points = [
                pad_connect,
                [
                    pad_connect[0] + straight_line_length * cos(pad_connect_rot),
                    pad_connect[1] + straight_line_length * sin(pad_connect_rot),
                ],
            ]
            straight_line_path = gdspy.FlexPath(
                straight_line_path_points,
                straight_line_width,
                gdsii_path=True,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            mask_builder.Main.add(straight_line_path)

            # adding aluminium pad over the connection
            al_pad_center = straight_line_path_points[1]
            al_pad_rect = gdspy.Rectangle(
                [
                    al_pad_center[0] + (al_pad_width / 2) * sin(pad_connect_rot),
                    al_pad_center[1] + (al_pad_length / 2) * sin(pad_connect_rot),
                ],
                [
                    al_pad_center[0] - (al_pad_width / 2) * sin(pad_connect_rot),
                    al_pad_center[1] - (al_pad_length / 2) * sin(pad_connect_rot),
                ],
                layer=mask_builder.layers.Aluminium.number,
                datatype=mask_builder.layers.Aluminium.datatype,
            )
            mask_builder.Main.add(al_pad_rect)

            al_straight_line_path = gdspy.FlexPath(
                [
                    al_pad_center,
                    [
                        al_pad_center[0] + al_straight_line_length * cos(pad_connect_rot),
                        al_pad_center[1] + al_straight_line_length * sin(pad_connect_rot),
                    ],
                ],
                al_straight_line_width,
                gdsii_path=True,
                layer=mask_builder.layers.Aluminium.number,
                datatype=mask_builder.layers.Aluminium.datatype,
            )
            al_straight_line_polys = mbu.get_polys_from_flexpath(al_straight_line_path)
            for i in range(len(al_straight_line_polys)):
                mask_builder.Main.add(
                    gdspy.Polygon(
                        al_straight_line_polys[i],
                        layer=mask_builder.layers.Aluminium.number,
                        datatype=mask_builder.layers.Aluminium.datatype,
                    )
                )

            mask_builder.Main.add(al_straight_line_path)

            pad_poly = gdspy.Polygon(
                pad_points,
                layer=mask_builder.layers.Nb_Antenna.number,
                datatype=mask_builder.layers.Nb_Antenna.datatype,
            )
            pad_cutout_poly = gdspy.Polygon(
                pad_cutout_points, layer=mask_builder.layers.SiN_dep.number, datatype=mask_builder.layers.SiN_dep.datatype
            )

            mask_builder.Main.add(pad_poly)
            mask_builder.silicon_nitride_cutouts.add(pad_cutout_poly)

    final_groundplane_cutout = gdspy.boolean(working_groundplane_cutout, main_square, "and")
    mask_builder.ground_plane_cutouts.add(final_groundplane_cutout)

    return


def make_DC_structure_pads_and_meander(
    mask_builder,
    x: float | int,
    y: float | int,
    rot: float | int,
    DC_structure_material: Layer,
) -> None:
    """This Makes two pads and a meander connecting them. Pads are 100x the
    size of the linewidth connecting them.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test DC structure.

    rot: float, int
        The angle (**in degrees**) the pad and meander structure should be
        rotated.

    DC_structure_material: Layer
        This is an instance of Layer. see maskpy.layers.Layer.
        Usually this is within the SoukMaskBuilder.layers.xxx.
        e.g. `self.layers.Aluminium`
    """
    pad_size = 1000
    meander_linewidth = 10

    meander_bend_height_bot = 50
    meander_bend_height_top = 50

    meander_length_mid = 1500
    meander_length_bot = 1700
    meander_length_top = 1700

    total_width_meander = meander_length_bot - meander_length_mid + meander_length_top

    pad_center_x_right = (total_width_meander / 2) + (pad_size / 2)
    pad_center_x_left = -pad_center_x_right

    pad_left = gdspy.Rectangle([pad_center_x_left - (pad_size / 2), -(pad_size / 2)], [pad_center_x_left + (pad_size / 2), +(pad_size / 2)])

    pad_right = gdspy.Rectangle(
        [pad_center_x_right - (pad_size / 2), -(pad_size / 2)], [pad_center_x_right + (pad_size / 2), +(pad_size / 2)]
    )

    meander_points = [
        [-(total_width_meander / 2), -meander_bend_height_bot],
        [-(total_width_meander / 2) + meander_length_bot, -meander_bend_height_bot],
        [-(total_width_meander / 2) + meander_length_bot, 0],
        [(total_width_meander / 2) - meander_length_top, 0],
        [(total_width_meander / 2) - meander_length_top, meander_bend_height_top],
        [(total_width_meander / 2), meander_bend_height_top],
    ]

    meander = gdspy.FlexPath(meander_points, meander_linewidth)

    pad_and_meander_structure = gdspy.Cell("pad_and_meander")

    pad_and_meander_structure.add(pad_left)
    pad_and_meander_structure.add(pad_right)
    pad_and_meander_structure.add(meander)

    pad_and_meander_polys = pad_and_meander_structure.get_polygons()
    for poly_points in pad_and_meander_polys:
        new_poly_points = mbu.rotate_and_move_points_list(poly_points, rot, x, y)
        poly = gdspy.Polygon(
            new_poly_points,
            layer=DC_structure_material.number,
            datatype=DC_structure_material.datatype,
        )
        mask_builder.Main.add(poly)

    return


def add_test_DC_structure_box_section(
    mask_builder,
    x: float | int,
    y: float | int,
    DC_structure_material: Layer,
    BL_text: str = "",
) -> None:
    """Adds a box of size 8000 centered on the x,y given that contains a series
    of 5 test DC lines and pads. Pads are 100x the size of the linewidth
    connecting them.

    Parameters
    ----------
    x, y: float, int
        The x, y coordinate to center the test DC structure.

    DC_structure_material: Layer
        This is an instance of Layer. see maskpy.layers.Layer.
        Usually this is within the SoukMaskBuilder.layers.xxx.
        e.g. `self.layers.Aluminium`

    KwArgs
    ------
    BL_text = ""
        This is the string that will be written in the bottom left of the
        structure. By default no string will be written.
    """
    if not isinstance(DC_structure_material, Layer):
        raise TypeError(f"DC_structure_material should be of type Layer, current type is {type(DC_structure_material)}")

    box_w_h = 8000

    main_square = gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2])

    mask_builder.ground_plane_cutouts.add(main_square)

    outer_square_side_length = box_w_h
    outer_lw = 100

    outer_square_points = [
        [x - outer_square_side_length / 2 - outer_lw / 2, y + outer_square_side_length / 2 + outer_lw],
        [x - outer_square_side_length / 2 - outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
        [x + outer_square_side_length / 2 + outer_lw / 2, y + outer_square_side_length / 2 + outer_lw / 2],
        [x - outer_square_side_length / 2 - outer_lw, y + outer_square_side_length / 2 + outer_lw / 2],
    ]
    outer_square = gdspy.FlexPath(
        outer_square_points,
        outer_lw,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    mask_builder.Main.add(outer_square)

    pad_and_meander_offset = 2800

    centers_and_rots: list[list[float]] = [
        [0, pad_and_meander_offset, 0],
        [-pad_and_meander_offset, 0, pi / 2],
        [0, -pad_and_meander_offset, 0],
        [pad_and_meander_offset, 0, pi / 2],
        [0, 0, pi / 4],
    ]

    for dx, dy, rot in centers_and_rots:
        mask_builder.make_DC_structure_pads_and_meander(x + dx, y + dy, rot, DC_structure_material)

    if BL_text != "":
        text_size = 300
        BL_text_x_offset = 50
        BL_text_y_offset = 50
        text_xy = [x - box_w_h / 2 + BL_text_x_offset, y - box_w_h / 2 + BL_text_y_offset]

        text = gdspy.Text(
            BL_text,
            text_size,
            position=text_xy,
            layer=DC_structure_material.number,
            datatype=DC_structure_material.datatype,
        )
        mask_builder.Main.add(text)

    return
