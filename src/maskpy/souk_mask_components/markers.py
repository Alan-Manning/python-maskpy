from collections.abc import Sequence

import gdspy
from phidl import geometry as phgeom

from .. import mask_builder_utils as mbu
from ..layers import Layer


def add_MLA_marker(
    mask_builder,
    x: float | int,
    y: float | int,
    materials: Sequence[Layer] | Layer,
    inner_lw: float | int = 5.0,
    outer_lw: float | int = 10.0,
) -> None:
    """Adds a singe or series of square markers with a cross inside centered on
    the x,y coord given. Each subsequent material will make a square larger
    than the previous by 1x the outer linewidth such that they are staggered.
    Inner crosses reamin the same and do not get staggered.

    Parameters
    ----------
    x, y: float | int
        x, y coordinate for the center of the MLA marker.

    materials: Sequence[Layer] | Layer
        Sequence of materials for the marker out of. This should be a list of
        Layer types or a single Layer type, e.g.
        >>> materials = [
        >>>     mask_builder.layers.Aluminium,
        >>>     mask_builder.layers.Nb_Antenna,
        >>>     mask_builder.layers.SiN_dep,
        >>> ]
        >>> #OR
        >>> materials = mask_builder.layers.Aluminium
        >>> #OR
        >>> materials = Layer("Ali", 1, 0)


    KwArgs
    ------
    inner_lw: float, int
        linewidth of the inner cross section of the marker.
        The default value = 5.

    outer_lw: float, int
        linewidth of the outer square section of the marker.
        The default value = 10.
    """

    outer_square_side_length = 1000
    inner_cross_length = 300

    vert_box_points = [[x - inner_lw / 2, y - inner_cross_length / 2], [x + inner_lw / 2, y + inner_cross_length / 2]]

    horz_box_points = [[x - inner_cross_length / 2, y - inner_lw / 2], [x + inner_cross_length / 2, y + inner_lw / 2]]

    # cutout_around_marker = gdspy.Rectangle([-175.0, 300.0], [325.0, -300.0])
    # cutout_around_marker.translate(x, y)
    # self.ground_plane_cutouts.add([cutout_around_marker])
    if not isinstance(materials, list):
        materials = [materials]

    for i, material in enumerate(materials):
        if not isinstance(material, Layer):
            raise TypeError(f"materials should be of type Layer. Found material {material} of type {type(material)}")

        square_offset = i * outer_lw

        square_points = [
            [
                x - outer_square_side_length / 2 - outer_lw / 2 - square_offset,
                y + outer_square_side_length / 2 + outer_lw + square_offset,
            ],
            [
                x - outer_square_side_length / 2 - outer_lw / 2 - square_offset,
                y - outer_square_side_length / 2 - outer_lw / 2 - square_offset,
            ],
            [
                x + outer_square_side_length / 2 + outer_lw / 2 + square_offset,
                y - outer_square_side_length / 2 - outer_lw / 2 - square_offset,
            ],
            [
                x + outer_square_side_length / 2 + outer_lw / 2 + square_offset,
                y + outer_square_side_length / 2 + outer_lw / 2 + square_offset,
            ],
            [
                x - outer_square_side_length / 2 - outer_lw - square_offset,
                y + outer_square_side_length / 2 + outer_lw / 2 + square_offset,
            ],
        ]

        square_path = gdspy.FlexPath(
            square_points,
            outer_lw,
            layer=material.number,
            datatype=material.datatype,
        )
        square_path_polys = mbu.get_polys_from_flexpath(square_path)
        for i in range(len(square_path_polys)):
            mask_builder.Main.add(
                gdspy.Polygon(
                    square_path_polys[i],
                    layer=material.number,
                    datatype=material.datatype,
                )
            )

        vert_box = gdspy.Rectangle(
            vert_box_points[0],
            vert_box_points[1],
            layer=material.number,
            datatype=material.datatype,
        )
        horz_box = gdspy.Rectangle(
            horz_box_points[0],
            horz_box_points[1],
            layer=material.number,
            datatype=material.datatype,
        )

        mask_builder.Main.add(vert_box)
        mask_builder.Main.add(horz_box)

    cutout_around_marker = gdspy.Rectangle(
        [
            x - outer_square_side_length / 2 - outer_lw / 2 - square_offset - outer_lw,
            y + outer_square_side_length / 2 + outer_lw + square_offset + outer_lw,
        ],
        [
            x + outer_square_side_length / 2 + outer_lw / 2 + square_offset + outer_lw,
            y - outer_square_side_length / 2 - outer_lw / 2 - square_offset - outer_lw,
        ],
    )
    mask_builder.ground_plane_cutouts.add([cutout_around_marker])
    mask_builder.silicon_nitride_cutouts.add([cutout_around_marker])
    mask_builder.silicon_nitride_membrane_cutouts.add([cutout_around_marker])
    mask_builder.silicon_oxide_cutouts.add([cutout_around_marker])

    return


def add_initial_alignment_markers(
    mask_builder,
    x: float | int,
    y: float | int,
    cross_length: float | int = 300,
    linewidth: float | int = 5,
    cutout_square_size: float | int = 1000,
) -> None:
    """Adds a singe MLA marker cross centered on the x,y coord given. Also adds
    a cutout window in layers that will be deposited ontop of this (namely
    Nb_groundplane).

    Parameters
    ----------
    x, y: float, int
        x, y coordinate for the center of the MLA marker.

    KwArgs
    ------
    cross_length: float, int
        Height of the marker cross in the center.
        Default = 300.

    linewidth: float, int
        Linewidth of the inner cross.
        Default = 5.

    cutout_square_size: float, int
        The width and heigh of the outer square cutout window.
        Default = 1000.
    """

    vert_box = gdspy.Rectangle(
        [x - (linewidth / 2), y - (cross_length / 2)], [x + (linewidth / 2), y + (cross_length / 2)], layer=0, datatype=0
    )

    horz_box = gdspy.Rectangle(
        [x - (cross_length / 2), y - (linewidth / 2)], [x + (cross_length / 2), y + (linewidth / 2)], layer=0, datatype=0
    )
    mask_builder.Main.add(vert_box)
    mask_builder.Main.add(horz_box)

    cutout_around_marker_rect = gdspy.Rectangle(
        [x - (cutout_square_size / 2), y - (cutout_square_size / 2)], [x + (cutout_square_size / 2), y + (cutout_square_size / 2)]
    )
    mask_builder.ground_plane_cutouts.add(cutout_around_marker_rect)

    return


def add_caliper_alignment_markers(
    mask_builder,
    x: float | int,
    y: float | int,
    rot: float | int,
    layer1: Layer,
    layer2: Layer,
    layer1_text: str,
    layer2_text: str,
) -> None:
    """Adds Alignment Markers to the mask in a cutout centered on the x,y
    given. This consists of 2 main corsses and a series of calipers made from
    the two materials given.

    Parameters
    ----------
    x, y: float, int
        The x,y coordinate to place the marker.

    rot: float, int
        The angle (**in degrees**) the marker should be rotated.

    layer1, layer2: Layer
        These are instances of Layer. see maskpy.layers.Layer.
        Usually this is within the SoukMaskBuilder.layers.xxx.
        e.g. `self.layers.Aluminium`

    layer1_text, layer2_text: str
        These should be strings that are the text to be placed above and
        below in the marker.
    """

    # adding the Text
    bot_text = gdspy.Text(
        layer1_text,
        100,
        (x - 100, y - 250),
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    top_text = gdspy.Text(
        layer2_text,
        100,
        (x - 100, y + 150),
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    mask_builder.Main.add(bot_text)
    mask_builder.Main.add(top_text)

    # adding the cutout around the aligment markers
    cutout_around_marker = gdspy.Rectangle([-175.0, 300.0], [325.0, -300.0])
    cutout_around_marker.translate(x, y)
    mask_builder.ground_plane_cutouts.add(cutout_around_marker)

    # adding the main crosses
    main_cross1_poly_points = [
        (-15.0, -100.0),
        (-15.0, -15.0),
        (-100.0, -15.0),
        (-100.0, 15.0),
        (-15.0, 15.0),
        (-15.0, 100.0),
        (15.0, 100.0),
        (15.0, 15.0),
        (100.0, 15.0),
        (100.0, -15.0),
        (15.0, -15.0),
        (15.0, -100.0),
    ]
    main_cross1_poly = gdspy.Polygon(
        main_cross1_poly_points,
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    main_cross1_poly.translate(x, y)
    mask_builder.Main.add(main_cross1_poly)

    main_cross2_poly_points = [
        (-11.0, -96.0),
        (-11.0, -11.0),
        (-96.0, -11.0),
        (-96.0, 11.0),
        (-11.0, 11.0),
        (-11.0, 96.0),
        (11.0, 96.0),
        (11.0, 11.0),
        (96.0, 11.0),
        (96.0, -11.0),
        (11.0, -11.0),
        (11.0, -96.0),
    ]
    main_cross2_poly = gdspy.Polygon(
        main_cross2_poly_points,
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    main_cross2_poly.translate(x, y)
    mask_builder.Main.add(main_cross2_poly)

    # adding the cross and squares in the top left
    TL_sq = gdspy.Rectangle(
        [x - 100, y + 100],
        [x - 75, y + 75],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    BL_sq = gdspy.Rectangle(
        [x - 100, y + 55],
        [x - 75, y + 30],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    TR_sq = gdspy.Rectangle(
        [x - 55, y + 100],
        [x - 30, y + 75],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    BR_sq = gdspy.Rectangle(
        [x - 55, y + 55],
        [x - 30, y + 30],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    mask_builder.Main.add(TL_sq)
    mask_builder.Main.add(BL_sq)
    mask_builder.Main.add(TR_sq)
    mask_builder.Main.add(BR_sq)

    top_cross_poly_points = [
        (-73, 100),
        (-57, 100),
        (-57, 73),
        (-30, 73),
        (-30, 57),
        (-57, 57),
        (-57, 30),
        (-73, 30),
        (-73, 57),
        (-100, 57),
        (-100, 73),
        (-73, 73),
    ]
    top_cross_poly = gdspy.Polygon(
        top_cross_poly_points,
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    top_cross_poly.translate(x, y)
    mask_builder.Main.add(top_cross_poly)

    # adding the bottom left rectangle markers
    bot_Rect1 = gdspy.Rectangle(
        [x - 100, y - 100],
        [x - 90, y - 50],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    bot_Rect2 = gdspy.Rectangle(
        [x - 80, y - 100],
        [x - 70, y - 50],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    bot_Rect3 = gdspy.Rectangle(
        [x - 60, y - 100],
        [x - 50, y - 50],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    bot_Rect4 = gdspy.Rectangle(
        [x - 40, y - 100],
        [x - 30, y - 50],
        layer=layer1.number,
        datatype=layer1.datatype,
    )

    bot_Rect5 = gdspy.Rectangle(
        [x - 100, y - 70],
        [x - 90, y - 30],
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    bot_Rect6 = gdspy.Rectangle(
        [x - 80, y - 70],
        [x - 70, y - 30],
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    bot_Rect7 = gdspy.Rectangle(
        [x - 60, y - 70],
        [x - 50, y - 30],
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    bot_Rect8 = gdspy.Rectangle(
        [x - 40, y - 70],
        [x - 30, y - 30],
        layer=layer2.number,
        datatype=layer2.datatype,
    )

    mask_builder.Main.add(bot_Rect1)
    mask_builder.Main.add(bot_Rect2)
    mask_builder.Main.add(bot_Rect3)
    mask_builder.Main.add(bot_Rect4)
    mask_builder.Main.add(bot_Rect5)
    mask_builder.Main.add(bot_Rect6)
    mask_builder.Main.add(bot_Rect7)
    mask_builder.Main.add(bot_Rect8)

    # adding the bottom right rectangle markers
    left_Rect1 = gdspy.Rectangle(
        [x + 100, y - 100],
        [x + 50, y - 90],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    left_Rect2 = gdspy.Rectangle(
        [x + 100, y - 80],
        [x + 50, y - 70],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    left_Rect3 = gdspy.Rectangle(
        [x + 100, y - 60],
        [x + 50, y - 50],
        layer=layer1.number,
        datatype=layer1.datatype,
    )
    left_Rect4 = gdspy.Rectangle(
        [x + 100, y - 40],
        [x + 50, y - 30],
        layer=layer1.number,
        datatype=layer1.datatype,
    )

    left_Rect5 = gdspy.Rectangle(
        [x + 70, y - 100],
        [x + 30, y - 90],
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    left_Rect6 = gdspy.Rectangle(
        [x + 70, y - 80],
        [x + 30, y - 70],
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    left_Rect7 = gdspy.Rectangle(
        [x + 70, y - 60],
        [x + 30, y - 50],
        layer=layer2.number,
        datatype=layer2.datatype,
    )
    left_Rect8 = gdspy.Rectangle(
        [x + 70, y - 40],
        [x + 30, y - 30],
        layer=layer2.number,
        datatype=layer2.datatype,
    )

    mask_builder.Main.add(left_Rect1)
    mask_builder.Main.add(left_Rect2)
    mask_builder.Main.add(left_Rect3)
    mask_builder.Main.add(left_Rect4)
    mask_builder.Main.add(left_Rect5)
    mask_builder.Main.add(left_Rect6)
    mask_builder.Main.add(left_Rect7)
    mask_builder.Main.add(left_Rect8)

    # adding the callipers
    caliper = phgeom.litho_calipers(
        notch_size=[6, 40],
        notch_spacing=6,
        num_notches=5,
        offset_per_notch=0.5,
        row_spacing=-20.0,
        layer1=layer1.number,
        layer2=layer2.number,
    )
    caliper.move([x + 60, y + 60])

    caliper_polys = caliper.get_polygons()
    for i, poly_points in enumerate(caliper_polys):
        if i % 2 == 0:
            mask_builder.Main.add(
                gdspy.Polygon(
                    poly_points,
                    layer=layer1.number,
                    datatype=layer1.datatype,
                )
            )
        else:
            mask_builder.Main.add(
                gdspy.Polygon(
                    poly_points,
                    layer=layer2.number,
                    datatype=layer2.datatype,
                )
            )

    rotcaliper = phgeom.litho_calipers(
        notch_size=[6, 40],
        notch_spacing=6,
        num_notches=5,
        offset_per_notch=0.5,
        row_spacing=-20.0,
        layer1=layer1.number,
        layer2=layer2.number,
    )
    rotcaliper.move([x + 170, y + 25])
    rotcaliper.rotate(-90.0, (x + 170, y + 25))

    rotcaliper_polys = rotcaliper.get_polygons()
    for i, poly_points in enumerate(rotcaliper_polys):
        if i % 2 == 0:
            mask_builder.Main.add(
                gdspy.Polygon(
                    poly_points,
                    layer=layer1.number,
                    datatype=layer1.datatype,
                )
            )
        else:
            mask_builder.Main.add(
                gdspy.Polygon(
                    poly_points,
                    layer=layer2.number,
                    datatype=layer2.datatype,
                )
            )

    return
