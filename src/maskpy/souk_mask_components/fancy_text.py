from collections.abc import Sequence
from typing import Literal

import gdspy
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

from ..layers import Layer
from ..logging import styled_type_error


def add_fancy_text(
    mask_builder,
    text: str,
    x: float | int,
    y: float | int,
    size: float,
    layer: Layer | Sequence[Layer],
    rotation: float = 0.0,
    horizontal_align: Literal["start", "center", "end"] = "start",
    vertical_align: Literal["baseline", "above", "center", "below"] = "baseline",
    font_properties: FontProperties | None = None,
    usetex: bool = False,
    return_polyset: bool = False,
    bb_cutout_in_grnd: bool = False,
    bb_cutout_in_sin_dep: bool = False,
) -> None | list[np.ndarray]:
    """Add fancy text to the mask. This is text with a custon font face,
    custom placement or rotation.

    Parameters
    ----------
    text: str
        The text to be added. This can include tex but will only be rendered
        if the usetex argument is set to True.

    x: float
        The x position for the text to be drawn at.

    y: float
        The y position for the text to be drawn at.

    size: float
        The size of the text. This is the extent from the lowest descender
        to the tallest ascender, so actual size of characters may differ
        depending on the font face or stylings.

    layer: Layer | Sequence[Layer]
        This is an instance of Layer of a Sequence of many layers. The text
        will be drawn on every layer given. see maskpy.layers.Layer.
        Usually this is within the SoukMaskBuilder.layers.xxx.
        e.g. `self.layers.Aluminium`

    KwArgs
    ------
    rotation: float = 0.0
        The angle (**in radians**) the text should be rotated. positive
        angles are rotations going counterclockwise. This rotaion is
        applied about the x,y coordinate given and is done after any
        translation of the text (i.e. different horizontal_align and
        vertical_align than start and baseline respectively).

    horizontal_align: Literal["start", "center", "end"] = "start",
        The horizontal aligment of the text. This takes values "start",
        "center", "end". Any other string will be interpreted as "start".
        "start" is such that the leftmost end of the text aligns inline
        with the x coordinate given. "end" aligns the rightmost end of the
        text inline with the x coordinate given. "center" aligns the middle
        of the leftmost and rightmost ends inline with the x coordinate
        given.

    vertical_align: str = "baseline",
        The vertical aligment of the text. This takes values "baseline",
        "above", "center", "below". Any other string will be interpreted as
        "baseline". "baseline" will align the bottom of characters without
        any descender inline with the y coordinate given. "above" will
        align the bottom of the lowest descender inline with the y
        coordinate given. "below" will align the top of the tallest
        ascender inline with the y coordinate given. "center" will align
        the middle of the lowest descender and the tallest ascender inline
        with the y coordinate given.

    font_properties: FontProperties | None = None,
        This is any options for the format of the text, family, style etc.
        See matplotlib.font_manager.FontProperties for all the options.
        If None provided, this will use a default FontProperties,
        font_properties=FontProperties(family="monospace", style="normal").

    usetex: bool = False,
        Whether to use tex rendering on the input text, defaults to False.

    return_polyset: bool = False,
        Whether to disable the drawing and return the polygon points used
        to create the text. Defaults to False. When True this not draw
        anything to the mask. Instead the polygon points that would have
        been used to draw that text will be returned, this is a list of
        NDArrays containing x,y coordinates.

    bb_cutout_in_grnd: bool = False
        Wether to make a cutout for the bounding box around the text in the
        NbGroundplane layer. Default if False.

    bb_cutout_in_sin_dep: bool = False
        Wether to make a cutout for the bounding box around the text in the
        SiN_dep layer. Default if False.
    """

    if font_properties is not None:
        if isinstance(font_properties, FontProperties):
            font_props = font_properties
        else:
            raise TypeError("font_properties argument should be of type FontProperties")
    else:
        font_props = FontProperties(family="monospace", style="normal")

    text_layers: list[Layer] = []

    if isinstance(layer, Layer):
        text_layers.append(layer)
    elif isinstance(layer, Sequence):
        if all(isinstance(lay, Layer) for lay in layer):
            for lay in layer:
                text_layers.append(lay)
        else:
            styled_type_error(layer, "layer", Sequence[Layer])
    else:
        styled_type_error(
            layer,
            "layer",
            Layer | Sequence[Layer],
        )

    full_text_path = TextPath((x, y), text, size=size, prop=font_props, usetex=usetex)
    list_of_text_polygon_sets: list[gdspy.PolygonSet] = []

    x_mins: list[float] = []
    y_mins: list[float] = []
    x_maxs: list[float] = []
    y_maxs: list[float] = []
    for points, path_code in full_text_path.iter_segments():
        match path_code:
            case full_text_path.MOVETO:
                curve = gdspy.Curve(*points)
            case full_text_path.LINETO:
                curve.L(*points)
            case full_text_path.CURVE3:
                curve.Q(*points)
            case full_text_path.CURVE4:
                curve.C(*points)
            case full_text_path.CLOSEPOLY:
                current_part_of_letter_polygon_points = curve.get_points()
                for text_polygon_points in list_of_text_polygon_sets:
                    # if the first point of the current letter part is in
                    # inside the previous letter part or the first point of
                    # the previous letter part is inside the current letter
                    # part then there is an ovelap and an xor needs to be
                    # completed on these letter parts to cutout the shape.
                    # i.e. cutout holes in a letter,for example 8, B, D, e.
                    first_point_of_current_letter_part = current_part_of_letter_polygon_points[:1]
                    first_point_of_previous_letter_part = text_polygon_points[:1]
                    all_of_current_inside_previous = gdspy.inside(
                        first_point_of_current_letter_part,
                        [text_polygon_points],
                    )[0]
                    all_of_previous_inside_current = gdspy.inside(
                        first_point_of_previous_letter_part,
                        [current_part_of_letter_polygon_points],
                    )[0]

                    if all_of_current_inside_previous or all_of_previous_inside_current:
                        previous_part_of_letter_polygon_points = list_of_text_polygon_sets.pop(-1)
                        current_part_of_letter_polygon_points = gdspy.boolean(
                            [previous_part_of_letter_polygon_points],
                            [current_part_of_letter_polygon_points],
                            "xor",
                            max_points=0,
                        ).polygons[0]

                list_of_text_polygon_sets.append(current_part_of_letter_polygon_points)

                x_mins.append(min([p[0] for p in current_part_of_letter_polygon_points]))
                x_maxs.append(max([p[0] for p in current_part_of_letter_polygon_points]))
                y_mins.append(min([p[1] for p in current_part_of_letter_polygon_points]))
                y_maxs.append(max([p[1] for p in current_part_of_letter_polygon_points]))

    # horizontal aligment
    match horizontal_align:
        case "center":
            delta_x = x - min(x_mins) - ((max(x_maxs) - min(x_mins)) / 2)
        case "end":
            delta_x = x - max(x_maxs)
        case _:  # "start". or anything else is start
            delta_x = x - min(x_mins)

    match vertical_align:
        case "above":
            delta_y = y - min(y_mins)
        case "center":
            delta_y = y - min(y_mins) - ((max(y_maxs) - min(y_mins)) / 2)
        case "below":
            delta_y = y - max(y_maxs)
        case _:  # "baseline". or anything else is baseline
            delta_y = 0.0

    text_polygons: list[gdspy.PolygonSet] = []

    for text_layer in text_layers:
        text_polygon = gdspy.PolygonSet(list_of_text_polygon_sets, layer=text_layer.number, datatype=text_layer.datatype)
        text_polygon.translate(delta_x, delta_y)
        if rotation != 0.0:
            text_polygon.rotate(rotation, (x, y))
        text_polygons.append(text_polygon)

    if bb_cutout_in_grnd:
        bb_cutout_rect = gdspy.Rectangle(
            (min(x_mins), min(y_mins)),
            (max(x_maxs), max(y_maxs)),
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        bb_cutout_rect.translate(delta_x, delta_y)
        if rotation != 0.0:
            bb_cutout_rect.rotate(rotation, (x, y))

        mask_builder.ground_plane_cutouts.add(bb_cutout_rect)

    if bb_cutout_in_sin_dep:
        bb_cutout_rect = gdspy.Rectangle(
            (min(x_mins), min(y_mins)),
            (max(x_maxs), max(y_maxs)),
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        bb_cutout_rect.translate(delta_x, delta_y)
        if rotation != 0.0:
            bb_cutout_rect.rotate(rotation, (x, y))

        mask_builder.silicon_nitride_cutouts.add(bb_cutout_rect)

    if return_polyset:
        return text_polygons[0].polygons

    for text_polygon in text_polygons:
        mask_builder.Main.add(text_polygon)

    return None
