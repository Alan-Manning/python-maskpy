from typing import Any, Literal, Sequence

import gdspy
import numpy as np
from numpy import cos, pi, sin

from .. import mask_builder_utils as mbu
from ..layers import Layer
from ..logging import TextColor, styled_text, styled_type_error


def _get_path_points_by_corner_type(
    cpw_points: Sequence[Sequence[float | int]],
    path_width: float | int | None,
    corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"],
    bend_radius: float | int | None,
    path_type: Literal["cutout", "dielectric", ""] = "",
) -> Sequence[Sequence[float | int]]:

    if (path_type == "cutout") or (path_type == "dielectric"):
        path_type_string = path_type + "_"
    else:
        path_type_string = ""

    MIN_NUMBER_OF_POINTS_IN_BEND: int = int(50)

    match corners_type:
        case "natural":
            path_points = cpw_points
        case "miter":
            if bend_radius is None:
                raise ValueError(f'`{path_type_string}bend_radius` should be defined for `{path_type_string}corners_type` = "miter"')

            path_points = mbu.get_rounded_path_points(
                cpw_points,
                bend_radius,
                no_points_in_bend=2,
            )
        case "bevel":
            raise NotImplementedError()
        case "round" | "circular bend":
            if bend_radius is None:
                raise ValueError(
                    f'`{path_type_string}bend_radius` should be defined for `{path_type_string}corners_type` = "circular bend" | "round"'
                )

            if path_width is not None:
                no_points_in_bend = max(int(path_width / bend_radius), MIN_NUMBER_OF_POINTS_IN_BEND)
            else:
                no_points_in_bend = MIN_NUMBER_OF_POINTS_IN_BEND

            path_points = mbu.get_rounded_path_points(
                cpw_points,
                bend_radius,
                no_points_in_bend=no_points_in_bend,
            )
        case "smooth":
            raise NotImplementedError()
        case _:
            styled_type_error(corners_type, "corners_type", Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"])
    return path_points


def _add_bridges(
    mask_builder,
    bridge_layer: Layer,
    path_points: Sequence[Sequence[float | int]],
    bridge_gap: float | int | Sequence[float | int],
    bridge_width: float | int | Sequence[float | int],
    bridge_height: float | int | Sequence[float | int],
    path_type: Literal["cutout", "dielectric", ""] = "",
):
    """."""
    if path_type == "cutout" or path_type == "dielectric":
        path_type_string = path_type + "_"
    else:
        path_type_string = ""

    test_path_points_path = gdspy.FlexPath(
        path_points,
        0.2,
        layer=mask_builder.layers.General_labeling.number,
        datatype=mask_builder.layers.General_labeling.datatype,
    )
    mask_builder.make_flexpath_into_polygons_and_add_to_main(test_path_points_path, mask_builder.layers.General_labeling)

    dls: list[float] = []
    for i in range(len(path_points) - 1):
        current_xy = path_points[i]
        next_xy = path_points[i + 1]
        dl = np.sqrt((current_xy[0] - next_xy[0]) ** 2 + (current_xy[1] - next_xy[1]) ** 2)
        dls.append(dl)

    total_length = np.sum(dls)

    if isinstance(bridge_gap, float | int):
        lengths_to_add_bridges = np.arange(
            bridge_gap,
            total_length,
            bridge_gap,
        )
    elif isinstance(bridge_gap, Sequence):
        lengths_to_add_bridges = np.array(bridge_gap)
    else:
        styled_type_error(
            bridge_gap,
            "bridge_gap",
            float | int | Sequence[float | int],
        )

    if isinstance(bridge_width, float | int):
        bridge_widths = np.ones_like(lengths_to_add_bridges) * bridge_width
    elif isinstance(bridge_width, Sequence):
        if isinstance(bridge_gap, Sequence):
            if len(bridge_width) != len(bridge_gap):
                raise ValueError(
                    f"Sequence of `{path_type_string}bridge_width` should be the same length as Sequence of `{path_type_string}bridge_gap`."
                )
        bridge_widths = np.array(bridge_width)
    else:
        styled_type_error(
            bridge_width,
            "bridge_width",
            float | int | Sequence[float | int],
        )

    if isinstance(bridge_height, float | int):
        bridge_heights = np.ones_like(lengths_to_add_bridges) * bridge_height
    elif isinstance(bridge_height, Sequence):
        if isinstance(bridge_gap, Sequence):
            if len(bridge_height) != len(bridge_gap):
                raise ValueError(
                    f"Sequence of `{path_type_string}bridge_height` should be the same length as Sequence of `{path_type_string}bridge_gap`."
                )
        bridge_heights = np.array(bridge_height)
    else:
        styled_type_error(
            bridge_height,
            "bridge_height",
            float | int | Sequence[float | int],
        )

    for length_to_add_bridge, bridge_w, bridge_h in zip(lengths_to_add_bridges, bridge_widths, bridge_heights):
        cumulative_dl = 0.0
        for index, dl in enumerate(dls):
            cumulative_dl += dl
            if length_to_add_bridge >= cumulative_dl:
                continue

            percentage_into_section = (length_to_add_bridge - (cumulative_dl - dl)) / dl

            prev_point = path_points[index]
            next_point = path_points[index + 1]

            dx = next_point[0] - prev_point[0]
            dy = next_point[1] - prev_point[1]

            bridge_center_x = prev_point[0] + (dx * percentage_into_section)
            bridge_center_y = prev_point[1] + (dy * percentage_into_section)

            angle_of_last_section = np.arctan2((next_point[1] - prev_point[1]), (next_point[0] - prev_point[0]))
            angle_for_bridge = angle_of_last_section

            bridge_rect = gdspy.Rectangle(
                [bridge_center_x - (bridge_w / 2), bridge_center_y - (bridge_h / 2)],
                [bridge_center_x + (bridge_w / 2), bridge_center_y + (bridge_h / 2)],
                layer=bridge_layer.number,
                datatype=bridge_layer.datatype,
            )
            bridge_rect.rotate(angle_for_bridge, [bridge_center_x, bridge_center_y])
            mask_builder.Main.add(bridge_rect)
            break
        else:  # NOTE This is a "for else" -> i.e. run if no break
            raise RuntimeError(f"length to add {path_type_string}bridge is greater than the total length of the path.")


def _get_cutout_path(
    path_points: Sequence[Sequence[float | int]],
    start_end_offset: Sequence[float | int | None],
) -> Sequence[Sequence[float | int]]:
    """."""

    start_offset = start_end_offset[0]
    end_offset = start_end_offset[1]

    if (start_offset is None) and (end_offset is None):
        return path_points

    dls: list[float] = []
    for i in range(len(path_points) - 1):
        current_xy = path_points[i]
        next_xy = path_points[i + 1]
        dl = np.sqrt((current_xy[0] - next_xy[0]) ** 2 + (current_xy[1] - next_xy[1]) ** 2)
        dls.append(dl)

    if start_offset is None:
        start_of_rounded_path_index = 0
        new_start_point = None
    else:
        cumulative_dl = 0.0
        for index, dl in enumerate(dls):
            cumulative_dl += dl
            if start_offset >= cumulative_dl:
                continue

            start_of_rounded_path_index = index + 1
            percentage_into_section = (start_offset - (cumulative_dl - dl)) / dl

            prev_point = path_points[index]
            next_point = path_points[index + 1]

            dx = next_point[0] - prev_point[0]
            dy = next_point[1] - prev_point[1]

            new_start_x = prev_point[0] + (dx * percentage_into_section)
            new_start_y = prev_point[1] + (dy * percentage_into_section)
            new_start_point = [new_start_x, new_start_y]
            break
        else:  # NOTE This is a "for else" -> i.e. run if no break
            raise RuntimeError(f"start offset is greater than the total length of the path.")

    if end_offset is None:
        end_of_rounded_path_index = -1
        new_end_point = None
    else:
        if end_offset < 0:
            # if less than zero then append to end of cutout path rather then cut into path.
            end_of_rounded_path_index = -1
            last_point = path_points[-1]
            penultimate_point = path_points[-2]
            angle_of_last_section = np.arctan2((last_point[1] - penultimate_point[1]), (last_point[0] - penultimate_point[0]))
            new_end_x = last_point[0] - (end_offset * cos(angle_of_last_section))
            new_end_y = last_point[1] - (end_offset * sin(angle_of_last_section))
            new_end_point = [new_end_x, new_end_y]

        else:
            total_length = np.sum(dls)
            cumulative_dl = 0.0
            for index, dl in enumerate(dls):
                cumulative_dl += dl
                if (total_length - end_offset) >= cumulative_dl:
                    continue

                end_of_rounded_path_index = index - 1
                percentage_into_section = ((total_length - end_offset) - (cumulative_dl - dl)) / dl

                prev_point = path_points[index]
                next_point = path_points[index + 1]

                dx = next_point[0] - prev_point[0]
                dy = next_point[1] - prev_point[1]

                new_end_x = prev_point[0] + (dx * percentage_into_section)
                new_end_y = prev_point[1] + (dy * percentage_into_section)
                new_end_point = [new_end_x, new_end_y]
                break
            else:  # NOTE This is a "for else" -> i.e. run if no break
                raise RuntimeError(f"end offset is greater than the total length of the path.")

    new_rounded_path_points = path_points[start_of_rounded_path_index:end_of_rounded_path_index]

    if new_start_point is not None:
        new_rounded_path_points.insert(0, new_start_point)
    if new_end_point is not None:
        new_rounded_path_points.append(new_end_point)

    return new_rounded_path_points


def add_cpw(
    # mask_builder: SoukMaskBuilder,
    mask_builder,
    cpw_points: Sequence[Sequence[float | int]],
    center_width: float | int,
    cutout_width: float | int,
    #
    center_material: Layer | Sequence[Layer] | None = None,
    dielectric_in_cutout_or_positive: Literal["positive", "cutout"] = "positive",
    dielectric_width: float | int | None = None,
    #
    center_corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] = "natural",
    center_bend_radius: float | int | None = None,
    cutout_corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] = "natural",
    cutout_bend_radius: float | int | None = None,
    dielectric_corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] = "natural",
    dielectric_bend_radius: float | int | None = None,
    #
    cutout_start_end_offset: Sequence[float | int | None] | None = None,
    dielectric_start_end_offset: Sequence[float | int | None] | None = None,
    #
    cutout_bridge_gap: float | int | Sequence[float | int] | None = None,
    cutout_bridge_width: float | int | Sequence[float | int] | None = None,
    cutout_bridge_height: float | int | Sequence[float | int] | None = None,
    #
    dielectric_bridge_gap: float | int | Sequence[float | int] | None = None,
    dielectric_bridge_width: float | int | Sequence[float | int] | None = None,
    dielectric_bridge_height: float | int | Sequence[float | int] | None = None,
) -> None:
    """Add a coplanar waveguide with a center line in defined material and a
    cutout in the groundplane layer. Can also add a dielectric with the cpw if
    its width is defined.

    Parameters
    ----------
    mask_builder: SoukMaskBuilder
       The mask builder to add the CPW to.

    cpw_points: Sequence[Sequence[float | int]]
        list of [x,y] lists which are the coordinates for the cpw to pass
        through.

    center_width: float | int
        The width of the center line of the CPW.

    cutout_width: float | int,
        The width of the ground plane cutout for the CPW centered around the
        center line.

    bend_radius: float, int
        The radius for the bends in the CPW to be made.

    KwArgs
    ------
    center_material: Layer | Sequence[Layer] | None = None,
        Default None. This is the Layer that the center line of the CPW should
        be drawn on. This is an instance of Layer or a Sequence of Layer's.
        When this isnt specified the center line will be drawn on the
        mask_builder.layers.Nb_Antenna Layer. When a Sequence of Layer's is
        given the center line will be drawn in all those layers. Usually this
        is within the SoukMaskBuilder.layers.xxx. e.g.
        `mask_builder.layers.Aluminium`. For more see maskpy.layers.Layer.

    dielectric_in_cutout_or_positive: Literal["positive", "cutout"] = "positive",
        This will by default draw the dielectric path to be positive.
        Alternatively is this is set to "cutout", then the dielectric path
        will be added as a cutout from the surrounding dielectric on the mask.

    dielectric_width: float | int | None = None,
        Default None. This is the width of the dielectric. When this isnt
        specified the dielectric will not be drawn else the width will take
        the value given.

    center_corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] = "natural"
        This is the type of corners for the center line of the CPW. This
        follows the coices provided by gdspy.FlexPath baring the custom user
        defined Callable.

    center_bend_radius: float | int | None = None
        The bend radius for the center line corners if they require it.

    cutout_corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] = "natural"
        This is the type of corners for the groundplane cutout of the CPW. This
        follows the coices provided by gdspy.FlexPath baring the custom user
        defined Callable.

    cutout_bend_radius: float | int | None = None
        The bend radius for the groundplane cutout corners if they require it.

    dielectric_corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] = "natural"
        This is the type of corners for the dielectric of the CPW. This
        follows the coices provided by gdspy.FlexPath baring the custom user
        defined Callable.

    dielectric_bend_radius: float | int | None = None
        The bend radius for the dielectric corners if they require it.

    cutout_start_end_offset: Sequence[float | int | None] | None = None
        Default None. When defined this will cause the groundplane cutout to be
        offset from the start and end respecively. For only offseting the start
        or end None can be passed in the Sequence. The values taken can be
        positive or negative to either shorted or lengthen the groundplane
        cutout respecively.

    dielectric_start_end_offset: Sequence[float | int | None] | None = None
        Default None. When defined this will cause the dielectric to be
        offset from the start and end respecively. For only offseting the start
        or end None can be passed in the Sequence. The values taken can be
        positive or negative to either shorted or lengthen the dielectric
        respecively.

    cutout_bridge_gap: float | int | Sequence[float | int] | None = None
        Default None will not add any bridges. Bridging over the groundplane
        cutout will only occur when this and `cutout_bridge_width` are defined.
        When this is defined it can either be a number which will space bridges
        evenly at that distance appart (center-to-center). If a Sequence is
        given then bridges will be placed centered at those lengths into the
        CPW path.

    cutout_bridge_width: float | int | Sequence[float | int] | None = None
        Default None will not add any bridges. Bridging over the groundplane
        cutout will only occur when this and `cutout_bridge_gap` are defined.
        When this is defined it can either be a number which will make all the
        bridges the given width. If a Sequence is given then bridges will each
        take thier width from the values in the Sequence, i.e. the first bridge
        will have width = cutout_bridge_width[0], the second will be [1] and so
        on. This limits the amount of bridges to the length of the Sequence
        given. *Note* If `cutout_bridge_gap` is also a Sequence then length of
        these Sequence's should match.

    cutout_bridge_height: float | int | Sequence[float | int] | None = None
        Default None will make the bridge's height (if they are to be drawn)
        equal to the `cutout_width` plus a little extra to ensure correct
        bridging on corners. When this is defined it can either be a number
        which will make all the bridges the given height. If a Sequence is
        given then bridges will each take thier height from the values in the
        Sequence, i.e. the first bridge's height = cutout_bridge_height[0],
        the second will be [1] and so on. This limits the amount of bridges to
        the length of the Sequence given. *Note* If `cutout_bridge_gap` is also
        a Sequence then length of these Sequence's should match.

    dielectric_bridge_gap: float | int | Sequence[float | int] | None = None
        Default None will not add any bridges. Dielectric bridging will only
        occur when this and `dielectric_bridge_width` are defined. When this is
        defined it can either be a number which will space bridges
        evenly at that distance appart (center-to-center). If a Sequence is
        given then bridges will be placed centered at those lengths into the
        CPW path.

    dielectric_bridge_width: float | int | Sequence[float | int] | None = None
        Default None will not add any bridges. Dielectric bridging will only
        occur when this and `dielectric_bridge_gap` are defined. When this is
        defined it can either be a number which will make all the bridges the
        given width. If a Sequence is given then bridges will each take thier
        width from the values in the Sequence, i.e. the first bridge will have
        width = dielectric_bridge_width[0], the second will be [1] and so
        on. This limits the amount of bridges to the length of the Sequence
        given. *Note* If `dielectric_bridge_gap` is also a Sequence then length
        of these Sequence's should match.

    dielectric_bridge_height: float | int | Sequence[float | int] | None = None
        Default None will make the bridge's height (if they are to be drawn)
        equal to the `dielectric_width` if defined else the `cutout_width` plus
        a little extra to ensure correct bridging on corners. When this is
        defined it can either be a number which will make all the bridges the
        given height. If a Sequence is given then bridges will each take thier
        height from the values in the Sequence, i.e. the first bridge's
        height = dielectric_bridge_height[0], the second will be [1] and so on.
        This limits the amount of bridges to the length of the Sequence given.
        *Note* If `dielectric_bridge_gap` is also a Sequence then length of
        these Sequence's should match.
    """

    # Getting center layers list.
    center_layers: list[Layer] = []

    if center_material is None:
        center_layers.append(mask_builder.layers.Nb_Antenna)
    elif isinstance(center_material, Layer):
        center_layers.append(center_material)
    elif isinstance(center_material, Sequence):
        if all(isinstance(mat, Layer) for mat in center_material):
            for material in center_material:
                center_layers.append(material)
        else:
            raise TypeError(
                "`center_material` elements should be of type Layer.\n"
                + styled_text(
                    f"current type is {type(center_material)}",
                    color=TextColor.ERROR,
                )
            )

    else:
        styled_type_error(
            center_material,
            "center_material",
            Layer | list[Layer] | None,
        )

    ###########################################################################
    # Center path
    ###########################################################################
    for layer in center_layers:
        center_path = gdspy.FlexPath(
            cpw_points,
            center_width,
            corners=center_corners_type,
            bend_radius=center_bend_radius,
            gdsii_path=True,
            layer=layer.number,
            datatype=layer.datatype,
        )
        mask_builder.make_flexpath_into_polygons_and_add_to_main(
            center_path,
            layer,
        )

    ###########################################################################
    # Groundplane cutout path
    ###########################################################################
    if cutout_start_end_offset is None:
        cutout_path_points = cpw_points
    else:
        path_points = _get_path_points_by_corner_type(
            cpw_points,
            cutout_width,
            cutout_corners_type,
            cutout_bend_radius,
            path_type="cutout",
        )

        cutout_path_points = _get_cutout_path(
            path_points,
            cutout_start_end_offset,
        )

    ground_cutout_path = gdspy.FlexPath(
        cutout_path_points,
        cutout_width,
        corners=cutout_corners_type,
        bend_radius=cutout_bend_radius,
        gdsii_path=True,
        layer=mask_builder.layers.Nb_Groundplane.number,
        datatype=mask_builder.layers.Nb_Groundplane.datatype,
    )
    for poly_points in mbu.get_polys_from_flexpath(ground_cutout_path):
        mask_builder.ground_plane_cutouts.add(
            gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
        )

    ###########################################################################
    # Dielectric path
    ###########################################################################
    if dielectric_in_cutout_or_positive == "positive":
        dielectric_cell = mask_builder.silicon_nitride_positives

    elif dielectric_in_cutout_or_positive == "cutout":
        dielectric_cell = mask_builder.silicon_nitride_cutouts
    else:
        allowed_cutout_types = ["positive", "cutout"]
        raise ValueError(f"`dielectric_in_cutout_or_positive` should only be of values { allowed_cutout_types }")

    if dielectric_width is not None:
        if dielectric_start_end_offset is None:
            dielectric_path_points = cpw_points
        else:
            path_points = _get_path_points_by_corner_type(
                cpw_points,
                dielectric_width,
                dielectric_corners_type,
                dielectric_bend_radius,
                path_type="dielectric",
            )

            dielectric_path_points = _get_cutout_path(
                path_points,
                dielectric_start_end_offset,
            )
        dielectric_path = gdspy.FlexPath(
            dielectric_path_points,
            dielectric_width,
            corners=dielectric_corners_type,
            bend_radius=dielectric_bend_radius,
            gdsii_path=True,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )

        for poly_points in mbu.get_polys_from_flexpath(dielectric_path):
            dielectric_path_segment_poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.SiN_dep.number,
                datatype=mask_builder.layers.SiN_dep.datatype,
            )
            dielectric_cell.add(dielectric_path_segment_poly)

    ###########################################################################
    # Bridging
    ###########################################################################
    DEFAULT_CUTOUT_BRIDGE_HEIGHT_MULTIPIER: float = 1.2
    DEFAULT_DIELECTRIC_BRIDGE_HEIGHT_MULTIPIER: float = 1.1

    ###########################################################################
    # Groundplane Bridging
    ###########################################################################
    if (cutout_bridge_gap is not None) or (cutout_bridge_width is not None):
        if cutout_bridge_gap is None:
            raise ValueError(f"`cutout_bridge_gap` needs to be defined in order to add groundplane bridges.")
        if cutout_bridge_width is None:
            raise ValueError(f"`cutout_bridge_width` needs to be defined in order to add groundplane bridges.")

        if cutout_bridge_height is None:
            cutout_bridge_height = DEFAULT_CUTOUT_BRIDGE_HEIGHT_MULTIPIER * cutout_width

        path_points = _get_path_points_by_corner_type(
            cpw_points,
            cutout_width,
            cutout_corners_type,
            cutout_bend_radius,
            path_type="cutout",
        )

        _add_bridges(
            mask_builder,
            mask_builder.layers.Nb_Groundplane,
            path_points,
            cutout_bridge_gap,
            cutout_bridge_width,
            cutout_bridge_height,
        )

    ###########################################################################
    # Dielectric Bridging
    ###########################################################################
    if (dielectric_bridge_gap is not None) or (dielectric_bridge_width is not None):
        if dielectric_bridge_gap is None:
            raise ValueError(f"`dielectric_bridge_gap` needs to be defined in order to add dielectric bridges.")
        if dielectric_bridge_width is None:
            raise ValueError(f"`dielectric_bridge_width` needs to be defined in order to add dielectric bridges.")

        if dielectric_bridge_height is None:
            if dielectric_width is None:
                dielectric_bridge_height = DEFAULT_DIELECTRIC_BRIDGE_HEIGHT_MULTIPIER * cutout_width
            else:
                dielectric_bridge_height = DEFAULT_DIELECTRIC_BRIDGE_HEIGHT_MULTIPIER * dielectric_width

        path_points = _get_path_points_by_corner_type(
            cpw_points,
            dielectric_width,
            dielectric_corners_type,
            dielectric_bend_radius,
            path_type="dielectric",
        )

        _add_bridges(
            mask_builder,
            mask_builder.layers.SiN_dep,
            path_points,
            dielectric_bridge_gap,
            dielectric_bridge_width,
            dielectric_bridge_height,
        )
