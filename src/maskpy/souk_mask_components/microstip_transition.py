from collections.abc import Sequence

import gdspy
import numpy as np
from numpy import cos, sin

from .. import mask_builder_utils as mbu
from ..layers import Layer
from ..logging import TextColor, styled_text
from ..souk_mask_configs import SoukMaskConfig, get_mask_default_config


def add_microstrip_to_cpw_transition(
    mask_builder,
    microstrip_linewidth: float | int,
    microstrip_end_x_y_rot: Sequence[float | int],
    cpw_center_width: float | int,
    cpw_cutout_width: float | int,
    cpw_end_x_y_rot: Sequence[float | int] | None = None,
    microstrip_layer: Layer | None = None,
    antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None,
) -> None:
    """Adds a co-planar waveguide feedline that transitons to a microstrip
    that. The microstrip will start at the microstrip_end_x_y_rot given. If the
    cpw_end_x_y_rot arg is given then it will try to connect the end of the
    transiton to that coordinate with the desired rotation, this could fail
    however if there isn't the space required. If the cpw_end_x_y_rot is not
    given the transiton will be drawn in a straight line out from the values in
    microstrip_end_x_y_rot. The parameters of the transiton geometry is within
    the config.

    Parameters
    ----------
    microstrip_linewidth: float | int,
        The line width for the microstrip end of the transiton.

    microstrip_end_x_y_rot: Sequence[float | int],
        The [x, y, rotation] for the microstrip end of the transiton. The
        rotation should be given in radians and zero angle is defined as
        pointing in the positive x direction.

    cpw_center_width: float | int,
        The center line width of the cpw end of the transiton.

    cpw_cutout_width: float | int,
        The cutout width of the cpw end of the transiton.

    KwArgs
    ------
    cpw_end_x_y_rot: Sequence[float | int] | None = None
        This is the [x, y, rotation] for the cpw end of the transiton. Default
        value None will create the transiton in a straight line out from the
        microstrip_end_x_y_rot. When this is defined it will try to add to
        connect the end of the transiton to that coordinate with the desired
        rotation in a smooth circular path. The rotation should be given in
        radians and zero angle is defined as pointing in the positive x
        direction.

    microstrip_layer: Layer | None = None
        The material layer for the center line of the transiton. Default value
        is None which will add it to mask_builder.layers.Nb_Antenna.

    antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None,
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.
    """
    cpw_trans_config = get_mask_default_config(
        SoukMaskConfig.ANTENNA_CPW_MICROSTRIP_TRANS,
        config_override=antenna_cpw_microstrip_trans_config_override,
    )

    if microstrip_layer is None:
        microstrip_layer = mask_builder.layers.Nb_Antenna

    elif not isinstance(microstrip_layer, Layer):
        raise TypeError(
            styled_text(
                f"microstrip_layer should be of type Layer, current type is {type(microstrip_layer)}",
                color=TextColor.WARNING,
            )
        )

    cpw_tans_lens = [
        cpw_trans_config["CPW_tans_len1"],
        cpw_trans_config["CPW_tans_len2"],
        cpw_trans_config["CPW_tans_len3"],
        cpw_trans_config["CPW_tans_len4"],
        cpw_trans_config["CPW_tans_len5"],
    ]
    gaps = [
        cpw_trans_config["gap1"],
        cpw_trans_config["gap2"],
        cpw_trans_config["gap3"],
        cpw_trans_config["gap4"],
        cpw_trans_config["gap5"],
    ]

    total_transition_length = np.sum(cpw_tans_lens) + np.sum(gaps)

    # the cutout rectangles
    trans_cutout_rect_1 = gdspy.Rectangle(
        [gaps[0], -cpw_cutout_width / 2],
        [gaps[0] + cpw_tans_lens[0], cpw_cutout_width / 2],
    )
    trans_cutout_rect_2 = gdspy.Rectangle(
        [np.sum(gaps[0:1]) + np.sum(cpw_tans_lens[0:1]) + gaps[1], -cpw_cutout_width / 2],
        [np.sum(gaps[0:1]) + np.sum(cpw_tans_lens[0:1]) + gaps[1] + cpw_tans_lens[1], cpw_cutout_width / 2],
    )
    trans_cutout_rect_3 = gdspy.Rectangle(
        [np.sum(gaps[0:2]) + np.sum(cpw_tans_lens[0:2]) + gaps[2], -cpw_cutout_width / 2],
        [np.sum(gaps[0:2]) + np.sum(cpw_tans_lens[0:2]) + gaps[2] + cpw_tans_lens[2], cpw_cutout_width / 2],
    )
    trans_cutout_rect_4 = gdspy.Rectangle(
        [np.sum(gaps[0:3]) + np.sum(cpw_tans_lens[0:3]) + gaps[3], -cpw_cutout_width / 2],
        [np.sum(gaps[0:3]) + np.sum(cpw_tans_lens[0:3]) + gaps[3] + cpw_tans_lens[3], cpw_cutout_width / 2],
    )
    trans_cutout_rect_5 = gdspy.Rectangle(
        [np.sum(gaps[0:4]) + np.sum(cpw_tans_lens[0:4]) + gaps[4], -cpw_cutout_width / 2],
        [np.sum(gaps[0:4]) + np.sum(cpw_tans_lens[0:4]) + gaps[4] + cpw_tans_lens[4], cpw_cutout_width / 2],
    )

    transitions_cell = gdspy.Cell("CPWTransitionForAntennas")
    transitions_cell.add(
        [
            trans_cutout_rect_1,
            trans_cutout_rect_2,
            trans_cutout_rect_3,
            trans_cutout_rect_4,
            trans_cutout_rect_5,
        ]
    )

    # the points for the center line of the transiton
    top_points = []
    bot_points = []
    for i in range(len(gaps)):
        top_points.append([np.sum(gaps[0:i]) + np.sum(cpw_tans_lens[0:i]), microstrip_linewidth / 2])
        top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(cpw_tans_lens[0:i]), microstrip_linewidth / 2])
        top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(cpw_tans_lens[0:i]), cpw_center_width / 2])
        top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(cpw_tans_lens[0 : i + 1]), cpw_center_width / 2])

        bot_points.append([np.sum(gaps[0:i]) + np.sum(cpw_tans_lens[0:i]), -microstrip_linewidth / 2])
        bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(cpw_tans_lens[0:i]), -microstrip_linewidth / 2])
        bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(cpw_tans_lens[0:i]), -cpw_center_width / 2])
        bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(cpw_tans_lens[0 : i + 1]), -cpw_center_width / 2])

    cpw_transition_poly_points = top_points
    cpw_transition_poly_points.extend(bot_points[::-1])

    end_of_transition = [
        microstrip_end_x_y_rot[0] + total_transition_length * cos(microstrip_end_x_y_rot[2]),
        microstrip_end_x_y_rot[1] + total_transition_length * sin(microstrip_end_x_y_rot[2]),
    ]

    transitions_cutouts = gdspy.CellReference(
        transitions_cell,
        (
            microstrip_end_x_y_rot[0],
            microstrip_end_x_y_rot[1],
        ),
        rotation=mbu.rad_to_deg(microstrip_end_x_y_rot[2]),
    )
    mask_builder.ground_plane_cutouts.add(transitions_cutouts)

    cpw_transition_center_points = mbu.rotate_and_move_points_list(
        cpw_transition_poly_points,
        microstrip_end_x_y_rot[2],
        microstrip_end_x_y_rot[0],
        microstrip_end_x_y_rot[1],
    )
    cpw_transition_center = gdspy.Polygon(
        cpw_transition_center_points,
        layer=microstrip_layer.number,
        datatype=microstrip_layer.datatype,
    )
    mask_builder.Main.add(cpw_transition_center)
    # straight line out from microstrip_end_x_y_rot.
    if cpw_end_x_y_rot is None:
        return
    # else try to connect up the cpw_end_x_y_rot.
    distance_between_microstrip_end_and_cpw_end = np.sqrt(
        (cpw_end_x_y_rot[0] - microstrip_end_x_y_rot[0]) ** 2 + (cpw_end_x_y_rot[1] - microstrip_end_x_y_rot[1]) ** 2
    )
    if distance_between_microstrip_end_and_cpw_end < total_transition_length:
        raise ValueError(
            styled_text(
                f"The distance between microstrip_end and cpw_end is shorter than {total_transition_length}um required by the total transition length.\nmicrostrip_end_x_y_rot = {microstrip_end_x_y_rot}\ncpw_end_x_y_rot = {cpw_end_x_y_rot}",
                TextColor.ERROR,
            )
        )

    mid_point, radius = mbu.get_intersection_point_and_curve_radius(
        end_of_transition[0],
        end_of_transition[1],
        microstrip_end_x_y_rot[2],
        cpw_end_x_y_rot[0],
        cpw_end_x_y_rot[1],
        cpw_end_x_y_rot[2],
    )

    # print(f"end_of_transition[0] = {end_of_transition[0]}")
    # print(f"end_of_transition[1] = {end_of_transition[1]}")
    # print(f"microstrip_end_x_y_rot[2] = {microstrip_end_x_y_rot[2]}")
    # print(f"cpw_end_x_y_rot[0] = {cpw_end_x_y_rot[0]}")
    # print(f"cpw_end_x_y_rot[1] = {cpw_end_x_y_rot[1]}")
    # print(f"cpw_end_x_y_rot[2] = {cpw_end_x_y_rot[2]}")
    #
    # print(f"mid_point: {mid_point}")
    # print(f"radius: {radius}")

    path_points = [
        end_of_transition,
        mid_point,
        [cpw_end_x_y_rot[0], cpw_end_x_y_rot[1]],
    ]
    path_from_transition_end_to_cpw_end = gdspy.FlexPath(
        path_points,
        microstrip_linewidth,
        corners="circular bend",
        bend_radius=radius,
        layer=microstrip_layer.number,
        datatype=microstrip_layer.datatype,
    )
    mask_builder.make_flexpath_into_polygons_and_add_to_main(
        path_from_transition_end_to_cpw_end,
        microstrip_layer,
    )

    # Testing
    mask_builder.make_flexpath_into_polygons_and_add_to_main(
        gdspy.FlexPath(
            path_points,
            cpw_center_width,
            layer=mask_builder.layers.General_labeling.number,
            datatype=mask_builder.layers.General_labeling.datatype,
        ),
        mask_builder.layers.General_labeling,
    )
    # end Testing
