import gdspy
from numpy import cos, sin

from .. import mask_builder_utils as mbu
from ..souk_mask_configs import SoukMaskConfig, get_mask_default_config


def add_antenna(
    mask_builder,
    x: float | int,
    y: float | int,
    rot: float | int,
    antenna_config_override: dict[str, float | int] | None = None,
    add_grnd_cutout: bool = True,
    add_SiN_dep_cutout: bool = True,
    add_backside_check: bool = True,
    return_configurator_points: bool = False,
):
    """Adds the antenna geometries to the Main cell. These are the 4
    antennas in at the center of each horn block. Optionally adds a ground
    plane cutout and Silicon Nitride deoposition cutout as a circle around
    the antenna structure. These geometries are defined by the dimensions
    within the Main_config_file_dict.

    Parameters
    ----------
    x,y: float, int
        The x,y coordinates about which to center the antenna structure.

    rot: float, int
        The angle (**in degrees**) the antenna structure should be rotated.

    KwArgs
    ------
    antenna_config_override: dict[str, float | int] | None = None
        This is an optional override dictionary containing key value pairs for
        variable name and that variable's value respectively. Any keys required
        that do not exist in this dict will be got from the default config. If
        extra keys that are not expected are provided a warnimg will be printed
        but nothing is done with those.

    add_grnd_cutout=True
        Whether or not to add a circular cutout in the Nb_Groundplane layer
        around the antenna structure.

    add_SiN_dep_cutout=True
        Whether or not to add a circular cutout in the SiN depositon layer
        around the antenna structure.

    add_backside_check=True
        Whether or not to add a circle in the backside check layer over the
        antenna structure.

    return_configurator_points=False
        return a the points for use in the configurator.
    """
    antenna_config = get_mask_default_config(SoukMaskConfig.ANTENNA, config_override=antenna_config_override)

    distance_from_center = antenna_config["distance_from_center"]
    base_width = antenna_config["base_width"]
    top_conect_width = antenna_config["top_conect_width"]
    straight_height = antenna_config["straight_height"]
    taper_height = antenna_config["taper_height"]
    backside_check_circle_radius = antenna_config["backside_check_circle_radius"]

    sin_dep_cutout_circle_radius = antenna_config["sin_dep_cutout_circle_radius"]
    grnd_cutout_circle_radius = antenna_config["grnd_cutout_circle_radius"]

    # NO LONGER HAS BACKWARD COMPATABILITY - provide an override for this.
    # default_cicle_cutout_radius = distance_from_center + straight_height + taper_height
    # # Checking for cutout circle radii in config to ensure backward compatability with old config files. Reverts to the old default value if not existing
    # if "sin_dep_cutout_circle_radius" in antenna_config:
    #     sin_dep_cutout_circle_radius = antenna_config["sin_dep_cutout_circle_radius"]
    # else:
    #     sin_dep_cutout_circle_radius = default_cicle_cutout_radius
    #
    # if "grnd_cutout_circle_radius" in antenna_config:
    #     grnd_cutout_circle_radius = antenna_config["grnd_cutout_circle_radius"]
    # else:
    #     grnd_cutout_circle_radius = default_cicle_cutout_radius

    antena_geometry = [
        [(top_conect_width / 2) + x, -(taper_height + straight_height + distance_from_center) + y],
        [-(top_conect_width / 2) + x, -(taper_height + straight_height + distance_from_center) + y],
        [-(base_width / 2) + x, -(straight_height + distance_from_center) + y],
        [-(base_width / 2) + x, -(distance_from_center) + y],
        [(base_width / 2) + x, -(distance_from_center) + y],
        [(base_width / 2) + x, -(straight_height + distance_from_center) + y],
    ]  # this is the shape of the antenna with the point at the bottom

    ant_bot = gdspy.Polygon(
        antena_geometry,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )  # this defines the bottom, right, top and left antennas as the same bottom antenna shape
    ant_right = gdspy.Polygon(
        antena_geometry,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    ant_top = gdspy.Polygon(
        antena_geometry,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )
    ant_left = gdspy.Polygon(
        antena_geometry,
        layer=mask_builder.layers.Nb_Antenna.number,
        datatype=mask_builder.layers.Nb_Antenna.datatype,
    )

    ant_bot.rotate(mbu.deg_to_rad(rot), (x, y))  # this then roatates the antenna by the roation passed into the method
    ant_top.rotate(
        mbu.deg_to_rad(rot + 180), (x, y)
    )  # this also rotates the top, left and right antennas which are all the same as the bottom to form 4 antennas orthogonal to each other
    ant_left.rotate(mbu.deg_to_rad(rot + 270), (x, y))
    ant_right.rotate(mbu.deg_to_rad(rot + 90), (x, y))

    mask_builder.Main.add(ant_bot)  # adding the antennas to the main cell
    mask_builder.Main.add(ant_top)
    mask_builder.Main.add(ant_left)
    mask_builder.Main.add(ant_right)

    # Adding the circular antenna cutout to the grnd plane
    if add_grnd_cutout:
        grnd_cutout_circle = gdspy.Round(
            [x, y],
            grnd_cutout_circle_radius,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        mask_builder.ground_plane_cutouts.add(grnd_cutout_circle)
    if add_SiN_dep_cutout:
        SiN_dep_cutout_circle = gdspy.Round(
            [x, y],
            sin_dep_cutout_circle_radius,
            layer=mask_builder.layers.SiN_dep.number,
            datatype=mask_builder.layers.SiN_dep.datatype,
        )
        mask_builder.silicon_nitride_cutouts.add(SiN_dep_cutout_circle)
    if add_backside_check:
        backside_check_circle = gdspy.Round(
            [x, y],
            backside_check_circle_radius,
            layer=mask_builder.layers.Backside_Check.number,
            datatype=mask_builder.layers.Backside_Check.datatype,
        )
        mask_builder.Main.add(backside_check_circle)

        mirrored_y_backside_check_circle = gdspy.Round(
            [-x, y],
            backside_check_circle_radius,
            layer=mask_builder.layers.Backside_Check.number,
            datatype=mask_builder.layers.Backside_Check.datatype,
        )
        mask_builder.MainBackside.add(mirrored_y_backside_check_circle)

    if not return_configurator_points:
        return

    configurator_points = {}

    ant_top_points = ant_top.polygons[0]
    dx = (ant_top_points[3][0] - ant_top_points[4][0]) / 2
    dy = (ant_top_points[3][1] - ant_top_points[4][1]) / 2

    configurator_points["distance_from_center"] = {
        "text": "distance_from_center",
        "start": [x, y],
        "end": [ant_top_points[3][0] - dx, ant_top_points[3][1] - dy],
    }

    configurator_points["base_width"] = {
        "text": "base_width",
        "start": [ant_top_points[3][0], ant_top_points[3][1]],
        "end": [ant_top_points[4][0], ant_top_points[4][1]],
    }

    configurator_points["top_conect_width"] = {
        "text": "top_conect_width",
        "start": [ant_top_points[0][0], ant_top_points[0][1]],
        "end": [ant_top_points[1][0], ant_top_points[1][1]],
    }

    configurator_points["straight_height"] = {
        "text": "straight_height",
        "start": [ant_top_points[2][0], ant_top_points[2][1]],
        "end": [ant_top_points[3][0], ant_top_points[3][1]],
    }

    mid_dx = (ant_top_points[2][0] - ant_top_points[5][0]) / 2
    mid_dy = (ant_top_points[2][1] - ant_top_points[5][1]) / 2
    top_dx = (ant_top_points[1][0] - ant_top_points[0][0]) / 2
    top_dy = (ant_top_points[1][1] - ant_top_points[0][1]) / 2

    configurator_points["taper_height"] = {
        "text": "taper_height",
        "start": [ant_top_points[2][0] - mid_dx, ant_top_points[2][1] - mid_dy],
        "end": [ant_top_points[1][0] - top_dx, ant_top_points[1][1] - top_dy],
    }

    configurator_points["antena_rotation"] = {
        "text": "antena_rotation",
        "start": [x, ant_top_points[1][1] - top_dy],
        "end": [ant_top_points[1][0] - top_dx, ant_top_points[1][1] - top_dy],
    }

    annotate_ang = 0.0
    configurator_points["grnd_cutout_circle_radius"] = {
        "text": "grnd_cutout_circle_radius",
        "start": [x, y],
        "end": [
            x + grnd_cutout_circle_radius * cos(annotate_ang),
            y + grnd_cutout_circle_radius * sin(annotate_ang),
        ],
    }

    annotate_ang -= 0.2
    configurator_points["backside_check_circle_radius"] = {
        "text": "backside_check_circle_radius",
        "start": [x, y],
        "end": [
            x + backside_check_circle_radius * cos(annotate_ang),
            y + backside_check_circle_radius * sin(annotate_ang),
        ],
    }

    annotate_ang -= 0.2
    configurator_points["sin_dep_cutout_circle_radius"] = {
        "text": "sin_dep_cutout_circle_radius",
        "start": [x, y],
        "end": [
            x + sin_dep_cutout_circle_radius * cos(annotate_ang),
            y + sin_dep_cutout_circle_radius * sin(annotate_ang),
        ],
    }

    return configurator_points
