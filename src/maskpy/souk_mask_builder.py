# import threading
# import concurrent.futures
import datetime
import time
from collections.abc import Callable, Sequence
from typing import Any, Literal

import gdspy
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from numpy import cos, pi, sin, tan
from numpy.typing import NDArray
from shapely import difference as shapely_difference
from shapely import geometry as shapely_geom
from tqdm import tqdm
from typing_extensions import runtime

from . import amber_resonators, layers
from . import mask_builder_utils as mbu
from . import measurements
from . import souk_mask_components as smc
from . import souk_resonators
from .amber_resonators import AmberResonatorType
from .layers import ANLMaskLayerSet, Layer, SoukMaskLayerSet
from .logging import TextColor, pretty_print, styled_type_error
from .souk_mask_configs import SoukMaskConfig, get_mask_default_config
from .souk_resonators import SoukResonator, SoukResonatorType


class SoukMaskBuilder:
    """SoukMaskBuilder is the base for making elements for an SOUK gds mask.

    KwArgs
    ------
    wafer_radius: float = measurements.six_inch,
        The wafer radius, the default is a six inch wafer.

    wafer_flat_length: float = measurements.six_inch_wafer_flat_length,
        The primary flat length for the wafer, The default is for a six
        inch wafer.

    main_cell_name: str = "MAIN_CELL",
        The name for the Main Cell.

    layers: SoukMaskLayerSet | None = None,
        The SOUK layer set for the mask. If not specified then a default
        SoukMaskLayerSet is created. See maskpy.layers.SoukMaskLayerSet.

    layers_ANL: ANLMaskLayerSet | None = None,
        The ANL layer set for the mask. If not specified then a default
        ANLMaskLayerSet is created. See maskpy.layers.ANLMaskLayerSet.


    Cell Attributes
    ---------------
    - Main
    - MainBackside
    - ground_plane_cutouts
    - ground_plane_positives
    - silicon_nitride_cutouts
    - silicon_nitride_positives
    - silicon_oxide_cutouts
    - silicon_oxide_positives
    - silicon_nitride_membrane_cutouts
    - silicon_nitride_membrane_positives

    layers Attribute
    ----------------
    The name of all the layers used in the functions to generate souk mask.
    See maskpy.layers.SoukMaskLayerSet and maskpy.layers.ANLMaskLayerSet.
    """

    def __init__(
        self,
        wafer_radius: float = measurements.three_inch,
        wafer_flat_length: float = measurements.six_inch_wafer_flat_length,
        main_cell_name: str = "MAIN_CELL",
        layers: SoukMaskLayerSet | None = None,
        layers_ANL: ANLMaskLayerSet | None = None,
    ):
        """Create the mask builder which builds the mask.

        KwArgs
        ------
        wafer_radius: float = measurements.six_inch,
            The wafer radius, the default is a six inch wafer.

        wafer_flat_length: float = measurements.six_inch_wafer_flat_length,
            The primary flat length for the wafer, The default is for a six
            inch wafer.

        main_cell_name: str = "MAIN_CELL",
            The name for the Main Cell.

        layers: SoukMaskLayerSet | None = None,
            The SOUK layer set for the mask. If not specified then a default
            SoukMaskLayerSet is created. See maskpy.layers.SoukMaskLayerSet.

        layers_ANL: ANLMaskLayerSet | None = None,
            The ANL layer set for the mask. If not specified then a default
            ANLMaskLayerSet is created. See maskpy.layers.ANLMaskLayerSet.
        """
        # Creating the main cell
        gdspy.current_library = gdspy.GdsLibrary()
        self.Main = gdspy.Cell(main_cell_name)
        self.MainBackside = gdspy.Cell(main_cell_name + "_BACKSIDE")

        # Creating cells for all the layers where boolean operations are required.
        self.ground_plane_cutouts = gdspy.Cell("GROUND_PLANE_CUTOUTS")
        self.ground_plane_positives = gdspy.Cell("GROUND_PLANE_POSITIVES")

        self.silicon_nitride_cutouts = gdspy.Cell("SILICON_NITRIDE_CUTOUTS")
        self.silicon_nitride_positives = gdspy.Cell("SILICON_NITRIDE_POSITIVES")

        self.silicon_oxide_cutouts = gdspy.Cell("SILICON_OXIDE_CUTOUTS")
        self.silicon_oxide_positives = gdspy.Cell("SILICON_OXIDE_POSITIVES")

        self.silicon_nitride_membrane_cutouts = gdspy.Cell("SILICON_NITRIDE_MEMBRANE_CUTOUTS")
        self.silicon_nitride_membrane_positives = gdspy.Cell("SILICON_NITRIDE_MEMBRANE_POSITIVES")

        self.aluminium_etch_cutouts = gdspy.Cell("ALUMINIUM_ETCH_CUTOUTS")
        self.aluminium_etch_positives = gdspy.Cell("ALUMINIUM_ETCH_POSITIVES")

        self.nb_etch_cutouts = gdspy.Cell("NIOBIUM_ETCH_CUTOUTS")
        self.nb_etch_positives = gdspy.Cell("NIOBIUM_ETCH_POSITIVES")

        self.nb_patch_cutouts = gdspy.Cell("NIOBIUM_PATCH_CUTOUTS")
        self.nb_patch_positives = gdspy.Cell("NIOBIUM_PATCH_POSITIVES")

        self.wafer_radius = wafer_radius
        self.wafer_flat_length = wafer_flat_length

        self.resonators_on_mask: list[SoukResonator] = []

        # prevents doing bool operations twice because its expensive.
        self.boolean_operations_completed_flag = False

        if layers is None:
            self.layers = SoukMaskLayerSet()
        elif isinstance(layers, SoukMaskLayerSet):
            self.layers = layers
        else:
            raise TypeError("Error creating SoukMaskBuilder, layers should be of type `SoukMaskLayerSet` or `None` for defualt values.")

        if layers_ANL is None:
            self.layers_ANL = ANLMaskLayerSet()
        elif isinstance(layers_ANL, ANLMaskLayerSet):
            self.layers_ANL = layers_ANL
        else:
            raise TypeError("Error creating SoukMaskBuilder, layers_ANl should be of type `ANLMaskLayerSet` or `None` for defualt values.")

        return

    def add_gds_file_to_mask(
        self,
        file_name: str,
        xy: list[float | int] | tuple[float | int, float | int],
        layer_map: dict[tuple[int, int] | list[int], Layer],
        file_path: str | None = None,
        placement: Literal[
            "center",
            "bot_left",
            "bot_center",
            "bot_right",
            "center_left",
            "center_right",
            "top_left",
            "top_center",
            "top_right",
        ] = "center",
    ) -> None:
        """Add a gds mask file to the mask.

        Parameters
        ----------
        file_name: str,
            The filename for the gds mask file. If this does not include the
            `.gds` file extention then it will be automatically inserted.

        xy: list[float | int] | tuple[float | int, float | int]
            list containing the [x, y] coordinates to add the mask file at.
            This means that if a gds file with a shape centered on (0,0) is
            added, that shpae is not centerd at the xy given.

        layer_map: dict[tuple[int, int] | list[int] (len(2)), Layer]
            This is a map from layer number in the mask file to add and the
            Layer to add that to in the mask builder. This is a dictionary
            with keys with tuple [layer_number, datatype] that exist in the
            mask file to be added. The values for those keys are Layer types
            which are the layer to add the key layer to. e.g.
            >>> layer_map = {
            >>>     (1,0): Layer("Aluminium", 1, 0),
            >>>     (2,0): Layer("Niobium", 2, 0),
            >>> }
            >>> # alternatively
            >>> layer_map = {
            >>>     [1,0]: mask_builder.layers.Aluminium,
            >>>     (2,0): mask_builder.layers.Niobium,
            >>> }

        KwArgs
        ------
        file_path: str | None = None
            The path for the file. When not specified it will look in the same
            directory as the python script for the file_name specified. When
            specified this should be an absolute path for the gds file.

        placement: Literal[
            "center",
            "bot_left",
            "bot_center",
            "bot_right",
            "center_left",
            "center_right",
            "top_left",
            "top_center",
            "top_right",
        ] = "center",
            The placement point to place the mask at. Default is center which
            will place the (0,0) point of the gds file at the xy given. When
            any of the other allowed values are given, e.g. "bot_left", the
            bottom left corner of the entire files bbox will be placed at the
            xy given. For "mid_right", the middle y of the right side of the
            entire files bbox is placed at the xy given. The allowed values
            for placement are:
            >>> "center"
            >>> "bot_left"
            >>> "bot_center"
            >>> "bot_right"
            >>> "center_left"
            >>> "center_right"
            >>> "top_left"
            >>> "top_center"
            >>> "top_right"
        """
        mbu.add_gds_file_to_mask(
            self,
            file_name,
            xy,
            layer_map,
            file_path=file_path,
            placement=placement,
        )

    def generate_hex_pack_circlepack(self, hexagon_radius, hexagon_center, x_pitch, y_pitch):
        print(x_pitch, y_pitch)
        hex_grid = []
        inside_horiz = True
        inside_vert = True
        horz_count = 0
        vert_count = 0

        bot_left_hex = [hexagon_center[0] + hexagon_radius * cos(4 * (pi / 3)), hexagon_center[1] + hexagon_radius * sin(4 * (pi / 3))]

        init_y = bot_left_hex[1] + (x_pitch / 2)

        while inside_vert:
            gy = init_y + vert_count * y_pitch
            if gy <= hexagon_center[1]:
                init_x = (
                    -hexagon_radius
                    - ((gy - (x_pitch / 2) - hexagon_center[0]) / np.sqrt(3))
                    + ((x_pitch / 2) / tan(pi / 3))
                    + hexagon_center[0]
                )
            elif gy > hexagon_center[1]:
                init_x += x_pitch / 2

            if mbu.inside_hexagon([init_x, gy + (y_pitch / 2)], hexagon_radius, hexagon_center):
                inside_horiz = True
                while inside_horiz:
                    gx = init_x + horz_count * x_pitch
                    horz_count += 1

                    if mbu.inside_hexagon([gx + (x_pitch / 2), gy], hexagon_radius, hexagon_center):
                        hex_grid.append([gx, gy])
                    else:
                        inside_horiz = False
                        horz_count = 0
                        vert_count += 1
            else:
                inside_vert = False

        return hex_grid

    def generate_hex_pack(self, hexagon_radius, hexagon_center, x_pitch, y_pitch, middle_gap_distance=0, adjust_center_hex_grid=True):
        hex_grid = []
        inside_horiz = True
        inside_vert = True
        horz_count = 0
        vert_count = 0

        bot_left_hex = [hexagon_center[0] + hexagon_radius * cos(4 * (pi / 3)), hexagon_center[1] + hexagon_radius * sin(4 * (pi / 3))]

        init_y = bot_left_hex[1] + (y_pitch / 2)

        while inside_vert:
            gy = init_y + vert_count * y_pitch
            if gy <= hexagon_center[1]:
                init_x = -hexagon_radius - ((gy - (y_pitch / 2) - hexagon_center[0]) / np.sqrt(3)) + (x_pitch / 2) + hexagon_center[0]
            elif gy > hexagon_center[1]:
                gy += middle_gap_distance
                # init_x += (x_pitch/2)
                init_x = -hexagon_radius + ((gy + (y_pitch / 2) - hexagon_center[0]) / np.sqrt(3)) + (x_pitch / 2) + hexagon_center[0]

            if mbu.inside_hexagon([init_x, gy + (y_pitch / 2)], hexagon_radius, hexagon_center):
                inside_horiz = True
                while inside_horiz:
                    gx = init_x + horz_count * x_pitch
                    horz_count += 1

                    if mbu.inside_hexagon([gx + (x_pitch / 2), gy], hexagon_radius, hexagon_center):
                        hex_grid.append([gx, gy])
                    else:
                        inside_horiz = False
                        horz_count = 0
                        vert_count += 1
            else:
                inside_vert = False

        if adjust_center_hex_grid:
            adjusted_center_hex_grid = self.center_up_hex_grid(hex_grid, hexagon_center)
            return adjusted_center_hex_grid

        return hex_grid

    def generate_half_split_hex_pack(
        self,
        hexagon_radius: float | int,
        hexagon_center: list[float | int],
        x_pitch: float | int,
        y_pitch: float | int,
        middle_gap_distance: float | int = 0,
        adjust_center_hex_grid: bool = False,
        return_two_grids: bool = False,
        return_split_point: bool = False,
    ) -> list[list[float | int]] | tuple[list[list[float | int]], int] | list[list[list[float | int]]]:
        """Generates a split hex pack grid within a hexagon. This is a hex pack
        that has a horizontal split along the center.

        Parameters
        ----------
        hexagon_radius: float, int
            The radius of the hexagon in which to hex pack within. Radius refers to
            the distance from the center of the hexagon to any hex point.

        hexagon_center: list[float | int]
            List of [x,y] coordinates of the center of the hexagon to generate the
            hex pack within.

        x_pitch, y_pitch: float, int
            The horizontal, vertical distance between the hex pack grid points.
            Each row in the grid is offset by half the horizontal pitch.

        KwArgs
        ------
        middle_gap_distance: float | int = 0
            The middle gap between the two halfs of the hex pack grid.

        adjust_center_hex_grid: bool = False
            Determines whether to adjust the halfs of the hex pack such they they
            sit evenly horizontally about the center of the hexagon they are
            confined within.

        return_two_grids: bool = False
            Determines whether or not to return two seperate grids. If False as by
            default, it will pack the two halfs of the grid into one long list
            containing [x,y] points defining the hex pack grid points. If True this
            function will return a list containing the two halfs of the grid,
            [bot_hex_pack, top_hex_pack], where each is a list of [x,y] points.

        return_split_point: bool = False
            Determines whether or not to return the point at which the hex pack
            transitions from the bot_hex_pack to the top_hex_pack. This is simply
            the length of the bot_hex_pack list. This has no effect if the
            return_two_grids argument is True.

        Returns
        -------
        DEFAULT
            hex_grid: lsit[list[float | int]]
            list containing [x,y] lists which define the hex pack grid points.

        OR if return_split_point=True
            RETURNS hex_grid, len(bot_hex_grid): tuple[list[list[float | int]], int]
            returns the hex grid and the point at which the hex pack transitions
            from the bot_hex_pack to the top_hex_pack.

        OR if return_two_grids=True
            RETURNS [bot_hex_grid, top_hex_grid]
                Returns a list containing the two halfs, [bot_hex_grid, top_hex_grid],
                where each is a list of [x,y] lists which define the hex pack grid
                points below and about the center of the bounding hexagon.
        """
        hex_grid = []

        bot_hex_grid = []
        top_hex_grid = []

        inside_horiz = True
        horz_count = 0

        bot_left_hex = [hexagon_center[0] + hexagon_radius * cos(2 * (pi / 3)), hexagon_center[1] + hexagon_radius * sin(4 * (pi / 3))]
        top_left_hex = [hexagon_center[0] + hexagon_radius * cos(2 * (pi / 3)), hexagon_center[1] + hexagon_radius * sin(2 * (pi / 3))]

        hex_height = top_left_hex[1] - bot_left_hex[1]

        number_of_full_rows = int(hex_height / y_pitch)
        number_of_rows_bot = int(number_of_full_rows / 2)

        init_middle_gap = hex_height - (number_of_full_rows * y_pitch)

        adjustment_to_middle_gap = (init_middle_gap - middle_gap_distance) / 2

        init_y_bot = bot_left_hex[1] + (y_pitch / 2) + adjustment_to_middle_gap
        init_y_top = top_left_hex[1] - (y_pitch / 2) - adjustment_to_middle_gap

        for row_number in range(number_of_full_rows):
            if row_number <= number_of_rows_bot:
                gy = init_y_bot + row_number * y_pitch
                init_x = -hexagon_radius - ((gy - (y_pitch / 2) - hexagon_center[0]) / np.sqrt(3)) + (x_pitch / 2) + hexagon_center[0]

            else:
                gy = init_y_top - (row_number - number_of_rows_bot - 1) * y_pitch
                init_x = -hexagon_radius + ((gy + (y_pitch / 2) - hexagon_center[0]) / np.sqrt(3)) + (x_pitch / 2) + hexagon_center[0]

            inside_horiz = True
            while inside_horiz:
                gx = init_x + horz_count * x_pitch
                horz_count += 1

                if mbu.inside_hexagon([gx + (x_pitch / 2), gy - (y_pitch / 2)], hexagon_radius, hexagon_center) and mbu.inside_hexagon(
                    [gx + (x_pitch / 2), gy + (y_pitch / 2)], hexagon_radius, hexagon_center
                ):
                    if row_number <= number_of_rows_bot:
                        bot_hex_grid.append([gx, gy])
                    else:
                        top_hex_grid.append([gx, gy])
                else:
                    inside_horiz = False
                    horz_count = 0

        if adjust_center_hex_grid:
            adjusted_center_bot_hex_grid = self.center_up_hex_grid(bot_hex_grid, hexagon_center, only_horizontal_center=True)
            adjusted_center_top_hex_grid = self.center_up_hex_grid(top_hex_grid, hexagon_center, only_horizontal_center=True)

            if return_two_grids:
                return [adjusted_center_bot_hex_grid, adjusted_center_top_hex_grid]

            for x, y in adjusted_center_bot_hex_grid:
                hex_grid.append([x, y])

            for x, y in adjusted_center_top_hex_grid:
                hex_grid.append([x, y])

            if return_split_point:
                return hex_grid, len(bot_hex_grid)

            return hex_grid

        if return_two_grids:
            return [bot_hex_grid, top_hex_grid]

        for x, y in bot_hex_grid:
            hex_grid.append([x, y])

        for x, y in top_hex_grid:
            hex_grid.append([x, y])

        if return_split_point:
            return hex_grid, len(bot_hex_grid)

        return hex_grid

    def center_up_hex_grid(self, hex_grid, wafer_center_xy, only_horizontal_center=False, only_vertical_center=False):
        """Takes in a list of [x,y] grid points and centers them about a
        central point.

        Parameters
        ----------
        hex_grid: list
            list of [x,y] lists that define the points in the hex pack grid.

        wafer_center_xy: list
            list containing the [x,y] coordinate of the center about which to
            center the grid.

        KwArgs
        ------
        only_horizontal_center=False
            Whether or not to only center up the grid in the horizontal direction.

        only_vertical_center=False
            Whether or not to only center up the grid in the vertical direction.

        Returns
        -------
        shifted_grid: list
            This is a list of the same form as the hex_grid passed in but each
            [x,y] has been adjusted accordingly by some dx, dy to center the grid.
        """
        sorted_hex_grid = sorted(hex_grid.copy(), key=lambda k: [k[1], k[0]])
        left = min(sorted_hex_grid)[0]
        right = max(sorted_hex_grid)[0]
        bot = sorted_hex_grid[0][1]
        top = sorted_hex_grid[-1][1]

        if only_horizontal_center and only_vertical_center:
            only_horizontal_center = False
            only_vertical_center = False

        dx = ((wafer_center_xy[0] - left) - (right - wafer_center_xy[0])) / 2 if not only_vertical_center else 0
        dy = ((wafer_center_xy[1] - bot) - (top - wafer_center_xy[1])) / 2 if not only_horizontal_center else 0

        shifted_grid = []
        for i in range(len(sorted_hex_grid)):
            shifted_grid.append([sorted_hex_grid[i][0] + dx, sorted_hex_grid[i][1] + dy])

        return shifted_grid

    def generate_octagonal_holder_port_dict(
        self,
        center_xy: list[float | int] | tuple[float | int, float | int],
        octagon_rotation: float | int,
        wafer_radius: float | int = 75000.0,
    ) -> dict[str, list[float]]:
        """Generates a dictionary of port x, y, and rotaions for the octagonal
        holder design. This contains 16 port locations where the two on each
        edge of the octagon are 10000um from the center point of that edge.
        This dictionary has keys "0" through "16" where "0" is the left most
        face upper, the numbering continues going counter-clockwise. each key
        contains the x,y and roation of the octagon port conection points.

        Parameters
        ----------
        center_xy: list
            list or tuple containing the [x,y] coordinate for the center point
            which to place the octagonal holder around.

        octagon_rotation: float | int
            the angle (**in Degrees**) for the rotaion of the holder around the
            center_xy argument.

        KwArgs
        ------
        wafer_radius: float | int
            This is the radius of the wafer where the holder positions will be.
            Default is 75000.0 which is a 3inch radius for a six inch wafer.

        Returns
        -------
        oct_ports_dict: dict[str, list[float]]
            dictionary containing the x,y,rot(degrees) for each port in the
            octagonal holder.
        """
        return smc.generate_octagonal_holder_port_dict(
            self,
            center_xy,
            octagon_rotation,
            wafer_radius=wafer_radius,
        )

    def make_Toms_6_inch_holder_and_get_ports(
        self,
        chip_center_xy: list[float | int] | tuple[float | int, float | int],
        wafer_flat_length: float | int,
        middle_flat_angle: float | int = (-pi / 2),
    ) -> dict[int, list[float]]:
        """Generates a dictionary of port x, y, and rotaions for the 6inch
        wafer holder design. This contains 18 port locations where the three on
        each edge of the center point and vertexes of a hexagon edge. This
        dictionary has keys "0" through "18" where "0" is botttom most middle
        right of the wafer where the flat is located at the bottom, the
        numbering continues going counter-clockwise. Each dict key contains the
        x,y and roation of the port conection points.

        Parameters
        ----------
        center_xy: list[float | int] | tuple[float | int, float | int]
            list or tuple containing the [x,y] coordinate for the center point
            which to place the octagonal holder around.

        wafer_flat_length: float | int
            The length of the flat on the silicon wafer.

        KwArgs
        ------
        middle_flat_angle: float | int
            Angle (**in radians**) the falt is located. Default is (-pi/2),
            ie the bottom middle side of the hexagon wafer. This angle is the
            angle the middle of the flat makes with the chip center.

        Returns
        -------
        ports_dict: dict[int, list[float]]
            dictionary containing the x,y,rot(**in radians**) for each port in
            the 6inch holder.
        """
        return smc.make_Toms_6_inch_holder_and_get_ports(
            self,
            chip_center_xy,
            wafer_flat_length,
            middle_flat_angle=middle_flat_angle,
        )

    def get_SOUK_boundary_hexagon(self) -> list[list[float | int]]:
        """Generates the boundary polygon points for the SOUK holder boundary.

        Returns
        -------
        polygon_points: list[list[float | int]]
            list of [x, y] lists that are the vertices of the boundary polygon.
        """
        middle_gap_width = 136.56e3

        # edge_straigh_length = 63.74e3
        corner_straight_length = 5.2e3

        corner_middle_inset = (corner_straight_length / 2) / tan(pi / 3)

        hex_corner_to_straight_length = np.sqrt(corner_middle_inset**2 + (corner_straight_length / 2) ** 2)

        hex_radius = (middle_gap_width + (2 * corner_middle_inset)) / 2

        polygon_points = []

        for i in range(6):
            angle_from_center = i * (pi / 3)

            hex_point_x = hex_radius * cos(angle_from_center)
            hex_point_y = hex_radius * sin(angle_from_center)

            right_side_point = [
                hex_point_x + hex_corner_to_straight_length * cos(angle_from_center - (2 * pi / 3)),
                hex_point_y + hex_corner_to_straight_length * sin(angle_from_center - (2 * pi / 3)),
            ]

            left_side_point = [
                hex_point_x + hex_corner_to_straight_length * cos(angle_from_center + (2 * pi / 3)),
                hex_point_y + hex_corner_to_straight_length * sin(angle_from_center + (2 * pi / 3)),
            ]

            polygon_points.append(right_side_point)
            polygon_points.append(left_side_point)

        return polygon_points

    def make_wafer_shape(
        self,
        layer: Layer,
        wafer_radius: float | int | None = None,
        wafer_flat_length: float | int | None = None,
        chip_center: list[float | int] = [0, 0],
        flat_angle: float = (-pi / 2),
    ) -> gdspy.PolygonSet:
        """Makes the shape of the silicon wafer with a flat on one side and
        returns this shape as a PolygonSet in layer 0 or the layer specified.

        Parameters
        ----------
        layer: Layer
            This is an instance of Layer. see maskpy.layers.Layer.
            Usually this is within the SoukMaskBuilder.layers.xxx.
            e.g. `self.layers.Aluminium`

        KwArgs
        ------
        wafer_radius: float | int | None = None
            The diameter of the wafer in microns.

        wafer_flat_length: float | int | None = None
            The length of the flat on the silicon wafer.

        chip_center: list
            List containing the [x,y] coordinate for the center of the chip.
            Default is [0,0].

        flat_angle: float, int
            Angle (**in radians**) the falt is located. Default is (-pi/2),
            ie the bottom middle side of the hexagon wafer. This angle is the
            angle the middle of the flat makes with the chip center.

        Returns
        -------
        final_wafer: gdspy PolygonSet
            The final wafer shape made with the cutout flat in it on the layer
            specified by the layer.
        """
        if not isinstance(layer, Layer):
            raise TypeError(f"layer should be of type Layer, current type is {type(layer)}")

        if wafer_radius is None:
            wafer_radius = self.wafer_radius

        if wafer_flat_length is None:
            wafer_flat_length = self.wafer_flat_length

        wafer_cirlce = gdspy.Round(chip_center, wafer_radius)

        center_to_flat_dist = (wafer_radius**2 - ((wafer_flat_length**2) / 4)) ** 0.5
        flat_to_edge_dist = wafer_radius - center_to_flat_dist

        middle_of_flat = [
            chip_center[0] + (center_to_flat_dist * cos(flat_angle)),
            chip_center[1] + (center_to_flat_dist * sin(flat_angle)),
        ]

        flat_cutout_rect = gdspy.Rectangle(
            [middle_of_flat[0] - (wafer_flat_length / 2), middle_of_flat[1]],
            [middle_of_flat[0] + (wafer_flat_length / 2), middle_of_flat[1] + flat_to_edge_dist],
        )

        flat_cutout_rect.rotate(flat_angle - (pi / 2), middle_of_flat)

        final_wafer = gdspy.boolean(
            wafer_cirlce,
            flat_cutout_rect,
            "not",
            layer=layer.number,
            datatype=layer.datatype,
        )

        if final_wafer is None:
            raise (RuntimeError("Error creating "))

        return final_wafer

    def add_fancy_text(
        self,
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

        return smc.add_fancy_text(
            self,
            text,
            x,
            y,
            size,
            layer,
            rotation=rotation,
            horizontal_align=horizontal_align,
            vertical_align=vertical_align,
            font_properties=font_properties,
            usetex=usetex,
            return_polyset=return_polyset,
            bb_cutout_in_grnd=bb_cutout_in_grnd,
            bb_cutout_in_sin_dep=bb_cutout_in_sin_dep,
        )

    def add_test_chip_quadrent_boundary_and_get_horn_positions(
        self,
        quadrent_center_xy: list[float] | tuple[float | int, float | int],
        test_chip_quad_config_override: dict[str, float | int] | None = None,
        bottom_left_text: str = "",
        top_right_text: str = "",
        top_left_text: str = "",
        time_stamp_position: Literal["BL", "BR", "ML", "MR", "TL", "TR"] | None = "TL",
        cardiff_logo: bool = True,
        souk_logo: bool = True,
        top_right_label_window: bool = True,
        return_outer_poly_points: bool = False,
        add_center_pin_cutout: bool = True,
        add_slotted_pin_cutout: bool = True,
        dice_line_tab_positions: list[str] = ["left"],
        add_groundplane_under_test_chip: bool = True,
        add_SiN_dep_under_test_chip: bool = False,
        SiN_dep_under_test_chip_edge_offset: float | int = 0,
        add_SiN_membrane_under_test_chip: bool = False,
        add_SiO_under_test_chip: bool = False,
        **kwargs,
    ) -> list[list[float]] | tuple[list[list[float]], list[list[float]]]:
        """Make the boundary for a test chip quad. Adds a centered pin hole and
        slotted pin hole.

        Parameters
        ----------
        quadrent_center_xy: list[flaot]
            list containing the [x,y] coordinates for the center of the test
            chip quad.

        KwArgs
        ------
        test_chip_quad_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        bottom_left_text: str = \"\"
            Text to add the bottom left of the chip. Default is a blank str
            which will not add any text. When not empty this text will be
            placed offset from the bot left corner at a textsize and offset
            specified in the config.

        top_right_text: str = \"\"
            Text to add the top right of the chip. Default is a blank str
            which will not add any text. When not empty this text will be
            placed offset from the top right corner at a textsize and offset
            specified in the config.

        top_left_text: str = \"\"
            Text to add the top left of the chip. Default is a blank str which
            will not add any text. When not empty this text will be placed
            offset from the top left corner at a textsize and offset
            specified in the config.

        time_stamp_position: Literal["BL", "BR", "ML", "MR", "TL", "TR"] | None = "ML",
            This will by default add a time stamp to the MidLeft side of the test
            chip. This can be disabled by passing None as a value. The other
            positions availible are ["BL", "BR", "ML", "MR", "TL", "TR"] where the
            first letter refer to the verical and horizontal position respectively.

        cardiff_logo: bool = True
            Default True will add the cardiff logo to the nb groundplane in the
            bottom right corner of the test chip quad. The size of this is
            specified in the config.

        souk_logo: bool = True
            Default True will add the souk logo to the nb groundplane in the
            bottom right corner of the test chip quad. The size of this is
            specified in the config.

        top_right_label_window: bool = True
            Default True will add a window cutout in the nb groundplane in the
            top right corner of the test chip quad. The size of this window
            is specified in the config.

        return_outer_poly_points: bool = False
            Default False. If true this will return the outer polygon points
            around the test chip quadrent with extra exclusion.

        add_center_pin_cutout: bool = True
            Default True. If false the center pin hole for the test chip will
            not be added to the mask.

        add_slotted_pin_cutout: bool = True
            Default True. If false the slotted pin hole for the test chip will
            not be added to the mask.

        dice_line_tab_positions: list[str]
            The edge positions that the dicing line should be tabbed. Any edges
            not specified will be solid dicing lines. This should be a list of
            strings. The default is a list of just ["left"], but can take str
            values "top", "left", "bot", "right", or "all".

        add_groundplane_under_test_chip: bool = True
            Default True. If false no groundplane will be added below the test
            chip boundary.

        add_SiN_dep_under_test_chip: bool = False
            Default False. If True an SiN dep will be added below the test
            chip boundary.

        SiN_dep_under_test_chip_edge_offset: float | int = 0,
            Default 0. When non-zero the SiN dep will be added below the test
            chip boundary with an offset from the boundary by the amount specified.
            This will only draw if `add_SiN_dep_under_test_chip` is set to `True`.
            Value here can be negaitve for an oversize.

        add_SiN_membrane_under_test_chip: bool = False,
            Default False. If True an SiN membrane will be added below the test
            chip boundary.

        add_SiO_under_test_chip: bool = False,
            Default False. If True an SiO layer will be added below the test
            chip boundary.

        Returns
        -------
        horn_centers: list
            list of [x,y] lists that define the coordinates for the horn
            centers for the 4 horns on the test chip quad. The coords are for
            the top left, top right, bot left, bot right horns respectively.

        outer_poly_points: list
            list of [x,y] lists defining the coordinates of the outer polygon.
            **Note This is only returned if return_outer_poly_points KwArg is
            True**.
        """
        return smc.add_test_chip_quadrent_boundary_and_get_horn_positions(
            self,
            quadrent_center_xy,
            test_chip_quad_config_override=test_chip_quad_config_override,
            bottom_left_text=bottom_left_text,
            top_right_text=top_right_text,
            top_left_text=top_left_text,
            time_stamp_position=time_stamp_position,
            cardiff_logo=cardiff_logo,
            souk_logo=souk_logo,
            top_right_label_window=top_right_label_window,
            return_outer_poly_points=return_outer_poly_points,
            add_center_pin_cutout=add_center_pin_cutout,
            add_slotted_pin_cutout=add_slotted_pin_cutout,
            dice_line_tab_positions=dice_line_tab_positions,
            add_groundplane_under_test_chip=add_groundplane_under_test_chip,
            add_SiN_dep_under_test_chip=add_SiN_dep_under_test_chip,
            SiN_dep_under_test_chip_edge_offset=SiN_dep_under_test_chip_edge_offset,
            add_SiN_membrane_under_test_chip=add_SiN_membrane_under_test_chip,
            add_SiO_under_test_chip=add_SiO_under_test_chip,
            **kwargs,
        )

    def add_cardiff_logo(
        self,
        logo_xy: Sequence[float | int],
        size: float | int = 100,
    ) -> None:
        """Adds the cardiff logo as a cutout in the ng groundplane layer.
        Requires the logo gds file. Logo will be placed centered on the logo_xy
        given.

        Parameters
        ----------
        logo_xy: Sequence[float | int]
            Sequence containing the [x,y] coordinate for the center of where
            the logo should be placed.

        KwArgs
        ------
        size: int, flaot = 100
            The size of the sides of the cardiff logo square.
            The default is 100.
        """
        return smc.add_cardiff_logo(
            self,
            logo_xy,
            size=size,
        )

    def add_so_logo(
        self,
        logo_xy: Sequence[float | int],
        size: float | int = 100,
        draw_all_in_one_layer: bool = True,
    ) -> None:
        """Adds the Simons Observatory logo as a positive in the ng groundplane
        layer by default. Requires the logo gds file. Logo will be placed
        centered on the logo_xy given.

        Parameters
        ----------
        logo_xy: Sequence[float | int]
            Sequence containing the [x,y] coordinate for the center of where the
            logo should be placed.

        KwArgs
        ------
        size: int, flaot = 100
            The size of the sides of the cardiff logo square.
            The default is 100.

        draw_all_in_one_layer: bool = True
            Default True which will draw the entire logo as a cutout in the
            NbGroundplane layer.
        """
        return smc.add_so_logo(
            self,
            logo_xy,
            size=size,
            draw_all_in_one_layer=draw_all_in_one_layer,
        )

    def add_souk_logo(
        self,
        logo_xy: Sequence[float | int],
        size: float | int = 100,
        draw_all_in_one_layer: bool = True,
        include_outer_ring_text: bool = True,
    ) -> None:
        """Adds the Simons Observatory United Kingdom logo as a positive in the
        ng groundplane layer by default. Requires the logo gds file. Logo will
        be placed centered on the logo_xy given.

        Parameters
        ----------
        logo_xy: Sequence[float | int]
            Sequence containing the [x,y] coordinate for the center of where the
            logo should be placed.

        KwArgs
        ------
        size: int | flaot = 100
            The size of the sides of the cardiff logo square.
            The default is 100.

        draw_all_in_one_layer: bool = True
            Default True which will draw the entire logo as a cutout in the
            NbGroundplane layer. When False this will draw the entire logo on
            seperate predefined layers.

        include_outer_ring_text: bool = True
            Default True which will draw the outer ring of text. When False
            everything bar that outer text ring will be drawn.
        """

        return smc.add_souk_logo(
            self,
            logo_xy,
            size=size,
            draw_all_in_one_layer=draw_all_in_one_layer,
            include_outer_ring_text=include_outer_ring_text,
        )

    def add_AM_signature(
        self,
        signature_xy: Sequence[float | int],
        size: float | int = 100,
        variant: Literal[0, 1, 2, 3] = 0,
        draw_as_cutout: bool = True,
    ) -> None:
        """Adds Alan Manning signature as a negative in the nb groundplane
        layer by default. Requires the signature gds file. Signiture will be
        placed centered on the signature_xy given.

        Parameters
        ----------
        signature_xy: Sequence[float | int]
            Sequence containing the [x,y] coordinate for the center of where the
            signature should be placed.

        KwArgs
        ------
        size: int, flaot
            The size of the largest extent of the signature. Default is 100.

        variant: Literal[0, 1, 2, 3] = 0,
            Default 0. Variants are 0 through 3. 0 draws both names side by
            side. 1 draws both names vertically stacked. 2 is jsut the
            forename. 3 is just the surname.

        draw_as_cutout = True
            Default True which will draw the signature as a cutout in the nb
            groundplane. When False the signature is drawn as a aluminium
            positve.
        """

        return smc.add_AM_signature(
            self,
            signature_xy,
            size=size,
            variant=variant,
            draw_as_cutout=draw_as_cutout,
        )

    def add_MLA_marker(
        self,
        x: float | int,
        y: float | int,
        materials: Sequence[Layer] | Layer,
        inner_lw: float | int = 5.0,
        outer_lw: float | int = 10.0,
    ) -> None:
        """Adds a singe or series of square markers with a cross inside
        centered on the x,y coord given. Each subsequent material will make a
        square larger than the previous by 1x the outer linewidth such that
        they are staggered. Inner crosses reamin the same and do not get
        staggered.

        Parameters
        ----------
        x, y: float | int
            x, y coordinate for the center of the MLA marker.

        materials: Sequence[Layer] | Layer
            Sequence of materials for the marker out of. This should be a
            Sequence of Layer types or a single Layer type, e.g.
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

        return smc.add_MLA_marker(
            self,
            x,
            y,
            materials,
            inner_lw=inner_lw,
            outer_lw=outer_lw,
        )

    def add_initial_alignment_markers(
        self,
        x: float | int,
        y: float | int,
        cross_length: float | int = 300,
        linewidth: float | int = 5,
        cutout_square_size: float | int = 1000,
    ) -> None:
        """Adds a singe MLA marker cross centered on the x,y coord given. Also
        adds a cutout window in layers that will be deposited ontop of this
        (namely Nb_groundplane).

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
        return smc.add_initial_alignment_markers(
            self,
            x,
            y,
            cross_length=cross_length,
            linewidth=linewidth,
            cutout_square_size=cutout_square_size,
        )

    def add_caliper_alignment_markers(
        self,
        x: float | int,
        y: float | int,
        rot: float | int,
        layer1: Layer,
        layer2: Layer,
        layer1_text: str,
        layer2_text: str,
    ) -> None:
        """Adds Alignment Markers to the mask in a cutout centered on the x,y
        given. This consists of 2 main corsses and a series of calipers made
        from the two materials given.

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
        return smc.add_caliper_alignment_markers(
            self,
            x,
            y,
            rot,
            layer1,
            layer2,
            layer1_text,
            layer2_text,
        )

    def add_test_H_pad_connected_box_section(
        self,
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
        """Adds a box of size 5500 centered on the x,y given that contains a
        series of 6 pads where 2 of those connect to a line connecting them.

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
        return smc.add_test_H_pad_connected_box_section(
            self,
            x,
            y,
            pad_materials,
            line_materials,
            rot=rot,
            linewidths=linewidths,
            TL_text=TL_text,
            add_text_for_linewidth=add_text_for_linewidth,
            add_dielectric_square=add_dielectric_square,
            add_groundplane_square=add_groundplane_square,
            add_dielectric_cutout_in_ports=add_dielectric_cutout_in_ports,
            add_dielectric_cutout_over_H=add_dielectric_cutout_over_H,
            add_dielectric_cutout_over_line=add_dielectric_cutout_over_line,
        )

    def add_test_linewidths_pad_connected_box_section(
        self,
        x: float | int,
        y: float | int,
        materials: Layer | Sequence[Layer],
        rot: float | int = 0,
        linewidths: Sequence[float | int] | None = None,
        TL_text: str = "",
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
        return smc.add_test_linewidths_pad_connected_box_section(
            self,
            x,
            y,
            materials,
            rot=rot,
            linewidths=linewidths,
            TL_text=TL_text,
            add_dielectric_square=add_dielectric_square,
            add_groundplane_square=add_groundplane_square,
            add_dielectric_cutout_in_ports=add_dielectric_cutout_in_ports,
        )

    def add_test_linewidth_structure_box_section(
        self,
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
        return smc.add_test_linewidth_structure_box_section(
            self,
            x,
            y,
            linewidths=linewidths,
        )

    def add_test_crossover_structure_box_section(
        self,
        x: float | int,
        y: float | int,
        filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None,
        rot: float | int = 0.0,
    ) -> None:
        """Adds a box of size 5500 centered on the x,y given that contains a
        series of 3 crossover sections that connect to pads.

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
        return smc.add_test_crossover_structure_box_section(
            self,
            x,
            y,
            filter_bank_ring_overlap_config_override=filter_bank_ring_overlap_config_override,
            rot=rot,
        )

    def add_test_straight_line_structure_box_section(
        self,
        x: float | int,
        y: float | int,
        rot: float | int = 0.0,
    ):
        """Adds a box of size 5500 centered on the x,y given that contains a
        series of 6 stright line sections that connect to pads.

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
        return smc.add_test_straight_line_structure_box_section(
            self,
            x,
            y,
            rot=rot,
        )

    def make_DC_structure_pads_and_meander(
        self,
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

        return smc.make_DC_structure_pads_and_meander(
            self,
            x,
            y,
            rot,
            DC_structure_material,
        )

    def add_test_DC_structure_box_section(
        self,
        x: float | int,
        y: float | int,
        DC_structure_material: Layer,
        BL_text: str = "",
    ) -> None:
        """Adds a box of size 8000 centered on the x,y given that contains a
        series of 5 test DC lines and pads. Pads are 100x the size of the
        linewidth connecting them.

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

        return smc.add_test_DC_structure_box_section(
            self,
            x,
            y,
            DC_structure_material,
            BL_text=BL_text,
        )

    def add_test_chip_quad_groundplane_edge_cuts(
        self,
        test_chip_quad_center_xy,
        test_chip_quad_width,
        test_chip_quad_height,
        top_mid_gap=0,
        bot_mid_gap=0,
        right_mid_gap=0,
        left_mid_gap=0,
    ):
        """Adds edge cutouts in the groundplane layer around the border of the
        test chip quad. These edge cuts sit within the defined border. There is
        by default edge cuts across every edge but there can be a gap in the
        middle of any edge specified by the #_mid_gap arguments. The size and
        spacing of these edge cuts is defined in config.

        Parameters
        ----------
        test_chip_quad_center_xy: list
            list containing the [x,y] coordinates for the center of the test
            chip quad.

        test_chip_quad_width: int, float
            The width of the test chip quad.

        test_chip_quad_height: int, float
            The height of the test chip quad

        KwArgs
        ------
        top_mid_gap, bot_mid_gap, right_mid_gap, left_mid_gap: float, int
            The gap in the edge cuts in the middle of the respective edge.
        """

        edge_cut_length = 250
        edge_cut_width = 100
        edge_cut_gap_length = 750

        # middle_gap_top = 5000
        # middle_gap_bot = 5000
        # middle_gap_left = 0
        # middle_gap_right = 0

        no_of_edge_cuts_horizontal = int(test_chip_quad_width / (edge_cut_length + edge_cut_gap_length))
        no_of_edge_cuts_vertical = int(test_chip_quad_height / (edge_cut_length + edge_cut_gap_length))

        # Adding the top and bot edge cuts
        for i in range(int(no_of_edge_cuts_horizontal / 2) + 1):
            edge_cut_offset_from_center = i * (edge_cut_length + edge_cut_gap_length)

            # top edge cuts
            if edge_cut_offset_from_center >= (top_mid_gap / 2):
                # right edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] + edge_cut_offset_from_center,
                    test_chip_quad_center_xy[1] + (test_chip_quad_height / 2) - (edge_cut_width / 2),
                ]

                top_edge_cut_rect_right = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_length / 2), edge_cut_center_xy[1] - (edge_cut_width / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_length / 2), edge_cut_center_xy[1] + (edge_cut_width / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(top_edge_cut_rect_right)

                # left edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] + (-1 * edge_cut_offset_from_center),
                    test_chip_quad_center_xy[1] + (test_chip_quad_height / 2) - (edge_cut_width / 2),
                ]

                top_edge_cut_rect_left = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_length / 2), edge_cut_center_xy[1] - (edge_cut_width / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_length / 2), edge_cut_center_xy[1] + (edge_cut_width / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(top_edge_cut_rect_left)

            # bot edge cuts
            if edge_cut_offset_from_center >= (bot_mid_gap / 2):
                # right edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] + edge_cut_offset_from_center,
                    test_chip_quad_center_xy[1] - (test_chip_quad_height / 2) + (edge_cut_width / 2),
                ]

                bot_edge_cut_rect_right = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_length / 2), edge_cut_center_xy[1] - (edge_cut_width / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_length / 2), edge_cut_center_xy[1] + (edge_cut_width / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(bot_edge_cut_rect_right)

                # left edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] + (-1 * edge_cut_offset_from_center),
                    test_chip_quad_center_xy[1] - (test_chip_quad_height / 2) + (edge_cut_width / 2),
                ]

                bot_edge_cut_rect_left = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_length / 2), edge_cut_center_xy[1] - (edge_cut_width / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_length / 2), edge_cut_center_xy[1] + (edge_cut_width / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(bot_edge_cut_rect_left)

        # Adding the left and right edge cuts
        for i in range(int(no_of_edge_cuts_vertical / 2) + 1):
            edge_cut_offset_from_center = i * (edge_cut_length + edge_cut_gap_length)

            # right edge cuts
            if edge_cut_offset_from_center >= (right_mid_gap / 2):
                # top edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] + (test_chip_quad_width / 2) - (edge_cut_width / 2),
                    test_chip_quad_center_xy[1] + edge_cut_offset_from_center,
                ]

                right_edge_cut_rect_top = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_width / 2), edge_cut_center_xy[1] - (edge_cut_length / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_width / 2), edge_cut_center_xy[1] + (edge_cut_length / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(right_edge_cut_rect_top)

                # bot edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] + (test_chip_quad_width / 2) - (edge_cut_width / 2),
                    test_chip_quad_center_xy[1] + (-1 * edge_cut_offset_from_center),
                ]

                right_edge_cut_rect_bot = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_width / 2), edge_cut_center_xy[1] - (edge_cut_length / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_width / 2), edge_cut_center_xy[1] + (edge_cut_length / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(right_edge_cut_rect_bot)

            # left edge cuts
            if edge_cut_offset_from_center >= (left_mid_gap / 2):
                # top edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] - (test_chip_quad_width / 2) + (edge_cut_width / 2),
                    test_chip_quad_center_xy[1] + edge_cut_offset_from_center,
                ]

                left_edge_cut_rect_top = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_width / 2), edge_cut_center_xy[1] - (edge_cut_length / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_width / 2), edge_cut_center_xy[1] + (edge_cut_length / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(left_edge_cut_rect_top)

                # bot edge cut
                edge_cut_center_xy = [
                    test_chip_quad_center_xy[0] - (test_chip_quad_width / 2) + (edge_cut_width / 2),
                    test_chip_quad_center_xy[1] + (-1 * edge_cut_offset_from_center),
                ]

                left_edge_cut_rect_bot = gdspy.Rectangle(
                    [edge_cut_center_xy[0] - (edge_cut_width / 2), edge_cut_center_xy[1] - (edge_cut_length / 2)],
                    [edge_cut_center_xy[0] + (edge_cut_width / 2), edge_cut_center_xy[1] + (edge_cut_length / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                self.ground_plane_cutouts.add(left_edge_cut_rect_bot)

        return

    def add_test_chip_quad_tabbed_dicing_line(
        self,
        test_chip_quad_center_xy: Sequence[float | int],
        test_chip_quad_width: float | int,
        test_chip_quad_height: float | int,
        tab_positions: Sequence[Literal["top", "left", "bot", "right", "all"]] = ["left"],
        corner_overrun_length: float | int = 200,
    ):
        """Adds a tabbed dicing line around the test chip quad. This is a
        tabbed dicing line on all tab_positions specified otherwise will be a
        solid dicing line. The tab positions is by default only the left side
        but can be any edge. There are corner markers by default that
        overextend by 200um but is variable.

        Parameters
        ----------
        test_chip_quad_center_xy: Sequence[float | int]
            list containing the [x,y] coordinates for the center of the test
            chip quad.

        test_chip_quad_width: float | int
            The width of the test chip quad.

        test_chip_quad_height: float | int
            The height of the test chip quad.

        KwArgs
        ------
        tab_positions: Sequence[Literal["top", "left", "bot", "right", "all"]]
            The edge positions that the dicing line should be tabbed. Any edges
            not specified will be solid dicing lines. This should be a list of
            strings. The default is a list of just ["left"], but can take str
            values "top", "left", "bot", "right", or "all".

        corner_overrun_length: int, float
            This is the length the corners of the tabbed dicing line should
            overrun at the corners of the test chip quad. The default is 200.
        """

        tab_length = 4000  # config
        tab_gap_length = 300  # config
        linewidth = 300  # config

        center_x = test_chip_quad_center_xy[0]
        center_y = test_chip_quad_center_xy[1]

        top_left_corner = [center_x - test_chip_quad_width / 2, center_y + test_chip_quad_height / 2]
        top_right_corner = [center_x + test_chip_quad_width / 2, center_y + test_chip_quad_height / 2]
        bot_left_corner = [center_x - test_chip_quad_width / 2, center_y - test_chip_quad_height / 2]
        bot_right_corner = [center_x + test_chip_quad_width / 2, center_y - test_chip_quad_height / 2]

        top_left_corner_poly_points = [
            top_left_corner,
            [top_left_corner[0] - linewidth - corner_overrun_length, top_left_corner[1]],
            [top_left_corner[0] - linewidth - corner_overrun_length, top_left_corner[1] + linewidth],
            [top_left_corner[0] - linewidth, top_left_corner[1] + linewidth],
            [top_left_corner[0] - linewidth, top_left_corner[1] + linewidth + corner_overrun_length],
            [top_left_corner[0], top_left_corner[1] + linewidth + corner_overrun_length],
        ]
        top_left_corner_poly = gdspy.Polygon(
            top_left_corner_poly_points,
            layer=self.layers.Tab_dicing_line.number,
            datatype=self.layers.Tab_dicing_line.datatype,
        )
        self.Main.add(top_left_corner_poly)

        top_right_corner_poly_points = [
            top_right_corner,
            [top_right_corner[0] + linewidth + corner_overrun_length, top_right_corner[1]],
            [top_right_corner[0] + linewidth + corner_overrun_length, top_right_corner[1] + linewidth],
            [top_right_corner[0] + linewidth, top_right_corner[1] + linewidth],
            [top_right_corner[0] + linewidth, top_right_corner[1] + linewidth + corner_overrun_length],
            [top_right_corner[0], top_right_corner[1] + linewidth + corner_overrun_length],
        ]
        top_right_corner_poly = gdspy.Polygon(
            top_right_corner_poly_points,
            layer=self.layers.Tab_dicing_line.number,
            datatype=self.layers.Tab_dicing_line.datatype,
        )
        self.Main.add(top_right_corner_poly)

        bot_left_corner_poly_points = [
            bot_left_corner,
            [bot_left_corner[0] - linewidth - corner_overrun_length, bot_left_corner[1]],
            [bot_left_corner[0] - linewidth - corner_overrun_length, bot_left_corner[1] - linewidth],
            [bot_left_corner[0] - linewidth, bot_left_corner[1] - linewidth],
            [bot_left_corner[0] - linewidth, bot_left_corner[1] - linewidth - corner_overrun_length],
            [bot_left_corner[0], bot_left_corner[1] - linewidth - corner_overrun_length],
        ]
        bot_left_corner_poly = gdspy.Polygon(
            bot_left_corner_poly_points,
            layer=self.layers.Tab_dicing_line.number,
            datatype=self.layers.Tab_dicing_line.datatype,
        )
        self.Main.add(bot_left_corner_poly)

        bot_right_corner_poly_points = [
            bot_right_corner,
            [bot_right_corner[0] + linewidth + corner_overrun_length, bot_right_corner[1]],
            [bot_right_corner[0] + linewidth + corner_overrun_length, bot_right_corner[1] - linewidth],
            [bot_right_corner[0] + linewidth, bot_right_corner[1] - linewidth],
            [bot_right_corner[0] + linewidth, bot_right_corner[1] - linewidth - corner_overrun_length],
            [bot_right_corner[0], bot_right_corner[1] - linewidth - corner_overrun_length],
        ]
        bot_right_corner_poly = gdspy.Polygon(
            bot_right_corner_poly_points,
            layer=self.layers.Tab_dicing_line.number,
            datatype=self.layers.Tab_dicing_line.datatype,
        )
        self.Main.add(bot_right_corner_poly)

        top_center = [center_x, center_y + test_chip_quad_height / 2]
        right_center = [center_x + test_chip_quad_width / 2, center_y]
        bot_center = [center_x, center_y - test_chip_quad_height / 2]
        left_center = [center_x - test_chip_quad_width / 2, center_y]

        sides = ["top", "right", "bot", "left"]
        side_centers = [top_center, right_center, bot_center, left_center]

        side_lengths = [test_chip_quad_width, test_chip_quad_height, test_chip_quad_width, test_chip_quad_height]
        side_rotations = [0, (3 / 2) * pi, pi, (1 / 2) * pi]

        # Drawing the tabbed lines, center line, then all the ones out from that.
        for side, side_center, side_length, side_rotation in zip(sides, side_centers, side_lengths, side_rotations, strict=False):
            half_length_of_side = side_length / 2

            # If this side is not in the tab_positions list then just draw a solid dice line and move to next side
            if side not in tab_positions and tab_positions != ["all"]:
                solid_dice_line = gdspy.FlexPath(
                    [
                        [side_center[0] - half_length_of_side, side_center[1] + (linewidth / 2)],
                        [side_center[0] + half_length_of_side, side_center[1] + (linewidth / 2)],
                    ],
                    linewidth,
                    ends="flush",
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                solid_dice_line.rotate(side_rotation, center=side_center)
                self.Main.add(solid_dice_line)

                continue

            len_of_first_half_mid_tab_and_post_gap = (tab_length / 2) + (linewidth / 2) + tab_gap_length
            len_of_full_tab_and_post_gap = (linewidth / 2) + tab_length + (linewidth / 2) + tab_gap_length
            no_of_tabs_after_middle = int((half_length_of_side - len_of_first_half_mid_tab_and_post_gap) / (len_of_full_tab_and_post_gap))

            center_tab = gdspy.FlexPath(
                [
                    [side_center[0] - (tab_length / 2), side_center[1] + (linewidth / 2)],
                    [side_center[0] + (tab_length / 2), side_center[1] + (linewidth / 2)],
                ],
                linewidth,
                ends="round",
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            center_tab.rotate(side_rotation, center=side_center)
            self.Main.add(center_tab)

            init_offset = len_of_first_half_mid_tab_and_post_gap + (linewidth / 2) + (tab_length / 2)

            # Loop for all the full tabs out from the center
            for i in range(no_of_tabs_after_middle):
                offset = init_offset + i * len_of_full_tab_and_post_gap
                left_side_tab = gdspy.FlexPath(
                    [
                        [side_center[0] - offset - (tab_length / 2), side_center[1] + (linewidth / 2)],
                        [side_center[0] - offset + (tab_length / 2), side_center[1] + (linewidth / 2)],
                    ],
                    linewidth,
                    ends="round",
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                left_side_tab.rotate(side_rotation, center=side_center)
                self.Main.add(left_side_tab)

                right_side_tab = gdspy.FlexPath(
                    [
                        [side_center[0] + offset - (tab_length / 2), side_center[1] + (linewidth / 2)],
                        [side_center[0] + offset + (tab_length / 2), side_center[1] + (linewidth / 2)],
                    ],
                    linewidth,
                    ends="round",
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                right_side_tab.rotate(side_rotation, center=side_center)
                self.Main.add(right_side_tab)

            # the finishing tabs to the corner that arent full length
            offset += len_of_full_tab_and_post_gap
            left_finishing_tab = gdspy.FlexPath(
                [
                    [side_center[0] - half_length_of_side, side_center[1] + (linewidth / 2)],
                    [side_center[0] - offset + (tab_length / 2), side_center[1] + (linewidth / 2)],
                ],
                linewidth,
                ends="round",
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            left_finishing_tab.rotate(side_rotation, center=side_center)
            self.Main.add(left_finishing_tab)

            right_finishing_tab = gdspy.FlexPath(
                [
                    [side_center[0] + half_length_of_side, side_center[1] + (linewidth / 2)],
                    [side_center[0] + offset - (tab_length / 2), side_center[1] + (linewidth / 2)],
                ],
                linewidth,
                ends="round",
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            right_finishing_tab.rotate(side_rotation, center=side_center)
            self.Main.add(right_finishing_tab)

        return

    def add_text_under_horn_in_test_chip_quad(
        self,
        text_string: str,
        horn_center_x: float | int,
        horn_center_y: float | int,
        font_props: FontProperties | None = None,
    ):
        """Adds a label under the horn on the test chip quad with the text
        string given.

        Parameters
        ----------
        text_string: str
            The text string to be placed under the horn.

        horn_center_x, horn_center_x: float, int
            The x, and y coordinate respectively for the center of the horn
            that needs the text underneath.

        KwArgs
        ------
        font_properties: FontProperties | None = None,
            This is any options for the format of the text, family, style etc.
            See matplotlib.font_manager.FontProperties for all the options.
            If None provided, this will use a default FontProperties,
            font_properties=FontProperties(family="monospace", style="normal").
        """
        x_offset_from_horn = 0
        y_offset_from_horn = -2800

        text_size = 300

        text_x = horn_center_x + x_offset_from_horn
        text_y = horn_center_y + y_offset_from_horn

        self.add_fancy_text(
            text_string,
            text_x,
            text_y,
            text_size,
            layer=self.layers.Aluminium,
            horizontal_align="center",
            vertical_align="below",
            font_properties=font_props,
            bb_cutout_in_grnd=True,
        )

        return

    def get_config_excel(self, file_name, sheet_name, path="config_files\\"):
        """Gets the config excel file containing the variable names and
        respective vales. different sheets contain different vairables.

        Parameters
        ----------
        file_name: str
            String including file extention, eg "myfile.xlsx"
        sheet_name: str
            String for the name of the sheet in the excel file.
        path: str
            Path to containing folder. Deafult is "config_files" folder in same
            directory as this file.

        Returns
        -------
        out: 'config_dict'
            Dictionary with keys of var names and values of var values.

        Raises
        ------
            FileNotFoundError: If files doesn't exist.
            Exception: For all other file opening issues.
        """
        raise NotImplementedError("Depreciated")

        # live_file = path + file_name
        #
        # try:
        #     df = pd.read_excel(live_file, sheet_name=sheet_name, header=None)
        # except FileNotFoundError:
        #     raise FileNotFoundError(f"Config file '{path + file_name}' not found.")
        #     print("Config file '" + path + file_name + "' not found.")
        #     return
        # except Exception as e:
        #     raise Exception(f"Issue with file '{path + file_name}'\n{e}")
        #
        # config_dict = {}
        # for i in range(len(df)):
        #     var_name = df[0][i]
        #     var_val = df[1][i]
        #     config_dict[var_name] = var_val
        #
        # return config_dict

    def make_flexpath_into_polygons_and_add_to_main(
        self,
        flexpath: gdspy.FlexPath,
        layer: Layer,
    ):
        """Gets the polygon points for the shape created by a gdspy Flexpath
        object. and then add that to Main as polygons.

        Parameters
        ----------
        flexpath: object
            an instanciated gdspy.Flexpath object.

        layer: Layer
            This is an instance of Layer. see maskpy.layers.Layer.
            Usually this is within the SoukMaskBuilder.layers.xxx.
            e.g. `self.layers.Aluminium`
        """
        if not isinstance(layer, Layer):
            raise TypeError(f"layer should be of type Layer, current type is {type(layer)}")

        poly_points_from_flexpath = mbu.get_polys_from_flexpath(flexpath)
        for i in range(len(poly_points_from_flexpath)):
            path_polygon = gdspy.Polygon(
                poly_points_from_flexpath[i],
                layer=layer.number,
                datatype=layer.datatype,
            )
            self.Main.add([path_polygon])

        return

    def add_antenna(
        self,
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
        return smc.add_antenna(
            self,
            x,
            y,
            rot,
            antenna_config_override=antenna_config_override,
            add_grnd_cutout=add_grnd_cutout,
            add_SiN_dep_cutout=add_SiN_dep_cutout,
            add_backside_check=add_backside_check,
            return_configurator_points=return_configurator_points,
        )

    def add_4_resonators_around_horn(
        self,
        x: float | int,
        y: float | int,
        rel_kid_positions: Sequence[Sequence[float | int]],
        KID_Nos: Sequence[int] | NDArray,
        f0s: Sequence[float | int] | NDArray,
        resonator_types: Sequence[SoukResonatorType],
        mux_func_overrides: Sequence[Callable | None] = [None, None, None, None],
        resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None],
        general_config_override: dict[str, float | int] | None = None,
        IDC_and_frame_materials: Sequence[str] = ("IDC_Nb", "IDC_Nb", "IDC_Nb", "IDC_Nb"),
        meander_materials: Sequence[str] = ("Al", "Al", "Al", "Al"),
        trim_lengths: Sequence[float | None] = (None, None, None, None),
        add_grnd_cutout: Sequence[bool] = (True, True, True, True),
        add_SiN_dep_dielectric_around: bool = True,
        add_SiN_dep_dielectric_cutout: Sequence[bool] = (True, True, True, True),
        add_SiO_cutout: Sequence[bool] = (True, True, True, True),
        add_SiN_membrane_cutout: Sequence[bool] = (True, True, True, True),
        add_backside_check: Sequence[bool] = (False, False, False, False),
        add_grnd_cutout_over_inductor: Sequence[bool] = (False, False, False, False),
        add_SiN_dep_dielectric_cutout_over_inductor: Sequence[bool] = (False, False, False, False),
        add_Aluminium_Patch_and_Etch: Sequence[bool] = (False, False, False, False),
        **kwargs,
    ):
        """Adds the 4 resonator geometries to the mask. The four resonators are
        placed around the center x, y point in order top right, top left, bot
        right, bot left (as should be the order of the rel_kid_positions list).
        The resonator geometries are defined by the dimensions within the
        Main_config_file_dict. By default it will, but optionally can choose
        not to, add all the neccessay cutouts for the structure.

        Parameters
        ----------
        x,y: float, int
            The x,y coordinates about which to center the antenna structure.

        rel_kid_positions: Sequence[Sequence[float | int]]
            Sequence of [x,y] Sequence's defining the points at which to make each KID
            relative to the center of the antenna. This KIDs are made from the very
            bottom center point of the inductive meander section. This list should
            be the positions, in order, of the top left, top right, bot left,
            bot right.

        KID_Nos: Sequence[int]
            Sequence of ints which are the numbers each KID should have drawn next to
            it. This Sequence of numbers for each KID should be the same order as the
            rel_kid_positions, TL, TR, BL, BR.

        f0s: Sequence[float, int]
            The resonant frequencies of the resonators. Should be in the same unit
            that the mux_funcs function takes.

        resonator_types: Sequence[SoukResonatorType]
            This is the type of resonators to be drawn. The values accepted
            here are members of the SoukResonatorType enum.
            The order of the values passed in will be attributed to each KID
            and should be the same order as the rel_kid_positions, TL, TR, BL,
            BR.

        KwArgs
        ------
        mux_func_overrides: Sequence[Callable | None] = [None, None, None, None]
            This is a Sequence of None or callable functions for getting the IDC
            and CC lengths from a given f0. When the default 4 long list of
            None is provided, The resonator's default muxing function will be
            used. The function should take a frequency as an arguments and
            should return an array-like (28 long) and a single float
            **in this order**, the array-like should contain all the lengths
            for each of the 28 arms for the IDC and the float should be the CC
            length. A simple example funtion:

            >>> def example_IDC_CC_func(
            >>>     f0: float | int
            >>>     )->tuple[list[float | int], float | int]:
            >>>     '''Example function for IDC and CC.
            >>>
            >>>     f0: float | int
            >>>         Resonant frequency in Hz.
            >>>     Returns
            >>>     -------
            >>>     IDC_lengths: list[float | int]
            >>>     CC_length: float | int
            >>>     '''
            >>>     if (f0 < 3.1e9):
            >>>         IDC_lengths = np.ones(28)*1900.0
            >>>         CC_length = 600.0
            >>>     else:
            >>>         IDC_lengths = np.ones(28)*1500.0
            >>>         CC_length = 300.0
            >>>
            >>>     return IDC_lengths, CC_length

        resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None]
            This is a Sequence of 4 optional override dictionarys containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        general_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        IDC_and_frame_materials=["IDC_Nb", "IDC_Nb", "IDC_Nb", "IDC_Nb"]
            The material to make each of the IDC and frame structures out of.
            By Default this is a tuple of "IDC_Nb" which will make all the KIDs
            out of the IDC_Nb material. When this is passed in it should be a
            length 4 list of strings that take any of the values "IDC_Nb",
            "Nb", or "Al", which will make the frame and IDC out of IDC_Nb,
            Nb_Antenna, or Aluminium layers respectively. The order of the
            values passed in will be attributed to each KID and should be the
            same order as the rel_kid_positions, TL, TR, BL, BR.

        meander_materials=["Al", "Al", "Al", "Al"]
            The material to make each of the inductive meander structure out of.
            By Default this is a list of "Al" which will make all the meanders
            out of the Aluminium material. When this is passed in it should be
            a length 4 list of strings that take any of the values "Al",
            "IDC_Nb", or "Nb", which will make the inductive meander out of
            Aluminium, IDC_Nb, or Nb_Antenna layers respectively. The order of
            the values passed in will be attributed to each meander and should
            be the same order as the rel_kid_positions, TL, TR, BL, BR.

        trim_lengths=[None,None,None,None]
            Whether not not to add trim boxes to the mask. When this is defined
            it should be a list of 4 ints of floats that define how long the
            trim arms should be on the mask. Trim boxes will be made to cover
            the trim arms to bring them down to the length specified.
            **More info can be found in the add_resonator_original method.**

        add_grnd_cutout=[True, True, True, True]
            Whether or not to add a cutout in the Nb_Groundplane layer in the
            neccary place for the KID structure. The order of the values passed
            in will be attributed to each KID and should be the same order as
            the rel_kid_positions, TL, TR, BL, BR.

        add_SiN_dep_dielectric_around=True
            Whether or not to add an SiN depositon layer around in neccary
            place for the KID structure, this is a box the of the pixel pitch
            around the center of the antenna structure.

        add_SiN_dep_dielectric_cutout=[True, True, True, True]
            Whether or not to add a cutout in the SiN depositon layer in the
            neccary place for the KID structure. The order of the values passed
            in will be attributed to each KID and should be the same order as
            the rel_kid_positions, TL, TR, BL, BR.

        add_SiO_cutout=[True, True, True, True]
            Whether or not to add a cutout in the Silicon Oxide layer in the
            neccary place for the KID structure. The order of the values passed
            in will be attributed to each KID and should be the same order as
            the rel_kid_positions, TL, TR, BL, BR.

        add_SiN_membrane_cutout=[True, True, True, True]
            Whether or not to add a cutout in the Silicon Nitride membrane layer
            in the neccary place for the KID structure. The order of the values
            passed in will be attributed to each KID and should be the same
            order as the rel_kid_positions, TL, TR, BL, BR.

        add_backside_check=[False, False, False, False]
            Whether or not to add a backside check cover in the neccary place
            for the KID structure. The order of the values passed in will be
            attributed to each KID and should be the same order as the
            rel_kid_positions, TL, TR, BL, BR.

        add_grnd_cutout_over_inductor=[False, False, False, False]
            Whether or not to add a groundplane cutout over the inductive
            meander. When true will create a cutout over the mander that is
            oversived by 30um in all directions relative to the center of the
            inductive meander. The order of the values passed in will be
            attributed to each KID and should be the same order as the
            rel_kid_positions, TL, TR, BL, BR.

        add_SiN_dep_dielectric_cutout_over_inductor=[False, False, False, False]
            Whether or not to add a SiN depositon cutout over the inductive
            meander. When true will create a cutout over the mander that is
            oversived by 20um in all directions relative to the center of the
            inductive meander. The order of the values passed in will be
            attributed to each KID and should be the same order as the
            rel_kid_positions, TL, TR, BL, BR.

        add_Aluminium_Patch_and_Etch=[False, False, False, False]
            Whether of not to add an Aluminium patch and etch around aluminium
            elements.
        """
        rot_angles = [pi / 2, -pi / 2, pi / 2, -pi / 2]
        to_mirror = [True, False, False, True]

        if not isinstance(trim_lengths, Sequence):
            styled_type_error(trim_lengths, "trim_lengths", Sequence[float | int | None])
        if not all(isinstance(x, (float | int | None)) for x in trim_lengths):
            raise TypeError("elements of list `trim_lengths` should be of type float | int | None.")

        if not isinstance(IDC_and_frame_materials, Sequence):
            styled_type_error(IDC_and_frame_materials, "IDC_and_frame_materials", Sequence[str])
        if not all(isinstance(x, (str)) for x in IDC_and_frame_materials):
            raise TypeError("elements of list `IDC_and_frame_materials` should be of type str.")

        if not isinstance(meander_materials, Sequence):
            styled_type_error(meander_materials, "meander_materials", Sequence[str])
        if not all(isinstance(x, (str)) for x in meander_materials):
            raise TypeError("elements of list `meander_materials` should be of type str.")

        for k, rel_kid_position in enumerate(rel_kid_positions):
            kid_x = x + rel_kid_position[0]
            kid_y = y + rel_kid_position[1]

            match resonator_types[k]:
                case (
                    SoukResonatorType.ORIGINAL_Q10K
                    | SoukResonatorType.ORIGINAL_Q20K
                    | SoukResonatorType.ORIGINAL_Q50K
                    | SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K
                    | SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K
                    | SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K
                ):
                    self.add_resonator_original(
                        resonator_types[k],
                        kid_x,
                        kid_y,
                        rot_angles[k],
                        f0s[k],
                        mux_func_override=mux_func_overrides[k],
                        resonator_config_override=resonator_config_overrides[k],
                        mirror=to_mirror[k],
                        IDC_and_frame_material=IDC_and_frame_materials[k],
                        meander_material=meander_materials[k],
                        trim_length=trim_lengths[k],
                        add_grnd_cutout=add_grnd_cutout[k],
                        add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout[k],
                        add_SiO_cutout=add_SiO_cutout[k],
                        add_SiN_membrane_cutout=add_SiN_membrane_cutout[k],
                        add_backside_check=add_backside_check[k],
                        add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor[k],
                        add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor[k],
                        add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch[k],
                    )
                case (
                    SoukResonatorType.HIGH_VOLUME_V1_Q20K
                    | SoukResonatorType.HIGH_VOLUME_V1_Q50K
                    | SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K
                    | SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K
                ):
                    self.add_resonator_high_volume_v1(
                        resonator_types[k],
                        kid_x,
                        kid_y,
                        rot_angles[k],
                        f0s[k],
                        mux_func_override=mux_func_overrides[k],
                        resonator_config_override=resonator_config_overrides[k],
                        mirror=to_mirror[k],
                        IDC_and_frame_material=IDC_and_frame_materials[k],
                        meander_material=meander_materials[k],
                        trim_length=trim_lengths[k],
                        add_grnd_cutout=add_grnd_cutout[k],
                        add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout[k],
                        add_SiO_cutout=add_SiO_cutout[k],
                        add_SiN_membrane_cutout=add_SiN_membrane_cutout[k],
                        add_backside_check=add_backside_check[k],
                        add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor[k],
                        add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor[k],
                    )
                case (
                    SoukResonatorType.HIGH_VOLUME_V2_Q20K
                    | SoukResonatorType.HIGH_VOLUME_V2_Q50K
                    | SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K
                    | SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K
                ):
                    self.add_resonator_high_volume_v2(
                        resonator_types[k],
                        kid_x,
                        kid_y,
                        rot_angles[k],
                        f0s[k],
                        mux_func_override=mux_func_overrides[k],
                        resonator_config_override=resonator_config_overrides[k],
                        mirror=to_mirror[k],
                        IDC_and_frame_material=IDC_and_frame_materials[k],
                        meander_material=meander_materials[k],
                        trim_length=trim_lengths[k],
                        add_grnd_cutout=add_grnd_cutout[k],
                        add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout[k],
                        add_SiO_cutout=add_SiO_cutout[k],
                        add_SiN_membrane_cutout=add_SiN_membrane_cutout[k],
                        add_backside_check=add_backside_check[k],
                        add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor[k],
                        add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor[k],
                    )
                case SoukResonatorType.CPW_COUPLED_V1:
                    self.add_resonator_cpw_coupled(
                        resonator_types[k],
                        kid_x,
                        kid_y,
                        rot_angles[k],
                        f0s[k],
                        mux_func_override=mux_func_overrides[k],
                        resonator_config_override=resonator_config_overrides[k],
                        mirror=to_mirror[k],
                        IDC_and_frame_material=IDC_and_frame_materials[k],
                        meander_material=meander_materials[k],
                        trim_length=trim_lengths[k],
                        add_grnd_cutout=add_grnd_cutout[k],
                        add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout[k],
                        add_SiO_cutout=add_SiO_cutout[k],
                        add_SiN_membrane_cutout=add_SiN_membrane_cutout[k],
                        add_backside_check=add_backside_check[k],
                        add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor[k],
                        add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor[k],
                    )

                case _:
                    raise NotImplementedError(f"resonator_type '{resonator_types[k]}' does not have an associated draw function.")

            # Adding the KID number Text
            resonator_config = souk_resonators.get_resonator_config(
                resonator_types[k], resonator_config_override=resonator_config_overrides[k]
            )

            text_size = resonator_config["text_size"]  # 90
            text_x_offset = resonator_config["text_x_offset"]  # 800
            text_y_offset = resonator_config["text_y_offset"]  # 900

            x_sign = 1 if k in [0, 2] else -1
            y_sign = 1 if k in [0, 1] else -1
            horizontal_align = "start" if x_sign == 1 else "end"
            vertical_align = "above" if y_sign == 1 else "below"

            font_props = FontProperties(family="monospace", style="normal")
            kid_no_text_x = kid_x + x_sign * text_x_offset
            kid_no_text_y = kid_y + y_sign * text_y_offset

            self.add_fancy_text(
                str(KID_Nos[k]),
                kid_no_text_x,
                kid_no_text_y,
                text_size,
                layer=self.layers.Aluminium,
                horizontal_align=horizontal_align,
                vertical_align=vertical_align,
                font_properties=font_props,
            )

        # Making the SiN dep dielectric to go around resonators
        if add_SiN_dep_dielectric_around:
            general_config = get_mask_default_config(SoukMaskConfig.GENERAL, config_override=general_config_override)

            dielectric_around_width = general_config["horizontal_pitch"]
            dielectric_around_height = general_config["vertical_pitch"]

            dielectric_around_box = gdspy.Rectangle(
                [x - dielectric_around_width / 2, y - dielectric_around_height / 2],
                [x + dielectric_around_width / 2, y + dielectric_around_height / 2],
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            self.silicon_nitride_positives.add(dielectric_around_box)

        return

    def add_amber_resonator(
        self,
        resonator_type: AmberResonatorType,
        x: float,
        y: float,
        rot_angle: float,
        f0: float | int,
        resonator_config_override: dict | None = None,
        mux_func_override: Callable | None = None,
        mirror=False,
        IDC_and_frame_material: Layer | None = None,
        meander_material: Layer | None = None,
        coupler_fork_material: Layer | None = None,
        add_grnd_cutout: bool = True,
        add_SiN_dep_dielectric_cutout: bool = True,
        add_SiO_cutout: bool = True,
        add_SiN_membrane_cutout: bool = True,
        add_backside_check: bool = False,
        add_inductor_cover: bool = False,
    ) -> None:
        """Test."""

        match resonator_type:
            case AmberResonatorType.ORIGINAL_AL_IND2:
                resonator = amber_resonators.OriginalAlInd2(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    coupler_fork_material=coupler_fork_material,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_inductor_cover=add_inductor_cover,
                )
                self.resonators_on_mask.append(resonator)
                return

            case AmberResonatorType.ORIGINAL_AL_IND4:
                resonator = amber_resonators.OriginalAlInd4(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    coupler_fork_material=coupler_fork_material,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_inductor_cover=add_inductor_cover,
                )
                self.resonators_on_mask.append(resonator)
                return

            case AmberResonatorType.ORIGINAL_NB_IND2:
                resonator = amber_resonators.OriginalNbInd2(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    coupler_fork_material=coupler_fork_material,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_inductor_cover=add_inductor_cover,
                )
                self.resonators_on_mask.append(resonator)
                return

            case AmberResonatorType.ORIGINAL_NB_IND4:
                resonator = amber_resonators.OriginalNbInd4(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    coupler_fork_material=coupler_fork_material,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_inductor_cover=add_inductor_cover,
                )
                self.resonators_on_mask.append(resonator)
                return
            case _:
                raise (ValueError(f"AmberResonatorType '{resonator_type}' does not have an associated draw function."))

    def add_resonator_cpw_coupled(
        self,
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
        add_Aluminium_Patch_and_Etch: bool = True,
        return_configurator_points: bool = False,
    ) -> None:
        """Adds the KID geometry to the Main cell at the x,y cooardinate given.
        The KID is placed where the base middle of the inductive meander is at
        this x,y. The KID geometry is defined by the dimensions within the
        Main_config_file_dict. By default it will, but optionally can choose
        not to, add all the neccessay cutouts for the structure.

        Parameters
        ----------
        x,y: float, int
            The x,y coordinates to place the KID. This is the very bottom center
            point of the inductive meander section.

        rot_angle: float, int
            The rotation angle (**in radians**) for the structure. Positive values
            are anti-clockwise, negative is clockwise. The default rotation is with
            the inductive meander at the bottom and coupler at the top with IDC
            arms running horizontally.

        f0: float, int
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
        match resonator_type:
            case SoukResonatorType.CPW_COUPLED_V1:
                resonator = souk_resonators.CpwCoupledV1(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
                self.resonators_on_mask.append(resonator)
                return
            case _:
                raise (ValueError(f"SoukResonatorType '{resonator_type}' does not have an associated draw function."))

    def add_resonator_original(
        self,
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
        add_Aluminium_Patch_and_Etch: bool = True,
        return_configurator_points: bool = False,
    ):
        """Adds the KID geometry to the Main cell at the x,y cooardinate given.
        The KID is placed where the base middle of the inductive meander is at
        this x,y. The KID geometry is defined by the dimensions within the
        Main_config_file_dict. By default it will, but optionally can choose
        not to, add all the neccessay cutouts for the structure.

        Parameters
        ----------
        resonator_type: SoukResonatorType
            This is the type of resonator to be drawn. The values accepted
            here are a subset of members of the SoukResonatorType enum:
            - SoukResonatorType.ORIGINAL_Q10K
            - SoukResonatorType.ORIGINAL_Q20K
            - SoukResonatorType.ORIGINAL_Q50K
            - SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K
            - SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K
            - SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K

        x,y: float, int
            The x,y coordinates to place the KID. This is the very bottom center
            point of the inductive meander section.

        rot_angle: float, int
            The rotation angle (**in radians**) for the structure. Positive values
            are anti-clockwise, negative is clockwise. The default rotation is with
            the inductive meander at the bottom and coupler at the top with IDC
            arms running horizontally.

        f0: float, int
            The resonant frequency of the resonator. Should be in the same unit
            that the mux_func function takes.

        KwArgs
        ------
        mux_func_override: Callable | None = None
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
            ...     f0: float | int
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

        mirror: bool = False
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

        add_Aluminium_Patch_and_Etch=True
            Whether of not to add an Aluminium patch and etch around aluminium
            elements.

        return_configurator_points=False
            return a the points for use in the configurator.
        """
        match resonator_type:
            case SoukResonatorType.ORIGINAL_Q10K:
                resonator = souk_resonators.OriginalQ10k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.ORIGINAL_Q20K:
                resonator = souk_resonators.OriginalQ20k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.ORIGINAL_Q50K:
                resonator = souk_resonators.OriginalQ50k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q10K:
                resonator = souk_resonators.OriginalLongTrunkQ10k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q20K:
                resonator = souk_resonators.OriginalLongTrunkQ20k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.ORIGINAL_LONG_TRUNK_Q50K:
                resonator = souk_resonators.OriginalLongTrunkQ50k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    add_Aluminium_Patch_and_Etch=add_Aluminium_Patch_and_Etch,
                    return_configurator_points=return_configurator_points,
                )
            case _:
                raise (ValueError(f"SoukResonatorType '{resonator_type}' does not have an associated draw function."))

        self.resonators_on_mask.append(resonator)

        return

    def add_resonator_high_volume_v1(
        self,
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
        add_Aluminium_Patch_and_Etch: bool = True,
        return_configurator_points: bool = False,
    ):
        """Adds the KID geometry to the Main cell athe the x,y cooardinate
        given. The KID is placed where the base middle of the inductive meander
        is at this x,y. The KID geometry is defined by the dimensions within
        the Main_config_file_dict. By default it will, but optionally can
        choose not to, add all the neccessay cutouts for the structure.

        Parameters
        ----------
        resonator_type: SoukResonatorType
            This is the type of resonator to be drawn. The values accepted
            here are a subset of members of the SoukResonatorType enum:
            - SoukResonatorType.HIGH_VOLUME_V1_Q20K
            - SoukResonatorType.HIGH_VOLUME_V1_Q50K
            - SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K
            - SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K

        x,y: float, int
            The x,y coordinates to place the KID. This is the very bottom center
            point of the inductive meander section.

        rot_angle: float, int
            The rotation angle (**in radians**) for the structure. Positive values
            are anti-clockwise, negative is clockwise. The default rotation is with
            the inductive meander at the bottom and coupler at the top with IDC
            arms running horizontally.

        f0: float, int
            The resonant frequency of the resonator. Should be in the same unit
            that the mux_func function takes.

        KwArgs
        ------
        mux_func_override: Callable | None = None
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
            ...     f0: float | int
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

        match resonator_type:
            case SoukResonatorType.HIGH_VOLUME_V1_Q20K:
                resonator = souk_resonators.HighVolumeV1Q20k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.HIGH_VOLUME_V1_Q50K:
                resonator = souk_resonators.HighVolumeV1Q50k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q20K:
                resonator = souk_resonators.HighVolumeV1LongTrunkQ20k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.HIGH_VOLUME_V1_LONG_TRUNK_Q50K:
                resonator = souk_resonators.HighVolumeV1LongTrunkQ50k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case _:
                raise (ValueError(f"SoukResonatorType '{resonator_type}' does not have an associated draw function."))

        self.resonators_on_mask.append(resonator)
        return

    def add_resonator_high_volume_v2(
        self,
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
        resonator_type: SoukResonatorType
            This is the type of resonator to be drawn. The values accepted
            here are a subset of members of the SoukResonatorType enum:
            - SoukResonatorType.HIGH_VOLUME_V2_Q20K
            - SoukResonatorType.HIGH_VOLUME_V2_Q50K
            - SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K
            - SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K

        x,y: float, int
            The x,y coordinates to place the KID. This is the very bottom center
            point of the inductive meander section.

        rot_angle: float, int
            The rotation angle (**in radians**) for the structure. Positive values
            are anti-clockwise, negative is clockwise. The default rotation is with
            the inductive meander at the bottom and coupler at the top with IDC
            arms running horizontally.

        f0: float, int
            The resonant frequency of the resonator. Should be in the same unit

        KwArgs
        ------
        mux_func_override: Callable | None = None
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
            ...     f0: float | int
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
        match resonator_type:
            case SoukResonatorType.HIGH_VOLUME_V2_Q20K:
                resonator = souk_resonators.HighVolumeV2Q20k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.HIGH_VOLUME_V2_Q50K:
                resonator = souk_resonators.HighVolumeV2Q50k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q20K:
                resonator = souk_resonators.HighVolumeV2LongTrunkQ20k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case SoukResonatorType.HIGH_VOLUME_V2_LONG_TRUNK_Q50K:
                resonator = souk_resonators.HighVolumeV2LongTrunkQ50k(
                    self,
                    resonator_type,
                    x,
                    y,
                    rot_angle,
                    f0,
                    mux_func_override=mux_func_override,
                    resonator_config_override=resonator_config_override,
                    mirror=mirror,
                    IDC_and_frame_material=IDC_and_frame_material,
                    meander_material=meander_material,
                    trim_length=trim_length,
                    add_grnd_cutout=add_grnd_cutout,
                    add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                    add_SiO_cutout=add_SiO_cutout,
                    add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                    add_backside_check=add_backside_check,
                    add_grnd_cutout_over_inductor=add_grnd_cutout_over_inductor,
                    add_SiN_dep_dielectric_cutout_over_inductor=add_SiN_dep_dielectric_cutout_over_inductor,
                    return_configurator_points=return_configurator_points,
                )
            case _:
                raise (ValueError(f"SoukResonatorType '{resonator_type}' does not have an associated draw function."))

        self.resonators_on_mask.append(resonator)
        return

    def add_Lo_pass_filters(
        self,
        x: float | int,
        y: float | int,
        inner_ring_line_width: float | int,
        inner_ring_radius: float | int,
        init_angle: float | int,
        direction: str,
        Lo_pass_filters_config_override: dict[str, float | int] | None = None,
        return_configurator_points: bool = False,
    ):
        """Adds the low pass filter arms around the inner ring of the filter
        bank.

        Parameters
        ----------
        x, y: float, int
            x, y coordinate of the center of the antenna to be places around.

        inner_ring_line_width: float, int
            Line width of the inner ring that the filters conect to.

        inner_ring_radius: float, int
            Radius of the inner ring that the filters conect to.

        init_angle: float, int
            Initial angle (*in radians*) to start at when making the filter arms.

        direction: str
            The direction around the inner ring in which to make the filters point.
            Allowed values are `clockwise` or `anti-clockwise`.

        KwArgs
        ------
        Lo_pass_filters_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        return_configurator_points=False
            return a the points for use in the configurator.
        """
        return smc.add_Lo_pass_filters(
            self,
            x,
            y,
            inner_ring_line_width,
            inner_ring_radius,
            init_angle,
            direction,
            Lo_pass_filters_config_override=Lo_pass_filters_config_override,
            return_configurator_points=return_configurator_points,
        )

    def add_Hi_pass_filters(
        self,
        x: float | int,
        y: float | int,
        inner_ring_line_width: float | int,
        inner_ring_radius: float | int,
        init_angle: float | int,
        direction: str,
        Hi_pass_filters_config_override: dict[str, float | int] | None = None,
        return_configurator_points: bool = False,
    ):
        """Adds the high pass filter arms around the inner ring of the filter
        bank.

        Parameters
        ----------
        x, y: float, int
            x, y coordinate of the center of the antenna to be places around.

        inner_ring_line_width: float, int
            Line width of the inner ring that the filters conect to.

        inner_ring_radius: float, int
            Radius of the inner ring that the filters conect to.

        init_angle: float, int
            Initial angle (*in radians*) to start at when making the filter arms.

        direction: str
            The direction around the inner ring in which to make the filters point.
            Allowed values are `clockwise` or `anti-clockwise`.

        KwArgs
        ------
        Hi_pass_filters_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        return_configurator_points=False
            return a the points for use in the configurator.
        """
        return smc.add_Hi_pass_filters(
            self,
            x,
            y,
            inner_ring_line_width,
            inner_ring_radius,
            init_angle,
            direction,
            Hi_pass_filters_config_override=Hi_pass_filters_config_override,
            return_configurator_points=return_configurator_points,
        )

    def add_combiner_section_and_get_conect_point(
        self,
        x: float | int,
        y: float | int,
        rot: float | int,
        outer_ring_conection_gap: float | int,
        outer_ring_linewidth: float | int,
        combiner_section_90ghz_config_override: dict[str, float | int] | None = None,
        combiner_section_150ghz_config_override: dict[str, float | int] | None = None,
        combiner_type: Literal["90GHZ", "150GHZ"] = "90GHZ",
        mirror_combiner: bool = False,
        return_configurator_points: bool = False,
    ) -> list[list[float | int]] | tuple[list[list[float | int]], dict | None]:
        """Adds the phase combiner to the outer ring of the filter bank. This
        is by default the 90GHz or optionally the 150GHz. The specific
        geometries dimensions is determined by the parameters in the config.
        This function will return the coordinate of the conection point.

        Parameters
        ----------
        x, y: float, int
            x, y coordinate for the base of the combiner structure, this is the
            very bottom middle where it connects to the outer ring of the filter
            bank. i.e. the middle of the gap in the outer ring.

        rot: float, int
            The angle (**in radians**) of the rotation of the whole combiner
            geometry.

        outer_ring_conection_gap: float, int
            The gap distance in the outer ring so the combiner can connect to
            either side.

        outer_ring_linewidth: float, int
            The line width of the outer ring of the filter bank this combiner
            attaches to.

        KwArgs
        ------
        combiner_section_90ghz_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        combiner_section_150ghz_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        combiner_type: str = "90GHZ"
            This defines what type of phase combiner to add. This can take string
            values of either "90GHZ" or "150GHZ".

        mirror_combiner: bool = False,
            This defines if the combiner should be mirrored around the line
            perpedicular to the outer ring of the filter bank. The default is
            False or can be True.

        return_configurator_points: bool = False
            return a the points for use in the configurator.

        Returns
        -------
        points_to_conect_kids_to: list[list[float | int]]
            list of [x,y] coordinate for the conection point which connects to a
            KID.
        """
        return smc.add_combiner_section_and_get_conect_point(
            self,
            x,
            y,
            rot,
            outer_ring_conection_gap,
            outer_ring_linewidth,
            combiner_section_90ghz_config_override=combiner_section_90ghz_config_override,
            combiner_section_150ghz_config_override=combiner_section_150ghz_config_override,
            combiner_type=combiner_type,
            mirror_combiner=mirror_combiner,
            return_configurator_points=return_configurator_points,
        )

    def add_filter_bank_ring_overlap_and_get_conections(
        self,
        x: float | int,
        y: float | int,
        ant_center_x: float | int,
        ant_center_y: float | int,
        overlap_no: int,
        rot: float | int,
        filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None,
        return_configurator_points: bool = False,
    ):
        """Adds a ring overlap conection to bridge the inner and outer rings
        over one another. This function will return the conection points where
        the inner and outer rings should connect to.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate of the center of the ring overlap.

        ant_center_x, ant_center_y: float, int
            The x, y coordinate of the center of the antenna structure,
            i.e. the center of the horn.

        overlap_no: int
            The number of the ring overlap. This determines where to draw it around
            the antenna. Starting at 0 for left middle and +1 for each subsequent
            overlap going anti-clockwise. Should not be more than 3, values more
            than this wrap back to 0 (left middle placement) because the overlap_no
            operates like it is modulo 4.

        rot: float | int
            The angle (**in radians**) which the overlap geometry should be rotated
            at. This rot angle defined as the anti-clockwise angle made with the
            positive x-axis.

        KwArgs
        ------
        filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        return_configurator_points: bool = False
            return a the points for use in the configurator.

        Returns
        -------
        conections_dict: dict
            This dictionary contains keys that map to an [x,y] list which are the
            coordinates defining the conection points where the inner and outer
            ring should connect to this overlap structure.

            This dict has keys: **'inner_conect_0', 'inner_conect_1',
            'outer_conect_0', 'outer_conect_1'**.
        """

        return smc.add_filter_bank_ring_overlap_and_get_conections(
            self,
            x,
            y,
            ant_center_x,
            ant_center_y,
            overlap_no,
            rot,
            filter_bank_ring_overlap_config_override=filter_bank_ring_overlap_config_override,
            return_configurator_points=return_configurator_points,
        )

    def add_filter_bank_and_get_conection_points(
        self,
        x: float | int,
        y: float | int,
        filter_bank_config_override: dict[str, float | int] | None = None,
        filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None,
        Hi_pass_filters_config_override: dict[str, float | int] | None = None,
        Lo_pass_filters_config_override: dict[str, float | int] | None = None,
        combiner_section_90ghz_config_override: dict[str, float | int] | None = None,
        combiner_section_150ghz_config_override: dict[str, float | int] | None = None,
        with_combiner: bool = True,
        with_crossover: bool = True,
        only_1_pol: bool = False,
        return_configurator_points: bool = False,
        return_configurator_points_for_Lo_pass: bool = False,
        return_configurator_points_for_Hi_pass: bool = False,
        return_configurator_points_for_combiner_150ghz: bool = False,
        return_configurator_points_for_combiner_90ghz: bool = False,
        **kwargs,
    ):
        """Adds the filter bank structure to the chip centered at the x,y
        coordinate given. By default this will be drawn with phase combiners
        and ring overlap crossovers and with both polorizations filtered.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate to center filter bank structure around.

        KwArgs
        ------
        filter_bank_config_override,
        Hi_pass_filters_config_override,
        Lo_pass_filters_config_override,
        combiner_section_90ghz_config_override,
        combiner_section_150ghz_config_override,
        filter_bank_ring_overlap_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        with_combiner = True
            Adds the combiner along with the filter bank rings by default. When
            Fasle the outer ring has no gap in it and the conection points will be
            the base of where the combiners would have gone, i.e. the middle of the
            outer ring arc sections.

        with_crossover = True
            Adds the crossover between inner and outer rings of the filter bank by
            default. When Fasle, the inner is only added for the first and 3rd
            antenna leaving the other two to conect to nothing. The outer ring is
            removed for the first and second antenna. The Hugh filters remain.

        only_1_pol = False
            If set True this will only draw one polorization of the filter bank.
            This will not draw the relevant inner and outer rings for one
            polorization and adds an inductive meander to absorb energy from the
            unused antenna polorization.

            **Note**:
            **There will still be some part of the structure that wont be used left
            on the mask**, namely the Hugh filters. These extra parts will have to
            be **removed manually**.

        return_configurator_points=False
            return a the points for use in the configurator.

        return_configurator_points_for_Lo_pass=False
            return a the points for use in the configurator.

        return_configurator_points_for_Hi_pass=False
            return a the points for use in the configurator.

        return_configurator_points_for_combiner_150ghz=False
            return a the points for use in the configurator.

        return_configurator_points_for_combiner_90ghz=False
            return a the points for use in the configurator.

        Returns
        -------
        points_to_conect_kids_dict: dict
            This dictionary contains keys that map to an [x,y] list which are the
            coordinates defining the conection points where the KIDs should connect
            to the filter bank.

            This dict has keys: **'TR', 'TL', 'BL', 'BR'**.
        """
        return smc.add_filter_bank_and_get_conection_points(
            self,
            x,
            y,
            filter_bank_config_override=filter_bank_config_override,
            filter_bank_ring_overlap_config_override=filter_bank_ring_overlap_config_override,
            Hi_pass_filters_config_override=Hi_pass_filters_config_override,
            Lo_pass_filters_config_override=Lo_pass_filters_config_override,
            combiner_section_90ghz_config_override=combiner_section_90ghz_config_override,
            combiner_section_150ghz_config_override=combiner_section_150ghz_config_override,
            with_combiner=with_combiner,
            with_crossover=with_crossover,
            only_1_pol=only_1_pol,
            return_configurator_points=return_configurator_points,
            return_configurator_points_for_Lo_pass=return_configurator_points_for_Lo_pass,
            return_configurator_points_for_Hi_pass=return_configurator_points_for_Hi_pass,
            return_configurator_points_for_combiner_150ghz=return_configurator_points_for_combiner_150ghz,
            return_configurator_points_for_combiner_90ghz=return_configurator_points_for_combiner_90ghz,
            **kwargs,
        )

    def connect_filter_bank_to_KIDs(
        self,
        x: float | int,
        y: float | int,
        absolute_filter_bank_points: dict[str, list[float | int]],
        relative_kid_positions: list[list[float | int]],
        only_1_pol_no_comb: bool = False,
        KID_conection_linewidth: float | int = 3,
        filter_conection_linewidth: float | int = 5,
        KID_extend_out_distance: float | int = 100.0,
        filter_extend_out_distance: float | int = 50.0,
        kid_side_bend_radius: float | int = 50,
        filt_side_bend_radius: float | int = 50,
    ):
        """Adds a microstrip feedline connecting the filter bank connection
        points (normally the end of the connection to the phase combiner) to
        the KIDs. There is a bend coming out of both the combiner on the filter
        bank and the KID and in the middle of the microstrip feedline is a
        taper if the two needed linewidths differ. The parameters for this
        microstrip line are defined at the top of this function.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate of the center of the antenna structure,
            i.e. the center of the horn.

        absolute_filter_bank_points: dict
            This dict should have keys: **'TR', 'TL', 'BL', 'BR'** that map to
            [x,y] lists which are the conection point coordinates for the filter
            bank.

        relative_kid_positions: list
            list of [x,y] lists defining the connection point coordinates for each
            of the KIDs. This is the very bottom center point of the inductive
            meander section. This list should be the coordinates, in order, of the
            TopLeft, TopRight, BotLeft, BotRight KIDs.

        KwArgs
        ------
        only_1_pol_no_comb: bool = False
            Default is False, when True this will alter the angles of the
            conections to the filter bank side of the conection to the KIDs.
            This is should be used when the function to add the filter bank,
            add_filter_bank_and_get_conection_points, has the arguments
            with_combiner=False, only_1_pol=True, with_crossover=False.

        KID_conection_linewidth: float | int = 3,
            The linewidth of the connection to the feedline.

        filter_conection_linewidth: float | int = 5,
            The linewidth of the connection to the filter bank.

        KID_extend_out_distance: float | int = 100.0,
            The distance over which the feedline is straight coming out of the
            KID conection.

        filter_extend_out_distance: float | int = 50.0,
            The distance over which the feedline is straight coming out of the
            filter bank conection.

        kid_side_bend_radius: float | int = 50,
            The bend radius of the path coming out of the KID conection.

        filt_side_bend_radius: float | int = 50,
            The bend radius of the path coming out of the filter bank conection.
        """

        # KID_conection_linewidth = 3
        # filter_conection_linewidth = 5
        # KID_extend_out_distance = 100.0
        # filter_extend_out_distance = 50.0
        # kid_side_bend_radius = 50
        # filt_side_bend_radius = 50

        Filter_TR = absolute_filter_bank_points["TR"]  # these are the coords of the end of the filters conection points
        Filter_TL = absolute_filter_bank_points["TL"]
        Filter_BL = absolute_filter_bank_points["BL"]
        Filter_BR = absolute_filter_bank_points["BR"]

        KID_TR = (
            relative_kid_positions[1][0] + x,
            relative_kid_positions[1][1] + y,
        )  # these are the coords of the base of each of the KIDS, TopLeft, TopRight, BotLeft, BotRight
        KID_TL = (relative_kid_positions[0][0] + x, relative_kid_positions[0][1] + y)
        KID_BL = (relative_kid_positions[2][0] + x, relative_kid_positions[2][1] + y)
        KID_BR = (relative_kid_positions[3][0] + x, relative_kid_positions[3][1] + y)

        FILT_coords = [Filter_TR, Filter_TL, Filter_BL, Filter_BR]
        KID_coords = [KID_TR, KID_TL, KID_BL, KID_BR]

        KID_extention_rots = [pi, 0, 0, pi]

        if not only_1_pol_no_comb:
            FILT_extention_rots = [-pi / 4, -3 * pi / 4, 3 * pi / 4, pi / 4]
        else:
            FILT_extention_rots = [0, pi / 2, 3 * pi / 4, pi / 4]

        for i in range(4):
            filt_coord = FILT_coords[i]
            kid_coord = KID_coords[i]

            filt_rot = FILT_extention_rots[i]
            kid_rot = KID_extention_rots[i]

            points_KID_to_taper = [
                kid_coord,
                (kid_coord[0] + KID_extend_out_distance * cos(kid_rot), kid_coord[1] + KID_extend_out_distance * sin(kid_rot)),
            ]  # defining the first points for the flexpath from the kid to the start of the taper

            points_taper_to_FILT = [
                [filt_coord[0] + filter_extend_out_distance * cos(filt_rot), filt_coord[1] + filter_extend_out_distance * sin(filt_rot)],
                filt_coord,
            ]  # defining the points for the end of the taper to the flex path coming out of the filter

            # this is the angle of the line conecting the two paths
            angle = np.arctan2(
                (points_taper_to_FILT[0][1] - points_KID_to_taper[-1][1]), (points_taper_to_FILT[0][0] - points_KID_to_taper[-1][0])
            )

            taper_kid_side_to_filt_side = [
                [
                    points_KID_to_taper[-1][0] + 2 * kid_side_bend_radius * cos(angle),
                    points_KID_to_taper[-1][1] + 2 * kid_side_bend_radius * sin(angle),
                ],
                [
                    points_taper_to_FILT[0][0] + 2 * filt_side_bend_radius * cos(angle + pi),
                    points_taper_to_FILT[0][1] + 2 * filt_side_bend_radius * sin(angle + pi),
                ],
            ]  # this calculates the middle ends points of the taper

            taper_points = [
                [
                    taper_kid_side_to_filt_side[0][0] + (KID_conection_linewidth / 2) * cos(angle - pi / 2),
                    taper_kid_side_to_filt_side[0][1] + (KID_conection_linewidth / 2) * sin(angle - pi / 2),
                ],
                [
                    taper_kid_side_to_filt_side[0][0] - (KID_conection_linewidth / 2) * cos(angle - pi / 2),
                    taper_kid_side_to_filt_side[0][1] - (KID_conection_linewidth / 2) * sin(angle - pi / 2),
                ],
                [
                    taper_kid_side_to_filt_side[1][0] - (filter_conection_linewidth / 2) * cos(angle - pi / 2),
                    taper_kid_side_to_filt_side[1][1] - (filter_conection_linewidth / 2) * sin(angle - pi / 2),
                ],
                [
                    taper_kid_side_to_filt_side[1][0] + (filter_conection_linewidth / 2) * cos(angle - pi / 2),
                    taper_kid_side_to_filt_side[1][1] + (filter_conection_linewidth / 2) * sin(angle - pi / 2),
                ],
            ]  # this calculates the geometry points for the trapezoidal taper

            taper = gdspy.Polygon(
                taper_points,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )  # makes the taper as a polygon

            points_KID_to_taper.append(
                taper_kid_side_to_filt_side[0]
            )  # appends the start of the taper to the flex path points from the kid so the taper conects the two paths

            points_taper_to_FILT.insert(
                0, taper_kid_side_to_filt_side[1]
            )  # inserts the end of the taper at the start of the flex path points for the filter side to again ensure it conects up

            # creates the flex paths for each side
            KID_to_taper = gdspy.FlexPath(
                points_KID_to_taper,
                KID_conection_linewidth,
                corners="circular bend",
                bend_radius=kid_side_bend_radius,
                gdsii_path=True,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            taper_to_FILT = gdspy.FlexPath(
                points_taper_to_FILT,
                filter_conection_linewidth,
                corners="circular bend",
                bend_radius=filt_side_bend_radius,
                gdsii_path=True,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )

            # adds the taper and the two flex paths for each side to the main cell
            self.Main.add(taper)
            self.make_flexpath_into_polygons_and_add_to_main(KID_to_taper, self.layers.Nb_Antenna)
            self.make_flexpath_into_polygons_and_add_to_main(taper_to_FILT, self.layers.Nb_Antenna)
            # self.Main.add(KID_to_taper)
            # self.Main.add(taper_to_FILT)

        return

    def connect_ants_to_filter_bank(
        self,
        x: float | int,
        y: float | int,
        relative_antena_conect_positions: list[list[float | int]],
        antena_rot: float | int,
        filter_bank_config_override: dict[str, float | int] | None = None,
        antenna_config_override: dict[str, float | int] | None = None,
        antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None,
        terminate_ants: list[str] = [],
        add_dielectric_under_conections: bool = True,
        return_configurator_points: bool = False,
        **kwargs,
    ):
        """Adds a co-planar waveguide feedline that transitons to a microstrip
        that connects the antenna pads to the inner ring of the filter bank.
        The co-planar waveguide is on the antenna pad side and the microstrip
        is on the filter bank side. The parameters of this geometry is within
        the config.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate of the center of the antennas,
            i.e. the center of the horn.

        relative_antena_conect_positions: list
            This list contains [x,y] lists which are the coordinates of the
            conection points of the tip of the antenna pads. These coords should be
            relative to the center of the antennas. The order of the antenna
            connections should be [bot, top, left, right], where each element is
            an [x, y].

        antena_rot: float, int
            The rotation angle (**in degrees**) of the antennas about thier center.

        KwArgs
        ------
        filter_bank_config_override,
        antenna_config_override,
        antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None,
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.


        terminate_ants: list of strings,
            Default empty list. This is a list of all the antennas that should
            be terminated by an aluminium meander. This should take a list of
            any of the following strings "L", "R", "T", "B". as in left
            e.g. terminate_ants=["L", "R"].

        add_dielectric_under_conections = True
            Default is True which adds a dielectric strip under the connection
            from the antenna to the filter bank.

        return_configurator_points=False
            return a the points for use in the configurator.
        """
        config = get_mask_default_config(
            SoukMaskConfig.ANTENNA_CPW_MICROSTRIP_TRANS,
            config_override=antenna_cpw_microstrip_trans_config_override,
        )

        antenna_config = get_mask_default_config(
            SoukMaskConfig.ANTENNA,
            config_override=antenna_config_override,
        )
        filter_bank_config = get_mask_default_config(
            SoukMaskConfig.FILTER_BANK,
            config_override=filter_bank_config_override,
        )

        inner_ring_radius = filter_bank_config["inner_ring_radius"]
        ant_line_width = antenna_config["top_conect_width"]

        microstrip_lw = config["microstrip_lw"]
        CPW_width = config["CPW_width"]

        CPW_tans_lens = [
            config["CPW_tans_len1"],
            config["CPW_tans_len2"],
            config["CPW_tans_len3"],
            config["CPW_tans_len4"],
            config["CPW_tans_len5"],
        ]
        gaps = [config["gap1"], config["gap2"], config["gap3"], config["gap4"], config["gap5"]]

        tot_trans_length = np.sum(CPW_tans_lens) + np.sum(gaps)

        trans1 = gdspy.Rectangle(
            [gaps[0], -CPW_width / 2],
            [gaps[0] + CPW_tans_lens[0], CPW_width / 2],
        )
        trans2 = gdspy.Rectangle(
            [np.sum(gaps[0:1]) + np.sum(CPW_tans_lens[0:1]) + gaps[1], -CPW_width / 2],
            [np.sum(gaps[0:1]) + np.sum(CPW_tans_lens[0:1]) + gaps[1] + CPW_tans_lens[1], CPW_width / 2],
        )
        trans3 = gdspy.Rectangle(
            [np.sum(gaps[0:2]) + np.sum(CPW_tans_lens[0:2]) + gaps[2], -CPW_width / 2],
            [np.sum(gaps[0:2]) + np.sum(CPW_tans_lens[0:2]) + gaps[2] + CPW_tans_lens[2], CPW_width / 2],
        )
        trans4 = gdspy.Rectangle(
            [np.sum(gaps[0:3]) + np.sum(CPW_tans_lens[0:3]) + gaps[3], -CPW_width / 2],
            [np.sum(gaps[0:3]) + np.sum(CPW_tans_lens[0:3]) + gaps[3] + CPW_tans_lens[3], CPW_width / 2],
        )
        trans5 = gdspy.Rectangle(
            [np.sum(gaps[0:4]) + np.sum(CPW_tans_lens[0:4]) + gaps[4], -CPW_width / 2],
            [np.sum(gaps[0:4]) + np.sum(CPW_tans_lens[0:4]) + gaps[4] + CPW_tans_lens[4], CPW_width / 2],
        )

        transition = gdspy.Cell("CPWTransitionForAntennas")
        transition.add([trans1, trans2, trans3, trans4, trans5])

        top_points = []
        bot_points = []
        for i in range(len(gaps)):
            top_points.append([np.sum(gaps[0:i]) + np.sum(CPW_tans_lens[0:i]), microstrip_lw / 2])
            top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), microstrip_lw / 2])
            top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), ant_line_width / 2])
            top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0 : i + 1]), ant_line_width / 2])

            bot_points.append([np.sum(gaps[0:i]) + np.sum(CPW_tans_lens[0:i]), -microstrip_lw / 2])
            bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), -microstrip_lw / 2])
            bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), -ant_line_width / 2])
            bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0 : i + 1]), -ant_line_width / 2])

        cpw_transition_poly_points = top_points
        cpw_transition_poly_points.extend(bot_points[::-1])

        ANT_L = [
            relative_antena_conect_positions[2][0] + x,
            relative_antena_conect_positions[2][1] + y,
        ]  # these are the coords of the end of the antennas , Left, Top, Right, Bot
        ANT_T = [relative_antena_conect_positions[1][0] + x, relative_antena_conect_positions[1][1] + y]
        ANT_R = [relative_antena_conect_positions[3][0] + x, relative_antena_conect_positions[3][1] + y]
        ANT_B = [relative_antena_conect_positions[0][0] + x, relative_antena_conect_positions[0][1] + y]

        ANTS = [ANT_R, ANT_T, ANT_L, ANT_B]

        FILT_TR = [
            x + inner_ring_radius * cos(pi / 4 + 0 * (pi / 2)),
            y + inner_ring_radius * sin(pi / 4 + 0 * (pi / 2)),
        ]  # these are the coords of the filter bank conect positions, TopRight, TopLeft, BotLeft, BotRight
        FILT_TL = [x + inner_ring_radius * cos(pi / 4 + 1 * (pi / 2)), y + inner_ring_radius * sin(pi / 4 + 1 * (pi / 2))]
        FILT_BL = [x + inner_ring_radius * cos(pi / 4 + 2 * (pi / 2)), y + inner_ring_radius * sin(pi / 4 + 2 * (pi / 2))]
        FILT_BR = [x + inner_ring_radius * cos(pi / 4 + 3 * (pi / 2)), y + inner_ring_radius * sin(pi / 4 + 3 * (pi / 2))]

        FILTS = [FILT_BR, FILT_TR, FILT_TL, FILT_BL]

        extend_out_distance = config[
            "extend_out_distance"
        ]  # the distance ove which the feedline is straight coming out of the antenna and Filter Bank
        flex_path_bend_radius = config["flex_path_bend_radius"]

        add_termination_meander_indexes = []
        if terminate_ants != []:
            if "R" in terminate_ants:
                add_termination_meander_indexes.append(0)
            if "T" in terminate_ants:
                add_termination_meander_indexes.append(1)
            if "L" in terminate_ants:
                add_termination_meander_indexes.append(2)
            if "B" in terminate_ants:
                add_termination_meander_indexes.append(3)

        for i in range(4):
            ANT = ANTS[i]
            FILT = FILTS[i]

            extended_ANT = [
                ANT[0] + extend_out_distance * cos(mbu.deg_to_rad(antena_rot + (i * 90))),
                ANT[1] + extend_out_distance * sin(mbu.deg_to_rad(antena_rot + (i * 90))),
            ]

            extended_FILT = [
                FILT[0] + extend_out_distance * cos(3 * pi / 4 + i * (pi / 2)),
                FILT[1] + extend_out_distance * sin(3 * pi / 4 + i * (pi / 2)),
            ]

            middle_straight_dx = extended_FILT[0] - extended_ANT[0]
            middle_straight_dy = extended_FILT[1] - extended_ANT[1]

            middle_straight_angle = np.arctan2(middle_straight_dy, middle_straight_dx)

            tenth_pnt_straight = [extended_ANT[0] + (middle_straight_dx / 10), extended_ANT[1] + (middle_straight_dy / 10)]

            trans_end = [
                tenth_pnt_straight[0] + tot_trans_length * cos(middle_straight_angle),
                tenth_pnt_straight[1] + tot_trans_length * sin(middle_straight_angle),
            ]

            ANT_to_TRANS_path = gdspy.FlexPath(
                [ANT, extended_ANT, tenth_pnt_straight],
                ant_line_width,
                corners="circular bend",
                bend_radius=flex_path_bend_radius,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            self.make_flexpath_into_polygons_and_add_to_main(ANT_to_TRANS_path, self.layers.Nb_Antenna)

            cutout_around_ANT_to_TANS_polys = mbu.get_polys_from_flexpath(
                gdspy.FlexPath(
                    [ANT, extended_ANT, tenth_pnt_straight],
                    CPW_width,
                    corners="circular bend",
                    bend_radius=flex_path_bend_radius,
                    gdsii_path=True,
                )
            )
            for i in range(len(cutout_around_ANT_to_TANS_polys)):
                self.ground_plane_cutouts.add(gdspy.Polygon(cutout_around_ANT_to_TANS_polys[i]))

            TRANS_to_FILT_path = gdspy.FlexPath(
                [trans_end, extended_FILT, FILT],
                microstrip_lw,
                corners="circular bend",
                bend_radius=flex_path_bend_radius,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            self.make_flexpath_into_polygons_and_add_to_main(TRANS_to_FILT_path, self.layers.Nb_Antenna)

            TRANS_cutouts = gdspy.CellReference(
                transition, tenth_pnt_straight, rotation=(mbu.rad_to_deg(middle_straight_angle))
            )  # adds the transition geometry between CPW and microstrip antenna a fith of the way down the middle straight
            self.ground_plane_cutouts.add(TRANS_cutouts)

            cpw_transition_center_points = mbu.rotate_and_move_points_list(
                cpw_transition_poly_points, middle_straight_angle, tenth_pnt_straight[0], tenth_pnt_straight[1]
            )
            cpw_transition_center = gdspy.Polygon(
                cpw_transition_center_points,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            self.Main.add(cpw_transition_center)

            if add_dielectric_under_conections:
                dielectric_path_points = [ANT, extended_ANT, extended_FILT, FILT]
                dielectric_under_connection_width = 5 * CPW_width
                dielectric_under_connection_path = gdspy.FlexPath(
                    dielectric_path_points,
                    dielectric_under_connection_width,
                    corners="circular bend",
                    bend_radius=flex_path_bend_radius,
                    layer=self.layers.SiN_dep.number,
                    datatype=self.layers.SiN_dep.datatype,
                )
                self.make_flexpath_into_polygons_and_add_to_main(dielectric_under_connection_path, self.layers.SiN_dep)

        mander_length = 2514
        meander_width = 3

        overlap_pad_width = 10
        overlap_pad_height = 10

        meander_bend_height = 30

        mander_bend_length = pi * meander_bend_height / 2

        for i in add_termination_meander_indexes:
            start_ang = -pi / 4 + i * (pi / 2)
            start_xy = [x + inner_ring_radius * cos(start_ang), y + inner_ring_radius * sin(start_ang)]

            arc_len_close = (mander_length - mander_bend_length) / 3
            arc_len_far = (mander_length - mander_bend_length) * 2 / 3

            arc_theta_close = arc_len_close / inner_ring_radius
            arc_theta_far = arc_len_far / (inner_ring_radius + meander_bend_height)

            bend_center_xy = [
                x + (inner_ring_radius + (meander_bend_height / 2)) * cos(start_ang - arc_theta_close),
                y + (inner_ring_radius + (meander_bend_height / 2)) * sin(start_ang - arc_theta_close),
            ]

            arc_close = gdspy.Round(
                [x, y],
                radius=(inner_ring_radius + (meander_width / 2)),
                inner_radius=(inner_ring_radius - (meander_width / 2)),
                initial_angle=(start_ang),
                final_angle=(start_ang - arc_theta_close),
                layer=self.layers.Aluminium_Direct.number,
                datatype=self.layers.Aluminium_Direct.datatype,
            )

            arc_bend = gdspy.Round(
                bend_center_xy,
                radius=((meander_bend_height / 2) + (meander_width / 2)),
                inner_radius=((meander_bend_height / 2) - (meander_width / 2)),
                initial_angle=(start_ang - arc_theta_close - pi),
                final_angle=(start_ang - arc_theta_close),
                layer=self.layers.Aluminium_Direct.number,
                datatype=self.layers.Aluminium_Direct.datatype,
            )

            arc_far = gdspy.Round(
                [x, y],
                radius=(inner_ring_radius + meander_bend_height + (meander_width / 2)),
                inner_radius=(inner_ring_radius + meander_bend_height - (meander_width / 2)),
                initial_angle=(start_ang - arc_theta_close + arc_theta_far),
                final_angle=(start_ang - arc_theta_close),
                layer=self.layers.Aluminium_Direct.number,
                datatype=self.layers.Aluminium_Direct.datatype,
            )

            overlap_ang = (overlap_pad_width / 2) / inner_ring_radius
            conect_overlap = gdspy.Round(
                [x, y],
                radius=(inner_ring_radius + (overlap_pad_height / 2)),
                inner_radius=(inner_ring_radius - (overlap_pad_height / 2)),
                initial_angle=(start_ang - overlap_ang),
                final_angle=(start_ang + overlap_ang),
                layer=self.layers.Aluminium_Direct.number,
                datatype=self.layers.Aluminium_Direct.datatype,
            )

            self.Main.add(arc_close)
            self.Main.add(arc_bend)
            self.Main.add(arc_far)
            self.Main.add(conect_overlap)

        if not return_configurator_points:
            return

        configurator_points = {}

        configurator_points["gap1"] = {
            "text": "gap1",
            "start": [cpw_transition_center_points[0][0], cpw_transition_center_points[0][1]],
            "end": [cpw_transition_center_points[1][0], cpw_transition_center_points[1][1]],
        }

        configurator_points["CPW_tans_len1"] = {
            "text": "CPW_tans_len1",
            "start": [cpw_transition_center_points[2][0], cpw_transition_center_points[2][1]],
            "end": [cpw_transition_center_points[3][0], cpw_transition_center_points[3][1]],
        }

        configurator_points["gap2"] = {
            "text": "gap2",
            "start": [cpw_transition_center_points[4][0], cpw_transition_center_points[4][1]],
            "end": [cpw_transition_center_points[5][0], cpw_transition_center_points[5][1]],
        }

        configurator_points["CPW_tans_len2"] = {
            "text": "CPW_tans_len2",
            "start": [cpw_transition_center_points[6][0], cpw_transition_center_points[6][1]],
            "end": [cpw_transition_center_points[7][0], cpw_transition_center_points[7][1]],
        }

        configurator_points["gap3"] = {
            "text": "gap3",
            "start": [cpw_transition_center_points[8][0], cpw_transition_center_points[8][1]],
            "end": [cpw_transition_center_points[9][0], cpw_transition_center_points[9][1]],
        }

        configurator_points["CPW_tans_len3"] = {
            "text": "CPW_tans_len3",
            "start": [cpw_transition_center_points[10][0], cpw_transition_center_points[10][1]],
            "end": [cpw_transition_center_points[11][0], cpw_transition_center_points[11][1]],
        }

        configurator_points["gap4"] = {
            "text": "gap4",
            "start": [cpw_transition_center_points[12][0], cpw_transition_center_points[12][1]],
            "end": [cpw_transition_center_points[13][0], cpw_transition_center_points[13][1]],
        }

        configurator_points["CPW_tans_len4"] = {
            "text": "CPW_tans_len4",
            "start": [cpw_transition_center_points[14][0], cpw_transition_center_points[14][1]],
            "end": [cpw_transition_center_points[15][0], cpw_transition_center_points[15][1]],
        }

        configurator_points["gap5"] = {
            "text": "gap5",
            "start": [cpw_transition_center_points[16][0], cpw_transition_center_points[16][1]],
            "end": [cpw_transition_center_points[17][0], cpw_transition_center_points[17][1]],
        }

        configurator_points["CPW_tans_len5"] = {
            "text": "CPW_tans_len5",
            "start": [cpw_transition_center_points[18][0], cpw_transition_center_points[18][1]],
            "end": [cpw_transition_center_points[19][0], cpw_transition_center_points[19][1]],
        }

        configurator_points["microstrip_lw"] = {
            "text": "microstrip_lw",
            "start": [cpw_transition_center_points[-1][0], cpw_transition_center_points[-1][1]],
            "end": [cpw_transition_center_points[0][0], cpw_transition_center_points[0][1]],
        }

        poins_for_last_trans_cutout = TRANS_cutouts.get_polygons()[0]
        configurator_points["CPW_width"] = {
            "text": "CPW_width",
            "start": [poins_for_last_trans_cutout[1][0], poins_for_last_trans_cutout[1][1]],
            "end": [poins_for_last_trans_cutout[0][0], poins_for_last_trans_cutout[0][1]],
        }

        configurator_points["extend_out_distance"] = {
            "text": "extend_out_distance",
            "start": [ANT[0], ANT[1]],
            "end": [extended_ANT[0], extended_ANT[1]],
        }

        extention_angle = mbu.deg_to_rad(antena_rot + 1 * 90) % (2 * pi)
        mid_of_angle = (middle_straight_angle % (2 * pi)) - (extention_angle)

        configurator_points["flex_path_bend_radius"] = {
            "text": "flex_path_bend_radius",
            "start": [
                extended_ANT[0] + flex_path_bend_radius * cos(mid_of_angle),
                extended_ANT[1] + flex_path_bend_radius * sin(mid_of_angle),
            ],
            "end": [
                extended_ANT[0] + CPW_width * cos(mid_of_angle),
                extended_ANT[1] + CPW_width * sin(mid_of_angle),
            ],
        }

        return configurator_points

    def connect_ants_to_KIDs(
        self,
        x: float | int,
        y: float | int,
        relative_antena_conect_positions: list[list[float | int]],
        relative_kid_positions: list[list[float | int]],
        antena_rot: float | int,
        antenna_config_override: dict[str, float | int] | None = None,
        antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None,
        add_dielectric_under_conections: bool = True,
        **kwargs,
    ) -> None:
        """Adds a co-planar waveguide feedline that transitons to a microstrip
        that connects the antenna pads to the meanders of the KIDs. The co-
        planar waveguide is on the antenna pad side and the microstrip is on
        the KID side. The parameters of this geometry is within the config.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate of the center of the antennas,
            i.e. the center of the horn.

        relative_antena_conect_positions: list
            This list contains [x,y] lists which are the coordinates of the
            conection points of the tip of the antenna pads. These coords should be
            relative to the center of the antennas. The order of the antenna
            connections should be [bot, top, left, right], where each element is
            an [x, y].

        relative_kid_positions: list
            list of [x,y] lists defining the connection point coordinates for each
            of the KIDs. This is the very bottom center point of the inductive
            meander section. This list should be the coordinates, in order, of the
            TopLeft, TopRight, BotLeft, BotRight KIDs.

        antena_rot: float, int
            The rotation angle (**in degrees**) of the antennas about thier center.

        KwArgs
        ------
        antenna_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        antenna_cpw_microstrip_trans_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        add_dielectric_under_conections = True
            Adds dielectric under the connection from the antenna to the KID
            meander.
        """
        antenna_config = get_mask_default_config(SoukMaskConfig.ANTENNA, config_override=antenna_config_override)
        antenna_cpw_microstrip_trans_config = get_mask_default_config(
            SoukMaskConfig.ANTENNA_CPW_MICROSTRIP_TRANS, config_override=antenna_cpw_microstrip_trans_config_override
        )

        ant_line_width = antenna_config["top_conect_width"]
        microstrip_lw = antenna_cpw_microstrip_trans_config["microstrip_lw"]
        CPW_width = antenna_cpw_microstrip_trans_config["CPW_width"]

        CPW_tans_lens = [
            antenna_cpw_microstrip_trans_config["CPW_tans_len1"],
            antenna_cpw_microstrip_trans_config["CPW_tans_len2"],
            antenna_cpw_microstrip_trans_config["CPW_tans_len3"],
            antenna_cpw_microstrip_trans_config["CPW_tans_len4"],
            antenna_cpw_microstrip_trans_config["CPW_tans_len5"],
        ]
        gaps = [
            antenna_cpw_microstrip_trans_config["gap1"],
            antenna_cpw_microstrip_trans_config["gap2"],
            antenna_cpw_microstrip_trans_config["gap3"],
            antenna_cpw_microstrip_trans_config["gap4"],
            antenna_cpw_microstrip_trans_config["gap5"],
        ]

        tot_trans_length = np.sum(CPW_tans_lens) + np.sum(gaps)

        trans1 = gdspy.Rectangle([gaps[0], -CPW_width / 2], [gaps[0] + CPW_tans_lens[0], CPW_width / 2])
        trans2 = gdspy.Rectangle(
            [np.sum(gaps[0:1]) + np.sum(CPW_tans_lens[0:1]) + gaps[1], -CPW_width / 2],
            [np.sum(gaps[0:1]) + np.sum(CPW_tans_lens[0:1]) + gaps[1] + CPW_tans_lens[1], CPW_width / 2],
        )
        trans3 = gdspy.Rectangle(
            [np.sum(gaps[0:2]) + np.sum(CPW_tans_lens[0:2]) + gaps[2], -CPW_width / 2],
            [np.sum(gaps[0:2]) + np.sum(CPW_tans_lens[0:2]) + gaps[2] + CPW_tans_lens[2], CPW_width / 2],
        )
        trans4 = gdspy.Rectangle(
            [np.sum(gaps[0:3]) + np.sum(CPW_tans_lens[0:3]) + gaps[3], -CPW_width / 2],
            [np.sum(gaps[0:3]) + np.sum(CPW_tans_lens[0:3]) + gaps[3] + CPW_tans_lens[3], CPW_width / 2],
        )
        trans5 = gdspy.Rectangle(
            [np.sum(gaps[0:4]) + np.sum(CPW_tans_lens[0:4]) + gaps[4], -CPW_width / 2],
            [np.sum(gaps[0:4]) + np.sum(CPW_tans_lens[0:4]) + gaps[4] + CPW_tans_lens[4], CPW_width / 2],
        )

        transition = gdspy.Cell("CPWTransitionForAntennas")
        transition.add([trans1, trans2, trans3, trans4, trans5])

        top_points = []
        bot_points = []
        for i in range(len(gaps)):
            top_points.append([np.sum(gaps[0:i]) + np.sum(CPW_tans_lens[0:i]), microstrip_lw / 2])
            top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), microstrip_lw / 2])
            top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), ant_line_width / 2])
            top_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0 : i + 1]), ant_line_width / 2])

            bot_points.append([np.sum(gaps[0:i]) + np.sum(CPW_tans_lens[0:i]), -microstrip_lw / 2])
            bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), -microstrip_lw / 2])
            bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0:i]), -ant_line_width / 2])
            bot_points.append([np.sum(gaps[0 : i + 1]) + np.sum(CPW_tans_lens[0 : i + 1]), -ant_line_width / 2])

        cpw_transition_poly_points = top_points
        cpw_transition_poly_points.extend(bot_points[::-1])

        ANT_L = [
            relative_antena_conect_positions[2][0] + x,
            relative_antena_conect_positions[2][1] + y,
        ]  # these are the coords of the end of the antennas , Left, Top, Right, Bot
        ANT_T = [relative_antena_conect_positions[1][0] + x, relative_antena_conect_positions[1][1] + y]
        ANT_R = [relative_antena_conect_positions[3][0] + x, relative_antena_conect_positions[3][1] + y]
        ANT_B = [relative_antena_conect_positions[0][0] + x, relative_antena_conect_positions[0][1] + y]

        ANTS = [ANT_R, ANT_T, ANT_L, ANT_B]

        KID_TL = [x + relative_kid_positions[0][0], y + relative_kid_positions[0][1]]
        KID_TR = [x + relative_kid_positions[1][0], y + relative_kid_positions[1][1]]
        KID_BL = [x + relative_kid_positions[2][0], y + relative_kid_positions[2][1]]
        KID_BR = [x + relative_kid_positions[3][0], y + relative_kid_positions[3][1]]

        KIDS = [KID_BR, KID_TR, KID_TL, KID_BL]
        KIDS_ROTS = [pi, pi, 0, 0]

        extend_out_distance = antenna_cpw_microstrip_trans_config[
            "extend_out_distance"
        ]  # the distance ove which the feedline is straight coming out of the antenna and Filter Bank
        flex_path_bend_radius = antenna_cpw_microstrip_trans_config["flex_path_bend_radius"]

        for i in range(4):
            ANT = ANTS[i]
            KID = KIDS[i]

            extended_ANT = [
                ANT[0] + extend_out_distance * cos(mbu.deg_to_rad(antena_rot + (i * 90))),
                ANT[1] + extend_out_distance * sin(mbu.deg_to_rad(antena_rot + (i * 90))),
            ]

            extended_KID = [KID[0] + extend_out_distance * cos(KIDS_ROTS[i]), KID[1] + extend_out_distance * sin(KIDS_ROTS[i])]

            middle_straight_dx = extended_KID[0] - extended_ANT[0]
            middle_straight_dy = extended_KID[1] - extended_ANT[1]

            middle_straight_angle = np.arctan2(middle_straight_dy, middle_straight_dx)

            tenth_pnt_straight = [extended_ANT[0] + (middle_straight_dx / 10), extended_ANT[1] + (middle_straight_dy / 10)]

            trans_end = [
                tenth_pnt_straight[0] + tot_trans_length * cos(middle_straight_angle),
                tenth_pnt_straight[1] + tot_trans_length * sin(middle_straight_angle),
            ]

            ANT_to_TRANS_path = gdspy.FlexPath(
                [ANT, extended_ANT, tenth_pnt_straight],
                ant_line_width,
                corners="circular bend",
                bend_radius=flex_path_bend_radius,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            # self.Main.add(ANT_to_TRANS_path)
            self.make_flexpath_into_polygons_and_add_to_main(ANT_to_TRANS_path, self.layers.Nb_Antenna)

            cutout_around_ANT_to_TANS_polys = mbu.get_polys_from_flexpath(
                gdspy.FlexPath(
                    [ANT, extended_ANT, tenth_pnt_straight],
                    CPW_width,
                    corners="circular bend",
                    bend_radius=flex_path_bend_radius,
                    gdsii_path=True,
                )
            )
            for i in range(len(cutout_around_ANT_to_TANS_polys)):
                self.ground_plane_cutouts.add(gdspy.Polygon(cutout_around_ANT_to_TANS_polys[i]))

            TRANS_to_KID_path = gdspy.FlexPath(
                [trans_end, extended_KID, KID],
                microstrip_lw,
                corners="circular bend",
                bend_radius=flex_path_bend_radius,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            # self.Main.add(TRANS_to_KID_path)
            self.make_flexpath_into_polygons_and_add_to_main(TRANS_to_KID_path, self.layers.Nb_Antenna)

            TRANS_cutouts = gdspy.CellReference(
                transition, tenth_pnt_straight, rotation=(mbu.rad_to_deg(middle_straight_angle))
            )  # adds the transition geometry between CPW and microstrip antenna a fith of the way down the middle straight
            self.ground_plane_cutouts.add(TRANS_cutouts)

            points = mbu.rotate_and_move_points_list(
                cpw_transition_poly_points, middle_straight_angle, tenth_pnt_straight[0], tenth_pnt_straight[1]
            )
            cpw_transition_center = gdspy.Polygon(
                points,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            self.Main.add(cpw_transition_center)

            if add_dielectric_under_conections:
                dielectric_path_points = [ANT, extended_ANT, extended_KID, KID]
                dielectric_under_connection_width = 5 * CPW_width
                dielectric_under_connection_path = gdspy.FlexPath(
                    dielectric_path_points,
                    dielectric_under_connection_width,
                    corners="circular bend",
                    bend_radius=flex_path_bend_radius,
                    layer=self.layers.SiN_dep.number,
                    datatype=self.layers.SiN_dep.datatype,
                )
                self.make_flexpath_into_polygons_and_add_to_main(dielectric_under_connection_path, self.layers.SiN_dep)

        return

    def get_rounded_path_from_passthough_points(self, feedline_pass_through_points, bend_radius, bend_points_gap=10):
        """Generates a list of [x,y] points of a path with rounded corners from
        an input list of [x,y] points of a path. Rounds corners with given bend
        radius.

        The bend is put before the corner point. e.g.
        bending the path [[0,1], [1,1], [1,0]] with bend radius of 1 will form a
        circular arc from [0,1] to [1,0] like the first quadrent of a circle with
        radius 1.

        Parameters
        ----------
        feedline_pass_through_points: list
            list of [x,y] lists which are the coordinates that the feedline passes
            through. **Note**, these points will not be neccessarily be included in
            the final returned rounded feedline.

        bend_radius: float, int
            radius for the bends to be created smoothing the feedline out.

        KwArgs
        ------
        bend_points_gap = 10
            This is the distance (in microns) between points in the rounded section
            of the line. This value can be made smaller or larger to make the
            rounded section smoother or more polygonal respectively. The default
            10um value however is smooth enough without creating a huge number of
            points.

        Returns
        -------
        rounded_feedline_points: list
            list of [x,y] lists that define coordinates that the rounded feedline
            passes through. Similar in form to the input list
            ("feedline_pass_through_points") however this new list will contain a
            lot more points. The
        """
        rounded_feedline_points = []
        rounded_feedline_points.append(feedline_pass_through_points[0])

        R = bend_radius

        for i in range(len(feedline_pass_through_points) - 2):
            x1 = feedline_pass_through_points[i][0]
            y1 = feedline_pass_through_points[i][1]

            x2 = feedline_pass_through_points[i + 1][0]
            y2 = feedline_pass_through_points[i + 1][1]

            x3 = feedline_pass_through_points[i + 2][0]
            y3 = feedline_pass_through_points[i + 2][1]

            len1 = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            len2 = np.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)

            v1x = (x2 - x1) / len1
            v1y = (y2 - y1) / len1

            v2x = (x3 - x2) / len2
            v2y = (y3 - y2) / len2

            if v1x * v2y - v1y * v2x > 0.0:
                bend = "left"
            else:
                bend = "right"

            if bend == "left":
                px1 = x1 - v1y * R
                py1 = y1 + v1x * R
                px2 = x2 - v2y * R
                py2 = y2 + v2x * R

            if bend == "right":
                px1 = x1 + v1y * R
                py1 = y1 - v1x * R
                px2 = x2 + v2y * R
                py2 = y2 - v2x * R

            den = v1x * v2y - v2x * v1y

            k1 = (v2y * (px2 - px1) - v2x * (py2 - py1)) / den
            k2 = (v1y * (px2 - px1) - v1x * (py2 - py1)) / den

            cx = px1 + k1 * v1x
            cy = py1 + k1 * v1y

            tx1 = x1 + k1 * v1x
            ty1 = y1 + k1 * v1y
            tx2 = x2 + k2 * v2x
            ty2 = y2 + k2 * v2y

            ang_inc = bend_points_gap / R

            init_angle = np.arctan2((ty1 - cy), (tx1 - cx))

            init_ang_unit_vec = [(tx1 - cx), (ty1 - cy)] / np.linalg.norm([(tx1 - cx), (ty1 - cy)])
            final_ang_unit_vec = [(tx2 - cx), (ty2 - cy)] / np.linalg.norm([(tx2 - cx), (ty2 - cy)])

            ang_diff = np.arccos(np.dot(init_ang_unit_vec, final_ang_unit_vec))

            init_angle = np.arctan2(init_ang_unit_vec[1], init_ang_unit_vec[0])

            loop_ang = init_angle

            for j in range(int(ang_diff / ang_inc)):
                bx = cx + R * cos(loop_ang)
                by = cy + R * sin(loop_ang)
                rounded_feedline_points.append([bx, by])

                if bend == "left":
                    loop_ang += ang_inc
                else:
                    loop_ang -= ang_inc

        rounded_feedline_points.append([feedline_pass_through_points[-1][0], feedline_pass_through_points[-1][1]])

        return rounded_feedline_points

    def add_feedline_and_dielectric(
        self,
        feed_points: list[list[float | int]],
        cpw_feedline_config_override: dict[str, float | int] | None = None,
        center_material: Layer | list[Layer] | None = None,
        add_bridges: bool = True,
        points_to_avoid_bridging: list[list[float | int]] | None = None,
        **kwargs,
    ):
        """Adds a Co-Planar Waveguide feedline through the points given along
        with the dielectric layer and adds brdiges across every so many units
        by default. The parameters determining the feeline geometry are the
        config. The CPW consists of the center line, cut out in the ground
        plane around this center line and a dielectric layer covering that.

        Parameters
        ----------
        feed_points: list[list[float | int]]
            list of [x,y] lists which are the coordinates for the feedline to pass
            through.

        KwArgs
        ----------
        cpw_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        center_material: Layer | list[Layer] | None = None
            This is an instance of Layer or list of two. see maskpy.layers.Layer.
            Usually this is within the SoukMaskBuilder.layers.xxx.
            e.g. `self.layers.Aluminium`.
            By default the value is None which will set the center material to
            Nb_Antenna.

        add_bridges: bool = True
            Adds bridges across the CPW by deafult, if set to False it will just
            add the feedline.

        points_to_avoid_bridging: list[list[float | int]] | None = None
            If not defined then the bridges will be added at the set distance
            if the add_bridges arg is True as by default.
            When this is defined this should be a list of [x,y] lists defining
            all the points in which bridges should not be placed. This
            typically would be a list of all the coordinates where the KIDs
            connect to the feedline. If a bridge would be placed within a set
            distance of one of these points it will try and shift the bridge in
            the positive x direction a distance defined in a variable below
            until it no longer is considered clashing.
        """
        config = get_mask_default_config(SoukMaskConfig.CPW_FEEDLINE, config_override=cpw_feedline_config_override)

        bend_radius = config["bend_radius"]
        feedline_width = config["feedline_width"]
        cutout_around_feedline_width = config["cutout_around_feedline_width"]
        dielectric_under_feedline_width = config["dielectric_under_feedline_width"]

        if center_material is None:
            center_material = self.layers.Nb_Antenna

        # makes and adds the center of the cpw flex path in the material specified
        if isinstance(center_material, Layer):
            center_path = gdspy.FlexPath(
                feed_points,
                feedline_width,
                corners="circular bend",
                bend_radius=bend_radius,
                gdsii_path=True,
                layer=center_material.number,
                datatype=center_material.datatype,
            )
            self.make_flexpath_into_polygons_and_add_to_main(center_path, self.layers.Nb_Antenna)
        elif isinstance(center_material, list):
            if all(isinstance(mat, Layer) for mat in center_material):
                for material in center_material:
                    center_path = gdspy.FlexPath(
                        feed_points,
                        feedline_width,
                        corners="circular bend",
                        bend_radius=bend_radius,
                        gdsii_path=True,
                        layer=material.number,
                        datatype=material.datatype,
                    )
                    self.make_flexpath_into_polygons_and_add_to_main(center_path, self.layers.Nb_Antenna)
            else:
                raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")
        else:
            raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")

        # # makes and adds the center of the cpw flex path in the material specified
        # if center_material == "Nb":
        #     center_path = gdspy.FlexPath(
        #         feed_points,
        #         feedline_width,
        #         corners="circular bend",
        #         bend_radius=bend_radius,
        #         gdsii_path=True,
        #         layer=self.layers.Nb_Antenna.number,
        #         datatype=self.layers.Nb_Antenna.datatype,
        #     )
        #     self.make_flexpath_into_polygons_and_add_to_main(center_path, self.layers.Nb_Antenna)
        #
        # if center_material == "Al":
        #     center_path = gdspy.FlexPath(
        #         feed_points,
        #         feedline_width,
        #         corners="circular bend",
        #         bend_radius=bend_radius,
        #         gdsii_path=True,
        #         layer=self.layers.Aluminium.number,
        #         datatype=self.layers.Aluminium.datatype,
        #     )
        #     self.make_flexpath_into_polygons_and_add_to_main(center_path, self.layers.Aluminium)
        #
        # if center_material == "Both":
        #     center_path1 = gdspy.FlexPath(
        #         feed_points,
        #         feedline_width,
        #         corners="circular bend",
        #         bend_radius=bend_radius,
        #         gdsii_path=True,
        #         layer=self.layers.Aluminium.number,
        #         datatype=self.layers.Aluminium.datatype,
        #     )
        #     center_path2 = gdspy.FlexPath(
        #         feed_points,
        #         feedline_width,
        #         corners="circular bend",
        #         bend_radius=bend_radius,
        #         gdsii_path=True,
        #         layer=self.layers.Nb_Antenna.number,
        #         datatype=self.layers.Nb_Antenna.datatype,
        #     )
        #     self.make_flexpath_into_polygons_and_add_to_main(center_path1, self.layers.Aluminium)
        #     self.make_flexpath_into_polygons_and_add_to_main(center_path2, self.layers.Nb_Antenna)

        ground_plane_cuttout_outer_path = gdspy.FlexPath(
            feed_points, cutout_around_feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True
        )  # makes the outer of the cpw flex path
        poly_points_flexpath = mbu.get_polys_from_flexpath(
            ground_plane_cuttout_outer_path
        )  # gets the polygon points of the outer cpw flex path via function and

        for i in range(len(poly_points_flexpath)):
            port_to_fdln_outer_path_polygon = gdspy.Polygon(
                poly_points_flexpath[i],
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            self.ground_plane_cutouts.add(
                [port_to_fdln_outer_path_polygon]
            )  # makes a polygon with that so it can be cutout of the ground plane (flex path cannot be added directly as boolean operation cant cut that out)

        dielectric_under_feedline_path = gdspy.FlexPath(
            feed_points,
            dielectric_under_feedline_width,
            corners="circular bend",
            bend_radius=bend_radius,
            gdsii_path=True,
            layer=self.layers.SiN_dep.number,
            datatype=self.layers.SiN_dep.datatype,
        )
        poly_points_flexpath = mbu.get_polys_from_flexpath(dielectric_under_feedline_path)

        for i in range(len(poly_points_flexpath)):
            port_to_fdln_outer_path_polygon = gdspy.Polygon(
                poly_points_flexpath[i],
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            self.silicon_nitride_positives.add([port_to_fdln_outer_path_polygon])

        if not add_bridges:
            return

        if points_to_avoid_bridging is not None:
            all_avoid_points = np.array(points_to_avoid_bridging)

        rounded_feedline_points = self.get_rounded_path_from_passthough_points(feed_points, bend_radius)

        # self.Main.add(gdspy.FlexPath(rounded_feedline_points, 20, layer=600))

        lengths_of_path_sections = np.zeros(len(rounded_feedline_points) - 1)
        angles_of_path_sections = np.zeros(len(rounded_feedline_points) - 1)
        for i in range(len(lengths_of_path_sections)):
            lengths_of_path_sections[i] = (
                (rounded_feedline_points[i][0] - rounded_feedline_points[i + 1][0]) ** 2
                + (rounded_feedline_points[i][1] - rounded_feedline_points[i + 1][1]) ** 2
            ) ** 0.5
            angles_of_path_sections[i] = np.arctan2(
                (rounded_feedline_points[i + 1][1] - rounded_feedline_points[i][1]),
                (rounded_feedline_points[i + 1][0] - rounded_feedline_points[i][0]),
            )

        total_conect_length = np.sum(lengths_of_path_sections)
        length_cumsum = np.cumsum(lengths_of_path_sections)

        bridge_gap = config["bridge_gap"]  # 5300/2
        bridge_width = config["bridge_width"]  # 6
        bridge_height = cutout_around_feedline_width + 5  # 70

        no_of_bridges = int(np.floor(total_conect_length / bridge_gap))

        bridge_offset_in_result_of_crash = 10
        distance_considered_clashing = 80

        for i in range(no_of_bridges):
            arg = np.argmin(np.abs(length_cumsum - (bridge_gap * (i + 1))))

            length_cumsum_added = np.sort(np.append(length_cumsum, bridge_gap * (i + 1)))
            arg = np.argwhere(length_cumsum_added == bridge_gap * (i + 1))[0][0]
            if arg == 0:
                length_from_prev_point = bridge_gap * (i + 1)
                length_perc_from_prev_point = length_from_prev_point / length_cumsum[0]
            else:
                length_from_prev_point = bridge_gap * (i + 1) - length_cumsum_added[arg - 1]
                length_perc_from_prev_point = length_from_prev_point / (length_cumsum[arg] - length_cumsum[arg - 1])

            x_bridge_pos = (
                rounded_feedline_points[arg][0]
                + (rounded_feedline_points[arg + 1][0] - rounded_feedline_points[arg][0]) * length_perc_from_prev_point
            )
            y_bridge_pos = (
                rounded_feedline_points[arg][1]
                + (rounded_feedline_points[arg + 1][1] - rounded_feedline_points[arg][1]) * length_perc_from_prev_point
            )

            if points_to_avoid_bridging == None:
                bridge = gdspy.Rectangle(
                    [x_bridge_pos - (bridge_width / 2), y_bridge_pos - (bridge_height / 2)],
                    [x_bridge_pos + (bridge_width / 2), y_bridge_pos + (bridge_height / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                bridge.rotate(angles_of_path_sections[arg], [x_bridge_pos, y_bridge_pos])
                self.Main.add(bridge)

            else:
                clash_testing = True
                while clash_testing:
                    x_dist_to_clash = all_avoid_points[:, 0] - x_bridge_pos
                    y_dist_to_clash = all_avoid_points[:, 1] - y_bridge_pos

                    if np.min((x_dist_to_clash**2 + y_dist_to_clash**2) ** 0.5) > distance_considered_clashing:
                        clash_testing = False
                        break
                    else:
                        x_bridge_pos += bridge_offset_in_result_of_crash

                bridge = gdspy.Rectangle(
                    [x_bridge_pos - (bridge_width / 2), y_bridge_pos - (bridge_height / 2)],
                    [x_bridge_pos + (bridge_width / 2), y_bridge_pos + (bridge_height / 2)],
                    layer=self.layers.Nb_Groundplane.number,
                    datatype=self.layers.Nb_Groundplane.datatype,
                )
                bridge.rotate(angles_of_path_sections[arg], [x_bridge_pos, y_bridge_pos])
                self.Main.add(bridge)

        return

    def add_arbitrary_CPW(
        self,
        feed_points: Sequence[Sequence[float | int]],
        bend_radius: float | int,
        feedline_width: float | int,
        cutout_around_feedline_width: float | int,
        corners_type: Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] | None = None,
        dielectric_under_feedline_width: float | int = 0.0,
        center_material: Layer | list[Layer] | None = None,
        add_bridges: bool = False,
        bridge_gap: float | int = 2650,
        bridge_width: float | int = 6,
        cutout_dilectric_around_end_distance: float | int = 0.0,
        cutout_groundplane_around_end_distance: float | int = 0.0,
        cutout_center_around_end_distance: float | int = 0.5,
        return_outer_poly_points: bool = False,
    ):
        """Adds a Co-Planar Waveguide (CPW) feedline through the points given
        along with the dielectric layer. The parameters determining the feeline
        geometry are passed as agurments. The CPW consists of the center line,
        cut out in the ground plane around this center line and a dielectric
        layer covering that.

        Parameters
        ----------
        feed_points: list
            list of [x,y] lists which are the coordinates for the feedline to pass
            through.

        bend_radius: float, int
            The radius for the bends in the CPW to be made.

        feedline_width: float, int
            The width of the center line of the CPW.

        cutout_around_feedline_width: float, int
            The width of the ground plane cutout centered around the
            center line.

        KwArgs
        ------
        corners_type : Literal["natural", "miter", "bevel", "round", "smooth", "circular bend"] | None = None
            This is the type of corners to add for the FlexPath. See gdspy
            FlexPath for more details. Default value is None which will set the
            corners to "circular bend" or "round" if the value for argument
            cutout_around_feedline_width is greather than or equal than
            (2 * bend_radius). This ensures no errors in generation.

        dielectric_under_feedline_width = 0
            This is a float or int and is the width of the dielectric that
            covers the center line and the cutout around that. This is only
            added when the value is non zero.

        center_material = "Nb"
            What material to make the center line of the CPW out of. Only takes
            string values of "Nb", "Al", "Both" or "Nb_Grnd", where the "Both"
            option will put both Niobium_Antenna and Aluminium.

        add_bridges = False
            Adds bridges across the CPW if set to True else it will just add
            the feedline withough bridging. Bridging will stop before the end
            of the dielectric, ie if the cutout_dielectric arg != 0, the
            bridging will also get removed back the same.

        cutout_dilectric_around_end_distance = 0
            If not zero this will create a cutout in the SiN dep layer at the
            end of the feedline of the length specified.

        cutout_center_around_end_distance = 0
            If not zero this will create a gap between the edge of the mask and
            center line of the feedline of the length specified.

        cutout_groundplane_around_end_distance = 0.5
            If not zero this will create a ground plane cutout at the end of
            the feedline that overextends the end of the feedline. The length
            of that overextention in the groundplane cutout will be this given
            value multiplied by the cutout_around_feedline_width.

        return_outer_poly_points = False
            Will return the outer polygon points if set to true.
        """

        if corners_type is not None:
            center_corners_type = corners_type
            cutout_corners_type = corners_type
            dielectirc_corners_type = corners_type
        else:
            if feedline_width >= 2 * bend_radius:
                center_corners_type = "round"
            else:
                center_corners_type = "circular bend"

            if cutout_around_feedline_width >= 2 * bend_radius:
                cutout_corners_type = "round"
            else:
                cutout_corners_type = "circular bend"

            if dielectric_under_feedline_width >= 2 * bend_radius:
                dielectirc_corners_type = "round"
            else:
                dielectirc_corners_type = "circular bend"

        if cutout_center_around_end_distance == 0:
            if isinstance(center_material, Layer):
                center_path = gdspy.FlexPath(
                    feed_points,
                    feedline_width,
                    corners=center_corners_type,
                    bend_radius=bend_radius,
                    gdsii_path=True,
                    layer=center_material.number,
                    datatype=center_material.datatype,
                )
                self.make_flexpath_into_polygons_and_add_to_main(center_path, self.layers.Nb_Antenna)
            elif isinstance(center_material, list):
                if all(isinstance(mat, Layer) for mat in center_material):
                    for material in center_material:
                        center_path = gdspy.FlexPath(
                            feed_points,
                            feedline_width,
                            corners=center_corners_type,
                            bend_radius=bend_radius,
                            gdsii_path=True,
                            layer=material.number,
                            datatype=material.datatype,
                        )
                        self.make_flexpath_into_polygons_and_add_to_main(center_path, self.layers.Nb_Antenna)
                else:
                    raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")
            else:
                raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")

        else:
            feedline_center_cell_cutouts = gdspy.Cell("feed_center_cutouts")

            angle_of_last_section = np.arctan2((feed_points[-2][1] - feed_points[-1][1]), (feed_points[-2][0] - feed_points[-1][0]))

            path_points = [
                [
                    feed_points[-1][0] + cutout_center_around_end_distance * cos(angle_of_last_section),
                    feed_points[-1][1] + cutout_center_around_end_distance * sin(angle_of_last_section),
                ],
                feed_points[-1],
            ]

            center_cutout_at_end_of_path = gdspy.FlexPath(
                path_points,
                feedline_width * 1.5,
                corners=center_corners_type,
                bend_radius=bend_radius,
                gdsii_path=True,
                layer=5000,
            )

            self.make_flexpath_into_polygons_and_add_to_main(center_cutout_at_end_of_path, self.layers.General_labeling)
            poly_points_flexpath = mbu.get_polys_from_flexpath(center_cutout_at_end_of_path)

            for i in range(len(poly_points_flexpath)):
                path_polygon = gdspy.Polygon(poly_points_flexpath[i], layer=5000)
                feedline_center_cell_cutouts.add([path_polygon])

            if isinstance(center_material, Layer):
                center_path = gdspy.FlexPath(
                    feed_points,
                    feedline_width,
                    corners=center_corners_type,
                    bend_radius=bend_radius,
                    gdsii_path=True,
                    layer=center_material.number,
                    datatype=center_material.datatype,
                )
                new_center_path = gdspy.boolean(
                    center_path,
                    feedline_center_cell_cutouts.get_polygons(),
                    "not",
                    layer=center_material.number,
                    datatype=center_material.datatype,
                )
                self.Main.add(new_center_path)
            elif isinstance(center_material, list):
                if all(isinstance(mat, Layer) for mat in center_material):
                    for material in center_material:
                        center_path = gdspy.FlexPath(
                            feed_points,
                            feedline_width,
                            corners=center_corners_type,
                            bend_radius=bend_radius,
                            gdsii_path=True,
                            layer=material.number,
                            datatype=material.datatype,
                        )
                        new_center_path = gdspy.boolean(
                            center_path,
                            feedline_center_cell_cutouts.get_polygons(),
                            "not",
                            layer=material.number,
                            datatype=material.datatype,
                        )
                        self.Main.add(new_center_path)
                else:
                    raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")
            else:
                raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")

        angle_of_last_section = np.arctan2((feed_points[-2][1] - feed_points[-1][1]), (feed_points[-2][0] - feed_points[-1][0]))

        if cutout_groundplane_around_end_distance == 0:
            ground_plane_cuttout_outer_path_points = feed_points
        else:
            ground_plane_cuttout_outer_path_points = feed_points.append(
                [
                    feed_points[-1][0] + 0.5 * cutout_around_feedline_width * cos(angle_of_last_section + pi),
                    feed_points[-1][1] + 0.5 * cutout_around_feedline_width * sin(angle_of_last_section + pi),
                ]
            )
        ground_plane_cuttout_outer_path = gdspy.FlexPath(
            ground_plane_cuttout_outer_path_points,
            cutout_around_feedline_width,
            corners=cutout_corners_type,
            bend_radius=bend_radius,
            gdsii_path=True,
        )

        for poly in mbu.get_polys_from_flexpath(ground_plane_cuttout_outer_path):
            port_to_fdln_outer_path_polygon = gdspy.Polygon(
                poly,
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            self.ground_plane_cutouts.add(port_to_fdln_outer_path_polygon)

        exclusion_around_feedline_extra_width = 3000
        exclusion_around_feedline_path = gdspy.FlexPath(
            feed_points,
            cutout_around_feedline_width + exclusion_around_feedline_extra_width,
            corners=cutout_corners_type,
            bend_radius=bend_radius,
            gdsii_path=True,
            layer=self.layers.SiN_dep.number,
            datatype=self.layers.SiN_dep.datatype,
        )
        outer_poly_points = mbu.get_polys_from_flexpath(exclusion_around_feedline_path)

        # Adding the dielectric under the feedlne.
        if dielectric_under_feedline_width != 0:
            dielectric_under_feedline_path = gdspy.FlexPath(
                feed_points,
                dielectric_under_feedline_width,
                corners=dielectirc_corners_type,
                bend_radius=bend_radius,
                gdsii_path=True,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            poly_points_flexpath = mbu.get_polys_from_flexpath(dielectric_under_feedline_path)

            for poly in mbu.get_polys_from_flexpath(dielectric_under_feedline_path):
                port_to_fdln_outer_path_polygon = gdspy.Polygon(
                    poly,
                    layer=self.layers.SiN_dep.number,
                    datatype=self.layers.SiN_dep.datatype,
                )
                self.silicon_nitride_positives.add(port_to_fdln_outer_path_polygon)

            exclusion_around_feedline_path = gdspy.FlexPath(
                feed_points,
                dielectric_under_feedline_width + exclusion_around_feedline_extra_width,
                corners=dielectirc_corners_type,
                bend_radius=bend_radius,
                gdsii_path=True,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            outer_poly_points = mbu.get_polys_from_flexpath(exclusion_around_feedline_path)

            if cutout_dilectric_around_end_distance != 0:
                angle_of_last_section = np.arctan2((feed_points[-2][1] - feed_points[-1][1]), (feed_points[-2][0] - feed_points[-1][0]))

                path_points = [
                    [
                        feed_points[-1][0] + cutout_dilectric_around_end_distance * cos(angle_of_last_section),
                        feed_points[-1][1] + cutout_dilectric_around_end_distance * sin(angle_of_last_section),
                    ],
                    [
                        feed_points[-1][0] + cutout_dilectric_around_end_distance * cos(angle_of_last_section + pi),
                        feed_points[-1][1] + cutout_dilectric_around_end_distance * sin(angle_of_last_section + pi),
                    ],
                ]

                dielectric_cutout_at_end_of_path = gdspy.FlexPath(
                    path_points,
                    dielectric_under_feedline_width * 1.5,
                    corners=dielectirc_corners_type,
                    bend_radius=bend_radius,
                    gdsii_path=True,
                    layer=self.layers.SiN_dep.number,
                    datatype=self.layers.SiN_dep.datatype,
                )

                for poly in mbu.get_polys_from_flexpath(dielectric_cutout_at_end_of_path):
                    port_to_fdln_outer_path_polygon = gdspy.Polygon(
                        poly,
                        layer=self.layers.SiN_dep.number,
                        datatype=self.layers.SiN_dep.datatype,
                    )
                    self.silicon_nitride_cutouts.add(port_to_fdln_outer_path_polygon)

        if not add_bridges:
            if not return_outer_poly_points:
                return
            else:
                return outer_poly_points

        rounded_feedline_points = self.get_rounded_path_from_passthough_points(feed_points, bend_radius)

        # self.Main.add(gdspy.FlexPath(rounded_feedline_points, 20, layer=600))

        lengths_of_path_sections = np.zeros(len(rounded_feedline_points) - 1)
        angles_of_path_sections = np.zeros(len(rounded_feedline_points) - 1)
        for i in range(len(lengths_of_path_sections)):
            lengths_of_path_sections[i] = (
                (rounded_feedline_points[i][0] - rounded_feedline_points[i + 1][0]) ** 2
                + (rounded_feedline_points[i][1] - rounded_feedline_points[i + 1][1]) ** 2
            ) ** 0.5
            angles_of_path_sections[i] = np.arctan2(
                (rounded_feedline_points[i + 1][1] - rounded_feedline_points[i][1]),
                (rounded_feedline_points[i + 1][0] - rounded_feedline_points[i][0]),
            )

        total_conect_length = np.sum(lengths_of_path_sections)
        length_cumsum = np.cumsum(lengths_of_path_sections)

        # bridge_gap = 5300 / 2
        # bridge_width = 6
        bridge_height = cutout_around_feedline_width + 10

        no_of_bridges = int(np.floor((total_conect_length - cutout_dilectric_around_end_distance) / bridge_gap))

        for i in range(no_of_bridges):
            arg = np.argmin(np.abs(length_cumsum - (bridge_gap * (i + 1))))

            length_cumsum_added = np.sort(np.append(length_cumsum, bridge_gap * (i + 1)))
            arg = np.argwhere(length_cumsum_added == bridge_gap * (i + 1))[0][0]
            if arg == 0:
                length_from_prev_point = bridge_gap * (i + 1)
                length_perc_from_prev_point = length_from_prev_point / length_cumsum[0]
            else:
                length_from_prev_point = bridge_gap * (i + 1) - length_cumsum_added[arg - 1]
                length_perc_from_prev_point = length_from_prev_point / (length_cumsum[arg] - length_cumsum[arg - 1])

            x_bridge_pos = (
                rounded_feedline_points[arg][0]
                + (rounded_feedline_points[arg + 1][0] - rounded_feedline_points[arg][0]) * length_perc_from_prev_point
            )
            y_bridge_pos = (
                rounded_feedline_points[arg][1]
                + (rounded_feedline_points[arg + 1][1] - rounded_feedline_points[arg][1]) * length_perc_from_prev_point
            )

            bridge = gdspy.Rectangle(
                [x_bridge_pos - (bridge_width / 2), y_bridge_pos - (bridge_height / 2)],
                [x_bridge_pos + (bridge_width / 2), y_bridge_pos + (bridge_height / 2)],
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            bridge.rotate(angles_of_path_sections[arg], [x_bridge_pos, y_bridge_pos])
            # self.ground_plane_positives.add(gdspy.Polygon(bridge.polygons[0], **self.Nb_Groundplane)) # TODO
            self.Main.add(bridge)

        if return_outer_poly_points:
            return outer_poly_points
        else:
            return

    def add_port_and_get_connection_point(
        self,
        x: float | int,
        y: float | int,
        rotation: float | int,
        cpw_feedline_config_override: dict[str, float | int] | None = None,
        port_config_override: dict[str, float | int] | None = None,
        center_material: Layer | None = None,
        add_extra_squares: bool = True,
        dielectric_cutout_in_port: bool = True,
        cutout_groundplane_around_end_distance: float | int = 0,
    ) -> tuple[float, float]:
        """Adds a port to the mask consisting of a tapered section and a back
        straight. The port will be added where the middle base of the back
        straight section sits at the x, y given.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate for the middle back of the port to be placed.

        rotation: float, int
            The rotation angle (**in radians**) of the port.

        KwArgs
        ------
        cpw_feedline_config_override,
        port_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        center_material: Layer | None = None
            This is an instance of Layer. see maskpy.layers.Layer.
            Usually this is within the SoukMaskBuilder.layers.xxx.
            e.g. `self.layers.Aluminium`.
            By default the value is None which will set the center material to
            Nb_Antenna.

        add_extra_squares: bool = True
            By default adds a series of 3 squares either side of the port in Nb.

        dielectric_cutout_in_port: bool = True
            By default adds a cutout in the dielectric over the center line of
            the port.

        cutout_groundplane_around_end_distance: float | int = 0
            Cutout the Nb groundplane behind the back of the port distance.
            This will be the same width as the rest of the port body and no
            dielectric will cover this region.
        """
        if center_material is None:
            center_material = self.layers.Nb_Antenna

        if not isinstance(center_material, Layer):
            raise TypeError(f"center_material should be of type Layer, current type is {type(center_material)}")

        cpw_feedline_config = get_mask_default_config(SoukMaskConfig.CPW_FEEDLINE, config_override=cpw_feedline_config_override)
        port_config = get_mask_default_config(SoukMaskConfig.PORT, config_override=port_config_override)

        outer_feedline_width = port_config["outer_feedline_width"]  # 500
        outer_cutout_around_feedline_width = port_config["outer_cutout_around_feedline_width"]  # 1125
        outer_dielectric_under_feedline_width = port_config["outer_dielectric_under_feedline_width"]  # 1245

        outer_back_length = port_config["outer_back_length"]  # 1000
        taper_length = port_config["taper_length"]  # 1500

        dielectric_cutout_in_port_width = port_config["dielectric_cutout_in_port_width"]  # 500
        dielectric_cutout_in_port_length = port_config["dielectric_cutout_in_port_length"]  # 1000

        inner_feedline_width = cpw_feedline_config["feedline_width"]
        inner_cutout_around_feedline_width = cpw_feedline_config["cutout_around_feedline_width"]
        inner_dielectric_under_feedline_width = cpw_feedline_config["dielectric_under_feedline_width"]

        port_center_points = [
            [x, y + (outer_feedline_width / 2)],
            [x + outer_back_length, y + (outer_feedline_width / 2)],
            [x + outer_back_length + taper_length, y + (inner_feedline_width / 2)],
            [x + outer_back_length + taper_length, y - (inner_feedline_width / 2)],
            [x + outer_back_length, y - (outer_feedline_width / 2)],
            [x, y - (outer_feedline_width / 2)],
        ]

        port_center = gdspy.Polygon(
            port_center_points,
            layer=center_material.number,
            datatype=center_material.datatype,
        )
        port_center.rotate(rotation, [x, y])
        self.Main.add(port_center)

        port_cutout_around_points = [
            [x - cutout_groundplane_around_end_distance, y + (outer_cutout_around_feedline_width / 2)],
            [x + outer_back_length, y + (outer_cutout_around_feedline_width / 2)],
            [x + outer_back_length + taper_length, y + (inner_cutout_around_feedline_width / 2)],
            [x + outer_back_length + taper_length, y - (inner_cutout_around_feedline_width / 2)],
            [x + outer_back_length, y - (outer_cutout_around_feedline_width / 2)],
            [x - cutout_groundplane_around_end_distance, y - (outer_cutout_around_feedline_width / 2)],
        ]

        port_cutout_around = gdspy.Polygon(
            port_cutout_around_points,
            layer=self.layers.Nb_Groundplane.number,
            datatype=self.layers.Nb_Groundplane.datatype,
        )
        port_cutout_around.rotate(rotation, [x, y])
        self.ground_plane_cutouts.add(port_cutout_around)

        port_dielectric_points = [
            [x, y + (outer_dielectric_under_feedline_width / 2)],
            [x + outer_back_length, y + (outer_dielectric_under_feedline_width / 2)],
            [x + outer_back_length + taper_length, y + (inner_dielectric_under_feedline_width / 2)],
            [x + outer_back_length + taper_length, y - (inner_dielectric_under_feedline_width / 2)],
            [x + outer_back_length, y - (outer_dielectric_under_feedline_width / 2)],
            [x, y - (outer_dielectric_under_feedline_width / 2)],
        ]

        port_dielectric = gdspy.Polygon(
            port_dielectric_points,
            layer=self.layers.SiN_dep.number,
            datatype=self.layers.SiN_dep.datatype,
        )
        port_dielectric.rotate(rotation, [x, y])
        self.silicon_nitride_positives.add(port_dielectric)

        if dielectric_cutout_in_port:
            port_dielectric_cutout_points = [
                [x, y + (dielectric_cutout_in_port_width / 2)],
                [x + dielectric_cutout_in_port_length, y + (dielectric_cutout_in_port_width / 2)],
                [x + dielectric_cutout_in_port_length, y - (dielectric_cutout_in_port_width / 2)],
                [x, y - (dielectric_cutout_in_port_width / 2)],
            ]
            port_dielectric_cutout = gdspy.Polygon(
                port_dielectric_cutout_points,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            port_dielectric_cutout.rotate(rotation, [x, y])
            self.silicon_nitride_cutouts.add(port_dielectric_cutout)

        if add_extra_squares:
            square_offset_from_edge = 200
            square_side_len = 400
            no_of_squares = 3

            for i in range(no_of_squares):
                offset_from_center = (outer_dielectric_under_feedline_width / 2) + square_side_len + (i) * 2 * square_side_len

                square_top = gdspy.Rectangle(
                    [x + square_offset_from_edge, y + offset_from_center],
                    [x + square_offset_from_edge + square_side_len, y + offset_from_center + square_side_len],
                    layer=self.layers.Nb_Antenna.number,
                    datatype=self.layers.Nb_Antenna.datatype,
                )
                square_top.rotate(rotation, [x, y])
                self.Main.add(square_top)

                square_bot = gdspy.Rectangle(
                    [x + square_offset_from_edge, y - offset_from_center],
                    [x + square_offset_from_edge + square_side_len, y - offset_from_center - square_side_len],
                    layer=self.layers.Nb_Antenna.number,
                    datatype=self.layers.Nb_Antenna.datatype,
                )
                square_bot.rotate(rotation, [x, y])
                self.Main.add(square_bot)

        feedline_connection_point = mbu.rotate(x, y, x + outer_back_length + taper_length, y, rotation)

        return feedline_connection_point

    def get_feedline_pass_through_points(
        self,
        x: float | int,
        y: float | int,
        relative_kid_positions: tuple[list[float | int], list[float | int], list[float | int], list[float | int]] | list[list[float | int]],
        resonator_types: tuple[SoukResonatorType, SoukResonatorType, SoukResonatorType, SoukResonatorType] | list[SoukResonatorType],
        cpw_feedline_config_override: dict[str, float | int] | None = None,
        resonator_config_overrides: list[dict[str, float | int] | None] = [None, None, None, None],
        **kwargs,
    ) -> list[list[float]]:
        """This will get a list of coordinates where the feedline should
        connect to the KIDs. This is the end of the coupler.

        Parameters
        ----------
        x, y: float, int
            The x, y coordinate of the center of the antenna structure,
            i.e. the center of the horn.

        relative_kid_positions: list
            list of [x,y] lists defining the connection point coordinates for each
            of the KIDs. This is the very bottom center point of the inductive
            meander section. This list should be the coordinates, in order, of the
            TopLeft, TopRight, BotLeft, BotRight KIDs.

        resonator_types: list[SoukResonatorType]
            This is the type of resonators drawn. The values accepted
            here are members of the SoukResonatorType enum.
            The order of the values passed in will be attributed to each KID
            and should be the same order as the rel_kid_positions, TL, TR, BL,
            BR.

        KwArgs
        ------
        cpw_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        resonator_config_overrides: list[dict[str, float | int] | None] = [None, None, None, None]
            This is a list of 4 optional override dictionarys containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        feed_pass_points: list
            list of [x,y] lists defining the coordinates where to feedline should
            connect to the KIDs. This list has the coonection coordinates for the
            KIDs in order top left, top right, bot left, bot right.
        """
        cpw_feedline_config = get_mask_default_config(SoukMaskConfig.CPW_FEEDLINE, config_override=cpw_feedline_config_override)

        half_feedline_width = cpw_feedline_config["feedline_width"] / 2

        if not isinstance(resonator_types, (list, tuple)):
            raise ValueError(f"resonator_types should be a list of 4 SoukResonatorTypes. Not '{resonator_types}'")

        feed_pass_points: list[list[float]] = []

        xy_signs = [
            [-1, +1],
            [+1, +1],
            [-1, -1],
            [+1, -1],
        ]

        for relative_kid_position, resonator_type, resonator_config_override, (x_sign, y_sign) in zip(
            relative_kid_positions, resonator_types, resonator_config_overrides, xy_signs, strict=False
        ):
            # note that vertical and horizontal are confused here. This is
            # because the souk_resonators functions refer to the refer to the
            # resonaotrs with zero rotation, i.e. the inductor at the bottom.
            # The vertical feedline offset is with the resonators rotated
            # by 90 degrees in the horn block.
            vertical_feedline_offset_from_rel_kid_position = souk_resonators.get_horizontal_coupler_end_to_meander_base_distance(
                resonator_type,
                resonator_config_override=resonator_config_override,
            )
            horizontal_feedline_offset_from_rel_kid_position = souk_resonators.get_vertical_coupler_center_to_meander_base_distance(
                resonator_type,
                resonator_config_override=resonator_config_override,
            )
            feed_pass_points.append(
                [
                    x + relative_kid_position[0] + (x_sign * horizontal_feedline_offset_from_rel_kid_position),
                    y + relative_kid_position[1] + (y_sign * (vertical_feedline_offset_from_rel_kid_position + half_feedline_width)),
                ]
            )

        return feed_pass_points

    def get_total_height_of_resonator(
        self,
        resonator_type: SoukResonatorType,
        resonator_config_override: dict[str, float | int] | None = None,
    ) -> float:
        """This will get the total height of the resonator from the base of the
        inductive meander to the end of the ground plane cutout at the top of
        the structure.

        Parameters
        ----------
        resonator_type: SoukResonatorType
            This is the type of resonator to get the total height of. The
            value accepted here are members of the SoukResonatorType enum.

        KwArgs
        ------
        resonator_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        KID_height: float
            The total height of the KID calculated from the config file.
        """

        return souk_resonators.get_total_height_of_resonator(resonator_type, resonator_config_override=resonator_config_override)

    def get_width_height_of_resonator_IDC_section(
        self,
        resonator_type: SoukResonatorType,
        resonator_config_override: dict[str, float | int] | None = None,
    ) -> tuple[float | int, float | int]:
        """This will get the total width and height of ground plane cutout
        around the IDC section of the KID.

        Parameters
        ----------
        resonator_type: SoukResonatorType
            This is the type of resonator to get the width and height of. The
            value accepted here are members of the SoukResonatorType enum.

        KwArgs
        ------
        resonator_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        [KID_width, KID_height]: list
            This is a list containing, in order, the KID_width and the KID_height
            of the KID's IDC section calculated from the config file.
        """
        resonator_width, resonator_height = souk_resonators.get_width_and_height_of_IDC_cutout_section(
            resonator_type,
            resonator_config_override=resonator_config_override,
        )
        return resonator_width, resonator_height

    def get_total_width_of_resonator(self, Main_config_file_dict):
        """This will get the total width of the resonator from the far left of
        the ground plane cutout to the far right of the ground plane cutout of
        the structure. This will always be the widest part of the KID.

        Parameters
        ----------
        Main_config_file_dict: dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        Returns
        -------
        KID_width: float
            The total width of the KID calculated from the config file.
        """
        required_key = "resonator"
        self.validate_if_config_dict_has_required_keys(Main_config_file_dict, required_key)
        res_config = Main_config_file_dict["resonator"]

        left_side_width = (
            (res_config["meander_bot_width"] / 2)
            - (res_config["meander_lw"] / 2)
            + (res_config["meander_left_width_1"] - res_config["meander_left_width_2"])
            + res_config["frame_bot_left_width"]
            + res_config["cutout_left_offset"]
        )

        right_side_width = (
            (res_config["meander_bot_width"] / 2)
            - (res_config["meander_lw"] / 2)
            + (res_config["meander_right_width_1"] - res_config["meander_right_width_2"])
            + res_config["frame_bot_right_width"]
            + res_config["cutout_right_offset"]
        )

        KID_width = right_side_width + left_side_width
        return KID_width

    def get_feedline_center_to_meander_base_distance(
        self,
        resonator_type: SoukResonatorType,
        resonator_config_override: dict[str, float | int] | None = None,
        cpw_feedline_config_override: dict[str, float | int] | None = None,
    ) -> float:
        """This will calculate the the vertical distance from the center of the
        feeline to where the base of the KIDs inductive meander will sit
        assuming the coupler arm buts up against the edge of the CPW feedline
        center strip. This distance is calculated based on the dimensions
        within the config file.

        Parameters
        ----------
        resonator_type: SoukResonatorType
            This is the type of resonator to get the feedline center to
            meander base distance for. The value accepted here are members of
            the SoukResonatorType enum.

        KwArgs
        ------
        resonator_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        cpw_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        feedline_to_meander_base_distance: float
            The vertical distance from the CPW feedline center to the middle of the
            base of the KIDs inductive meander calculated from the config file.
            This distance is similar to, but not the same as, the value
            "right_side_width" from within the "get_total_width_of_resonator"
            function but here this also calculates the distance to feedline
            dimension.
        """
        cpw_feedline_config = get_mask_default_config(SoukMaskConfig.CPW_FEEDLINE, config_override=cpw_feedline_config_override)

        coupler_end_to_meander_base_distance = souk_resonators.get_horizontal_coupler_end_to_meander_base_distance(
            resonator_type, resonator_config_override=resonator_config_override
        )
        feedline_width = cpw_feedline_config["feedline_width"]
        return (feedline_width / 2) + coupler_end_to_meander_base_distance

    def get_relative_kid_positions(
        self,
        resonator_types: list[SoukResonatorType] | tuple[SoukResonatorType, SoukResonatorType, SoukResonatorType, SoukResonatorType],
        resonator_config_overrides: list[dict[str, float | int] | None] = [None, None, None, None],
        general_config_override: dict[str, float | int] | None = None,
        **kwargs,
    ) -> tuple[list[float | int], list[float | int], list[float | int], list[float | int]]:
        """Gets the positions of the base of the KID meanders relative to the
        center of the antennas in the middle. The relative KID positions is
        returned as a list is the order, [top_left, top_right, bot_left,
        bot_right], where each element is an [x, y] list.

        Parameters
        ----------
        resonator_types: list[SoukResonatorType]
            This is the type of resonators to get the rel_kid_positions for.
            The values accepted here are members of the SoukResonatorType enum.
            The order of the resonator_types passed in should be [top_left,
            top_right, bot_left, bot_right]

        KwArgs
        ------
        resonator_config_overrides: list[dict[str, float | int] | None] = [None, None, None, None]
            This is a list of 4 optional override dictionarys containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        general_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        rel_kid_positions: tuple[list[float]]
            tuple of [x,y] lists that are the coordinates of the base of the KIDs
            meanders. These KID positons are in order top_left, top_right,
            bot_left, bot_right.
        """

        general_config = get_mask_default_config(SoukMaskConfig.GENERAL, config_override=general_config_override)

        rel_kid_positions = []

        xy_signs = [
            [-1, +1],
            [+1, +1],
            [-1, -1],
            [+1, -1],
        ]

        for i in range(4):
            resonator_config = souk_resonators.get_resonator_config(resonator_types[i], resonator_config_overrides[i])

            tot_kid_height = souk_resonators.get_total_height_of_resonator(
                resonator_types[i], resonator_config_override=resonator_config_overrides[i]
            )

            ground_gap_between_resonators = resonator_config["grndpln_gap_between_adjacent_resonators"]

            horizontal_offset_from_center_of_antenna = (
                (general_config["horizontal_pitch"] / 2) - tot_kid_height - (ground_gap_between_resonators / 2)
            )

            feed_to_meander_base_distance = self.get_feedline_center_to_meander_base_distance(
                resonator_types[i],
                resonator_config_override=resonator_config_overrides[i],
            )

            vertical_offset_from_center_of_antenna = (general_config["vertical_pitch"] / 2) - feed_to_meander_base_distance

            rel_kid_positions.append(
                [
                    xy_signs[i][0] * horizontal_offset_from_center_of_antenna,
                    xy_signs[i][1] * vertical_offset_from_center_of_antenna,
                ]
            )

        return tuple(rel_kid_positions)

    def get_relative_antenna_conect_positions(
        self,
        rotation: float | int,
        antenna_config_override: dict[str, float | int] | None = None,
    ):
        r"""Gets the connect positions of the end of the antennas relative to
        the center of the antennas. The relative antenna conect positions is
        returned as a list is the order, [bot, top, left, right], where each
        element is an [x, y] list. If the rotation setting in config is non-
        zero  this should be included in the rotation arg. The relative
        positions are only rotated by this rotation amount.

        Parameters
        ----------
        rotation: float
            The angle (**in radians**) of rotation of the antennas about the center
            of the antennas. A negative rotation angle is clockwise, positive is
            anti-clockwise.

        KwArgs
        ------
        antenna_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Return
        ------
        rel_ant_conect_positions: list
            list of [x,y] lists that describe the coordinates for the end of the
            antennas wher they will feed into a KID or filer bank. These antenna
            connection positons are in order bot, top, left, right.
            **Note** These refer to the locations of the end of the antennas before
            any rotation has been done. i.e. if rotated $\pi$ radians the top
            refers to the antenna at the bottom now, left will refer to right
            antenna and so on.
        """
        antenna_config = get_mask_default_config(SoukMaskConfig.ANTENNA, config_override=antenna_config_override)

        ant_connect_offset_from_center = (
            antenna_config["distance_from_center"] + antenna_config["straight_height"] + antenna_config["taper_height"]
        )

        rot = mbu.deg_to_rad(rotation)

        angles = [(-pi / 2 + rot), (pi / 2 + rot), (pi + rot), (0 + rot)]

        rel_ant_conect_positions = []

        for angle in angles:
            rel_ant_conect_positions.append([ant_connect_offset_from_center * cos(angle), ant_connect_offset_from_center * sin(angle)])

        return rel_ant_conect_positions

    def add_sma_connector_and_launcher_and_get_connection(
        self,
        x: float | int,
        y: float | int,
        rot: float | int,
        sma_config_override: dict[str, float | int] | None = None,
        cpw_feedline_config_override: dict[str, float | int] | None = None,
        bend: str = "none",
        cutout_dielectric_over_taper: bool = True,
        return_configurator_points: bool = False,
    ) -> list[float | int] | tuple[list[float | int], dict[str, dict[str, Any]]]:
        """Creates an SMA connector where the center pin is located at the x,y
        given and returns the point where a feedline should connect.

        Parameters
        ----------
        x, y: float | int
            The x, y coordinate for the center of the SMA conector pin.

        rot: float | int
            Angle (**in radians**), the rotation of the whole assembly around the
            center x, y given. positive is anti-clockwise, negative is clockwise.

        KwArgs
        ------
        sma_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        cpw_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        bend: str, {"none", "left", "right"}
            If not "none" it will bend the cpw coming out of the sma connector
            around to either the left or right 90 degrees and then calculate the
            taper to the edge of the connector. Anything other than "left" will be regarded as "none"

        cutout_dielectric_over_taper: bool = True
            Default True. This will cutout the dielectirc over the taper
            section and some of the section behind that also. This is overall
            5000um of cutout. When False this will not be cutout.

        return_configurator_points: bool = False
            return a the points for use in the configurator. Default False.

        Returns
        -------
        conection: list[float | int]
            The [x, y] connection point at the end of the taper for the SMA
            conector where a feedline should connect.
        """
        sma_config = get_mask_default_config(SoukMaskConfig.SMA_CONNECTOR, config_override=sma_config_override)
        feed_config = get_mask_default_config(SoukMaskConfig.CPW_FEEDLINE, config_override=cpw_feedline_config_override)

        sma_square_offset_left = sma_config["sma_square_offset_left"]  # 0
        # sma_square_offset_right = 7000  # sma_config["sma_square_offset_right"]#5000
        sma_square_offset_right = sma_config["sma_square_offset_right"]  # 5000
        sma_square_offset_top = sma_config["sma_square_offset_top"]  # 1250
        sma_square_offset_bot = sma_config["sma_square_offset_bot"]  # 1250

        sma_circle_radius = sma_config["sma_circle_radius"]  # 2160

        sma_square_width = sma_config["sma_square_width"]  # 5000
        sma_square_height = sma_config["sma_square_height"]  # 5000

        taper_length = sma_config["taper_length"]  # 3000

        central_linewidth = sma_config["central_linewidth"]  # 500
        gap = sma_config["gap"]  # 500
        dielectric_overlap_distance = sma_config["dielectric_overlap_distance"]  # 300

        feedline_width = feed_config["feedline_width"]
        cutout_around_feedline_width = feed_config["cutout_around_feedline_width"]
        dielectric_under_feedline_width = feed_config["dielectric_under_feedline_width"]

        # showing where conenctor is
        connector_circle = gdspy.Round([x, y], sma_circle_radius, layer=400)
        self.Main.add(connector_circle)

        big_box = gdspy.Rectangle(
            [x - (sma_square_width / 2) - sma_square_offset_left, y - (sma_square_height / 2) - sma_square_offset_bot],
            [x + (sma_square_width / 2) + sma_square_offset_right, y + (sma_square_height / 2) + sma_square_offset_top],
            layer=401,
        )
        big_box.rotate(rot, [x, y])
        self.Main.add(big_box)

        sma_box = gdspy.Rectangle(
            [x - (sma_square_width / 2), y - (sma_square_height / 2)], [x + (sma_square_width / 2), y + (sma_square_height / 2)], layer=401
        )
        sma_box.rotate(rot, [x, y])
        self.Main.add(sma_box)

        # TODO Make this a propper implimentation for bends in path also
        # making the dielectric cutout in the wide section if the launcher
        if cutout_dielectric_over_taper:
            dielectirc_cutout_length = 5000
            dielectirc_cutout_height = 5000
            cutout_box = gdspy.Rectangle(
                [x + (sma_square_width / 2) + sma_square_offset_right - dielectirc_cutout_length, y - (dielectirc_cutout_height / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right, y + (dielectirc_cutout_height / 2)],
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            cutout_box.rotate(rot, [x, y])
            self.silicon_nitride_cutouts.add([cutout_box])

        configurator_points = {}

        # making cpw and taper
        if (bend != "left") and (bend != "right"):
            central_section_poly_points = [
                [x - (sma_square_width / 2) - sma_square_offset_left, y + (central_linewidth / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right - taper_length, y + (central_linewidth / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right, y + (feedline_width / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right, y - (feedline_width / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right - taper_length, y - (central_linewidth / 2)],
                [x - (sma_square_width / 2) - sma_square_offset_left, y - (central_linewidth / 2)],
            ]

            grnd_cutout_poly_points = [
                [x - (sma_square_width / 2) - sma_square_offset_left, y + (central_linewidth / 2) + gap],
                [x + (sma_square_width / 2) + sma_square_offset_right - taper_length, y + (central_linewidth / 2) + gap],
                [x + (sma_square_width / 2) + sma_square_offset_right, y + (cutout_around_feedline_width / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right, y - (cutout_around_feedline_width / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right - taper_length, y - (central_linewidth / 2) - gap],
                [x - (sma_square_width / 2) - sma_square_offset_left, y - (central_linewidth / 2) - gap],
            ]

            dielectric_poly_points = [
                [x - (sma_square_width / 2) - sma_square_offset_left, y + (central_linewidth / 2) + gap + dielectric_overlap_distance],
                [
                    x + (sma_square_width / 2) + sma_square_offset_right - taper_length,
                    y + (central_linewidth / 2) + gap + dielectric_overlap_distance,
                ],
                [x + (sma_square_width / 2) + sma_square_offset_right, y + (dielectric_under_feedline_width / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right, y - (dielectric_under_feedline_width / 2)],
                [
                    x + (sma_square_width / 2) + sma_square_offset_right - taper_length,
                    y - (central_linewidth / 2) - gap - dielectric_overlap_distance,
                ],
                [x - (sma_square_width / 2) - sma_square_offset_left, y - (central_linewidth / 2) - gap - dielectric_overlap_distance],
            ]

            # adding the central path
            central_section_poly = gdspy.Polygon(
                central_section_poly_points,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            central_section_poly.rotate(rot, [x, y])
            self.Main.add(central_section_poly)

            # adding the ground plane cutout around the central path
            grnd_cutout_poly = gdspy.Polygon(
                grnd_cutout_poly_points,
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            grnd_cutout_poly.rotate(rot, [x, y])
            self.ground_plane_cutouts.add(grnd_cutout_poly)

            # adding the dielectric around the ground plane cutout
            dielectric_poly = gdspy.Polygon(
                dielectric_poly_points,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            dielectric_poly.rotate(rot, [x, y])
            self.silicon_nitride_positives.add(dielectric_poly)

            # getting the conection point
            connection_point_xy = mbu.rotate(x, y, x + (sma_square_width / 2) + sma_square_offset_right, y, rot)

            configurator_points["taper_length"] = {
                "text": "taper_length",
                "start": [grnd_cutout_poly_points[1][0], grnd_cutout_poly_points[1][1]],
                "end": [grnd_cutout_poly_points[2][0], grnd_cutout_poly_points[1][1]],
            }

            configurator_points["dielectric_overlap_distance_top"] = {
                "text": "dielectric_overlap_distance",
                "start": [x, y + (central_linewidth / 2) + gap],
                "end": [x, y + (central_linewidth / 2) + gap + dielectric_overlap_distance],
            }

            configurator_points["dielectric_overlap_distance_bot"] = {
                "text": "dielectric_overlap_distance",
                "start": [x, y - (central_linewidth / 2) - gap],
                "end": [x, y - (central_linewidth / 2) - gap - dielectric_overlap_distance],
            }

        else:
            if bend == "left":
                final_bend_angle = -pi / 2
                sign = +1
            else:
                final_bend_angle = +pi / 2
                sign = -1

            bend_offset_from_end = sma_config["bend_offset_from_end"]  # 200
            extra_bend_center = sma_config["extra_bend_center"]  # 50

            taper_length = (
                (sma_square_height / 2)
                + sma_square_offset_top
                - (central_linewidth / 2)
                - gap
                - dielectric_overlap_distance
                - extra_bend_center
            )

            central_section_poly_points = [
                [x - (sma_square_width / 2) - sma_square_offset_left, y + (central_linewidth / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end, y + (central_linewidth / 2)],
                [x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end, y - (central_linewidth / 2)],
                [x - (sma_square_width / 2) - sma_square_offset_left, y - (central_linewidth / 2)],
            ]

            grnd_cutout_poly_points = [
                [x - (sma_square_width / 2) - sma_square_offset_left, y + (central_linewidth / 2) + gap],
                [
                    x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end,
                    y + (central_linewidth / 2) + gap,
                ],
                [
                    x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end,
                    y - (central_linewidth / 2) - gap,
                ],
                [x - (sma_square_width / 2) - sma_square_offset_left, y - (central_linewidth / 2) - gap],
            ]

            dielectric_poly_points = [
                [x - (sma_square_width / 2) - sma_square_offset_left, y + (central_linewidth / 2) + gap + dielectric_overlap_distance],
                [
                    x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end,
                    y + (central_linewidth / 2) + gap + dielectric_overlap_distance,
                ],
                [
                    x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end,
                    y - (central_linewidth / 2) - gap - dielectric_overlap_distance,
                ],
                [x - (sma_square_width / 2) - sma_square_offset_left, y - (central_linewidth / 2) - gap - dielectric_overlap_distance],
            ]

            central_section_poly = gdspy.Polygon(
                central_section_poly_points,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            central_section_poly.rotate(rot, [x, y])
            self.Main.add(central_section_poly)

            # adding the ground plane cutout around the central path
            grnd_cutout_poly = gdspy.Polygon(
                grnd_cutout_poly_points,
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            grnd_cutout_poly.rotate(rot, [x, y])
            self.ground_plane_cutouts.add(grnd_cutout_poly)

            # adding the dielectric around the ground plane cutout
            dielectric_poly = gdspy.Polygon(
                dielectric_poly_points,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            dielectric_poly.rotate(rot, [x, y])
            self.Main.add(dielectric_poly)

            # making the bend for each section
            bend_center_xy = [
                x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end,
                y + sign * ((central_linewidth / 2) + gap + dielectric_overlap_distance + extra_bend_center),
            ]

            dielectric_section_inner_bend_radius = extra_bend_center
            dielectric_section_outer_bend_radius = (
                dielectric_section_inner_bend_radius + (2 * dielectric_overlap_distance) + (2 * gap) + central_linewidth
            )

            dielectric_section_bend = gdspy.Round(
                bend_center_xy,
                dielectric_section_outer_bend_radius,
                inner_radius=dielectric_section_inner_bend_radius,
                initial_angle=0,
                final_angle=final_bend_angle,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            dielectric_section_bend.rotate(rot, [x, y])
            self.Main.add(dielectric_section_bend)

            grnd_cutout_inner_bend_radius = extra_bend_center + dielectric_overlap_distance
            grnd_cutout_outer_bend_radius = grnd_cutout_inner_bend_radius + (2 * gap) + central_linewidth

            grnd_cutout_bend = gdspy.Round(
                bend_center_xy,
                grnd_cutout_outer_bend_radius,
                inner_radius=grnd_cutout_inner_bend_radius,
                initial_angle=0,
                final_angle=final_bend_angle,
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            grnd_cutout_bend.rotate(rot, [x, y])
            self.ground_plane_cutouts.add(grnd_cutout_bend)

            central_section_inner_bend_radius = extra_bend_center + dielectric_overlap_distance + gap
            central_section_outer_bend_radius = central_section_inner_bend_radius + central_linewidth

            central_section_bend = gdspy.Round(
                bend_center_xy,
                central_section_outer_bend_radius,
                inner_radius=central_section_inner_bend_radius,
                initial_angle=0,
                final_angle=final_bend_angle,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            central_section_bend.rotate(rot, [x, y])
            self.Main.add(central_section_bend)

            end_of_round_x = bend_center_xy[0] + extra_bend_center + dielectric_overlap_distance + gap + (central_linewidth / 2)
            end_of_round_y = bend_center_xy[1]

            central_taper_poly_points = [
                [end_of_round_x - (central_linewidth / 2), end_of_round_y],
                [end_of_round_x - (feedline_width / 2), end_of_round_y + sign * taper_length],
                [end_of_round_x + (feedline_width / 2), end_of_round_y + sign * taper_length],
                [end_of_round_x + (central_linewidth / 2), end_of_round_y],
            ]

            grnd_cutout_taper_poly_points = [
                [end_of_round_x - (central_linewidth / 2) - gap, end_of_round_y],
                [end_of_round_x - (cutout_around_feedline_width / 2), end_of_round_y + sign * taper_length],
                [end_of_round_x + (cutout_around_feedline_width / 2), end_of_round_y + sign * taper_length],
                [end_of_round_x + (central_linewidth / 2) + gap, end_of_round_y],
            ]

            dielectric_taper_poly_points = [
                [end_of_round_x - (central_linewidth / 2) - gap - dielectric_overlap_distance, end_of_round_y],
                [end_of_round_x - (dielectric_under_feedline_width / 2), end_of_round_y + sign * taper_length],
                [end_of_round_x + (dielectric_under_feedline_width / 2), end_of_round_y + sign * taper_length],
                [end_of_round_x + (central_linewidth / 2) + gap + dielectric_overlap_distance, end_of_round_y],
            ]

            # adding the central taper
            central_taper_poly = gdspy.Polygon(
                central_taper_poly_points,
                layer=self.layers.Nb_Antenna.number,
                datatype=self.layers.Nb_Antenna.datatype,
            )
            central_taper_poly.rotate(rot, [x, y])
            self.Main.add(central_taper_poly)

            # adding the ground plane cutout taper the central path
            grnd_cutout_taper_poly = gdspy.Polygon(
                grnd_cutout_taper_poly_points,
                layer=self.layers.Nb_Groundplane.number,
                datatype=self.layers.Nb_Groundplane.datatype,
            )
            grnd_cutout_taper_poly.rotate(rot, [x, y])
            self.ground_plane_cutouts.add(grnd_cutout_taper_poly)

            # adding the dielectric taper around the ground plane cutout
            dielectric_taper_poly = gdspy.Polygon(
                dielectric_taper_poly_points,
                layer=self.layers.SiN_dep.number,
                datatype=self.layers.SiN_dep.datatype,
            )
            dielectric_taper_poly.rotate(rot, [x, y])
            self.Main.add(dielectric_taper_poly)

            # getting the conection point
            connection_point_xy = mbu.rotate(x, y, end_of_round_x, end_of_round_y + sign * taper_length, rot)

        if not return_configurator_points:
            return [connection_point_xy[0], connection_point_xy[1]]

        configurator_points["sma_square_offset_left"] = {
            "text": "sma_square_offset_left",
            "start": [x - (sma_square_width / 2), y],
            "end": [x - (sma_square_width / 2) - sma_square_offset_left, y],
        }

        configurator_points["sma_square_offset_right"] = {
            "text": "sma_square_offset_right",
            "start": [x + (sma_square_width / 2), y],
            "end": [x + (sma_square_width / 2) + sma_square_offset_right, y],
        }

        configurator_points["sma_square_offset_top"] = {
            "text": "sma_square_offset_top",
            "start": [x, y + (sma_square_width / 2)],
            "end": [x, y + (sma_square_width / 2) + sma_square_offset_top],
        }

        configurator_points["sma_square_offset_bot"] = {
            "text": "sma_square_offset_bot",
            "start": [x, y - (sma_square_width / 2)],
            "end": [x, y - (sma_square_width / 2) - sma_square_offset_bot],
        }

        configurator_points["sma_circle_radius"] = {
            "text": "sma_circle_radius",
            "start": [x, y],
            "end": [x + sma_circle_radius * cos(pi / 4), y + sma_circle_radius * sin(pi / 4)],
        }

        configurator_points["sma_square_width"] = {
            "text": "sma_square_width",
            "start": [x - (sma_square_width / 2), y - (sma_square_height / 4)],
            "end": [x + (sma_square_width / 2), y - (sma_square_height / 4)],
        }

        configurator_points["sma_square_height"] = {
            "text": "sma_square_height",
            "start": [x - (sma_square_width / 4), y - (sma_square_height / 2)],
            "end": [x - (sma_square_width / 4), y + (sma_square_height / 2)],
        }

        configurator_points["central_linewidth"] = {
            "text": "central_linewidth",
            "start": [x, y - (central_linewidth / 2)],
            "end": [x, y + (central_linewidth / 2)],
        }

        configurator_points["gap_top"] = {
            "text": "gap",
            "start": [x, y + (central_linewidth / 2)],
            "end": [x, y + (central_linewidth / 2) + gap],
        }

        configurator_points["gap_bot"] = {
            "text": "gap",
            "start": [x, y - (central_linewidth / 2)],
            "end": [x, y - (central_linewidth / 2) - gap],
        }

        return [connection_point_xy[0], connection_point_xy[1]], configurator_points

    def get_sma_connector_and_laucher_bounding_box(
        self,
        x: float | int,
        y: float | int,
        rot: float | int,
        sma_config_override: dict[str, float | int] | None = None,
    ) -> list[list[float | int]]:
        """Gets the outer bounding box around an SMA connector where the center
        pin is located at the x,y given.

        Parameters
        ----------
        x, y: float | int
            The x, y coordinate for the center of the SMA conector pin.

        rot: float | int
            Angle (**in radians**), the rotation of the whole assembly around the
            center x, y given. positive is anti-clockwise, negative is clockwise.

        KwArgs
        ------
        sma_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        bounding_box_points: list[list[float | int]]
            Returns the bounding box points of the SMA connector. This is a list
            of [x,y] lists which are coordinates for the top left, top right,
            bot right, bot left points in this order.
        """
        sma_config = get_mask_default_config(SoukMaskConfig.SMA_CONNECTOR, config_override=sma_config_override)

        sma_square_offset_left = sma_config["sma_square_offset_left"]  # 0
        sma_square_offset_right = sma_config["sma_square_offset_right"]  # 5000
        sma_square_offset_top = sma_config["sma_square_offset_top"]  # 1250
        sma_square_offset_bot = sma_config["sma_square_offset_bot"]  # 1250

        sma_square_width = sma_config["sma_square_width"]  # 5000
        sma_square_height = sma_config["sma_square_height"]  # 5000

        # top left, top right, bot right, bot left
        bounding_box_points = [
            [x - (sma_square_width / 2) - sma_square_offset_left, y + (sma_square_height / 2) + sma_square_offset_top],
            [x + (sma_square_width / 2) + sma_square_offset_right, y + (sma_square_height / 2) + sma_square_offset_top],
            [x + (sma_square_width / 2) + sma_square_offset_right, y - (sma_square_height / 2) - sma_square_offset_bot],
            [x - (sma_square_width / 2) - sma_square_offset_left, y - (sma_square_height / 2) - sma_square_offset_bot],
        ]

        bounding_box_points = mbu.rotate_and_move_points_list(bounding_box_points, rot, 0, 0, ox=x, oy=y)

        return bounding_box_points

    def remove_exclusions_from_hex_pack(
        self,
        hex_pack: list[list[float | int]],
        exclusions_list: list[list[float | int]],
        general_config_override: dict[str, float | int] | None = None,
    ) -> list[list[float | int]]:
        """Removes hex pack points where any part of the KID block would
        intersect with another feature in the exclusions list.

        Parameters
        ----------
        hex_pack: list
            list containing individual [x,y] lists defining the center of all the
            hex pack grid points.

        exclusions_list: list
            list containing one or many lists which each are a list of [x, y] lists
            defining the boundary points for the exclusion polygon shape.

        KwArgs
        ------
        general_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        new_hex_pack: list[list[float | int]]
            Returns a list containing [x, y] lists defining the center of the hex
            points (similar in form to the hex_pack argument) that no longer
            clashes with any features from the exclusion list.
        """
        general_config = get_mask_default_config(SoukMaskConfig.GENERAL, config_override=general_config_override)

        kid_block_horizontal_size = general_config["horizontal_pitch"]
        kid_block_vertical_size = general_config["vertical_pitch"]

        new_hex_pack = []

        for hex_pack_point_xy in hex_pack:
            x = hex_pack_point_xy[0]
            y = hex_pack_point_xy[1]

            polygon_surrounding_kid_block = shapely_geom.Polygon(
                [
                    [x - (kid_block_horizontal_size / 2), y + (kid_block_vertical_size / 2)],
                    [x + (kid_block_horizontal_size / 2), y + (kid_block_vertical_size / 2)],
                    [x + (kid_block_horizontal_size / 2), y - (kid_block_vertical_size / 2)],
                    [x - (kid_block_horizontal_size / 2), y - (kid_block_vertical_size / 2)],
                ]
            )
            collision = False

            for exclusion_poly_points in exclusions_list:
                exclusion_polygon = shapely_geom.Polygon(exclusion_poly_points)

                if exclusion_polygon.intersects(polygon_surrounding_kid_block):
                    collision = True
                    break

            if not collision:
                new_hex_pack.append(hex_pack_point_xy)

        return new_hex_pack

    def remove_exclusions_from_split_hex_pack(
        self,
        split_hex_pack: list[list[list[float | int]]],
        exclusions_list: list[list[float | int]],
        general_config_override: dict[str, float | int] | None = None,
    ) -> list[list[list[float | int]]]:
        """Removes hex pack points where any part of the KID block would
        intersect with another feature in the exclusions list.

        Parameters
        ----------
        split_hex_pack: list
            list containing two seperate lists of the hex pack grid. Each of those
            lists should be a list of [x, y] lists defining the center of the hex
            point.

        exclusions_list: list
            list containing one or many lists which each are a list of [x, y] lists
            defining the boundary points for the exclusion polygon shape.

        KwArgs
        ------
        general_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        Returns
        -------
        [new_bot_hex_pack, new_top_hex_pack]: list[list[list[float | int]]]
            Returns a list (similar in form to the split_hex_pack argument) of two
            lists containing [x, y] lists defining the center of the hex points
            that no longer clashes with any features from the exclusion list.
        """
        general_config = get_mask_default_config(SoukMaskConfig.GENERAL, config_override=general_config_override)

        kid_block_horizontal_size = general_config["horizontal_pitch"]
        kid_block_vertical_size = general_config["vertical_pitch"]

        new_bot_hex_pack = []
        new_top_hex_pack = []

        bot_hex_pack = split_hex_pack[0]
        top_hex_pack = split_hex_pack[1]

        for hex_pack_point_xy in bot_hex_pack:
            x = hex_pack_point_xy[0]
            y = hex_pack_point_xy[1]

            polygon_surrounding_kid_block = shapely_geom.Polygon(
                [
                    [x - (kid_block_horizontal_size / 2), y + (kid_block_vertical_size / 2)],
                    [x + (kid_block_horizontal_size / 2), y + (kid_block_vertical_size / 2)],
                    [x + (kid_block_horizontal_size / 2), y - (kid_block_vertical_size / 2)],
                    [x - (kid_block_horizontal_size / 2), y - (kid_block_vertical_size / 2)],
                ]
            )
            collision = False

            for exclusion_poly_points in exclusions_list:
                exclusion_polygon = shapely_geom.Polygon(exclusion_poly_points)

                if exclusion_polygon.intersects(polygon_surrounding_kid_block):
                    collision = True
                    break

            if not collision:
                new_bot_hex_pack.append(hex_pack_point_xy)

        for hex_pack_point_xy in top_hex_pack:
            x = hex_pack_point_xy[0]
            y = hex_pack_point_xy[1]

            polygon_surrounding_kid_block = shapely_geom.Polygon(
                [
                    [x - (kid_block_horizontal_size / 2), y + (kid_block_vertical_size / 2)],
                    [x + (kid_block_horizontal_size / 2), y + (kid_block_vertical_size / 2)],
                    [x + (kid_block_horizontal_size / 2), y - (kid_block_vertical_size / 2)],
                    [x - (kid_block_horizontal_size / 2), y - (kid_block_vertical_size / 2)],
                ]
            )
            collision = False

            for exclusion_poly_points in exclusions_list:
                exclusion_polygon = shapely_geom.Polygon(exclusion_poly_points)

                if exclusion_polygon.intersects(polygon_surrounding_kid_block):
                    collision = True
                    break

            if not collision:
                new_top_hex_pack.append(hex_pack_point_xy)

        return [new_bot_hex_pack, new_top_hex_pack]

    def combine_split_hex_pack_grid(self, split_hex_pack_grid):
        """Combines a split hex pack grid into one long list instead of two
        seperate lists.

        Parameters
        ----------
        split_hex_pack_grid: list
            list containing two seperate lists of the hex pack grid. Each of those
            lists should be a list of [x, y] lists defining the center of the hex
            point.


        Returns
        -------

        combined_hex_pack_grid: list
            list containing individual [x,y] lists defining the center of all the
            hex pack grid points.
        """
        combined_hex_pack_grid = []

        for xy in split_hex_pack_grid[0]:
            combined_hex_pack_grid.append(xy)

        for xy in split_hex_pack_grid[1]:
            combined_hex_pack_grid.append(xy)

        return combined_hex_pack_grid

    def get_feedline_running_list(
        self,
        feed_pass_points: Sequence[Sequence[float | int]],
        cpw_feedline_config_override: dict[str, float | int] | None = None,
        init_direction: str = "right",
        **kwargs,
    ) -> list[list[float | int]]:
        """This will take in a Sequence of all the points that a feedline
        should pass through and then organizes them into a line running left to
        right on one row then right to left on the next row above and so on
        untill all points have been mapped like this. This can start right to
        left or left to right depending upon the init_direction. This is used
        to then be able to draw a feedline that meanders nicely through the
        mask passing over all the KIDs.

        Parameters
        ----------
        feed_pass_points: Sequence[Sequence[float | int]]
            Sequence of [x,y] Sequences that define the coordinates of all the points the
            feedline should pass through.

        KwArgs
        ------
        cpw_feedline_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        init_direction: str = "right"
            The initial direction to order the fitst row of points. Only takes
            strings "right" or "left". Default is "right" which will take the
            lowest row and start from the left most point and go right.

        Returns
        -------
        running_list: list[list[float | int]]
            list of [x,y] lists for all the points the feedline should pass through
            ordered in rows running right to left on one row and then left to right
            on the row above and so on untill all the points in the original
            feed_pass_points list have been ordered.
        """
        cpw_feedline_config = get_mask_default_config(SoukMaskConfig.CPW_FEEDLINE, config_override=cpw_feedline_config_override)

        sorted_list = sorted(feed_pass_points.copy(), key=lambda k: [k[1], k[0]])
        xs_sorted = np.array([a[0] for a in sorted_list])
        ys_sorted = [a[1] for a in sorted_list]
        indexes = [index for index, _ in enumerate(ys_sorted) if ys_sorted[index] != ys_sorted[index - 1]]
        indexes.append(len(ys_sorted))

        running_list = []

        bend_rad = cpw_feedline_config["bend_radius"]
        extra_straight_length = cpw_feedline_config["extra_straight_length"]

        bend_extend_length = bend_rad * tan(pi / 3)

        extend_length = extra_straight_length + bend_extend_length

        for i in range(len(indexes) - 1):
            y_val = ys_sorted[indexes[i]]
            x_vals = xs_sorted[indexes[i] : indexes[i + 1]]

            left_side_point = [x_vals[0] - extend_length, y_val]
            right_side_point = [x_vals[-1] + extend_length, y_val]

            if init_direction == "right":
                direction = +1 if (i) % 2 == 1 else -1  # + is right, - is left
            else:
                direction = -1 if (i) % 2 == 1 else +1  # + is right, - is left

            if direction == 1:
                # going right
                running_list.append(right_side_point)
                running_list.append(left_side_point)
            else:
                # going left
                running_list.append(left_side_point)
                running_list.append(right_side_point)

        for count, xy in enumerate(feed_pass_points):
            marker_bl_side = gdspy.Rectangle(
                (xy[0] - 5, xy[1] - 5),
                (xy[0], xy[1]),
                layer=self.layers.General_labeling.number,
                datatype=self.layers.General_labeling.datatype,
            )
            marker_tr_side = gdspy.Rectangle(
                (xy[0], xy[1]),
                (xy[0] + 5, xy[1] + 5),
                layer=self.layers.General_labeling.number,
                datatype=self.layers.General_labeling.datatype,
            )
            self.add_fancy_text(
                f"FPP{count}",
                xy[0],
                xy[1],
                5,
                self.layers.General_labeling,
                horizontal_align="end",
                vertical_align="above",
            )
            self.Main.add(marker_bl_side)
            self.Main.add(marker_tr_side)

        for count, xy in enumerate(running_list):
            marker_bl_side = gdspy.Rectangle(
                (xy[0] - 5, xy[1] - 5),
                (xy[0], xy[1]),
                layer=self.layers.General_labeling.number,
                datatype=self.layers.General_labeling.datatype,
            )
            marker_tr_side = gdspy.Rectangle(
                (xy[0], xy[1]),
                (xy[0] + 5, xy[1] + 5),
                layer=self.layers.General_labeling.number,
                datatype=self.layers.General_labeling.datatype,
            )
            self.add_fancy_text(
                f"RL{count}",
                xy[0],
                xy[1],
                5,
                self.layers.General_labeling,
                horizontal_align="end",
                vertical_align="above",
            )
            self.Main.add(marker_bl_side)
            self.Main.add(marker_tr_side)

        return running_list

    def add_center_pin_and_get_bounding_points(
        self,
        center_pin_xy: list[float | int],
        center_pin_radius: float | int,
    ) -> list[list[float | int]]:
        """Adds a center pin hole cutout to all layers and adds a pin hole
        positive circle to the pin hole positve later. This also returns a list
        of [x, y] points defining the boundary points for the pin.

        Adds a center pin cutout to all layers and adds the pin positive to a pin
        positives layer. Also returns a list of [x, y] points defining the boundary
        points for the pin.

        Parameters
        ----------
        center_pin_xy: list[float | int]
            list containing the [x,y] coordinate for the center of the pin.

        center_pin_radius: float, int
            Radius of the center pin.

        Returns
        -------
        bounding_points: list[list[float | int]]
            Returns a list containing [x, y] lists defining the boundary points of
            the center pin.
        """

        ground_plane_cutout_for_center_pin = gdspy.Round(
            center_pin_xy,
            center_pin_radius,
            tolerance=0.001,
            layer=self.layers.Nb_Groundplane.number,
            datatype=self.layers.Nb_Groundplane.datatype,
        )
        self.ground_plane_cutouts.add(ground_plane_cutout_for_center_pin)

        silicon_nitride_cutout_for_center_pin = gdspy.Round(
            center_pin_xy,
            center_pin_radius,
            tolerance=0.001,
            layer=self.layers.SiN_dep.number,
            datatype=self.layers.SiN_dep.datatype,
        )
        self.silicon_nitride_cutouts.add(silicon_nitride_cutout_for_center_pin)

        silicon_oxide_cutout_for_center_pin = gdspy.Round(
            center_pin_xy,
            center_pin_radius,
            tolerance=0.001,
            layer=self.layers.SiO.number,
            datatype=self.layers.SiO.datatype,
        )
        self.silicon_oxide_cutouts.add(silicon_oxide_cutout_for_center_pin)

        silicon_nitride_membrane_cutout_for_center_pin = gdspy.Round(
            center_pin_xy,
            center_pin_radius,
            tolerance=0.001,
            layer=self.layers.SiN_Membrane.number,
            datatype=self.layers.SiN_Membrane.datatype,
        )
        self.silicon_nitride_membrane_cutouts.add(silicon_nitride_membrane_cutout_for_center_pin)

        center_pin_positive = gdspy.Round(
            center_pin_xy,
            center_pin_radius,
            tolerance=0.001,
            layer=self.layers.Pin_hole_positives.number,
            datatype=self.layers.Pin_hole_positives.datatype,
        )
        self.Main.add(center_pin_positive)

        bounding_points = []
        for i in range(200):
            bounding_points.append(
                [center_pin_xy[0] + center_pin_radius * cos(i * center_pin_radius), center_pin_xy[1]]
                + center_pin_radius * sin(i * center_pin_radius)
            )

        return bounding_points

    def add_slotted_pin_and_get_bounding_points(
        self,
        slotted_pin_center: list[float | int],
        slotted_pin_length: float | int,
        slotted_pin_radius: float | int,
        center_pin_xy: list[float | int],
        points_in_curve: int = 300,
    ) -> list[list[float | int]]:
        """Adds a slotted pin cutout to all layers and adds the pin positive to
        a pin positives layer. Also returns a list of [x, y] points defining
        the boundary points for the pin. The slotted pin will always point
        towards the center of the chip where the center pin hole would be.

        Parameters
        ----------
        slotted_pin_center: list[float | int]
            list containing the [x, y] coordinate for the center of the slotted
            pin hole.

        slotted_pin_length: float, int
            The maximal length of the slot length.
                i.e. (2*slotted_pin_radius + central straight length).

        slotted_pin_radius: float, int
            Radius of the corners of the slotted pin.

        center_pin_xy: list[float | int]
            list containing the [x, y] coordinate for the center pin on the chip.
            This is used to ensure the slotted pin is pointing in the direction of
            the center pin hole.

        KwArgs
        ------
        points_in_curve: int = 300
            The number of points in each curved section of the slotted pin.

        Returns
        -------
        bounding_points: list[list[float | int]]
            Returns a list containing [x, y] lists defining the boundary points of
            the slotted pin.
        """

        slotted_pin_angle = np.arctan2((slotted_pin_center[1] - center_pin_xy[1]), (slotted_pin_center[0] - center_pin_xy[0]))

        mbu.rad_to_deg(slotted_pin_angle)

        end_points_slot = [
            [
                slotted_pin_center[0] + slotted_pin_length * cos(slotted_pin_angle),
                slotted_pin_center[1] + slotted_pin_length * sin(slotted_pin_angle),
            ],
            [
                slotted_pin_center[0] + slotted_pin_length * cos(slotted_pin_angle + pi),
                slotted_pin_center[1] + slotted_pin_length * sin(slotted_pin_angle + pi),
            ],
        ]

        bounding_points = []

        # ang_diff = pi / 100
        # loop_angle = slotted_pin_angle + pi / 2

        start_loop_angle_end_0 = slotted_pin_angle + pi / 2
        end_loop_angle_end_0 = start_loop_angle_end_0 - pi

        start_loop_angle_end_1 = end_loop_angle_end_0
        end_loop_angle_end_1 = start_loop_angle_end_1 - pi

        loop_angles_end_0 = np.linspace(
            start_loop_angle_end_0,
            end_loop_angle_end_0,
            points_in_curve,
        )
        loop_angles_end_1 = np.linspace(
            start_loop_angle_end_1,
            end_loop_angle_end_1,
            points_in_curve,
        )

        for loop_angle in loop_angles_end_0:
            bounding_points.append(
                [
                    end_points_slot[0][0] + slotted_pin_radius * cos(loop_angle),
                    end_points_slot[0][1] + slotted_pin_radius * sin(loop_angle),
                ]
            )

        for loop_angle in loop_angles_end_1:
            bounding_points.append(
                [
                    end_points_slot[1][0] + slotted_pin_radius * cos(loop_angle),
                    end_points_slot[1][1] + slotted_pin_radius * sin(loop_angle),
                ]
            )

        ground_plane_cutout_for_slotted_pin = gdspy.Polygon(
            bounding_points,
            layer=self.layers.Nb_Groundplane.number,
            datatype=self.layers.Nb_Groundplane.datatype,
        )
        self.ground_plane_cutouts.add(ground_plane_cutout_for_slotted_pin)

        silicon_nitride_cutout_for_slotted_pin = gdspy.Polygon(
            bounding_points,
            layer=self.layers.SiN_dep.number,
            datatype=self.layers.SiN_dep.datatype,
        )
        self.silicon_nitride_cutouts.add(silicon_nitride_cutout_for_slotted_pin)

        silicon_oxide_cutout_for_slotted_pin = gdspy.Polygon(
            bounding_points,
            layer=self.layers.SiO.number,
            datatype=self.layers.SiO.datatype,
        )
        self.silicon_oxide_cutouts.add(silicon_oxide_cutout_for_slotted_pin)

        silicon_nitride_membrane_cutout_for_slotted_pin = gdspy.Polygon(
            bounding_points,
            layer=self.layers.SiN_Membrane.number,
            datatype=self.layers.SiN_Membrane.datatype,
        )
        self.silicon_nitride_membrane_cutouts.add(silicon_nitride_membrane_cutout_for_slotted_pin)

        slotted_pin_positive = gdspy.Polygon(
            bounding_points,
            layer=self.layers.Pin_hole_positives.number,
            datatype=self.layers.Pin_hole_positives.datatype,
        )
        self.Main.add(slotted_pin_positive)

        return bounding_points

    def add_top_choke_features(
        self,
        x: float | int,
        y: float | int,
        top_choke_config_override: dict[str, float | int] | None = None,
        **kwargs,
    ) -> None:
        """Adds the top choke anulus and waveguide hole to the mask at the xy
        given. These are added on seperate layers. The parameters controling
        the dimesions of these are defined in the config.

        Parameters
        ----------
        x, y: float, int
            the x, y position to place the top choke feature

        KwArgs
        ------
        top_choke_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.
        """
        top_choke_config = get_mask_default_config(SoukMaskConfig.TOP_CHOKE, config_override=top_choke_config_override)

        waveguide_hole_radius = top_choke_config["waveguide_hole_radius"]  # 1200
        anulus_width = top_choke_config["anulus_width"]  # 600

        aunulus_outer_radius = waveguide_hole_radius + anulus_width
        aunulus_inner_radius = waveguide_hole_radius

        # Making the top choke waveguide hole
        top_choke_waveguide_hole = gdspy.Round(
            [x, y],
            waveguide_hole_radius,
            layer=self.layers.Top_choke_waveguide_hole.number,
            datatype=self.layers.Top_choke_waveguide_hole.datatype,
        )
        self.Main.add(top_choke_waveguide_hole)

        # Making the top choke anulus
        top_choke_anulus = gdspy.Round(
            [x, y],
            aunulus_outer_radius,
            inner_radius=aunulus_inner_radius,
            layer=self.layers.Top_choke_anulus.number,
            datatype=self.layers.Top_choke_anulus.datatype,
        )
        self.Main.add(top_choke_anulus)

        return

    def add_bottom_choke_features(
        self,
        x: float | int,
        y: float | int,
        rel_kid_positions: Sequence[Sequence[float | int]],
        resonator_types: Sequence[SoukResonatorType],
        bottom_choke_config_override: dict[str, float | int] | None = None,
        resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None],
        **kwargs,
    ) -> None:
        """At the xy given, adds the bottom choke waveguide hole and creates
        IDC cutout holes. Also adds pads to support the wafer. These are added
        on seperate layers. The parameters controling the dimesions of these
        are defined in the config.

        Parameters
        ----------
        x: float | int
            The x position to place the bottom choke features.

        y: float | int
            The y position to place the bottom choke features.

        rel_kid_positions: Sequence[Sequence[float | int]]
            Sequence of [x, y] lists that describe where the base of each KIDs meander
            is placed relative to the center of the antenna.
            Expected order is top_left, top_right, bot_left, bot_right

        resonator_types: Sequence[SoukResonatorType]
            This is the type of resonators drawn. The values accepted here are
            members of the SoukResonatorType enum. The order of the values
            passed in will be attributed to each KID and should be the same
            order as the rel_kid_positions, TL, TR, BL,
            BR.

        KwArgs
        ------
        bottom_choke_config_override: dict[str, float | int] | None = None
            This is an optional override dictionary containing key value pairs for
            variable name and that variable's value respectively. Any keys required
            that do not exist in this dict will be got from the default config. If
            extra keys that are not expected are provided a warnimg will be printed
            but nothing is done with those.

        resonator_config_overrides: Sequence[dict[str, float | int] | None] = [None, None, None, None]
            This is a Sequence of 4 optional override dictionarys containing key
            value pairs for variable name and that variable's value respectively.
            Any keys required that do not exist in this dict will be got from the
            default config. If extra keys that are not expected are provided a
            warning will be printed but nothing is done with those.
        """
        smc.add_bottom_choke_wave_guide_hole(
            self,
            x,
            y,
            bottom_choke_config_override=bottom_choke_config_override,
        )
        smc.add_bottom_choke_pads(
            self,
            x,
            y,
            bottom_choke_config_override=bottom_choke_config_override,
        )
        smc.add_bottom_choke_IDC_holes(
            self,
            x,
            y,
            rel_kid_positions,
            resonator_types,
            bottom_choke_config_override=bottom_choke_config_override,
            resonator_config_overrides=resonator_config_overrides,
        )
        smc.add_bottom_choke_backshort_hole(
            self,
            x,
            y,
            rel_kid_positions,
            resonator_types,
            bottom_choke_config_override=bottom_choke_config_override,
            resonator_config_overrides=resonator_config_overrides,
        )

    def add_hexagonal_tabbed_dicing_line(
        self,
        hexagon_points: list[float | int],
        hex_rad: float | int,
        chip_center_xy: list[float | int],
        all_sma_conectors_xy_rot: list[list[float | int]],
    ) -> None:
        """Adds a tabbed dicing line around the SOUK boundary.

        Parameters
        ----------
        hexagon_points: list[float | int]
            list of [x,y] lists that define the corner coordinates of the
            hexagonal boundary to create dicing line for.

        hex_rad: float, int
            The radius of the hexagon.

        chip_center_xy: list[float | int]
            list of [x,y] coordinates of the center of the chip.

        all_sma_conectors_xy_rot: list[list[float | int]]
            list of [x, y, rot] lists for all the sma connectors on the mask.
        """

        diceline_linewidth = 300
        diceline_end_extend_length = 500

        diceline_tab_length = 4000
        diceline_tab_gap = 700

        # diceline_tab_gap+=diceline_linewidth

        dice_x_offset_from_center = 4500
        dice_y_offset_from_center = 2500

        test_cell = gdspy.Cell("smatest")

        for sma_xy_rot in all_sma_conectors_xy_rot:
            x = sma_xy_rot[0]
            y = sma_xy_rot[1]
            rot = sma_xy_rot[2]

            main_dice_line_center = [x + dice_x_offset_from_center * cos(rot), y + dice_x_offset_from_center * sin(rot)]

            main_dice_line_start_end = [
                [
                    main_dice_line_center[0] + dice_y_offset_from_center * cos(rot + pi / 2),
                    main_dice_line_center[1] + dice_y_offset_from_center * sin(rot + pi / 2),
                ],
                [
                    main_dice_line_center[0] + dice_y_offset_from_center * cos(rot - pi / 2),
                    main_dice_line_center[1] + dice_y_offset_from_center * sin(rot - pi / 2),
                ],
            ]

            dice_line = gdspy.FlexPath(
                main_dice_line_start_end,
                diceline_linewidth,
                offset=(diceline_linewidth / 2),
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            # ####self.Main.add(dice_line)

            test_box = gdspy.Rectangle(
                [x + dice_x_offset_from_center, y + dice_y_offset_from_center],
                [x - 4 * dice_x_offset_from_center, y - dice_y_offset_from_center],
                layer=5000,
            )
            test_box.rotate(rot, [x, y])
            # self.Main.add(test_box)
            test_cell.add(test_box)
        hexagon_poly = gdspy.Polygon(hexagon_points, layer=5000)
        new_outline = gdspy.boolean(hexagon_poly, test_cell.get_polygons([5000, 0]), "not")
        self.Main.add(new_outline)

        new_outline_points = np.roll(new_outline.polygons[0], 2, axis=0)

        solid_hex_line_path = gdspy.FlexPath(
            new_outline_points,
            diceline_linewidth,
            offset=-diceline_linewidth / 2,
            ends="flush",
            layer=self.layers.Tab_dicing_line.number,
            datatype=self.layers.Tab_dicing_line.datatype,
        )
        solid_hex_line_polys = mbu.get_polys_from_flexpath(solid_hex_line_path)
        for poly_points in solid_hex_line_polys:
            solid_hex_line_poly = gdspy.Polygon(
                poly_points,
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            self.Main.add(solid_hex_line_poly)

        left_end = [
            new_outline_points[-1][0] + (diceline_linewidth / 2) * cos(-5 * pi / 6),
            new_outline_points[-1][1] + (diceline_linewidth / 2) * sin(-5 * pi / 6),
        ]

        right_end = [
            new_outline_points[0][0] + (diceline_linewidth / 2) * cos(-5 * pi / 6),
            new_outline_points[0][1] + (diceline_linewidth / 2) * sin(-5 * pi / 6),
        ]

        length_between_ends = ((left_end[0] - right_end[0]) ** 2 + (left_end[1] - right_end[1]) ** 2) ** 0.5

        center_point = [(left_end[0] + right_end[0]) / 2, (left_end[1] + right_end[1]) / 2]

        side_angle = -pi / 3

        middle_tab = gdspy.FlexPath(
            [[center_point[0] - diceline_tab_length / 2, center_point[1]], [center_point[0] + diceline_tab_length / 2, center_point[1]]],
            diceline_linewidth,
            ends="round",
            layer=self.layers.Tab_dicing_line.number,
            datatype=self.layers.Tab_dicing_line.datatype,
        )

        middle_tab.rotate(side_angle, center_point)
        self.Main.add(middle_tab)

        for i in range(int(length_between_ends / (diceline_tab_length + diceline_tab_gap + diceline_linewidth))):
            if (diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap)) > length_between_ends / 2:
                tab_right = gdspy.FlexPath(
                    [
                        [
                            center_point[0] - diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                            center_point[1],
                        ],
                        [center_point[0] + length_between_ends / 2 + diceline_end_extend_length, center_point[1]],
                    ],
                    diceline_linewidth,
                    ends="flush",
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                tab_right.rotate(side_angle, center_point)
                self.Main.add(tab_right)

                cap_right = gdspy.Round(
                    [
                        center_point[0] - diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                    diceline_linewidth / 2,
                    initial_angle=pi / 2,
                    final_angle=3 * pi / 2,
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                cap_right.rotate(side_angle, center_point)
                self.Main.add(cap_right)

                tab_left = gdspy.FlexPath(
                    [
                        [center_point[0] - length_between_ends / 2 - diceline_end_extend_length, center_point[1]],
                        [
                            center_point[0] + diceline_tab_length / 2 - i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                            center_point[1],
                        ],
                    ],
                    diceline_linewidth,
                    ends="flush",
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                tab_left.rotate(side_angle, center_point)
                self.Main.add(tab_left)

                cap_left = gdspy.Round(
                    [
                        center_point[0] + diceline_tab_length / 2 - i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                    diceline_linewidth / 2,
                    initial_angle=pi / 2,
                    final_angle=-pi / 2,
                    layer=self.layers.Tab_dicing_line.number,
                    datatype=self.layers.Tab_dicing_line.datatype,
                )
                cap_left.rotate(side_angle, center_point)
                self.Main.add(cap_left)

                break

            tab_right = gdspy.FlexPath(
                [
                    [
                        center_point[0] - diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                    [
                        center_point[0] + diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                ],
                diceline_linewidth,
                ends="round",
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            tab_right.rotate(side_angle, center_point)
            self.Main.add(tab_right)

            tab_left = gdspy.FlexPath(
                [
                    [
                        center_point[0] - diceline_tab_length / 2 - i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                    [
                        center_point[0] + diceline_tab_length / 2 - i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                ],
                diceline_linewidth,
                ends="round",
                layer=self.layers.Tab_dicing_line.number,
                datatype=self.layers.Tab_dicing_line.datatype,
            )
            tab_left.rotate(side_angle, center_point)
            self.Main.add(tab_left)

        return

    def add_hexagonal_bottom_choke_tabbed_dicing_line(
        self,
        hexagon_points: list[list[float | int]],
        hex_rad: float | int,
        chip_center_xy: list[float | int],
        all_sma_conectors_xy_rot: list[float | int],
        bend_rad_corner: float | int,
    ) -> None:
        """Adds a Bottom choke tabbed dicing line around the SOUK boundary.

        Parameters
        ----------
        hexagon_points: list
            list of [x,y] lists that define the corner coordinates of the
            hexagonal boundary to create dicing line for.

        hex_rad: float, int
            The radius of the hexagon.

        chip_center_xy: list
            list of [x,y] coordinates of the center of the chip.

        all_sma_conectors_xy_rot: list
            list of [x,y,rot] lists for all the sma connectors on the mask.

        bend_rad_corner: float, int
            The bend radius for the outer corners
        """

        diceline_linewidth = 300
        diceline_end_extend_length = 500

        diceline_tab_length = 4000
        diceline_tab_gap = 700

        # diceline_tab_gap+=diceline_linewidth

        # dice_x_offset_from_center = 8800
        dice_x_offset_from_center = 9500
        dice_y_offset_from_center = 2500

        test_cell = gdspy.Cell("smatest")

        for sma_xy_rot in all_sma_conectors_xy_rot:
            x = sma_xy_rot[0]
            y = sma_xy_rot[1]
            rot = sma_xy_rot[2]

            main_dice_line_center = [x + dice_x_offset_from_center * cos(rot), y + dice_x_offset_from_center * sin(rot)]

            main_dice_line_start_end = [
                [
                    main_dice_line_center[0] + dice_y_offset_from_center * cos(rot + pi / 2),
                    main_dice_line_center[1] + dice_y_offset_from_center * sin(rot + pi / 2),
                ],
                [
                    main_dice_line_center[0] + dice_y_offset_from_center * cos(rot - pi / 2),
                    main_dice_line_center[1] + dice_y_offset_from_center * sin(rot - pi / 2),
                ],
            ]

            # dice_line = gdspy.FlexPath(
            #     main_dice_line_start_end,
            #     diceline_linewidth,
            #     offset=diceline_linewidth / 2,
            #     layer=self.Bottom_choke_Tab_dicing_line["layer"],
            #     datatype=self.Bottom_choke_Tab_dicing_line["datatype"],
            # )
            # ####self.Main.add(dice_line)

            test_box = gdspy.Rectangle(
                [x + dice_x_offset_from_center, y + dice_y_offset_from_center],
                [x - 4 * dice_x_offset_from_center, y - dice_y_offset_from_center],
                layer=5000,
            )
            test_box.rotate(rot, [x, y])
            # self.Main.add(test_box)
            test_cell.add(test_box)
        hexagon_poly = gdspy.Polygon(hexagon_points, layer=5000)
        new_outline = gdspy.boolean(hexagon_poly, test_cell.get_polygons([5000, 0]), "not")
        self.Main.add(new_outline)

        new_outline_points = np.roll(new_outline.polygons[0], 2, axis=0)

        hex_line_part1_points = new_outline_points[1:11]
        hex_line_part1 = gdspy.FlexPath(
            hex_line_part1_points,
            diceline_linewidth,
            offset=-diceline_linewidth / 2,
            ends="flush",
            corners="natural",
            layer=self.layers.Bottom_choke_Tab_dicing_line.number,
            datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
        )
        self.Main.add(hex_line_part1)

        hex_line_part2_points = new_outline_points[10:13]
        hex_line_part2 = gdspy.FlexPath(
            hex_line_part2_points,
            diceline_linewidth,
            offset=-diceline_linewidth / 2,
            ends="extended",
            corners="circular bend",
            bend_radius=bend_rad_corner,
            layer=self.layers.Bottom_choke_Tab_dicing_line.number,
            datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
        )
        self.Main.add(hex_line_part2)

        hex_line_part3_points = new_outline_points[12:15]
        hex_line_part3 = gdspy.FlexPath(
            hex_line_part3_points,
            diceline_linewidth,
            offset=-diceline_linewidth / 2,
            ends="flush",
            corners="natural",
            layer=self.layers.Bottom_choke_Tab_dicing_line.number,
            datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
        )
        self.Main.add(hex_line_part3)

        left_end = [
            new_outline_points[-1][0] + (diceline_linewidth / 2) * cos(-5 * pi / 6),
            new_outline_points[-1][1] + (diceline_linewidth / 2) * sin(-5 * pi / 6),
        ]

        right_end = [
            new_outline_points[0][0] + (diceline_linewidth / 2) * cos(-5 * pi / 6),
            new_outline_points[0][1] + (diceline_linewidth / 2) * sin(-5 * pi / 6),
        ]

        length_between_ends = ((left_end[0] - right_end[0]) ** 2 + (left_end[1] - right_end[1]) ** 2) ** 0.5

        center_point = [(left_end[0] + right_end[0]) / 2, (left_end[1] + right_end[1]) / 2]

        side_angle = -pi / 3

        middle_tab = gdspy.FlexPath(
            [[center_point[0] - diceline_tab_length / 2, center_point[1]], [center_point[0] + diceline_tab_length / 2, center_point[1]]],
            diceline_linewidth,
            ends="round",
            layer=self.layers.Bottom_choke_Tab_dicing_line.number,
            datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
        )

        middle_tab.rotate(side_angle, center_point)
        self.Main.add(middle_tab)

        for i in range(int(length_between_ends / (diceline_tab_length + diceline_tab_gap + diceline_linewidth))):
            if (diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap)) > length_between_ends / 2:
                tab_right = gdspy.FlexPath(
                    [
                        list(tab_right.points[0]),
                        [new_outline_points[0][0], new_outline_points[0][1] - (diceline_linewidth / 2)],
                        [new_outline_points[1][0], new_outline_points[1][1] - (diceline_linewidth / 2)],
                    ],
                    diceline_linewidth,
                    ends="flush",
                    corners="circular bend",
                    bend_radius=bend_rad_corner,
                    layer=self.layers.Bottom_choke_Tab_dicing_line.number,
                    datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
                )

                self.Main.add(tab_right)

                # cap_right = gdspy.Round([center_point[0] - diceline_tab_length/2 + i*(diceline_tab_length + diceline_tab_gap + diceline_linewidth), center_point[1]],
                #                         diceline_linewidth/2, initial_angle=pi/2, final_angle=3*pi/2, **self.Bottom_choke_Tab_dicing_line)
                # cap_right.rotate(side_angle, center_point)
                # self.Main.add(cap_right)

                tab_left = gdspy.FlexPath(
                    [
                        list(tab_left.points[0]),
                        [new_outline_points[-2][0], new_outline_points[-2][1] + (diceline_linewidth / 2)],
                        [new_outline_points[-3][0], new_outline_points[-3][1] + (diceline_linewidth / 2)],
                    ],
                    diceline_linewidth,
                    ends="flush",
                    corners="circular bend",
                    bend_radius=200,
                    layer=self.layers.Bottom_choke_Tab_dicing_line.number,
                    datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
                )

                self.Main.add(tab_left)

                # cap_left = gdspy.Round([center_point[0] + diceline_tab_length/2 - i*(diceline_tab_length + diceline_tab_gap + diceline_linewidth), center_point[1]],
                #                         diceline_linewidth/2, initial_angle=pi/2, final_angle=-pi/2, **self.Bottom_choke_Tab_dicing_line)
                # cap_left.rotate(side_angle, center_point)
                # self.Main.add(cap_left)

                break

            tab_right = gdspy.FlexPath(
                [
                    [
                        center_point[0] - diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                    [
                        center_point[0] + diceline_tab_length / 2 + i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                ],
                diceline_linewidth,
                ends="round",
                layer=self.layers.Bottom_choke_Tab_dicing_line.number,
                datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
            )
            tab_right.rotate(side_angle, center_point)
            self.Main.add(tab_right)

            tab_left = gdspy.FlexPath(
                [
                    [
                        center_point[0] - diceline_tab_length / 2 - i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                    [
                        center_point[0] + diceline_tab_length / 2 - i * (diceline_tab_length + diceline_tab_gap + diceline_linewidth),
                        center_point[1],
                    ],
                ],
                diceline_linewidth,
                ends="round",
                layer=self.layers.Bottom_choke_Tab_dicing_line.number,
                datatype=self.layers.Bottom_choke_Tab_dicing_line.datatype,
            )
            tab_left.rotate(side_angle, center_point)
            self.Main.add(tab_left)

        return

    def add_all_positive_layers_to_mask(self) -> None:
        """This adds all the positive cells to the mask under new layer names.

        This is the ground_plane, silicon_nitride, silicon_oxide and the
        silicon_nitride_membrane.
        """

        raise RuntimeError("This needs fixing.")

        positive_layers_dict = {
            "silicon_nitride_positives": {"layer": 2003, "datatype": 0},
            "silicon_oxide_positives": {"layer": 2005, "datatype": 0},
            "silicon_nitride_membrane_positives": {"layer": 2006, "datatype": 0},
            "ground_plane_positives": {"layer": 2030, "datatype": 0},
        }

        self.all_layers_name_lookup_from_number.update(
            {
                positive_layers_dict["silicon_nitride_positives"]["layer"]: "silicon_nitride_positives",
                positive_layers_dict["silicon_oxide_positives"]["layer"]: "silicon_oxide_positives",
                positive_layers_dict["silicon_nitride_membrane_positives"]["layer"]: "silicon_nitride_membrane_positives",
                positive_layers_dict["ground_plane_positives"]["layer"]: "ground_plane_positives",
            }
        )

        silicon_nitride_positive_polygons = self.silicon_nitride_positives.get_polygons()
        for poly_points in silicon_nitride_positive_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=positive_layers_dict["silicon_nitride_positives"]["layer"],
                datatype=positive_layers_dict["silicon_nitride_positives"]["datatype"],
            )
            self.Main.add(poly)

        silicon_oxide_positive_polygons = self.silicon_oxide_positives.get_polygons()
        for poly_points in silicon_oxide_positive_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=positive_layers_dict["silicon_oxide_positives"]["layer"],
                datatype=positive_layers_dict["silicon_oxide_positives"]["datatype"],
            )
            self.Main.add(poly)

        silicon_nitride_membrane_positive_polygons = self.silicon_nitride_membrane_positives.get_polygons()
        for poly_points in silicon_nitride_membrane_positive_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=positive_layers_dict["silicon_nitride_membrane_positives"]["layer"],
                datatype=positive_layers_dict["silicon_nitride_membrane_positives"]["datatype"],
            )
            self.Main.add(poly)

        ground_plane_positive_polygons = self.ground_plane_positives.get_polygons()
        for poly_points in ground_plane_positive_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=positive_layers_dict["ground_plane_positives"]["layer"],
                datatype=positive_layers_dict["ground_plane_positives"]["datatype"],
            )
            self.Main.add(poly)

        return

    def add_all_cutout_layers_to_mask(self) -> None:
        """This adds all the cutout cells to the mask under new layer names.

        This is the ground_plane, silicon_nitride, silicon_oxide and the
        silicon_nitride_membrane.
        """

        raise RuntimeError("This needs fixing.")

        cutout_layers_dict = {
            "silicon_nitride_cutouts": {"layer": 3003, "datatype": 0},
            "silicon_oxide_cutouts": {"layer": 3005, "datatype": 0},
            "silicon_nitride_membrane_cutouts": {"layer": 3006, "datatype": 0},
            "ground_plane_cutouts": {"layer": 3030, "datatype": 0},
        }

        self.all_layers_name_lookup_from_number.update(
            {
                cutout_layers_dict["silicon_nitride_cutouts"]["layer"]: "silicon_nitride_cutouts",
                cutout_layers_dict["silicon_oxide_cutouts"]["layer"]: "silicon_oxide_cutouts",
                cutout_layers_dict["silicon_nitride_membrane_cutouts"]["layer"]: "silicon_nitride_membrane_cutouts",
                cutout_layers_dict["ground_plane_cutouts"]["layer"]: "ground_plane_cutouts",
            }
        )

        silicon_nitride_cutout_polygons = self.silicon_nitride_cutouts.get_polygons()
        for poly_points in silicon_nitride_cutout_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=cutout_layers_dict["silicon_nitride_cutouts"]["layer"],
                datatype=cutout_layers_dict["silicon_nitride_cutouts"]["datatype"],
            )
            self.Main.add(poly)

        silicon_oxide_cutout_polygons = self.silicon_oxide_cutouts.get_polygons()
        for poly_points in silicon_oxide_cutout_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=cutout_layers_dict["silicon_oxide_cutouts"]["layer"],
                datatype=cutout_layers_dict["silicon_oxide_cutouts"]["datatype"],
            )
            self.Main.add(poly)

        silicon_nitride_membrane_cutout_polygons = self.silicon_nitride_membrane_cutouts.get_polygons()
        for poly_points in silicon_nitride_membrane_cutout_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=cutout_layers_dict["silicon_nitride_membrane_cutouts"]["layer"],
                datatype=cutout_layers_dict["silicon_nitride_membrane_cutouts"]["datatype"],
            )
            self.Main.add(poly)

        ground_plane_cutout_polygons = self.ground_plane_cutouts.get_polygons()
        for poly_points in ground_plane_cutout_polygons:
            poly = gdspy.Polygon(
                poly_points,
                layer=cutout_layers_dict["ground_plane_cutouts"]["layer"],
                datatype=cutout_layers_dict["ground_plane_cutouts"]["datatype"],
            )
            self.Main.add(poly)

        return

    def add_Aluminium_Bulk_around_all_Aluminium(self, region_to_avoid=[]):
        """Adds a bulk aluminium block bounding box around all the aluminium in
        the Main cell of the mask. This is such that the aluminium can be
        etched out of the bulk aluminium layer.

        KwArgs
        -------
        regoin_to_avoid: list
            This by default is an empty list. When defined this should be a list
            of [x,y] list points defining a polygon (similar to defining a gdspy
            polygon) that should be enclose a regoin where no bulk aluminium
            should be added. e.g. [[x1,y1], [x2,y2], [x3,y3], [x4,y4]].
            This can be multiple different regions as seperate lists. e.g.
            [ [[x1,y1], [x2,y2] .. [xN,yN]], [[x1,y1], [x2,y2] .. [xN,yN]] ].
        """
        # TODO add a region to only copy and not expand the bulk aluminium added.

        Al_polys = self.Main.get_polygons([self.layers.Aluminium.number, self.layers.Aluminium.datatype])

        bulk_aluminium_positves = gdspy.Cell("BULK_ALUMINIUM_POSITVES")
        bulk_aluminium_negatives = gdspy.Cell("BULK_ALUMINIUM_NEGATIVES")

        for poly in Al_polys:
            new_poly = gdspy.Polygon(
                poly,
                layer=self.layers.Aluminium_Bulk.number,
                datatype=self.layers.Aluminium_Bulk.datatype,
            )
            [BB_xmin, BB_ymin], [BB_xmax, BB_ymax] = new_poly.get_bounding_box()

            y_offset_top = 100
            y_offset_bot = 100
            x_offset_left = 10
            x_offset_right = 10

            oversized_rectangle = gdspy.Rectangle(
                [BB_xmin - x_offset_left, BB_ymin - y_offset_bot],
                [BB_xmax + x_offset_right, BB_ymax + y_offset_top],
                layer=self.layers.Aluminium_Bulk.number,
                datatype=self.layers.Aluminium_Bulk.datatype,
            )
            bulk_aluminium_positves.add(oversized_rectangle)

        if len(region_to_avoid) != 0:
            if isinstance(region_to_avoid[0][0], list):  # if this element is list then there are multiple regions to avoid.
                for region in region_to_avoid:
                    negative_poly = gdspy.Polygon(
                        region,
                        layer=self.layers.Aluminium_Bulk.number,
                        datatype=self.layers.Aluminium_Bulk.datatype,
                    )
                    bulk_aluminium_negatives.add(negative_poly)
            elif isinstance(region_to_avoid[0][0], (int, float)):
                negative_poly = gdspy.Polygon(
                    region_to_avoid,
                    layer=self.layers.Aluminium_Bulk.number,
                    datatype=self.layers.Aluminium_Bulk.datatype,
                )
                bulk_aluminium_negatives.add(negative_poly)
            else:
                raise TypeError("Incorrect format for region_to_avoid")

        bulk_aluminium = gdspy.boolean(
            bulk_aluminium_positves.get_polygons([self.layers.Aluminium_Bulk.number, self.layers.Aluminium_Bulk.datatype]),
            bulk_aluminium_negatives,
            "not",
            precision=0.01,
            layer=self.layers.Aluminium_Bulk.number,
            datatype=self.layers.Aluminium_Bulk.datatype,
        )
        self.Main.add(bulk_aluminium)

        return

    def add_Aluminium_Patch_and_Etch(
        self,
        patch_size_diff: float = 5.0,
        etch_size_diff: float = 8.0,
        oversize_aluminium_size_diff: float = 1.0,
        alternate_cell: gdspy.Cell | None = None,
        alternate_patch_layer: Layer | None = None,
        alternate_etch_layer: Layer | None = None,
        printing: bool = True,
        precision: float = 0.001,
        include_existing_patch_etch_polys: bool = True,
    ) -> None:
        """Take all the Aluminium on the mask and create a patch and etch layer
        for it. These layers are the aluminium polygons on the mask sized by
        the respective size_diff arguments given and added to the mask. The
        etch layer can also have the aluminium cutout oversized for CD
        tolerance.

        KwArgs
        ------
        patch_size_diff: float = 5.0
            The size diff for the patch.

        etch_size_diff: float = 8.0
            The size diff for the etch

        oversize_aluminium_size_diff: float = 1.0
            The oversize of the aluminium that will be cutout from the etch
            size.

        alternate_cell: gdspy.Cell | None = None
            This will be the cell that the geometry is added to in place of
            self.Main. If Not specified the geometry will get added to Main.

        alternate_patch_layer: Layer | None = None
            This is an instance of Layer. see maskpy.layers.Layer.
            Usually this is within the SoukMaskBuilder.layers.xxx.
            e.g. `self.layers.Aluminium`

        alternate_etch_layer: Layer | None = None
            This is an instance of Layer. see maskpy.layers.Layer.
            Usually this is within the SoukMaskBuilder.layers.xxx.
            e.g. `self.layers.Aluminium`

        include_existing_patch_etch_polys: bool = True
            This will combine the existing patch and etch polys already on the
            mask with the result generated from this method. Set this to False
            to ignore the patch and etch polygons that already exists in the
            mask_builder.layers.Aluminium_Patch and the
            mask_builder.layers.Aluminium_Etch layers.
        """

        if alternate_cell is None:
            cell = self.Main
        else:
            cell = alternate_cell

        if alternate_patch_layer is None:
            patch_layer = self.layers.Aluminium_Patch
        elif isinstance(alternate_patch_layer, Layer):
            patch_layer = alternate_patch_layer
        else:
            raise TypeError(f"alternate_patch_layer should be of type Layer, current type is {type(alternate_patch_layer)}")

        if alternate_etch_layer is None:
            etch_layer = self.layers.Aluminium_Etch
        elif isinstance(alternate_etch_layer, Layer):
            etch_layer = alternate_etch_layer
        else:
            raise TypeError(f"alternate_etch_layer should be of type Layer, current type is {type(alternate_etch_layer)}")

        all_aluminium_polys: list[list[list[float]]] = self.Main.get_polygons(
            [
                self.layers.Aluminium.number,
                self.layers.Aluminium.datatype,
            ]
        )

        all_nb_ant_polys: list[list[list[float]]] = self.Main.get_polygons(
            [
                self.layers.Nb_Antenna.number,
                self.layers.Nb_Antenna.datatype,
            ]
        )
        all_nb_idc_polys: list[list[list[float]]] = self.Main.get_polygons(
            [
                self.layers.IDC_Nb.number,
                self.layers.IDC_Nb.datatype,
            ]
        )

        # if no aluminium exists.
        if len(all_aluminium_polys) == 0:
            print("No Aluminium on the mask to add Patch and Etch for.")
            return

        if printing:
            print("Adding Aluminium Patch and Etch to the mask.")
            print("├─  Aluminium Patch")

        patch_positive_polys = gdspy.offset(all_aluminium_polys, patch_size_diff, join_first=True)

        if include_existing_patch_etch_polys:
            existing_patch_polys: list[list[list[float]]] = self.Main.get_polygons(
                [
                    self.layers.Aluminium_Patch.number,
                    self.layers.Aluminium_Patch.datatype,
                ]
            )
            patch_positive_polys = gdspy.boolean(
                patch_positive_polys,
                existing_patch_polys,
                "or",
                precision=precision,
                layer=patch_layer.number,
                datatype=patch_layer.datatype,
            )

        patch_negative_polys = None
        cell.add(
            gdspy.boolean(
                patch_positive_polys,
                patch_negative_polys,
                "not",
                precision=precision,
                layer=patch_layer.number,
                datatype=patch_layer.datatype,
            )
        )

        if printing:
            print("├─  Aluminium Etch")

        etch_positive_polys = gdspy.offset(all_aluminium_polys, etch_size_diff, join_first=True)

        # combine with aluminium etch positives if then exist
        if self.aluminium_etch_positives.get_polygons():
            etch_positive_polys = gdspy.boolean(
                etch_positive_polys,
                self.aluminium_etch_positives.get_polygons(),
                "or",
                precision=precision,
            )

        etch_negative_polys = gdspy.boolean(
            gdspy.boolean(all_nb_ant_polys, all_nb_idc_polys, "or", precision=precision),
            gdspy.offset(all_aluminium_polys, oversize_aluminium_size_diff, join_first=True),
            "or",
            precision=precision,
        )

        cell.add(
            gdspy.boolean(
                etch_positive_polys,
                etch_negative_polys,
                "not",
                precision=precision,
                layer=etch_layer.number,
                datatype=etch_layer.datatype,
            )
        )
        if printing:
            print("└─  Done\n")

    def add_Nb_Antenna_Patch_and_Etch(
        self,
        al_oversize_for_patch: float = 6.0,
        al_oversize_for_etch: float = 2.0,
        nb_oversize_for_patch: float = 1.0,
    ) -> None:
        """Take all the Nb_Antenna and Nb_IDC polys on the mask and create a
        patch and etch layer for it. These layers are the Nb_Antenna polygons
        on the mask sized by the respective size_diff arguments given and added
        to the mask. The etch layer can also have the Nb cutout oversized for
        CD tolerance.

        KwArgs
        ------
        al_oversize_for_patch: float = 6.0
            The size diff for the aluminium to make the patch

        al_oversize_for_etch: float = 2.0
            The size diff for the aluminium to make the etch

        nb_oversize_for_patch: float = 1.0
            The oversize of the niobium that will be cutout.
        """
        PRECISION = 0.001

        all_aluminium_polys: list[list[list[float]]] = self.Main.get_polygons(
            [self.layers.Aluminium.number, self.layers.Aluminium.datatype],
        )

        all_nb_ant_polys: list[list[list[float]]] = self.Main.get_polygons(
            [self.layers.Nb_Antenna.number, self.layers.Nb_Antenna.datatype],
        )
        all_nb_idc_polys: list[list[list[float]]] = self.Main.get_polygons(
            [self.layers.IDC_Nb.number, self.layers.IDC_Nb.datatype],
        )

        all_nb_polyset = gdspy.boolean(
            all_nb_ant_polys,
            all_nb_idc_polys,
            "or",
            precision=PRECISION,
        )

        print("Adding Niobium Patch and Etch to the mask.")
        print("├─  Nb_Patch")
        # oversize the AL for the patch
        # oversize the Nb on the mask for the patch
        # Cutout the oversized Nb on the mask from that oversized AL patch

        oversized_Al_for_patch = gdspy.offset(all_aluminium_polys, al_oversize_for_patch, join_first=True)
        oversized_Nb_for_patch = gdspy.offset(all_nb_polyset, nb_oversize_for_patch, join_first=True)

        patch_polyset_pre_positives = gdspy.boolean(
            oversized_Al_for_patch,
            oversized_Nb_for_patch,
            "not",
            precision=PRECISION,
            # layer=self.Nb_Patch["layer"],
            # datatype=self.Nb_Patch["datatype"],
        )
        patch_polyset = gdspy.boolean(
            patch_polyset_pre_positives,
            self.nb_patch_positives,
            "or",
            precision=PRECISION,
            layer=self.layers.Nb_Patch.number,
            datatype=self.layers.Nb_Patch.datatype,
        )

        self.Main.add(patch_polyset)

        print("├─  Nb_Etch")
        # oversize the AL for the etch (less oversize than the patch oversize)
        # etch is the Nb on the mask + oversized AL etch
        oversized_Al_for_etch = gdspy.offset(all_aluminium_polys, al_oversize_for_etch, join_first=True)

        etch_polyset = gdspy.boolean(
            oversized_Al_for_etch,
            all_nb_polyset,
            "or",
            precision=PRECISION,
            layer=self.layers.Nb_Etch.number,
            datatype=self.layers.Nb_Etch.datatype,
        )
        self.Main.add(etch_polyset)

        print("└─  Done\n")

    def expected_time_progress_bar(self, expected_time, progress_bar_title, steps=100):
        """Shows a progress bar based on the expected_time taken to complete an
        operation. This is intended to be used in a daemon thread and get
        alongside the other long running operation, e.g. the
        do_boolean_operations() method used this for each of the longer running
        boolean operation being completed.

        Parameters
        ----------
        expected_time: int
            The expected amount of time (**in seconds**) the operation will
            last and hence how long this progress bar will run for.

        progress_bar_title: str
            The title to put in the progress bar.

        KwArgs
        ------
        steps: int
            The number of steps to divide the progress into.
        """
        time_inc = expected_time / steps
        # with alive_bar(expected_time, manual=True, title=progress_bar_title, refresh_secs=math.ceil(time_inc), calibrate=50, stats=True, force_tty=True, length=20, max_cols=70) as bar:
        #     time_elapsed = 0
        #     while time_elapsed<expected_time:
        #         time.sleep(time_inc)
        #         time_elapsed += time_inc
        #         bar(time_elapsed/expected_time)

        for _ in tqdm(range(steps)):
            time.sleep(time_inc)

        return

    def _boolean_ground_plane(self, precision: float):
        start_time = time.time()
        groundplane = gdspy.boolean(
            self.ground_plane_positives.get_polygons(),
            self.ground_plane_cutouts,
            "not",
            precision=precision,
            layer=self.layers.Nb_Groundplane.number,
            datatype=self.layers.Nb_Groundplane.datatype,
        )

        if groundplane:
            self.Main.add(groundplane)

        time_taken = np.round(time.time() - start_time, 2)
        print("│   Ground Plane boolean operation DONE")
        print(f"│   time taken to do ground = {time_taken}\n│")
        return

    def _boolean_SiN(self, precision: float):
        start_time = time.time()
        dielectric_layer = gdspy.boolean(
            self.silicon_nitride_positives.get_polygons(),
            self.silicon_nitride_cutouts,
            "not",
            precision=precision,
            layer=self.layers.SiN_dep.number,
            datatype=self.layers.SiN_dep.datatype,
        )

        if dielectric_layer:
            self.Main.add(dielectric_layer)

        time_taken = np.round(time.time() - start_time, 2)
        print("│   dielectirc layer boolean operation DONE")
        print(f"│   time taken to do sin_dep = {time_taken}\n│")
        return

    def _boolean_SiO(self, precision: float):
        start_time = time.time()
        SiO_layer = gdspy.boolean(
            self.silicon_oxide_positives.get_polygons(),
            self.silicon_oxide_cutouts,
            "not",
            precision=precision,
            layer=self.layers.SiO.number,
            datatype=self.layers.SiO.datatype,
        )

        if SiO_layer:
            self.Main.add(SiO_layer)

        time_taken = np.round(time.time() - start_time, 2)
        print("│   Silicon DiOxide layer boolean operation DONE")
        print(f"│   time taken to do SiO = {time_taken}\n│")
        return

    def _boolean_SiN_mem(self, precision: float):
        start_time = time.time()
        SiN_membrane_layer = gdspy.boolean(
            self.silicon_nitride_membrane_positives.get_polygons(),
            self.silicon_nitride_membrane_cutouts,
            "not",
            precision=precision,
            layer=self.layers.SiN_Membrane.number,
            datatype=self.layers.SiN_Membrane.datatype,
        )

        if SiN_membrane_layer:
            self.Main.add(SiN_membrane_layer)

        time_taken = np.round(time.time() - start_time, 2)
        print("│   Silicon Nitride membrane layer boolean operation DONE")
        print(f"│   time taken to do SiN mem = {time_taken}\n│")

    def do_boolean_operations(self, precision: float = 0.001):
        """This does all the boolean operations on the the positive cells and
        cutout cells for the mask and adds these to the Main cell ready for
        writing the final gds file.

        This is the ground_plane, silicon_nitride, silicon_oxide,
        silicon_nitride_membrane.
        """
        if self.boolean_operations_completed_flag:
            return

        self.boolean_operations_completed_flag = True

        # TODO Make an alive loop for each of the steps in this process.
        # time is roughly proportioanl to number of polygons in the positives and cutouts combined.

        print("Doing the boolean operations to the mask.")
        very_first_time = time.time()

        # Groundplane
        combined_number_of_ground_plane_positives_cutouts = 0
        for li in self.ground_plane_positives.get_polygons():
            combined_number_of_ground_plane_positives_cutouts += len(li)

        for li in self.ground_plane_cutouts.get_polygons():
            combined_number_of_ground_plane_positives_cutouts += len(li)

        ground_plane_expected_time = combined_number_of_ground_plane_positives_cutouts * 4.94628e-04

        # SiNDep
        combined_number_of_SiN_dep_positives_cutouts = 0
        for li in self.silicon_nitride_positives.get_polygons():
            combined_number_of_SiN_dep_positives_cutouts += len(li)

        for li in self.silicon_nitride_cutouts.get_polygons():
            combined_number_of_SiN_dep_positives_cutouts += len(li)

        SiN_expected_time = combined_number_of_SiN_dep_positives_cutouts * 2.955534e-09

        # SiOLayer
        combined_number_of_SiO_positives_cutouts = 0
        for li in self.silicon_oxide_positives.get_polygons():
            combined_number_of_SiO_positives_cutouts += len(li)

        for li in self.silicon_oxide_cutouts.get_polygons():
            combined_number_of_SiO_positives_cutouts += len(li)

        SiO_expected_time = combined_number_of_SiO_positives_cutouts * 2.73806e-05

        # SiNmembrane
        combined_number_of_SiN_mem_positives_cutouts = 0
        for li in self.silicon_nitride_positives.get_polygons():
            combined_number_of_SiN_mem_positives_cutouts += len(li)

        for li in self.silicon_oxide_cutouts.get_polygons():
            combined_number_of_SiN_mem_positives_cutouts += len(li)

        SiN_mem_expected_time = combined_number_of_SiN_mem_positives_cutouts * 1.92823e-05

        # longest_expected_time = max([ground_plane_expected_time, SiN_mem_expected_time, SiO_expected_time, SiN_mem_expected_time])
        # expected_end = datetime.datetime.fromtimestamp( very_first_time + int(longest_expected_time))
        # print(f"    Expected time - {np.round(longest_expected_time, 2)} secs")
        # print(f"    Expected end  - {expected_end.hour}:{expected_end.minute}\n")

        sequential_time = ground_plane_expected_time + SiN_mem_expected_time + SiO_expected_time + SiN_mem_expected_time
        expected_end = datetime.datetime.fromtimestamp(very_first_time + int(sequential_time))
        print(f"│   Expected time taken: {np.round(sequential_time, 2)} secs")
        print(f"│   Expected end time  : {expected_end.hour}:{str(expected_end.minute).zfill(2)}\n│")

        print(f"├─  SiN ~ {np.round(SiN_mem_expected_time, 2)} secs")
        self._boolean_SiN(precision)

        print(f"├─  SiO ~ {np.round(SiO_expected_time, 2)} secs")
        self._boolean_SiO(precision)

        print(f"├─  SiN_mem ~ {np.round(SiN_expected_time, 2)} secs")
        self._boolean_SiN_mem(precision)

        print(f"├─  Grnd ~ {np.round(ground_plane_expected_time, 2)} secs")
        self._boolean_ground_plane(precision)

        # IN TESTING MAKING THE BOOL OPERATIONS CONCURRENT TO SPEED UP.

        # with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        #     bool_ground_done = executor.submit(
        #         self.boolean_ground_plane
        #     )  # , desc="Ground Plane Boolean Operation", total=ground_plane_expected_time
        #     bool_SiN_done = executor.submit(self.boolean_SiN)  # , desc="Silicon Nitride Boolean Operation", total=SiN_expected_time
        #     bool_SiO_done = executor.submit(self.boolean_SiO)  # , desc="Silicon DiOxide Boolean Operation", total=SiO_expected_time
        #     Boolean_SiN_mem_done = executor.submit(
        #         self.boolean_SiN_mem
        #     )  # , desc="Silicon Nitride membrane Boolean Operation", total=SiN_mem_expected_time

        very_last_time = time.time()
        print("└─  The whole time taken overall = " + str(very_last_time - very_first_time) + "\n\n\n")

        return

    def save_layer_props_ANL_version(self, filename: str, ANL_Main: gdspy.Cell) -> None:
        """Save a layer properties file for the mask under the name given. This
        should be done after do_boolean_operations() has been run to ensure
        that all the layers that should exist are included within the generated
        layer properies file.

        Parameters
        ----------
        filename: str
            This is the name for the layer_properties file. Save will
            automatically append `_ANL` to the filename. If this does not
            include a ".lyp" file extention it will be automatically added.

        ANL_Main: gdspy.Cell
            The Main Cell with the ANL mask on.

        Output
        ------
        xml_layer_prop .lyp file.
        """
        if filename.split(".")[-1] == "lyp":
            filename = str(filename.split(".")[:-1])
        filename += "_ANL.lyp"

        layer_numbers_in_mask = list(ANL_Main.get_layers())

        layers.generate_layer_properties(
            filename,
            layer_numbers_in_mask,
            layers.LayerPropType.ANL,
            self.layers_ANL,
            add_layer_number_to_name=True,
        )
        return

    def save_layer_props(self, filename: str) -> None:
        """Save a layer properties file for the mask under the name given. This
        should be done after do_boolean_operations() has been run to ensure
        that all the layers that should exist are included within the generated
        layer properies file.

        Parameters
        ----------
        filename: str
            This is the name for the layer_properties file. If this does not
            include a ".lyp" file extention it will be automatically added.

        Output
        ------
        xml_layer_prop .lyp file.
        """

        layer_numbers_in_mask = list(self.Main.get_layers())

        layers.generate_layer_properties(
            filename,
            layer_numbers_in_mask,
            layers.LayerPropType.SOUK,
            self.layers,
            add_layer_number_to_name=True,
        )

        return

    def save_resonator_details_file(self, filename: str) -> None:
        """Write and save a file detailing all the resonators added to the
        mask.

        Parameters
        ----------
        filename: str
            The name for the output file for resonator details.
        """
        if len(self.resonators_on_mask) == 0:
            print("No resonators on mask. Writing resonator detials will be skipped.")
            return

        data_for_dataframe: dict[str, list[Any]] = {
            "uid": [],
            "KID_type": [],
            "KID_No": [],
            "x_coord": [],
            "y_coord": [],
            "rot": [],
            "f0": [],
            "mux_IDC": [],
            "mux_CC": [],
            "trim": [],
            "mux_override": [],
            "config_override": [],
        }

        for uid, resonator in enumerate(self.resonators_on_mask):
            details = resonator.get_details()
            data_for_dataframe["uid"].append(uid)
            data_for_dataframe["KID_type"].append(details.get("KID_type"))
            data_for_dataframe["KID_No"].append(details.get("KID_No"))
            data_for_dataframe["x_coord"].append(details.get("x_coord"))
            data_for_dataframe["y_coord"].append(details.get("y_coord"))
            data_for_dataframe["rot"].append(details.get("rot"))
            data_for_dataframe["f0"].append(details.get("f0"))
            data_for_dataframe["mux_IDC"].append(details.get("mux_IDC"))
            data_for_dataframe["mux_CC"].append(details.get("mux_CC"))
            data_for_dataframe["trim"].append(details.get("trim"))
            data_for_dataframe["mux_override"].append(details.get("mux_override"))
            data_for_dataframe["config_override"].append(details.get("config_override"))

        df = pd.DataFrame.from_dict(data_for_dataframe)

        if filename[:-4] != ".csv":
            filename += ".csv"
        df.to_csv(filename, index=False)

        print(f"\nWritten resonator details file '{filename}'.\n")

        return

    def convert_main_to_ANL(self) -> gdspy.Cell:
        """Convert the Main cell to ANL Main.

        Returns
        -------
        ANL_Main: gdspy.Cell
            A gdspy Cell with the ANL mask conversion.
        """
        print("Converting Mask to ANL version.")
        ANL_Main = gdspy.Cell("ANL_MAIN_CELL")

        PRECISION = 0.0001

        wafer_polyset = self.make_wafer_shape(Layer("None", 0, 0))

        # nb_wiring
        print("├─  nb_wiring")
        # Merging all the Nb to one layer
        anl_nb_wiring = gdspy.boolean(
            self.Main.get_polygons(
                [self.layers.Nb_Antenna.number, self.layers.Nb_Antenna.datatype],
            ),
            self.Main.get_polygons(
                [self.layers.IDC_Nb.number, self.layers.IDC_Nb.datatype],
            ),
            "or",
            precision=PRECISION,
            layer=self.layers_ANL.nb_wiring.number,
            datatype=self.layers_ANL.nb_wiring.datatype,
        )
        ANL_Main.add(anl_nb_wiring)

        # aluminium
        print("├─  aluminium")
        anl_aluminium = gdspy.boolean(
            self.Main.get_polygons(
                [self.layers.Aluminium.number, self.layers.Aluminium.datatype],
            ),
            self.Main.get_polygons(
                [self.layers.Aluminium_Direct.number, self.layers.Aluminium_Direct.datatype],
            ),
            "or",
            precision=PRECISION,
            layer=self.layers_ANL.aluminium.number,
            datatype=self.layers_ANL.aluminium.datatype,
        )
        ANL_Main.add(anl_aluminium)

        # dielectric
        print("├─  dielectric")
        anl_dielectric = gdspy.boolean(
            self.Main.get_polygons(
                [self.layers.SiN_dep.number, self.layers.SiN_dep.datatype],
            ),
            None,
            "or",
            precision=PRECISION,
            layer=self.layers_ANL.dielectric.number,
            datatype=self.layers_ANL.dielectric.datatype,
        )
        ANL_Main.add(anl_dielectric)

        # groundplane
        print("├─  groundplane")
        nb_groundplane_polygons = self.Main.get_polygons(
            [self.layers.Nb_Groundplane.number, self.layers.Nb_Groundplane.datatype],
        )

        print(f"nb_groundplane_polygons:\nlen = {len(nb_groundplane_polygons)}")
        anl_groundplane = gdspy.boolean(
            wafer_polyset,
            nb_groundplane_polygons,
            "not",
            precision=PRECISION,
            layer=self.layers_ANL.groundplane.number,
            datatype=self.layers_ANL.groundplane.datatype,
        )
        ANL_Main.add(anl_groundplane)

        # oxide
        print("├─  oxide")
        anl_oxide = gdspy.boolean(
            self.silicon_oxide_cutouts.get_polygons(),
            None,
            "not",
            precision=PRECISION,
            layer=self.layers_ANL.oxide.number,
            datatype=self.layers_ANL.oxide.datatype,
        )
        ANL_Main.add(anl_oxide)

        # membrane
        print("├─  membrane")
        anl_membrane = gdspy.boolean(
            self.silicon_nitride_membrane_cutouts.get_polygons(),
            None,
            "not",
            precision=PRECISION,
            layer=self.layers_ANL.membrane.number,
            datatype=self.layers_ANL.membrane.datatype,
        )
        ANL_Main.add(anl_membrane)

        # al_patch
        # al_etch
        print("├─  al_patch")
        print("├─  al_etch")
        self.add_Aluminium_Patch_and_Etch(
            oversize_aluminium_size_diff=0.0,
            alternate_cell=ANL_Main,
            alternate_patch_layer=self.layers_ANL.al_patch,
            alternate_etch_layer=self.layers_ANL.al_etch,
            printing=False,
        )

        # nb_patch
        print("├─  nb_patch")
        nb_patch_polygons = self.Main.get_polygons(
            [self.layers.Nb_Patch.number, self.layers.Nb_Patch.datatype],
        )
        if len(nb_patch_polygons) == 0:
            print("├───   NO Nb_Patch Polygons found!")
        else:
            anl_nb_patch = gdspy.boolean(
                nb_patch_polygons,
                None,
                "or",
                precision=PRECISION,
                layer=self.layers_ANL.nb_patch.number,
                datatype=self.layers_ANL.nb_patch.datatype,
            )
            ANL_Main.add(anl_nb_patch)

        # nb_etch
        print("├─  nb_etch")
        nb_etch_polygons = self.Main.get_polygons(
            [self.layers.Nb_Etch.number, self.layers.Nb_Etch.datatype],
        )
        if len(nb_etch_polygons) == 0:
            print("├───   NO Nb_Etch Polygons found!")
        else:
            anl_nb_etch = gdspy.boolean(
                nb_etch_polygons,
                None,
                "or",
                precision=PRECISION,
                layer=self.layers_ANL.nb_etch.number,
                datatype=self.layers_ANL.nb_etch.datatype,
            )
            ANL_Main.add(anl_nb_etch)

        print("└─  Done")

        return ANL_Main

    def convert_main_to_ANL_OLD(self) -> gdspy.Cell:
        """Convert the Main cell to ANL Main.

        Returns
        -------
        ANL_Main: gdspy.Cell
            A gdspy Cell with the ANL mask conversion.
        """
        print("Converting Mask to ANL version.")
        ANL_Main = gdspy.Cell("ANL_MAIN_CELL")

        PRECISION = 0.0001

        wafer_polyset = self.make_wafer_shape(Layer("None", 0, 0))

        # nb_wiring
        print("├─  nb_wiring")
        # Merging all the Nb to one layer
        anl_nb_wiring = gdspy.boolean(
            self.Main.get_polygons(
                [self.layers.Nb_Antenna.number, self.layers.Nb_Antenna.datatype],
            ),
            self.Main.get_polygons(
                [self.layers.IDC_Nb.number, self.layers.IDC_Nb.datatype],
            ),
            "or",
            precision=PRECISION,
            layer=self.layers_ANL.nb_wiring.number,
            datatype=self.layers_ANL.nb_wiring.datatype,
        )
        ANL_Main.add(anl_nb_wiring)

        # aluminium
        print("├─  aluminium")
        anl_aluminium = gdspy.boolean(
            self.Main.get_polygons(
                [self.layers.Aluminium.number, self.layers.Aluminium.datatype],
            ),
            self.Main.get_polygons(
                [self.layers.Aluminium_Direct.number, self.layers.Aluminium_Direct.datatype],
            ),
            "or",
            precision=PRECISION,
            layer=self.layers_ANL.aluminium.number,
            datatype=self.layers_ANL.aluminium.datatype,
        )
        ANL_Main.add(anl_aluminium)

        # dielectric
        print("├─  dielectric")

        # dielectric_polygons = self.Main.get_polygons(
        #     [self.layers.SiN_dep.number, self.layers.SiN_dep.datatype],
        # )
        #
        # print(f"dielectric_polygons:\nlen = {len(dielectric_polygons)}")
        #
        # # Chunking the dielectric boolean operation.
        # chunk_size = 15000
        # no_of_chunks = 1 + (len(dielectric_polygons) // chunk_size)  # "//" is integer divide
        #
        # test_dielectric = wafer_polyset
        # for chunk_number in range(no_of_chunks):
        #     chunk_start = chunk_number * chunk_size
        #     chunk_end = (chunk_number + 1) * chunk_size
        #     print(f"├─────  chunk_number:{chunk_number}")
        #     print(f"├─────  [{chunk_start}:{chunk_end}]")
        #
        #     test_dielectric = gdspy.boolean(
        #         test_dielectric,
        #         dielectric_polygons[chunk_start:chunk_end],
        #         "not",
        #         # precision=PRECISION,
        #         precision=0.1,
        #         layer=self.layers_ANL.dielectric.number,
        #         datatype=self.layers_ANL.dielectric.datatype,
        #     )
        # ANL_Main.add(test_dielectric)
        # print("├─  ADDED dielectric")

        dielectric_polygons = self.Main.get_polygons(
            [self.layers.SiN_dep.number, self.layers.SiN_dep.datatype],
        )
        wafer_polygons = wafer_polyset.polygons

        print(f"len of dielectric_polygons:{len(dielectric_polygons)}")
        print(f"len of wafer_polygons:{len(wafer_polygons)}")

        s_diel_polys: list[shapely_geom.Polygon] = []
        for i, poly_points in enumerate(dielectric_polygons):
            s_poly = shapely_geom.Polygon(poly_points)
            if s_poly.is_valid:
                s_diel_polys.append(s_poly)
            else:
                cleaned = s_poly.buffer(0)
                if cleaned.is_valid:
                    s_diel_polys.append(cleaned)
                else:
                    print(f"error at poly number {i}")
                    print(f"poly_points = {poly_points}")
                    print(f"s_poly = {s_poly}")
                    print(f"cleanded.is_valid = {cleaned.is_valid}")
                    print(f"cleanded = {cleaned}")
                    raise RuntimeError
        s_diel_multipoly = shapely_geom.MultiPolygon(s_diel_polys)
        cleaned_s_diel_multipoly = s_diel_multipoly.buffer(0)

        s_wafer_polys: list[shapely_geom.Polygon] = []
        for i, poly_points in enumerate(wafer_polygons):
            s_poly = shapely_geom.Polygon(poly_points)
            if s_poly.is_valid:
                s_wafer_polys.append(s_poly)
            else:
                cleaned = s_poly.buffer(0)
                if cleaned.is_valid:
                    s_wafer_polys.append(cleaned)
                else:
                    print(f"error at poly number {i}")
                    print(f"poly_points = {poly_points}")
                    print(f"s_poly = {s_poly}")
                    print(f"cleanded.is_valid = {cleaned.is_valid}")
                    print(f"cleanded = {cleaned}")
                    raise RuntimeError
        s_wafer_multipoly = shapely_geom.MultiPolygon(s_wafer_polys)
        cleaned_s_wafer_multipoly = s_wafer_multipoly.buffer(0)

        test_diff = cleaned_s_diel_multipoly
        for polygon in cleaned_s_diel_multipoly.geoms:
            ints = polygon.interiors
            for interior in ints:
                ixs, iys = interior.xy
                int_points = np.vstack((ixs, iys)).T
                poly_int = gdspy.Polygon(
                    int_points,
                    layer=1235,
                    datatype=0,
                )
                ANL_Main.add(poly_int)

                s_poly = shapely_geom.Polygon(int_points)
                test_diff = test_diff.difference(s_poly)

            exts = polygon.exterior
            exs, eys = exts.xy
            ext_points = np.vstack((exs, eys)).T
            poly_ext = gdspy.Polygon(
                ext_points,
                layer=1234,
                datatype=0,
            )
            ANL_Main.add(poly_ext)

        # for polygon in cleaned_s_wafer_multipoly.geoms:
        #     xs, ys = polygon.exterior.xy
        #     points = np.vstack((xs, ys)).T
        #     anl_dielectric_poly = gdspy.Polygon(
        #         points,
        #         layer=1235,
        #         datatype=0,
        #     )
        #     ANL_Main.add(anl_dielectric_poly)

        # res = s_wafer_multipoly.difference(s_diel_multipoly)
        # diff = shapely_difference(cleaned_s_wafer_multipoly, cleaned_s_diel_multipoly)
        # diff = cleaned_s_diel_multipoly.difference(cleaned_s_wafer_multipoly)
        # diff = cleaned_s_wafer_multipoly.difference(cleaned_s_diel_multipoly)
        print("doing diff")
        diff = cleaned_s_wafer_multipoly.difference(test_diff)
        print(f"diff type = {type(diff)}")
        print(type(diff))

        for polygon in diff.geoms:
            exs, eys = polygon.exterior.xy
            epoints = np.vstack((exs, eys)).T
            ext_poly = gdspy.Polygon(
                epoints,
                layer=2000,
                datatype=0,
            )
            ANL_Main.add(ext_poly)

            for interior in polygon.interiors:
                ixs, iys = polygon.exterior.xy
                ipoints = np.vstack((ixs, iys)).T
                int_poly = gdspy.Polygon(
                    ipoints,
                    layer=2001,
                    datatype=0,
                )
                ANL_Main.add(int_poly)

        print("├─  ADDED dielectric")

        # groundplane
        print("├─  groundplane")
        nb_groundplane_polygons = self.Main.get_polygons(
            [self.layers.Nb_Groundplane.number, self.layers.Nb_Groundplane.datatype],
        )

        print(f"nb_groundplane_polygons:\nlen = {len(nb_groundplane_polygons)}")
        anl_groundplane = gdspy.boolean(
            wafer_polyset,
            nb_groundplane_polygons,
            "not",
            precision=PRECISION,
            layer=self.layers_ANL.groundplane.number,
            datatype=self.layers_ANL.groundplane.datatype,
        )
        ANL_Main.add(anl_groundplane)

        # oxide
        print("├─  oxide")
        anl_oxide = gdspy.boolean(
            self.silicon_oxide_cutouts.get_polygons(),
            None,
            "not",
            precision=PRECISION,
            layer=self.layers_ANL.oxide.number,
            datatype=self.layers_ANL.oxide.datatype,
        )
        ANL_Main.add(anl_oxide)

        # membrane
        print("├─  membrane")
        anl_membrane = gdspy.boolean(
            self.silicon_nitride_membrane_cutouts.get_polygons(),
            None,
            "not",
            precision=PRECISION,
            layer=self.layers_ANL.membrane.number,
            datatype=self.layers_ANL.membrane.datatype,
        )
        ANL_Main.add(anl_membrane)

        # al_patch
        # al_etch
        print("├─  al_patch")
        print("├─  al_etch")
        self.add_Aluminium_Patch_and_Etch(
            oversize_aluminium_size_diff=0.0,
            alternate_cell=ANL_Main,
            alternate_patch_layer=self.layers_ANL.al_patch,
            alternate_etch_layer=self.layers_ANL.al_etch,
            printing=False,
        )

        # nb_patch
        print("├─  nb_patch")
        nb_patch_polygons = self.Main.get_polygons(
            [self.layers.Nb_Patch.number, self.layers.Nb_Patch.datatype],
        )
        if len(nb_patch_polygons) == 0:
            print("├───   NO Nb_Patch Polygons found!")
        else:
            anl_nb_patch = gdspy.boolean(
                nb_patch_polygons,
                None,
                "or",
                precision=PRECISION,
                layer=self.layers_ANL.nb_patch.number,
                datatype=self.layers_ANL.nb_patch.datatype,
            )
            ANL_Main.add(anl_nb_patch)

        # nb_etch
        print("├─  nb_etch")
        nb_etch_polygons = self.Main.get_polygons(
            [self.layers.Nb_Etch.number, self.layers.Nb_Etch.datatype],
        )
        if len(nb_etch_polygons) == 0:
            print("├───   NO Nb_Etch Polygons found!")
        else:
            anl_nb_etch = gdspy.boolean(
                nb_etch_polygons,
                None,
                "or",
                precision=PRECISION,
                layer=self.layers_ANL.nb_etch.number,
                datatype=self.layers_ANL.nb_etch.datatype,
            )
            ANL_Main.add(anl_nb_etch)

        print("└─  Done")

        return ANL_Main

    def save_gds_file_ANL_version(
        self,
        filename: str,
        make_mirrored_along_x: bool = False,
        mirrored_x_filename: str = "",
        make_mirrored_along_y: bool = False,
        mirrored_y_filename: str = "",
        save_layer_props: bool = True,
    ) -> None:
        """Write and save the mask in the ANL format. This saves everything in
        the Main cell to a gds file under the given filename. Additionally can
        save mirrored versions.

        Parameters
        ----------
        filename: str
            This is the name for the gds file that will be written. Save will
            automatically append `_ANL` to the filename. If this filename does
            not include a ".gds" file extention it will be automatically added.

        KwArgs
        ------
        make_mirrored_along_x, make_mirrored_along_y: bool = False
            Whether or not a mirrored version of the mask (mirrored along the x
            or y axis respectively) should also be saved.
            See mirrored_filename for naming this.

        mirrored_x_filename, mirrored_x_filename: str = ""
            This is the name for the mirrored (in the x or y axis respectively)
            version of the mask. Like the filename arg, if any of these do not
            include a ".gds" file extention it will be automatically added. If
            this argument is not provided the mirrored version of the mask will
            take same name as the filename provided with "_MIRRORED_ACROSS_{X/Y}"
            ,depending upon the mirror axis, appended at the end.

        save_layer_props: bool = True
            Save a layer properties file for the mask under with the same
            naming convention as the filename argument.
            See save_layer_props_ANL_version for more info.

        Output
        ------
        xml_layer_prop .lyp file. Disabled by save_layer_props KwArg.
        """
        from phidl import Device

        if filename.split(".")[-1] == "gds":
            filename = str(filename.split(".")[:-1])
        filename += "_ANL.gds"

        pretty_print(f'Writing ANL version to GDS\n"{filename}"', color=TextColor.BLUE)

        self.do_boolean_operations()
        ANL_Main = self.convert_main_to_ANL()

        MainDevice = Device("MainDevice")
        MainDevice.add(ANL_Main)
        MainDevice.write_gds(filename)

        if make_mirrored_along_x:
            if mirrored_x_filename == "":
                mirrored_x_filename = filename[:-4] + "_MIRRORED_ACROSS_X.gds"
            elif mirrored_x_filename[-4:] != ".gds":
                mirrored_x_filename += ".gds"

            pretty_print(f'Writing Mask mirrored along x-axis to GDS\n"{mirrored_x_filename}"', color=TextColor.BLUE)
            Main_mirrored_across_x_axis = ANL_Main.copy("MainCell_xax_mirror", deep_copy=True, x_reflection=True)

            MainDevice_mirrored_x = Device("MainDevice_mirrored_x")
            MainDevice_mirrored_x.add(Main_mirrored_across_x_axis)
            MainDevice_mirrored_x.write_gds(mirrored_x_filename)

        if make_mirrored_along_y:
            # print("Making mirrored along y may produce unexpected results with FlexPaths!")  # Deep moral failings with phidl and gdspy
            if mirrored_y_filename == "":
                mirrored_y_filename = filename[:-4] + "_MIRRORED_ACROSS_Y.gds"
            elif mirrored_y_filename[-4:] != ".gds":
                mirrored_y_filename += ".gds"

            pretty_print(f'Writing Mask mirrored along y-axis to GDS\n"{mirrored_y_filename}"', color=TextColor.BLUE)
            Main_mirrored_across_x_axis_rot90 = ANL_Main.copy(
                "MainCell_xax_mirror_rot90", deep_copy=True, rotation=(pi / 2), x_reflection=True
            )
            Main_mirrored_across_y_axis = Main_mirrored_across_x_axis_rot90.copy("MainCell_yax_mirror", deep_copy=True, rotation=(pi / 2))

            MainDevice_mirrored_y = Device("MainDevice_mirrored_y")
            MainDevice_mirrored_y.add(Main_mirrored_across_y_axis)
            MainDevice_mirrored_y.write_gds(mirrored_y_filename)

        pretty_print("Finished writing files to disk.", color=TextColor.BLUE)

        if save_layer_props:
            self.save_layer_props_ANL_version(filename.split("_ANL.gds")[0], ANL_Main)

        return

    def save_gds_file(
        self,
        filename: str,
        try_to_make_backside: bool = False,
        make_mirrored_along_x: bool = False,
        mirrored_x_filename: str = "",
        make_mirrored_along_y: bool = False,
        mirrored_y_filename: str = "",
        save_layer_props: bool = True,
    ) -> None:
        """Write and save the mask. This saves everything in the Main cell to a
        gds file under the given filename. Additionally can save a mirrored
        version. If Anything exists in the MainBackside cell. this will also be
        saved with "_BACKSIDE" appended to the end.

        Parameters
        ----------
        filename: str
            This is the name for the gds file that will be written. If this
            filename does not include a ".gds" file extention it will be
            automatically added.

        KwArgs
        ------
        try_to_make_backside: bool
            Default False, This will try to make the backside of the mask file
            if anything exists in the MainBackside cell. If nothing exists in
            this cell, Nothing will get made.

        make_mirrored_along_x, make_mirrored_along_y: bool = False
            Whether or not a mirrored version of the mask (mirrored along the x
            or y axis respectively) should also be saved.
            See mirrored_filename for naming this.

        mirrored_x_filename, mirrored_x_filename: str = ""
            This is the name for the mirrored (in the x or y axis respectively)
            version of the mask. Like the filename arg, if any of these do not
            include a ".gds" file extention it will be automatically added. If
            this argument is not provided the mirrored version of the mask will
            take same name as the filename provided with "_MIRRORED_ACROSS_{X/Y}"
            ,depending upon the mirror axis, appended at the end.

        save_layer_props: bool = True
            Save a layer properties file for the mask under with the same
            naming convention as the filename argument.
            See save_layer_props for more info.

        Output
        ------
        xml_layer_prop .lyp file. Disabled by save_layer_props KwArg.
        """
        from phidl import Device

        if filename.split(".")[-1] != "gds":
            filename += ".gds"

        pretty_print(f'Writing to GDS\n"{filename}"', color=TextColor.BLUE)
        self.do_boolean_operations()
        MainDevice = Device("MainDevice")
        MainDevice.add(self.Main)
        MainDevice.write_gds(filename)

        if try_to_make_backside:
            if len(self.MainBackside.get_layers()) != 0:
                backside_filename = filename[:-4] + "_BACKSIDE.gds"
                print(f'Writing backside to GDS "{backside_filename}"')
                MainBacksideDevice = Device("MainBacksideDevice")
                MainBacksideDevice.add(self.MainBackside)
                MainBacksideDevice.write_gds(backside_filename)
            else:
                print("Nothing in MainBackside to procude backside mask.")

        if make_mirrored_along_x:
            if mirrored_x_filename == "":
                mirrored_x_filename = filename[:-4] + "_MIRRORED_ACROSS_X.gds"
            elif mirrored_x_filename[-4:] != ".gds":
                mirrored_x_filename += ".gds"

            pretty_print(f'Writing Mask mirrored along x-axis to GDS\n"{mirrored_x_filename}"', color=TextColor.BLUE)
            Main_mirrored_across_x_axis = self.Main.copy("MainCell_xax_mirror", deep_copy=True, x_reflection=True)

            MainDevice_mirrored_x = Device("MainDevice_mirrored_x")
            MainDevice_mirrored_x.add(Main_mirrored_across_x_axis)
            MainDevice_mirrored_x.write_gds(mirrored_x_filename)

        if make_mirrored_along_y:
            # print("Making mirrored along y may produce unexpected results with FlexPaths!")  # Deep moral failings with phidl and gdspy
            if mirrored_y_filename == "":
                mirrored_y_filename = filename[:-4] + "_MIRRORED_ACROSS_Y.gds"
            elif mirrored_y_filename[-4:] != ".gds":
                mirrored_y_filename += ".gds"

            pretty_print(f'Writing Mask mirrored along y-axis to GDS\n"{mirrored_y_filename}"', color=TextColor.BLUE)
            Main_mirrored_across_x_axis_rot90 = self.Main.copy(
                "MainCell_xax_mirror_rot90", deep_copy=True, rotation=(pi / 2), x_reflection=True
            )
            Main_mirrored_across_y_axis = Main_mirrored_across_x_axis_rot90.copy("MainCell_yax_mirror", deep_copy=True, rotation=(pi / 2))

            MainDevice_mirrored_y = Device("MainDevice_mirrored_y")
            MainDevice_mirrored_y.add(Main_mirrored_across_y_axis)
            MainDevice_mirrored_y.write_gds(mirrored_y_filename)

        pretty_print("Finished writing files to disk.", color=TextColor.BLUE)

        if save_layer_props:
            self.save_layer_props(filename.split(".gds")[0])
        return
