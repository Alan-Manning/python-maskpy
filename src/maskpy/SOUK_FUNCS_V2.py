"""
This is a package for making elements for a SOUK gds mask.

This package relys upon gdspy.
"""
# import threading
import concurrent.futures
import copy
import os
import time
import xml.etree.ElementTree as ET

import gdspy
import numpy as np
import pandas as pd
from numpy import cos as cos

# import math
from numpy import pi as pi
from numpy import sin as sin
from numpy import tan as tan
from phidl import geometry as phgeom
from shapely import geometry as shapely_geom

# from alive_progress import alive_bar
from tqdm import tqdm

"""
###############################################################################
Functions
###############################################################################
"""


class SOUK_Functions:
    """
    Functions for making elements for an SOUK gds mask.

    ...

    gdspy Cell Attributes
    ---------------------
    Main

    MainBackside

    ground_plane_cutouts

    ground_plane_positives

    silicon_nitride_cutouts

    silicon_nitride_positives

    silicon_oxide_cutouts

    silicon_oxide_positives

    silicon_nitride_membrane_cutouts

    silicon_nitride_membrane_positives


    Layer Attributes
    ----------------
    The name of all the layers used in the functions to generate souk mask.
    Below are all the default layer numbers and datatypes. The layer numbers
    amd datatypes can be changed when the class is instantiated but layernames
    need to be as stated below.

    Aluminium : {"layer": 1,   "datatype": 0}

    Nb_Antenna : {"layer": 2,   "datatype": 0}

    SiN_dep : {"layer": 3,   "datatype": 0}

    IDC_Nb : {"layer": 4,   "datatype": 0}

    SiO : {"layer": 5,   "datatype": 0}

    SiN_Membrane : {"layer": 6,   "datatype": 0}

    Aluminium_Bulk : {"layer": 7,   "datatype": 0}

    Backside_Check : {"layer": 9,   "datatype": 0}

    Backside_Flipped : {"layer": 29,  "datatype": 0}

    Nb_Groundplane : {"layer": 30,  "datatype": 0}

    Pin_hole_positives : {"layer": 31,  "datatype": 0}

    TrimLayer : {"layer": 90,  "datatype": 0}

    Top_choke_waveguide_hole : {"layer": 150, "datatype": 0}

    Top_choke_anulus : {"layer": 151, "datatype": 0}

    Bottom_choke_waveguide_hole : {"layer": 155, "datatype": 0}

    Bottom_choke_IDC_hole : {"layer": 156, "datatype": 0}

    Bottom_choke_pads : {"layer": 157, "datatype": 0}

    Tab_dicing_line : {"layer": 160, "datatype": 0}

    Chip_holder : {"layer": 175, "datatype": 0}

    General_labeling : {"layer": 65535, "datatype": 0}

    Methods
    -------


    """

    def __init__(self, main_cell_name="MAIN_CELL", layers="Default"):
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

        default_layers = {
            "Aluminium": {"layer": 1, "datatype": 0},
            "Nb_Antenna": {"layer": 2, "datatype": 0},
            "SiN_dep": {"layer": 3, "datatype": 0},
            "IDC_Nb": {"layer": 4, "datatype": 0},
            "SiO": {"layer": 5, "datatype": 0},
            "SiN_Membrane": {"layer": 6, "datatype": 0},
            "Aluminium_Bulk": {"layer": 7, "datatype": 0},
            "Backside_Check": {"layer": 9, "datatype": 0},
            "Backside_Flipped": {"layer": 29, "datatype": 0},
            "Nb_Groundplane": {"layer": 30, "datatype": 0},
            "Pin_hole_positives": {"layer": 31, "datatype": 0},
            "TrimLayer": {"layer": 90, "datatype": 0},
            "Top_choke_waveguide_hole": {"layer": 150, "datatype": 0},
            "Top_choke_anulus": {"layer": 151, "datatype": 0},
            "Bottom_choke_waveguide_hole": {"layer": 155, "datatype": 0},
            "Bottom_choke_IDC_hole": {"layer": 156, "datatype": 0},
            "Bottom_choke_pads": {"layer": 157, "datatype": 0},
            "Tab_dicing_line": {"layer": 160, "datatype": 0},
            "Bottom_choke_Tab_dicing_line": {"layer": 161, "datatype": 0},
            "Chip_holder": {"layer": 175, "datatype": 0},
            "General_labeling": {"layer": 65535, "datatype": 0},
        }

        if layers == "Default":
            self.Aluminium = default_layers["Aluminium"]
            self.Nb_Antenna = default_layers["Nb_Antenna"]
            self.SiN_dep = default_layers["SiN_dep"]
            self.IDC_Nb = default_layers["IDC_Nb"]
            self.SiO = default_layers["SiO"]
            self.SiN_Membrane = default_layers["SiN_Membrane"]
            self.Aluminium_Bulk = default_layers["Aluminium_Bulk"]
            self.Backside_Check = default_layers["Backside_Check"]
            self.Backside_Flipped = default_layers["Backside_Flipped"]
            self.Nb_Groundplane = default_layers["Nb_Groundplane"]
            self.Pin_hole_positives = default_layers["Pin_hole_positives"]
            self.TrimLayer = default_layers["TrimLayer"]
            self.Top_choke_waveguide_hole = default_layers["Top_choke_waveguide_hole"]
            self.Top_choke_anulus = default_layers["Top_choke_anulus"]
            self.Bottom_choke_waveguide_hole = default_layers["Bottom_choke_waveguide_hole"]
            self.Bottom_choke_IDC_hole = default_layers["Bottom_choke_IDC_hole"]
            self.Bottom_choke_pads = default_layers["Bottom_choke_pads"]
            self.Tab_dicing_line = default_layers["Tab_dicing_line"]
            self.Bottom_choke_Tab_dicing_line = default_layers["Bottom_choke_Tab_dicing_line"]
            self.Chip_holder = default_layers["Chip_holder"]
            self.General_labeling = default_layers["General_labeling"]

        elif type(layers) == dict:
            self.Aluminium = layers["Aluminium"] if "Aluminium" in layers else default_layers["Aluminium"]
            self.Nb_Antenna = layers["Nb_Antenna"] if "Nb_Antenna" in layers else default_layers["Nb_Antenna"]
            self.SiN_dep = layers["SiN_dep"] if "SiN_dep" in layers else default_layers["SiN_dep"]
            self.IDC_Nb = layers["IDC_Nb"] if "IDC_Nb" in layers else default_layers["IDC_Nb"]
            self.SiO = layers["SiO"] if "SiO" in layers else default_layers["SiO"]
            self.SiN_Membrane = layers["SiN_Membrane"] if "SiN_Membrane" in layers else default_layers["SiN_Membrane"]
            self.Aluminium_Bulk = layers["Aluminium_Bulk"] if "Aluminium_Bulk" in layers else default_layers["Aluminium_Bulk"]
            self.Backside_Check = layers["Backside_Check"] if "Backside_Check" in layers else default_layers["Backside_Check"]
            self.Backside_Flipped = layers["Backside_Flipped"] if "Backside_Flipped" in layers else default_layers["Backside_Flipped"]
            self.Nb_Groundplane = layers["Nb_Groundplane"] if "Nb_Groundplane" in layers else default_layers["Nb_Groundplane"]
            self.Pin_hole_positives = (
                layers["Pin_hole_positives"] if "Pin_hole_positives" in layers else default_layers["Pin_hole_positives"]
            )
            self.TrimLayer = layers["TrimLayer"] if "TrimLayer" in layers else default_layers["TrimLayer"]
            self.Top_choke_waveguide_hole = (
                layers["Top_choke_waveguide_hole"] if "Top_choke_waveguide_hole" in layers else default_layers["Top_choke_waveguide_hole"]
            )
            self.Top_choke_anulus = layers["Top_choke_anulus"] if "Top_choke_anulus" in layers else default_layers["Top_choke_anulus"]
            self.Bottom_choke_waveguide_hole = (
                layers["Bottom_choke_waveguide_hole"]
                if "Bottom_choke_waveguide_hole" in layers
                else default_layers["Bottom_choke_waveguide_hole"]
            )
            self.Bottom_choke_IDC_hole = (
                layers["Bottom_choke_IDC_hole"] if "Bottom_choke_IDC_hole" in layers else default_layers["Bottom_choke_IDC_hole"]
            )
            self.Bottom_choke_pads = layers["Bottom_choke_pads"] if "Bottom_choke_pads" in layers else default_layers["Bottom_choke_pads"]
            self.Tab_dicing_line = layers["Tab_dicing_line"] if "Tab_dicing_line" in layers else default_layers["Tab_dicing_line"]
            self.Bottom_choke_Tab_dicing_line = (
                layers["Bottom_choke_Tab_dicing_line"]
                if "Bottom_choke_Tab_dicing_line" in layers
                else default_layers["Bottom_choke_Tab_dicing_line"]
            )

        else:
            raise Exception(
                'layers not in correct format to read. should be dictionary with keys of layer names and values of dictionarys in the form "{"layer": 157, "datatype": 0}"'
            )

        self.all_layers_name_lookup_from_number = {
            self.Aluminium["layer"]: "Aluminium",
            self.Nb_Antenna["layer"]: "Nb_Antenna",
            self.SiN_dep["layer"]: "SiN_dep",
            self.IDC_Nb["layer"]: "IDC_Nb",
            self.SiO["layer"]: "SiO",
            self.SiN_Membrane["layer"]: "SiN_Membrane",
            self.Aluminium_Bulk["layer"]: "Aluminium_Bulk",
            self.Backside_Check["layer"]: "Backside_Check",
            self.Backside_Flipped["layer"]: "Backside_Flipped",
            self.Nb_Groundplane["layer"]: "Nb_Groundplane",
            self.Pin_hole_positives["layer"]: "Pin_hole_positives",
            self.TrimLayer["layer"]: "TrimLayer",
            self.Top_choke_waveguide_hole["layer"]: "Top_choke_waveguide_hole",
            self.Top_choke_anulus["layer"]: "Top_choke_anulus",
            self.Bottom_choke_waveguide_hole["layer"]: "Bottom_choke_waveguide_hole",
            self.Bottom_choke_IDC_hole["layer"]: "Bottom_choke_IDC_hole",
            self.Bottom_choke_pads["layer"]: "Bottom_choke_pads",
            self.Tab_dicing_line["layer"]: "Tab_dicing_line",
            self.Bottom_choke_Tab_dicing_line["layer"]: "Bottom_choke_Tab_dicing_line",
            self.Chip_holder["layer"]: "Chip_holder",
            self.General_labeling["layer"]: "General_labeling",
        }

    def inside_hexagon(self, xy, hex_rad, hex_center=[0, 0]):
        """
        Checks if an xy point sits within a hexagon with given radius.

        Parameters
        ----------
        xy : list
            list of [x,y] coordinates of a point to check if its inside the
            hexagon.

        hex_rad : float, int
            The radius of the hexagon to check if point sits inside.

        KwArgs
        ------
        hex_center = [0,0]
            list of the [x,y] coordinates for the center of the haxagon.

        Returns
        -------
        boolean :
            True if the xy point given sits within the hexagon.
            False if the xy point sits on the boundary or outside the hexagon.

        """

        R = hex_rad
        x, y = map(abs, [xy[0] - hex_center[0], xy[1] - hex_center[0]])
        cx, cy = hex_center
        return y < np.sqrt(3) * min((R - x), (R / 2))

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

            if self.inside_hexagon([init_x, gy + (y_pitch / 2)], hexagon_radius, hexagon_center):
                inside_horiz = True
                while inside_horiz:
                    gx = init_x + horz_count * x_pitch
                    horz_count += 1

                    if self.inside_hexagon([gx + (x_pitch / 2), gy], hexagon_radius, hexagon_center):
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

            if self.inside_hexagon([init_x, gy + (y_pitch / 2)], hexagon_radius, hexagon_center):
                inside_horiz = True
                while inside_horiz:
                    gx = init_x + horz_count * x_pitch
                    horz_count += 1

                    if self.inside_hexagon([gx + (x_pitch / 2), gy], hexagon_radius, hexagon_center):
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
        hexagon_radius,
        hexagon_center,
        x_pitch,
        y_pitch,
        middle_gap_distance=0,
        adjust_center_hex_grid=False,
        return_two_grids=False,
        return_split_point=False,
    ):
        """
        Generates a split hex pack grid within a hexagon. This is a hex pack that has a
        horizontal split along the center.

        Parameters
        ----------
        hexagon_radius : float, int
            The radius of the hexagon in which to hex pack within. Radius refers to
            the distance from the center of the hexagon to any hex point.

        hexagon_center : list
            List of [x,y] coordinates of the center of the hexagon to generate the
            hex pack within.

        x_pitch, y_pitch :  : float, int
            The horizontal, vertical distance between the hex pack grid points.
            Each row in the grid is offset by half the horizontal pitch.

        KwArgs
        ------

        middle_gap_distance = 0
            The middle gap between the two halfs of the hex pack grid.

        adjust_center_hex_grid = False
            Determines whether to adjust the halfs of the hex pack such they they
            sit evenly horizontally about the center of the hexagon they are
            confined within.

        return_two_grids = False
            Determines whether or not to return two seperate grids. If False as by
            default, it will pack the two halfs of the grid into one long list
            containing [x,y] points defining the hex pack grid points. If True this
            function will return a list containing the two halfs of the grid,
            [bot_hex_pack, top_hex_pack], where each is a list of [x,y] points.

        return_split_point=False
            Determines whether or not to return the point at which the hex pack
            transitions from the bot_hex_pack to the top_hex_pack. This is simply
            the length of the bot_hex_pack list. This has no effect if the
            return_two_grids argument is True.

        Returns
        -------
        DEFAULT
            hex_grid : lsit
            list containing [x,y] lists which define the hex pack grid points.

        OR if return_split_point=True
            RETURNS hex_grid, len(bot_hex_grid)
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

                if self.inside_hexagon([gx + (x_pitch / 2), gy - (y_pitch / 2)], hexagon_radius, hexagon_center) and self.inside_hexagon(
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
        """
        Takes in a list of [x,y] grid points and centers them about a central point.

        Parameters
        ----------
        hex_grid : list
            list of [x,y] lists that define the points in the hex pack grid.

        wafer_center_xy : list
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
        shifted_grid : list
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

    def generate_octagonal_holder_port_dict(self, center_xy, octagon_rotation, wafer_radius=75000.0):
        """
        Generates a dictionary of port x, y, and rotaions for the octagonal holder
        design. This contains 16 port locations where the two on each edge of the
        octagon are 10000um from the center point of that edge. This dictionary has
        keys "0" through "16" where "0" is the left most face upper, the numbering
        continues going counter-clockwise. each key contains the x,y and roation of
        the octagon port conection points.

        Parameters
        ----------
        center_xy : list
            list containing the [x,y] coordinate for the center point which to
            place the octagonal holder around.

        octagon_rotation : float, int
            the angle (**in Degrees**) for the rotaion of the holder around the
            center_xy argument.

        KwArgs
        ------
        wafer_radius : float, int
            This is the radius of the wafer where the holder positions will be.
            Default is 75000.0 which is a 3inch radius for a six inch wafer.

        Returns
        -------
        oct_ports_dict : dict
            dictionary containing the x,y,rot(degrees) for each port in the
            octagonal holder.

        """

        six_inch = 150000
        three_inch = six_inch / 2
        wafer_diameter = 2 * wafer_radius

        center_wafer_x, center_wafer_y = center_xy
        # the roation of the the octagonal holder for the ports to connect into
        oct_ports_dict = {}  # setting up a dictionary of {x,y,rotation} for each possible conection in the octagonal holder

        oct_ports_dict["0"] = list(
            self.rotate(
                center_wafer_x,
                center_wafer_y,
                (center_wafer_x - np.sqrt(wafer_radius**2 - 10000**2)),
                (center_wafer_y + 10000),
                self.deg_to_rad(octagon_rotation),
            )
        ) + [
            octagon_rotation
        ]  # sets the x,y and roation of the first two octagon port conection points
        oct_ports_dict["1"] = list(
            self.rotate(
                center_wafer_x,
                center_wafer_y,
                (center_wafer_x - np.sqrt(wafer_radius**2 - 10000**2)),
                (center_wafer_y - 10000),
                self.deg_to_rad(octagon_rotation),
            )
        ) + [octagon_rotation]

        for i in range(1, 8):  # makes the rest of the ports around the octagonal holder
            oct_ports_dict[f"{2 * i}"] = list(
                self.rotate(
                    center_wafer_x,
                    center_wafer_y,
                    oct_ports_dict[f"{2 * (i - 1)}"][0],
                    oct_ports_dict[f"{2 * (i - 1)}"][1],
                    pi / 4,
                )
            ) + [(oct_ports_dict[f"{2 * (i - 1)}"][2] + 45) % 360]
            oct_ports_dict[f"{(2 * i) + 1}"] = list(
                self.rotate(
                    center_wafer_x,
                    center_wafer_y,
                    oct_ports_dict[f"{2 * i - 1}"][0],
                    oct_ports_dict[f"{2 * i - 1}"][1],
                    pi / 4,
                )
            ) + [(oct_ports_dict[f"{2 * i - 1}"][2] + 45) % 360]

        for i in range(8):  # adds a rectangle around the outside of the wafer shoing where the ports would be+
            self.Main.add(
                gdspy.Rectangle(
                    (oct_ports_dict[f"{2 * i}"][0] - 1000, oct_ports_dict[f"{2 * i}"][1] - 250),
                    (oct_ports_dict[f"{2 * i}"][0], oct_ports_dict[f"{2 * i}"][1] + 250),
                ).rotate(
                    self.deg_to_rad(oct_ports_dict[f"{2 * i}"][2]),
                    (oct_ports_dict[f"{2 * i}"][0], oct_ports_dict[f"{2 * i}"][1]),
                )
            )
            self.Main.add(
                gdspy.Rectangle(
                    (oct_ports_dict[f"{(2 * i) + 1}"][0] - 1000, oct_ports_dict[f"{(2 * i) + 1}"][1] - 250),
                    (oct_ports_dict[f"{(2 * i) + 1}"][0], oct_ports_dict[f"{(2 * i) + 1}"][1] + 250),
                ).rotate(
                    self.deg_to_rad(oct_ports_dict[f"{(2 * i) + 1}"][2]),
                    (oct_ports_dict[f"{(2 * i) + 1}"][0], oct_ports_dict[f"{(2 * i) + 1}"][1]),
                )
            )

        return oct_ports_dict

    def make_Toms_6_inch_holder_and_get_ports(self, chip_center_xy, wafer_flat_length, middle_flat_angle=(-pi / 2)):
        """
        Generates a dictionary of port x, y, and rotaions for the 6inch wafer
        holder design. This contains 18 port locations where the three on each
        edge of the center point and vertexes of a hexagon edge. This
        dictionary has keys "0" through "18" where "0" is botttom most middle
        right of the wafer where the flat is located at the bottom, the
        numbering continues going counter-clockwise. Each dict key contains the
        x,y and roation of the port conection points.

        Parameters
        ----------
        center_xy : list
            list containing the [x,y] coordinate for the center point which to
            place the octagonal holder around.

        wafer_flat_length : float, int
            The length of the flat on the silicon wafer.

        KwArgs
        ------
        middle_flat_angle : float, int
            Angle (**in radians**) the falt is located. Default is (-pi/2),
            ie the bottom middle side of the hexagon wafer. This angle is the
            angle the middle of the flat makes with the chip center.

        Returns
        -------
        ports_dict : dict
            dictionary containing the x,y,rot(**in radians**) for each port in
            the 6inch holder.

        """

        six_inch = 150000
        three_inch = six_inch / 2

        wafer_diameter = six_inch
        wafer_radius = three_inch

        holder_outer_hex_radius = six_inch * 4 / 6
        holder_outer_port_width = 5000

        port_edge_pitch = 18000

        ports_dict = {}
        holder_outer_hex_points = []

        holder_positives = gdspy.Cell("HOLDER_POSITIVES")
        holder_cutouts = gdspy.Cell("HOLDER_CUTOUTS")

        # loop making the hexagonal outer and getting port locations and rotations
        for i in range(6):
            ang_cent_to_edge = i * (pi / 3)
            ang_edge_to_cent = ang_cent_to_edge - pi

            holder_outer_hex_points.append(
                [
                    chip_center_xy[0] + holder_outer_hex_radius * cos(ang_cent_to_edge - (pi / 6)),
                    chip_center_xy[1] + holder_outer_hex_radius * sin(ang_cent_to_edge - (pi / 6)),
                ]
            )

            offset_from_cent_angle = np.arcsin(port_edge_pitch / wafer_radius)

            top_edge_port_xy_rot = [
                chip_center_xy[0] + 1.0025 * wafer_radius * cos(ang_cent_to_edge + offset_from_cent_angle),
                chip_center_xy[1] + 1.0025 * wafer_radius * sin(ang_cent_to_edge + offset_from_cent_angle),
                ang_edge_to_cent,
            ]

            mid_edge_port_xy_rot = [
                chip_center_xy[0] + 1.0025 * wafer_radius * cos(ang_cent_to_edge),
                chip_center_xy[1] + 1.0025 * wafer_radius * sin(ang_cent_to_edge),
                ang_edge_to_cent,
            ]

            bot_edge_port_xy_rot = [
                chip_center_xy[0] + 1.0025 * wafer_radius * cos(ang_cent_to_edge - offset_from_cent_angle),
                chip_center_xy[1] + 1.0025 * wafer_radius * sin(ang_cent_to_edge - offset_from_cent_angle),
                ang_edge_to_cent,
            ]

            ports_dict[(3 * i) + 0] = bot_edge_port_xy_rot
            ports_dict[(3 * i) + 1] = mid_edge_port_xy_rot
            ports_dict[(3 * i) + 2] = top_edge_port_xy_rot

            top_port_rect = gdspy.Rectangle(
                [top_edge_port_xy_rot[0] - wafer_radius, top_edge_port_xy_rot[1] + (holder_outer_port_width / 2)],
                [top_edge_port_xy_rot[0] + wafer_radius, top_edge_port_xy_rot[1] - (holder_outer_port_width / 2)],
                **self.Chip_holder,
            )
            top_port_rect.rotate(top_edge_port_xy_rot[2], [top_edge_port_xy_rot[0], top_edge_port_xy_rot[1]])
            holder_cutouts.add(top_port_rect)

            mid_port_rect = gdspy.Rectangle(
                [mid_edge_port_xy_rot[0] - wafer_radius, mid_edge_port_xy_rot[1] + (holder_outer_port_width / 2)],
                [mid_edge_port_xy_rot[0] + wafer_radius, mid_edge_port_xy_rot[1] - (holder_outer_port_width / 2)],
                **self.Chip_holder,
            )
            mid_port_rect.rotate(mid_edge_port_xy_rot[2], [mid_edge_port_xy_rot[0], mid_edge_port_xy_rot[1]])
            holder_cutouts.add(mid_port_rect)

            bot_port_rect = gdspy.Rectangle(
                [bot_edge_port_xy_rot[0] - wafer_radius, bot_edge_port_xy_rot[1] + (holder_outer_port_width / 2)],
                [bot_edge_port_xy_rot[0] + wafer_radius, bot_edge_port_xy_rot[1] - (holder_outer_port_width / 2)],
                **self.Chip_holder,
            )
            bot_port_rect.rotate(bot_edge_port_xy_rot[2], [bot_edge_port_xy_rot[0], bot_edge_port_xy_rot[1]])
            holder_cutouts.add(bot_port_rect)

            # offset_from_cent_angle = np.arcsin((18000/wafer_radius))

            # top_edge_port_TEST = [chip_center_xy[0] + wafer_radius*cos(ang_cent_to_edge + offset_from_cent_angle),
            #                       chip_center_xy[1] + wafer_radius*sin(ang_cent_to_edge + offset_from_cent_angle),
            #                       ang_edge_to_cent
            #                       ]

            # top_port_rect_TEST = gdspy.Rectangle([top_edge_port_TEST[0], top_edge_port_TEST[1] + (10/2)],
            #                                      [top_edge_port_TEST[0] + wafer_radius, top_edge_port_TEST[1] - (10/2)],
            #                                       layer=0)
            # self.Main.add(top_port_rect_TEST)

        holder_outer_hex = gdspy.Polygon(holder_outer_hex_points, **self.Chip_holder)
        holder_positives.add(holder_outer_hex)

        wafer_circle_cutout = self.make_silicon_wafer(
            wafer_diameter, wafer_flat_length, self.Chip_holder, chip_center=chip_center_xy, flat_angle=middle_flat_angle
        )
        holder_cutouts.add(wafer_circle_cutout)

        holder = gdspy.boolean(holder_positives, holder_cutouts, "not", **self.Chip_holder)
        self.Main.add(holder)

        return ports_dict

    def get_SOUK_boundary_hexagon(self):
        """
        Generates the boundary polygon points for the SOUK holder boundary.

        Returns
        -------
        polygon_points : list
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

    def make_silicon_wafer(self, wafer_diameter, flat_length, layer, chip_center=[0, 0], flat_angle=(-pi / 2)):
        """
        Makes the shape of the silicon wafer with a flat on one side and
        returns this shape as a PolygonSet in layer 0 or the layer specified.

        Parameters
        ----------
        wafer_diameter : float, int
            The diameter of the wafer in microns.

        flat_length : float, int
            The length of the flat on the silicon wafer.

        layer : dict
            The layer number and datatype to make the wafer shape in.
            e.g. {"layer": 1, "datatype": 0}.

        KwArgs
        ------
        chip_center : list
            List containing the [x,y] coordinate for the center of the chip.
            Default is [0,0].

        flat_angle : float, int
            Angle (**in radians**) the falt is located. Default is (-pi/2),
            ie the bottom middle side of the hexagon wafer. This angle is the
            angle the middle of the flat makes with the chip center.

        Returns
        -------
        final_wafer : gdspy PolygonSet
            The final wafer shape made with the cutout flat in it on the layer
            specified by the layer.

        """

        wafer_radius = wafer_diameter / 2

        wafer_cirlce = gdspy.Round(chip_center, wafer_radius)

        center_to_flat_dist = (wafer_radius**2 - ((flat_length**2) / 4)) ** 0.5
        flat_to_edge_dist = wafer_radius - center_to_flat_dist

        middle_of_flat = [center_to_flat_dist * cos(flat_angle), center_to_flat_dist * sin(flat_angle)]

        flat_cutout_rect = gdspy.Rectangle(
            [middle_of_flat[0] - (flat_length / 2), middle_of_flat[1]],
            [middle_of_flat[0] + (flat_length / 2), middle_of_flat[1] + flat_to_edge_dist],
        )

        flat_cutout_rect.rotate(flat_angle - (pi / 2), middle_of_flat)

        final_wafer = gdspy.boolean(wafer_cirlce, flat_cutout_rect, "not", **layer)

        return final_wafer

    def make_test_chip_quadrent_boundary_and_get_horn_positions(
        self,
        quadrent_center_xy,
        Main_config_file_dict,
        bottom_left_text="",
        top_right_text="",
        top_left_text="",
        cardiff_logo=True,
        top_right_label_window=True,
        return_outer_poly_points=False,
    ):
        """
        Make the boundary for a test chip quad. Adds a centered pin hole and
        slotted pin hole.

        Parameters
        ----------
        quadrent_center_xy : list
            list containing the [x,y] coordinates for the center of the test
            chip quad.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "test_chip_quad".

        KwArgs
        ------
        bottom_left_text : str
            Text to add the bottom left of the chip. Default is a blank str
            which will not add any text. When not empty this text will be
            placed offset from the bot left corner at a textsize and offset
            specified in the config.

        top_right_text : str
            Text to add the top right of the chip. Default is a blank str
            which will not add any text. When not empty this text will be
            placed offset from the top right corner at a textsize and offset
            specified in the config.

        top_left_text : str
            Text to add the top left of the chip. Default is a blank str which
            will not add any text. When not empty this text will be placed
            offset from the top left corner at a textsize and offset
            specified in the config.

        cardiff_logo : Boolean
            Default True will add the cardiff logo to the nb groundplane in the
            bottom right corner of the test chip quad. The size of this is
            specified in the config.

        top_right_label_window : Boolean
            Default True will add a window cutout in the nb groundplane in the
            top right corner of the test chip quad. The size of this window
            is specified in the config.

        return_outer_poly_points : Boolean
            Default False. If true this will return the outer polygon points
            around the test chip quadrent with extra exclusion.

        Returns
        -------
        horn_centers : list
            list of [x,y] lists that define the coordinates for the horn
            centers for the 4 horns on the test chip quad. The coords are for
            the top left, top right, bot left, bot right horns respectively.

        outer_poly_points : list
            list of [x,y] lists defining the coordinates of the outer polygon.
            **Note This is only returned if return_outer_poly_points KwArg is
            True**.

        """

        config = Main_config_file_dict["test_chip_quad"]

        test_chip_quad_width = config["test_chip_quad_width"]  # 27000
        test_chip_quad_height = config["test_chip_quad_height"]  # 27000

        horn_offset_from_center_x = config["horn_offset_from_center_x"]  # 6000
        horn_offset_from_center_y = config["horn_offset_from_center_y"]  # 6000

        top_right_label_window_offset_x = config["top_right_label_window_offset_x"]  # 2000
        top_right_label_window_offset_y = config["top_right_label_window_offset_y"]  # 500
        top_right_label_window_width = config["top_right_label_window_width"]  # 4000
        top_right_label_window_height = config["top_right_label_window_height"]  # 1000

        top_right_chip_corner = [quadrent_center_xy[0] + (test_chip_quad_width / 2), quadrent_center_xy[1] + (test_chip_quad_height / 2)]

        top_left_chip_corner = [quadrent_center_xy[0] - (test_chip_quad_width / 2), quadrent_center_xy[1] + (test_chip_quad_height / 2)]

        bottom_left_chip_corner = [quadrent_center_xy[0] - (test_chip_quad_width / 2), quadrent_center_xy[1] - (test_chip_quad_height / 2)]

        bottom_right_chip_corner = [quadrent_center_xy[0] + (test_chip_quad_width / 2), quadrent_center_xy[1] - (test_chip_quad_height / 2)]

        # Adding the bottom_left_text if True
        if bottom_left_text != "":
            bottom_left_text_offset_x = config["bottom_left_text_offset_x"]  # 1000
            bottom_left_text_offset_y = config["bottom_left_text_offset_y"]  # 1000
            bottom_left_text_size = config["bottom_left_text_size"]  # 90

            BL_text = gdspy.Text(
                bottom_left_text,
                bottom_left_text_size,
                position=[bottom_left_chip_corner[0] + bottom_left_text_offset_x, bottom_left_chip_corner[1] + bottom_left_text_offset_y],
                **self.Aluminium,
            )

            self.Main.add(BL_text)
            self.ground_plane_cutouts.add(gdspy.Rectangle(BL_text.get_bounding_box()[0], BL_text.get_bounding_box()[1]))

        # Adding the top_left_text if True
        if top_left_text != "":
            top_left_text_offset_x = config["top_left_text_offset_x"]  # 1000
            top_left_text_offset_y = config["top_left_text_offset_y"]  # 1000
            top_left_text_size = config["top_left_text_size"]  # 90

            TL_text = gdspy.Text(
                top_left_text,
                top_left_text_size,
                position=[top_left_chip_corner[0] + top_left_text_offset_x, top_left_chip_corner[1] - top_left_text_offset_y],
                **self.Aluminium,
            )

            self.Main.add(TL_text)
            self.ground_plane_cutouts.add(gdspy.Rectangle(BL_text.get_bounding_box()[0], BL_text.get_bounding_box()[1]))

        # Adding the top_right_text if True
        if top_right_text != "":
            top_right_text_offset_x = config["top_right_text_offset_x"]  # 1000
            top_right_text_offset_y = config["top_right_text_offset_y"]  # 1000
            top_right_text_size = config["top_right_text_size"]  # 90
            no_of_chars_with_post_gap = len(top_right_text) - 1
            length_of_text = top_right_text_size * (no_of_chars_with_post_gap * (8 / 9) + (5 / 9))

            TR_text = gdspy.Text(
                top_right_text,
                top_right_text_size,
                position=[
                    top_right_chip_corner[0] - top_right_text_offset_x - length_of_text,
                    top_right_chip_corner[1] - top_right_text_offset_y,
                ],
                **self.Aluminium,
            )

            self.Main.add(TR_text)
            self.ground_plane_cutouts.add(gdspy.Rectangle(TR_text.get_bounding_box()[0], TR_text.get_bounding_box()[1]))

        # Adding the top_right_label_window if True
        if top_right_label_window:
            top_right_label_window_rect = gdspy.Rectangle(
                [top_right_chip_corner[0] - top_right_label_window_offset_x, top_right_chip_corner[1] - top_right_label_window_offset_y],
                [
                    top_right_chip_corner[0] - top_right_label_window_offset_x - top_right_label_window_width,
                    top_right_chip_corner[1] - top_right_label_window_offset_y - top_right_label_window_height,
                ],
            )
            self.ground_plane_cutouts.add(top_right_label_window_rect)

        # Adding the cardiff logo if True
        if cardiff_logo:
            cardiff_logo_size = config["cardiff_logo_size"]  # 2000
            cardiff_logo_offset_x = config["cardiff_logo_offset_x"]  # 1000
            cardiff_logo_offset_y = config["cardiff_logo_offset_y"]  # 1000

            logo_xy = [
                bottom_right_chip_corner[0] - cardiff_logo_offset_x - (cardiff_logo_size / 2),
                bottom_right_chip_corner[1] + cardiff_logo_offset_y + (cardiff_logo_size / 2),
            ]

            self.add_cardiff_logo(logo_xy, cardiff_logo_size)

        # Adding the dicing line in the groundplane and the tabbed dicing line
        self.add_test_chip_quad_tabbed_dicing_line(quadrent_center_xy, test_chip_quad_width, test_chip_quad_height)

        groundplane_edge_cuts_top_mid_gap = config["groundplane_edge_cuts_top_mid_gap"]  # 5500
        groundplane_edge_cuts_bot_mid_gap = config["groundplane_edge_cuts_bot_mid_gap"]  # 5500
        groundplane_edge_cuts_right_mid_gap = config["groundplane_edge_cuts_right_mid_gap"]  # 0
        groundplane_edge_cuts_left_mid_gap = config["groundplane_edge_cuts_left_mid_gap"]  # 0

        self.add_test_chip_quad_groundplane_edge_cuts(
            quadrent_center_xy,
            test_chip_quad_width,
            test_chip_quad_height,
            top_mid_gap=groundplane_edge_cuts_top_mid_gap,
            bot_mid_gap=groundplane_edge_cuts_bot_mid_gap,
            right_mid_gap=groundplane_edge_cuts_right_mid_gap,
            left_mid_gap=groundplane_edge_cuts_left_mid_gap,
        )

        # Adding the center pin and slotted pin
        center_pin_radius = config["center_pin_radius"]  # 2020/2

        center_pin_xy = quadrent_center_xy

        slotted_pin_radius = config["slotted_pin_radius"]  # 2020/2
        slotted_pin_length = config["slotted_pin_length"]  # 1000/2
        slotted_pin_offset_x_from_center_pin = config["slotted_pin_offset_x_from_center_pin"]  # 8000
        slotted_pin_offset_y_from_center_pin = config["slotted_pin_offset_y_from_center_pin"]  # 0

        slotted_pin_xy = [
            quadrent_center_xy[0] + slotted_pin_offset_x_from_center_pin,
            quadrent_center_xy[1] + slotted_pin_offset_y_from_center_pin,
        ]

        _ = self.add_center_pin_and_get_bounding_points(quadrent_center_xy, center_pin_radius)

        _ = self.add_slotted_pin_and_get_bounding_points(slotted_pin_xy, slotted_pin_length, slotted_pin_radius, center_pin_xy)

        horn_centers = [
            [quadrent_center_xy[0] - horn_offset_from_center_x, quadrent_center_xy[1] + horn_offset_from_center_y],
            [quadrent_center_xy[0] + horn_offset_from_center_x, quadrent_center_xy[1] + horn_offset_from_center_y],
            [quadrent_center_xy[0] - horn_offset_from_center_x, quadrent_center_xy[1] - horn_offset_from_center_y],
            [quadrent_center_xy[0] + horn_offset_from_center_x, quadrent_center_xy[1] - horn_offset_from_center_y],
        ]

        extra_exclusion_around_test_chip_quadrent = 2000  # config

        outer_poly_points = [
            [
                quadrent_center_xy[0] - (test_chip_quad_width / 2) - extra_exclusion_around_test_chip_quadrent,
                quadrent_center_xy[1] + (test_chip_quad_height / 2) + extra_exclusion_around_test_chip_quadrent,
            ],
            [
                quadrent_center_xy[0] + (test_chip_quad_width / 2) + extra_exclusion_around_test_chip_quadrent,
                quadrent_center_xy[1] + (test_chip_quad_height / 2) + extra_exclusion_around_test_chip_quadrent,
            ],
            [
                quadrent_center_xy[0] + (test_chip_quad_width / 2) + extra_exclusion_around_test_chip_quadrent,
                quadrent_center_xy[1] - (test_chip_quad_height / 2) - extra_exclusion_around_test_chip_quadrent,
            ],
            [
                quadrent_center_xy[0] - (test_chip_quad_width / 2) - extra_exclusion_around_test_chip_quadrent,
                quadrent_center_xy[1] - (test_chip_quad_height / 2) - extra_exclusion_around_test_chip_quadrent,
            ],
        ]

        if return_outer_poly_points:
            return horn_centers, outer_poly_points

        return horn_centers

    def add_cardiff_logo(self, logo_xy, size=100):
        """
        Adds the cardiff logo as a cutout in the ng groundplane layer.
        Requires the logo gds file. Logo will be placed centered on the logo_xy
        given.

        Parameters
        ----------
        logo_xy : list
            list containing the [x,y] coordinate for the center of where the
            logo should be placed.

        KwArgs
        ------
        size : int, flaot
            The size of the sides of the cardiff logo square.
            The default is 100.

        """
        x = logo_xy[0]
        y = logo_xy[1]

        default_size = 100
        scale_factor = size / default_size

        logo_file_name = "Cardiff_Uni_Logo.gds"
        logo_file_path = os.path.dirname(os.path.realpath(__file__))

        try:
            lib = gdspy.GdsLibrary(infile=(logo_file_path + "\\" + logo_file_name))
        except FileNotFoundError:
            print("Logo gds file '" + logo_file_name + "' not found.")
            return
        except Exception:
            print("Issue with file " + logo_file_name + ".")
            return

        main_cell = lib.top_level()[0]
        pol_dict = main_cell.get_polygons(by_spec=True)
        text_polygons = pol_dict[(0, 0)]
        bound_polygons = pol_dict[(2, 0)]
        box_polygons = pol_dict[(3, 0)]
        outline_frame_polygons = pol_dict[(11, 0)]

        for i in range(len(text_polygons)):
            poly = gdspy.Polygon(text_polygons[i], **self.Nb_Groundplane)
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            self.Main.add(poly)

        for i in range(len(box_polygons)):
            poly = gdspy.Polygon(box_polygons[i], **self.Nb_Groundplane)
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            self.ground_plane_cutouts.add(poly)

        return

    def add_MLA_marker(self, x, y, materials, inner_lw=5, outer_lw=10):
        """
        Adds a singe or series of square markers with a cross inside centered
        on the x,y coord given. Each subsequent material will make a square
        larger than the previous by 1x the outer linewidth such that they are
        staggered. Inner crosses reamin the same and do not get staggered.

        Parameters
        ----------
        x, y : float, int
            x, y coordinate for the center of the MLA marker.

        materials : list
            list of materials for the marker out of. This should be a list of
            layer dicts, e.g. [souk.Aluminium, souk.Nb_Antenna, souk.SiN_dep].

        KwArgs
        ------
        inner_lw : float, int
            linewidth of the inner cross section of the marker.
            The default value = 5.

        outer_lw : float, int
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

        for i, material in enumerate(materials):
            lay = material["layer"]
            dtype = material["datatype"]

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

            square_path = gdspy.FlexPath(square_points, outer_lw, layer=lay, datatype=dtype)
            square_path_polys = self.get_polys_from_flexpath(square_path)
            for i in range(len(square_path_polys)):
                self.Main.add(gdspy.Polygon(square_path_polys[i], layer=lay, datatype=dtype))

            vert_box = gdspy.Rectangle(vert_box_points[0], vert_box_points[1], layer=lay, datatype=dtype)
            horz_box = gdspy.Rectangle(horz_box_points[0], horz_box_points[1], layer=lay, datatype=dtype)

            self.Main.add(vert_box)
            self.Main.add(horz_box)

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
        self.ground_plane_cutouts.add([cutout_around_marker])
        self.silicon_nitride_cutouts.add([cutout_around_marker])
        self.silicon_nitride_membrane_cutouts.add([cutout_around_marker])
        self.silicon_oxide_cutouts.add([cutout_around_marker])

        return

    def add_initial_allignment_markers(self, x, y, cross_length=300, linewidth=5, cutout_square_size=1000):
        """
        Adds a singe MLA marker cross centered on the x,y coord given. Also
        adds a cutout window in layers that will be deposited ontop of this
        (namely Nb_groundplane).

        Parameters
        ----------
        x, y : float, int
            x, y coordinate for the center of the MLA marker.

        KwArgs
        ------
        cross_length : float, int
            Height of the marker cross in the center.
            Default = 300.

        linewidth : float, int
            Linewidth of the inner cross.
            Default = 5.

        cutout_square_size : float, int
            The width and heigh of the outer square cutout window.
            Default = 1000.

        """

        vert_box = gdspy.Rectangle(
            [x - (linewidth / 2), y - (cross_length / 2)], [x + (linewidth / 2), y + (cross_length / 2)], layer=0, datatype=0
        )

        horz_box = gdspy.Rectangle(
            [x - (cross_length / 2), y - (linewidth / 2)], [x + (cross_length / 2), y + (linewidth / 2)], layer=0, datatype=0
        )
        self.Main.add(vert_box)
        self.Main.add(horz_box)

        cutout_around_marker_rect = gdspy.Rectangle(
            [x - (cutout_square_size / 2), y - (cutout_square_size / 2)], [x + (cutout_square_size / 2), y + (cutout_square_size / 2)]
        )
        self.ground_plane_cutouts.add(cutout_around_marker_rect)

        return

    def add_caliper_alignment_markers(self, x, y, rot, layer1, layer2, layer1_text, layer2_text):
        """
        Adds Allignment Markers to the mask in a cutout centered on the x,y
        given. This consists of 2 main corsses and a series of calipers made
        from the two materials given.

        Parameters
        ----------
        x, y : float, int
            The x,y coordinate to place the marker.

        rot : float, int
            The angle (**in degrees**) the marker should be rotated.

        layer1, layer2 : dict
            These should be layer dicts that contain keys for 'layer' and
            'datatype' e.g. souk.Aluminium, souk.Nb_Antenna, or even custom
            {"layer": 1234, "datatype": 0}.

        layer1_text, layer2_text : str
            These should be strings that are the text to be placed above and
            below in the marker.

        """

        # adding the Text
        bot_text = gdspy.Text(layer1_text, 100, (x - 100, y - 250), **layer1)
        top_text = gdspy.Text(layer2_text, 100, (x - 100, y + 150), **layer2)
        self.Main.add(bot_text)
        self.Main.add(top_text)

        # adding the cutout around the alligment markers
        cutout_around_marker = gdspy.Rectangle([-175.0, 300.0], [325.0, -300.0])
        cutout_around_marker.translate(x, y)
        self.ground_plane_cutouts.add(cutout_around_marker)

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
        main_cross1_poly = gdspy.Polygon(main_cross1_poly_points, **layer1)
        main_cross1_poly.translate(x, y)
        self.Main.add(main_cross1_poly)

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
        main_cross2_poly = gdspy.Polygon(main_cross2_poly_points, **layer2)
        main_cross2_poly.translate(x, y)
        self.Main.add(main_cross2_poly)

        # adding the cross and squares in the top left
        TL_sq = gdspy.Rectangle([x - 100, y + 100], [x - 75, y + 75], **layer1)
        BL_sq = gdspy.Rectangle([x - 100, y + 55], [x - 75, y + 30], **layer1)
        TR_sq = gdspy.Rectangle([x - 55, y + 100], [x - 30, y + 75], **layer1)
        BR_sq = gdspy.Rectangle([x - 55, y + 55], [x - 30, y + 30], **layer1)
        self.Main.add(TL_sq)
        self.Main.add(BL_sq)
        self.Main.add(TR_sq)
        self.Main.add(BR_sq)

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
        top_cross_poly = gdspy.Polygon(top_cross_poly_points, **layer2)
        top_cross_poly.translate(x, y)
        self.Main.add(top_cross_poly)

        # adding the bottom left rectangle markers
        bot_Rect1 = gdspy.Rectangle([x - 100, y - 100], [x - 90, y - 50], **layer1)
        bot_Rect2 = gdspy.Rectangle([x - 80, y - 100], [x - 70, y - 50], **layer1)
        bot_Rect3 = gdspy.Rectangle([x - 60, y - 100], [x - 50, y - 50], **layer1)
        bot_Rect4 = gdspy.Rectangle([x - 40, y - 100], [x - 30, y - 50], **layer1)

        bot_Rect5 = gdspy.Rectangle([x - 100, y - 70], [x - 90, y - 30], **layer2)
        bot_Rect6 = gdspy.Rectangle([x - 80, y - 70], [x - 70, y - 30], **layer2)
        bot_Rect7 = gdspy.Rectangle([x - 60, y - 70], [x - 50, y - 30], **layer2)
        bot_Rect8 = gdspy.Rectangle([x - 40, y - 70], [x - 30, y - 30], **layer2)

        self.Main.add(bot_Rect1)
        self.Main.add(bot_Rect2)
        self.Main.add(bot_Rect3)
        self.Main.add(bot_Rect4)
        self.Main.add(bot_Rect5)
        self.Main.add(bot_Rect6)
        self.Main.add(bot_Rect7)
        self.Main.add(bot_Rect8)

        # adding the bottom right rectangle markers
        left_Rect1 = gdspy.Rectangle([x + 100, y - 100], [x + 50, y - 90], **layer1)
        left_Rect2 = gdspy.Rectangle([x + 100, y - 80], [x + 50, y - 70], **layer1)
        left_Rect3 = gdspy.Rectangle([x + 100, y - 60], [x + 50, y - 50], **layer1)
        left_Rect4 = gdspy.Rectangle([x + 100, y - 40], [x + 50, y - 30], **layer1)

        left_Rect5 = gdspy.Rectangle([x + 70, y - 100], [x + 30, y - 90], **layer2)
        left_Rect6 = gdspy.Rectangle([x + 70, y - 80], [x + 30, y - 70], **layer2)
        left_Rect7 = gdspy.Rectangle([x + 70, y - 60], [x + 30, y - 50], **layer2)
        left_Rect8 = gdspy.Rectangle([x + 70, y - 40], [x + 30, y - 30], **layer2)

        self.Main.add(left_Rect1)
        self.Main.add(left_Rect2)
        self.Main.add(left_Rect3)
        self.Main.add(left_Rect4)
        self.Main.add(left_Rect5)
        self.Main.add(left_Rect6)
        self.Main.add(left_Rect7)
        self.Main.add(left_Rect8)

        # adding the callipers
        caliper = phgeom.litho_calipers(
            notch_size=[6, 40],
            notch_spacing=6,
            num_notches=5,
            offset_per_notch=0.5,
            row_spacing=-20.0,
            layer1=layer1["layer"],
            layer2=layer2["layer"],
        )
        caliper.move([x + 60, y + 60])

        caliper_polys = caliper.get_polygons()
        for i, poly_points in enumerate(caliper_polys):
            if i % 2 == 0:
                self.Main.add(gdspy.Polygon(poly_points, **layer1))
            else:
                self.Main.add(gdspy.Polygon(poly_points, **layer2))

        rotcaliper = phgeom.litho_calipers(
            notch_size=[6, 40],
            notch_spacing=6,
            num_notches=5,
            offset_per_notch=0.5,
            row_spacing=-20.0,
            layer1=layer1["layer"],
            layer2=layer2["layer"],
        )
        rotcaliper.move([x + 170, y + 25])
        rotcaliper.rotate(-90.0, (x + 170, y + 25))

        rotcaliper_polys = rotcaliper.get_polygons()
        for i, poly_points in enumerate(rotcaliper_polys):
            if i % 2 == 0:
                self.Main.add(gdspy.Polygon(poly_points, **layer1))
            else:
                self.Main.add(gdspy.Polygon(poly_points, **layer2))

        return

    def add_test_linewidth_structure_box_section(self, x, y, linewidths=[0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10]):
        """
        Adds a box of size 5500 centered on the x,y given that contains a
        series of lines to test the linewidth. Test linewidths are [0.25, 0.5,
        0.75, 1, 1.5, 2, 3, 5, 10]. Materials are ["Aluminium", "Nb_Antenna",
        "Nb_grnd", "SiN_Dep","IDC_Nb"].

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate to center the test linewidth structure.

        KwArgs
        ------
        linewidths : list
            This is a list of float or ints that define the linewidths in the
            test structure. By default these linewiths are:
            >>> [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10].



        """

        box_w_h = 5500
        main_square = gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2], **self.Nb_Groundplane)

        self.ground_plane_cutouts.add(main_square)

        outer_square_side_length = box_w_h
        outer_lw = 100

        outer_square_points = [
            [x - outer_square_side_length / 2 - outer_lw / 2, y + outer_square_side_length / 2 + outer_lw],
            [x - outer_square_side_length / 2 - outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
            [x + outer_square_side_length / 2 + outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
            [x + outer_square_side_length / 2 + outer_lw / 2, y + outer_square_side_length / 2 + outer_lw / 2],
            [x - outer_square_side_length / 2 - outer_lw, y + outer_square_side_length / 2 + outer_lw / 2],
        ]
        outer_square = gdspy.FlexPath(outer_square_points, outer_lw, **self.Nb_Antenna)
        self.Main.add(outer_square)

        gap = 50
        big_gap = 200

        # linewidths = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10]

        height = 4000
        text_height = 36
        label_height = 5 * text_height

        materials = ["Aluminium", "Nb_Antenna", "Nb_grnd", "SiN_Dep", "IDC_Nb"]

        lays = [
            self.Aluminium["layer"],
            self.Nb_Antenna["layer"],
            self.Nb_Groundplane["layer"],
            self.SiN_dep["layer"],
            self.IDC_Nb["layer"],
        ]

        dtypes = [
            self.Aluminium["datatype"],
            self.Nb_Antenna["datatype"],
            self.Nb_Groundplane["datatype"],
            self.SiN_dep["datatype"],
            self.IDC_Nb["datatype"],
        ]

        y_pos = y
        init_x = (
            x
            - (
                len(materials) * (np.sum(linewidths) + (len(linewidths) - 1) * gap)
                + (len(materials) - 1) * (big_gap)
                + len(materials) * (9.5 / 9) * label_height
            )
            / 2
            + ((9.5 / 9) * label_height + 0.125)
        )

        for (
            j,
            material,
        ) in enumerate(materials):
            x_pos = init_x + j * (big_gap + (len(linewidths) - 1) * 50 + np.sum(linewidths) + (9.5 / 9) * label_height)

            for i, lw in enumerate(linewidths):
                if i == 0:
                    self.Main.add(
                        gdspy.Text(
                            material,
                            label_height,
                            [x_pos - gap - (lw / 2) + ((2 / 9) * label_height), y_pos - height / 2 + text_height],
                            angle=pi / 2,
                            layer=lays[j],
                            datatype=dtypes[j],
                        )
                    )

                self.Main.add(
                    gdspy.Rectangle(
                        [x_pos - lw / 2, y_pos - height / 2], [x_pos + lw / 2, y_pos + height / 2], layer=lays[j], datatype=dtypes[j]
                    )
                )

                self.Main.add(
                    gdspy.Text(
                        str(lw) + "um",
                        text_height,
                        [x_pos - ((5.5 / 9) * text_height), y_pos - height / 2 - ((1 / 9) * text_height)],
                        angle=3 * pi / 2,
                        layer=lays[j],
                        datatype=dtypes[j],
                    )
                )

                if i < len(linewidths) - 1:
                    x_pos += gap + lw / 2 + linewidths[i + 1] / 2

        return

    def add_test_crossover_structure_box_section(self, x, y, Main_config_file_dict, rot=0.0):
        """
        Adds a box of size 5500 centered on the x,y given that contains a
        series of 3 crossover sections that connect to pads.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate to center the test crossover structure.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "filter_bank_ring_overlap".

        KwArgs
        ------
        rot : float, int
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
        self.silicon_nitride_positives.add(
            gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2], **self.SiN_dep)
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
        outer_square = gdspy.FlexPath(outer_square_points, outer_lw, **self.Nb_Antenna)
        outer_square.rotate(rot, [x, y])
        self.Main.add(outer_square)

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
            cross_xy = self.rotate(x, y, cross_xy[0], cross_xy[1], rot)
            cross_connections = self.add_filter_bank_ring_overlap_and_get_conections(
                cross_xy[0], cross_xy[1], cross_xy[0] + 100, cross_xy[1], 0, cross_angle + rot, Main_config_file_dict
            )
            ground_under_cross = gdspy.Rectangle(
                [cross_xy[0] - (ground_under_crossover_w_h / 2), cross_xy[1] - (ground_under_crossover_w_h / 2)],
                [cross_xy[0] + (ground_under_crossover_w_h / 2), cross_xy[1] + (ground_under_crossover_w_h / 2)],
                **self.Nb_Groundplane,
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

            cross_connect_rots = [rot + 3 * pi / 4, rot + pi / 4, rot - 3 * pi / 4, rot - pi / 4]

            for i in range(4):
                pad_no += 1

                pad_points = pad_right_poly_points if pad_no % 2 == 0 else pad_left_poly_points
                pad_cutout_points = pad_right_cutout_poly_points if pad_no % 2 == 0 else pad_left_cutout_poly_points
                pad_connect = pad_right_connect if pad_no % 2 == 0 else pad_left_connect
                pad_connect_rot = rot + pi if pad_no % 2 == 0 else rot

                pad_points = self.move_points_list(pad_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
                pad_cutout_points = self.move_points_list(pad_cutout_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
                pad_connect = self.move_single_point(pad_connect, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])

                pad_points = self.rotate_points_list(pad_points, rot, x, y)
                pad_cutout_points = self.rotate_points_list(pad_cutout_points, rot, x, y)
                pad_connect = self.rotate(x, y, pad_connect[0], pad_connect[1], rot)

                pad_poly = gdspy.Polygon(pad_points, **self.Nb_Antenna)
                pad_cutout_poly = gdspy.Polygon(pad_cutout_points, self.SiN_dep)

                self.Main.add(pad_poly)
                self.silicon_nitride_cutouts.add(pad_cutout_poly)

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
                    **self.Nb_Antenna,
                )
                self.Main.add(pad1_path)

        final_groundplane_cutout = gdspy.boolean(working_groundplane_cutout, main_square, "and")
        self.ground_plane_cutouts.add(final_groundplane_cutout)

        return

    def add_test_straight_line_structure_box_section(self, x, y, rot=0.0):
        """
        Adds a box of size 5500 centered on the x,y given that contains a
        series of 6 stright line sections that connect to pads.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate to center the test crossover structure.

        KwArgs
        ------
        rot : float, int
            The angle (**in Radians**) the structure should be rotated.
            Default is 0 which has the pads on the left and right. This will
            only take the values in the range 0 to pi/2.

        """

        box_w_h = 5500

        main_square = gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2])
        main_square.rotate(rot, [x, y])

        working_groundplane_cutout = main_square
        # working_sin_dep_cutout = main_square
        self.silicon_nitride_positives.add(
            gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2], **self.SiN_dep)
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
        outer_square = gdspy.FlexPath(outer_square_points, outer_lw, **self.Nb_Antenna)
        outer_square.rotate(rot, [x, y])
        self.Main.add(outer_square)

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

                pad_points = self.move_points_list(pad_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
                pad_cutout_points = self.move_points_list(pad_cutout_points, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])
                pad_connect = self.move_single_point(pad_connect, pads_offset[i][0], set_y_offset[set_no] + pads_offset[i][1])

                pad_points = self.rotate_points_list(pad_points, rot, x, y)
                pad_cutout_points = self.rotate_points_list(pad_cutout_points, rot, x, y)
                pad_connect = self.rotate(x, y, pad_connect[0], pad_connect[1], rot)

                straight_line_path_points = [
                    pad_connect,
                    [
                        pad_connect[0] + straight_line_length * cos(pad_connect_rot),
                        pad_connect[1] + straight_line_length * sin(pad_connect_rot),
                    ],
                ]
                straight_line_path = gdspy.FlexPath(straight_line_path_points, straight_line_width, gdsii_path=True, **self.Nb_Antenna)
                self.Main.add(straight_line_path)

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
                    **self.Aluminium,
                )
                self.Main.add(al_pad_rect)

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
                    **self.Aluminium,
                )
                al_straight_line_polys = self.get_polys_from_flexpath(al_straight_line_path)
                for i in range(len(al_straight_line_polys)):
                    self.Main.add(gdspy.Polygon(al_straight_line_polys[i], **self.Aluminium))

                self.Main.add(al_straight_line_path)

                pad_poly = gdspy.Polygon(pad_points, **self.Nb_Antenna)
                pad_cutout_poly = gdspy.Polygon(pad_cutout_points, **self.SiN_dep)

                self.Main.add(pad_poly)
                self.silicon_nitride_cutouts.add(pad_cutout_poly)

        final_groundplane_cutout = gdspy.boolean(working_groundplane_cutout, main_square, "and")
        self.ground_plane_cutouts.add(final_groundplane_cutout)

        return

    def make_DC_structure_pads_and_meander(self, x, y, rot, DC_structure_material):
        """
        This Makes two pads and a meander connecting them. Pads are 100x the
        size of the linewidth connecting them.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate to center the test DC structure.

        rot : float, int
            The angle (**in degrees**) the pad and meander structure should be
            rotated.

        DC_structure_material : dict
            This is the layer the pad and meander structure should be made of.
            This should be a layer dict that contains keys for 'layer' and
            'datatype' e.g. souk.Aluminium, souk.Nb_Antenna, or even custom
            {"layer": 1234, "datatype": 0}.

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

        pad_left = gdspy.Rectangle(
            [pad_center_x_left - (pad_size / 2), -(pad_size / 2)], [pad_center_x_left + (pad_size / 2), +(pad_size / 2)]
        )

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
            new_poly_points = self.rotate_and_move_points_list(poly_points, rot, x, y)
            poly = gdspy.Polygon(new_poly_points, **DC_structure_material)
            self.Main.add(poly)

        return

    def add_test_DC_structure_box_section(self, x, y, DC_structure_material, BL_text=""):
        """
        Adds a box of size 8000 centered on the x,y given that contains a
        series of 5 test DC lines and pads. Pads are 100x the size of the
        linewidth connecting them.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate to center the test DC structure.

        DC_structure_material : dict
            This is the layer the pad and meander structure should be made of.
            This should be a layer dict that contains keys for 'layer' and
            'datatype' e.g. souk.Aluminium, souk.Nb_Antenna, or even custom
            {"layer": 1234, "datatype": 0}.

        KwArgs
        ------
        BL_text = ""
            This is the string that will be written in the bottom left of the
            structure. By default no string will be written.
        """

        box_w_h = 8000

        main_square = gdspy.Rectangle([x - box_w_h / 2, y - box_w_h / 2], [x + box_w_h / 2, y + box_w_h / 2])

        self.ground_plane_cutouts.add(main_square)

        outer_square_side_length = box_w_h
        outer_lw = 100

        outer_square_points = [
            [x - outer_square_side_length / 2 - outer_lw / 2, y + outer_square_side_length / 2 + outer_lw],
            [x - outer_square_side_length / 2 - outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
            [x + outer_square_side_length / 2 + outer_lw / 2, y - outer_square_side_length / 2 - outer_lw / 2],
            [x + outer_square_side_length / 2 + outer_lw / 2, y + outer_square_side_length / 2 + outer_lw / 2],
            [x - outer_square_side_length / 2 - outer_lw, y + outer_square_side_length / 2 + outer_lw / 2],
        ]
        outer_square = gdspy.FlexPath(outer_square_points, outer_lw, **self.Nb_Antenna)
        self.Main.add(outer_square)

        pad_and_meander_offset = 2800

        centers_and_rots = [
            [0, pad_and_meander_offset, 0],
            [-pad_and_meander_offset, 0, pi / 2],
            [0, -pad_and_meander_offset, 0],
            [pad_and_meander_offset, 0, pi / 2],
            [0, 0, pi / 4],
        ]

        for dx, dy, rot in centers_and_rots:
            self.make_DC_structure_pads_and_meander(x + dx, y + dy, rot, DC_structure_material)

        if BL_text != "":
            text_size = 300
            BL_text_x_offset = 50
            BL_text_y_offset = 50
            text_xy = [x - box_w_h / 2 + BL_text_x_offset, y - box_w_h / 2 + BL_text_y_offset]

            text = gdspy.Text(BL_text, text_size, position=text_xy, **DC_structure_material)
            self.Main.add(text)

        return

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
        """
        Adds edge cutouts in the groundplane layer around the border of the
        test chip quad. These edge cuts sit within the defined border. There
        is by default edge cuts across every edge but there can be a gap in the
        middle of any edge specified by the #_mid_gap arguments. The size and
        spacing of these edge cuts is defined in config.

        Parameters
        ----------
        test_chip_quad_center_xy : list
            list containing the [x,y] coordinates for the center of the test
            chip quad.

        test_chip_quad_width : int, float
            The width of the test chip quad.

        test_chip_quad_height : int, float
            The height of the test chip quad

        KwArgs
        ------
        top_mid_gap, bot_mid_gap, right_mid_gap, left_mid_gap : float, int
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
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
                    **self.Nb_Groundplane,
                )
                self.ground_plane_cutouts.add(left_edge_cut_rect_bot)

        return

    def add_test_chip_quad_tabbed_dicing_line(
        self, test_chip_quad_center_xy, test_chip_quad_width, test_chip_quad_height, tab_position="left", corner_overrun_length=200
    ):
        """
        Adds a tabbed dicing line around the test chip quad. This is a solid
        dicing line of three edges and a tabbed edge on the other. This tab
        position is by default on the left side but can be any edge. The dice
        lines by default overextend by 200um but this can be changed.

        Parameters
        ----------
        test_chip_quad_center_xy : list
            list containing the [x,y] coordinates for the center of the test
            chip quad.

        test_chip_quad_width : int, float
            The width of the test chip quad.

        test_chip_quad_height : int, float
            The height of the test chip quad.

        KwArgs
        ------
        tab_position : str
            The edge to position the tabbed edge of the dicing line on.
            The default is "left" but can take str values "top", "left", "bot",
            or "right".

        corner_overrun_length : int, float
            This is the length the corners of the tabbed dicing line should
            overrun at the corners of the test chip quad. The default is 200.

        """

        tab_length = 4000  # config
        tab_gap_length = 300  # config
        linewidth = 300  # config

        tab_rot_dict = {"top": (-pi / 2), "left": 0, "bot": (pi / 2), "right": pi}

        if tab_position in tab_rot_dict.keys():
            rotation = tab_rot_dict[tab_position]
        else:
            raise ValueError("tab_position arg should be a str that takes the value of 'right' 'left' 'top' or 'bot'")

        center_x = test_chip_quad_center_xy[0]
        center_y = test_chip_quad_center_xy[1]

        top_left_corner = [center_x - test_chip_quad_width / 2, center_y + test_chip_quad_height / 2]
        top_right_corner = [center_x + test_chip_quad_width / 2, center_y + test_chip_quad_height / 2]
        bot_left_corner = [center_x - test_chip_quad_width / 2, center_y - test_chip_quad_height / 2]
        bot_right_corner = [center_x + test_chip_quad_width / 2, center_y - test_chip_quad_height / 2]

        right_rect = gdspy.Rectangle(
            [bot_right_corner[0], bot_right_corner[1] - linewidth - corner_overrun_length],
            [top_right_corner[0] + linewidth, top_right_corner[1] + linewidth + corner_overrun_length],
            **self.Tab_dicing_line,
        )
        right_rect.rotate(rotation, [center_x, center_y])
        self.Main.add(right_rect)

        top_rect = gdspy.Rectangle(
            [top_left_corner[0] - linewidth - corner_overrun_length, top_left_corner[1]],
            [top_right_corner[0] + linewidth + corner_overrun_length, top_right_corner[1] + linewidth],
            **self.Tab_dicing_line,
        )
        top_rect.rotate(rotation, [center_x, center_y])
        self.Main.add(top_rect)

        bot_rect = gdspy.Rectangle(
            [bot_left_corner[0] - linewidth - corner_overrun_length, bot_left_corner[1] - linewidth],
            [bot_right_corner[0] + linewidth + corner_overrun_length, bot_right_corner[1]],
            **self.Tab_dicing_line,
        )
        bot_rect.rotate(rotation, [center_x, center_y])
        self.Main.add(bot_rect)

        middle_left_center = [center_x - test_chip_quad_width / 2 - linewidth / 2, center_y]

        path_middle = gdspy.FlexPath(
            [
                [middle_left_center[0], middle_left_center[1] - tab_length / 2],
                [middle_left_center[0], middle_left_center[1] + tab_length / 2],
            ],
            linewidth,
            ends="round",
            **self.Tab_dicing_line,
        )
        path_middle.rotate(rotation, [center_x, center_y])
        self.Main.add(path_middle)

        for i in range(int(np.round((test_chip_quad_height / 2) / (tab_length + tab_gap_length + linewidth), 0)) + 1):
            if (tab_length / 2 + i * (tab_length + 2 * tab_gap_length)) > test_chip_quad_height / 2:
                edge_tab = gdspy.FlexPath(
                    [
                        [middle_left_center[0], middle_left_center[1] - tab_length / 2 + i * (tab_length + 2 * tab_gap_length)],
                        [middle_left_center[0], top_left_corner[1] + linewidth + corner_overrun_length],
                    ],
                    linewidth,
                    **self.Tab_dicing_line,
                )
                edge_tab.rotate(rotation, [center_x, center_y])
                self.Main.add(edge_tab)

                cap = gdspy.Round(
                    [middle_left_center[0], middle_left_center[1] - tab_length / 2 + i * (tab_length + 2 * tab_gap_length)],
                    linewidth / 2,
                    initial_angle=0,
                    final_angle=-pi,
                    **self.Tab_dicing_line,
                )
                cap.rotate(rotation, [center_x, center_y])
                self.Main.add(cap)

                edge_tab = gdspy.FlexPath(
                    [
                        [middle_left_center[0], bot_left_corner[1] - linewidth - corner_overrun_length],
                        [middle_left_center[0], middle_left_center[1] + tab_length / 2 + -i * (tab_length + 2 * tab_gap_length)],
                    ],
                    linewidth,
                    **self.Tab_dicing_line,
                )
                edge_tab.rotate(rotation, [center_x, center_y])
                self.Main.add(edge_tab)

                cap = gdspy.Round(
                    [middle_left_center[0], middle_left_center[1] + tab_length / 2 + -i * (tab_length + 2 * tab_gap_length)],
                    linewidth / 2,
                    initial_angle=0,
                    final_angle=pi,
                    **self.Tab_dicing_line,
                )
                cap.rotate(rotation, [center_x, center_y])
                self.Main.add(cap)

                break

            edge_tab = gdspy.FlexPath(
                [
                    [middle_left_center[0], middle_left_center[1] - tab_length / 2 + i * (tab_length + 2 * tab_gap_length)],
                    [middle_left_center[0], middle_left_center[1] + tab_length / 2 + i * (tab_length + 2 * tab_gap_length)],
                ],
                linewidth,
                ends="round",
                **self.Tab_dicing_line,
            )
            edge_tab.rotate(rotation, [center_x, center_y])
            self.Main.add(edge_tab)

            edge_tab = gdspy.FlexPath(
                [
                    [middle_left_center[0], middle_left_center[1] - tab_length / 2 + -i * (tab_length + 2 * tab_gap_length)],
                    [middle_left_center[0], middle_left_center[1] + tab_length / 2 + -i * (tab_length + 2 * tab_gap_length)],
                ],
                linewidth,
                ends="round",
                **self.Tab_dicing_line,
            )
            edge_tab.rotate(rotation, [center_x, center_y])
            self.Main.add(edge_tab)

        return

    def add_text_under_horn_in_test_chip_quad(self, text_string, horn_center_x, horn_center_y):
        """
        Adds a label under the horn on the test chip quad with the text string
        given.

        Parameters
        ----------
        text_string : str
            The text string to be placed under the horn.

        horn_center_x, horn_center_x : float, int
            The x, and y coordinate respectively for the center of the horn
            that needs the text underneath.

        """
        x_offset_from_horn = 0
        y_offset_from_horn = -3100
        text_size = 300

        under_text = gdspy.Text(
            text_string, text_size, position=[horn_center_x + x_offset_from_horn, horn_center_y + y_offset_from_horn], **self.Aluminium
        )
        under_text.translate(-abs(under_text.get_bounding_box()[0][0] - under_text.get_bounding_box()[1][0]) / 2, 0)
        self.Main.add(under_text)
        self.ground_plane_cutouts.add(
            gdspy.Rectangle(under_text.get_bounding_box()[0], under_text.get_bounding_box()[1], **self.Nb_Groundplane)
        )

        return

    def get_config_excel(self, file_name, sheet_name, path="config_files\\"):
        """
        gets the config excel file containing the variable names and respective
        vales. different sheets contain different vairables.

        Parameters
        ----------
        file_name : str
            String including file extention, eg "myfile.xlsx"
        sheet_name : str
            String for the name of the sheet in the excel file.
        path : str
            Path to containing folder. Deafult is "config_files" folder in same
            directory as this file.

        Returns
        -------
        out : 'config_dict'
            Dictionary with keys of var names and values of var values.

        Raises
        ------
            FileNotFoundError: If files doesn't exist.
            Exception: For all other file opening issues.
        """

        live_file = path + file_name

        try:
            df = pd.read_excel(live_file, sheet_name=sheet_name, header=None)
        except FileNotFoundError:
            print("Config file '" + path + file_name + "' not found.")
            return
        except Exception:
            print("Issue with file " + path + file_name + ".")
            return

        config_dict = {}
        for i in range(len(df)):
            var_name = df[0][i]
            var_val = df[1][i]
            config_dict[var_name] = var_val

        return config_dict

    def deg_to_rad(self, deg):
        """
        Converts an angle given in degrees to radiands.
        """
        return (deg / 180.0) * pi

    def rad_to_deg(self, rad):
        """
        Converts an angle given in radiands to degrees.
        """
        return rad * (180.0 / pi)

    def rotate(self, ox, oy, px, py, angle):
        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """
        qx = ox + cos(angle) * (px - ox) - sin(angle) * (py - oy)
        qy = oy + sin(angle) * (px - ox) + cos(angle) * (py - oy)
        return qx, qy

    def rotate_points_list(self, points, rot_angle, ox=0, oy=0):
        """
        Rotates a list of points lists all by angle.

        Rotates counterclockwise by a given angle around a given origin.
        The angle should be given in radians.
        """
        points_copy = points.copy()
        new_points = []
        for i in range(len(points_copy)):  # rotating each element
            old = list(points_copy[i])
            px = old[0]
            py = old[1]
            qx = ox + cos(rot_angle) * (px - ox) - sin(rot_angle) * (py - oy)
            qy = oy + sin(rot_angle) * (px - ox) + cos(rot_angle) * (py - oy)
            new_points.append([qx, qy])

        return new_points

    def rotate_and_move_points_list(self, points, rot_angle, dx, dy, ox=0, oy=0):
        """
        Rotates and moves a list of points lists all by angle and dx,dy.

        Rotates counterclockwise by a given angle around a given origin.
        The angle should be given in radians.

        The rotation is done before the translation.
        """
        points_copy = points.copy()
        new_points = []
        for i in range(len(points_copy)):  # shifting each element
            old = list(points_copy[i])
            px = old[0]
            py = old[1]
            qx = ox + cos(rot_angle) * (px - ox) - sin(rot_angle) * (py - oy)
            qy = oy + sin(rot_angle) * (px - ox) + cos(rot_angle) * (py - oy)
            new_points.append([qx + dx, qy + dy])

        return new_points

    def move_points_list(self, points, dx, dy):
        """
        Moves a list of points lists all by dx,dy.

        Parameters
        ----------
        points : list
            list of [x,y] lists to be moved.

        dx, dy : float, int
            The delta x, delta y to move the points in the points list.

        Returns
        -------
        new_points : list
            Returns a new list of points that have been Moved. List is of the
            same form as input points list.

        """
        points_copy = points.copy()
        new_points = []
        for i in range(len(points_copy)):  # shifting each element
            old = list(points_copy[i])
            px = old[0]
            py = old[1]
            new_points.append([px + dx, py + dy])

        return new_points

    def move_single_point(self, point, dx, dy):
        """
        Moves an [x,y] point by dx,dy.

        Parameters
        ----------
        point : list
            Point [x,y] list to be moved.

        dx, dy : float, int
            The delta x, delta y to move the point by.

        Returns
        -------
        new_point : list
            Returns a new point list that has been Moved. List is of the
            same form as input point list.

        """
        point_copy = point.copy()
        px = point_copy[0]
        py = point_copy[1]
        new_point = [px + dx, py + dy]

        return new_point

    def rotate_and_move_single_point(self, point, rot_angle, dx, dy, ox=0, oy=0):
        """
        Rotates and moves a list of tupple points all by angle and dx,dy.

        Rotates counterclockwise by a given angle around a given origin.
        The angle should be given in radians.
        """

        old = point
        px = old[0]
        py = old[1]
        qx = ox + cos(rot_angle) * (px - ox) - sin(rot_angle) * (py - oy)
        qy = oy + sin(rot_angle) * (px - ox) + cos(rot_angle) * (py - oy)
        new_point = [qx + dx, qy + dy]

        return new_point

    def mirror_points_around_yaxis(self, points):
        """
        mirrors a set of points around the y axis.

        points : list
            list of [x,y] lists defining the x,y for each point to mirror.
        """
        points = points.copy()
        for i in range(len(points)):
            old = list(points[i])
            px = -old[0]
            py = old[1]
            points[i] = [px, py]

        return points

    def mirror_points_around_xaxis(self, points):
        """
        mirrors a set of points around the x axis.

        points : list
            list of [x,y] lists defining the x,y for each point to mirror.
        """
        points = points.copy()
        for i in range(len(points)):
            old = list(points[i])
            px = old[0]
            py = -old[1]
            points[i] = [px, py]

        return points

    def create_miter_join(self, p0, v0, p1, v1, p2, w):
        """
        creates a miter on the corner of the bends in a gdspy Flexpath object.
        takes 6 arguments (vertex and direction vector from both segments being
        joined, the center and width of the path) and return a list of vertices
        that make the join.

        Parameters
        ----------
        p0 : array-like[2]
            Vertex [x, y] describing the end of the first path segment.
        v0 : rray-like[2]
            Vector [dx, dy] describing the first path section.
        p1 : array-like[2]
            Vertex [x, y] describing the start of the second path segment.
        v1 : array-like[2]
            Vector [dx, dy] describing the second path section.
        p2 : array-like[2]
            Vertex [x, y] describing center of the joins between the path segments.
        w : int, deciaml
            The width of the path.

        Returns
        ----------
        out : list
            List of array-like[2] describing the vertices that make the corner join.
        """

        # Calculate intersection point p between lines defined by
        # p0 + u0 * v0 (for all u0) and p1 + u1 * v1 (for all u1)
        den = v1[1] * v0[0] - v1[0] * v0[1]
        lim = 1e-12 * (v0[0] ** 2 + v0[1] ** 2) * (v1[0] ** 2 + v1[1] ** 2)
        if den**2 < lim:
            # Lines are parallel: use mid-point
            u0 = u1 = 0
            p = 0.5 * (p0 + p1)
        else:
            dx = p1[0] - p0[0]
            dy = p1[1] - p0[1]
            u0 = (v1[1] * dx - v1[0] * dy) / den
            u1 = (v0[1] * dx - v0[0] * dy) / den
            p = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
        if u0 <= 0 and u1 >= 0:
            # Inner corner
            return [p]
        # Outer corner
        angle0 = np.arctan2(v0[1], v0[0])
        angle1 = np.arctan2(v1[1], v1[0])

        return [
            [p0[0] + (w / 2) * cos(angle0 + pi), p0[1] + (w / 2) * sin(angle0 + pi)],
            [p1[0] + (w / 2) * cos(angle1), p1[1] + (w / 2) * sin(angle1)],
        ]

    def get_polys_from_flexpath(self, flexpath):
        """
        Gets the polygon points for the shape created by a gdspy Flexpath object.

        Parameters
        ----------
        flexpath : object
            an instanciated gdspy.Flexpath object.

        Returns
        -------
        polys : list
            list of [x,y] lists defining the points in the polygon.

        """
        polys = []
        for i in range(len(flexpath.get_polygons())):
            poly_set = []
            for k in range(len(flexpath.get_polygons()[i])):
                poly_set.append([flexpath.get_polygons()[i][k][0], flexpath.get_polygons()[i][k][1]])

            polys.append(poly_set)

        return polys

    def add_antena(self, x, y, rot, Main_config_file_dict, add_grnd_cutout=True, add_SiN_dep_cutout=True, add_backside_check=True):
        """
        Adds the antenna geometries to the Main cell. These are the 4 antennas in
        at the center of each horn block. Optionally adds a ground plane cutout
        and Silicon Nitride deoposition cutout as a circle around the antenna
        structure. These geometries are defined by the dimensions within the
        Main_config_file_dict.

        Parameters
        ----------
        x,y : float, int
            The x,y coordinates about which to center the antenna structure.

        rot : float, int
            The angle (**in degrees**) the antenna structure should be rotated.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "antenna".

        KwArgs
        ------
        add_grnd_cutout=True
            Whether or not to add a circular cutout in the Nb_Groundplane layer
            around the antenna structure.

        add_SiN_dep_cutout=True
            Whether or not to add a circular cutout in the SiN depositon layer
            around the antenna structure.

        add_backside_check=True
            Whether or not to add a circle in the backside check layer over the
            antenna structure.

        """
        config = Main_config_file_dict["antenna"]

        distance_from_center = config["distance_from_center"]
        base_width = config["base_width"]
        top_conect_width = config["top_conect_width"]
        straight_height = config["straight_height"]
        taper_height = config["taper_height"]
        backside_check_circle_radius = config["backside_check_circle_radius"]

        default_cicle_cutout_radius = distance_from_center + straight_height + taper_height
        # Checking for cutout circle radii in config to ensure backward compatability with old config files. Reverts to the old default value if not existing
        if "sin_dep_cutout_circle_radius" in config:
            sin_dep_cutout_circle_radius = config["sin_dep_cutout_circle_radius"]
        else:
            sin_dep_cutout_circle_radius = default_cicle_cutout_radius

        if "grnd_cutout_circle_radius" in config:
            grnd_cutout_circle_radius = config["grnd_cutout_circle_radius"]
        else:
            grnd_cutout_circle_radius = default_cicle_cutout_radius

        antena_geometry = [
            [(top_conect_width / 2) + x, -(taper_height + straight_height + distance_from_center) + y],
            [-(top_conect_width / 2) + x, -(taper_height + straight_height + distance_from_center) + y],
            [-(base_width / 2) + x, -(straight_height + distance_from_center) + y],
            [-(base_width / 2) + x, -(distance_from_center) + y],
            [(base_width / 2) + x, -(distance_from_center) + y],
            [(base_width / 2) + x, -(straight_height + distance_from_center) + y],
        ]  # this is the shape of the antenna with the point at the bottom

        ant_bot = gdspy.Polygon(
            antena_geometry, **self.Nb_Antenna
        )  # this defines the bottom, right, top and left antennas as the same bottom antenna shape
        ant_right = gdspy.Polygon(antena_geometry, **self.Nb_Antenna)
        ant_top = gdspy.Polygon(antena_geometry, **self.Nb_Antenna)
        ant_left = gdspy.Polygon(antena_geometry, **self.Nb_Antenna)

        ant_bot.rotate(self.deg_to_rad(rot), (x, y))  # this then roatates the antenna by the roation passed into the method
        ant_top.rotate(
            self.deg_to_rad(rot + 180), (x, y)
        )  # this also rotates the top, left and right antennas which are all the same as the bottom to form 4 antennas orthogonal to each other
        ant_left.rotate(self.deg_to_rad(rot + 270), (x, y))
        ant_right.rotate(self.deg_to_rad(rot + 90), (x, y))

        self.Main.add(ant_bot)  # adding the antennas to the main cell
        self.Main.add(ant_top)
        self.Main.add(ant_left)
        self.Main.add(ant_right)

        # Adding the circular antenna cutout to the grnd plane
        if add_grnd_cutout:
            grnd_cutout_circle = gdspy.Round([x, y], grnd_cutout_circle_radius, **self.Nb_Groundplane)
            self.ground_plane_cutouts.add(grnd_cutout_circle)
        if add_SiN_dep_cutout:
            SiN_dep_cutout_circle = gdspy.Round([x, y], sin_dep_cutout_circle_radius, **self.SiN_dep)
            self.silicon_nitride_cutouts.add(SiN_dep_cutout_circle)
        if add_backside_check:
            backside_check_circle = gdspy.Round([x, y], backside_check_circle_radius, **self.Backside_Check)
            self.Main.add(backside_check_circle)

            mirrored_y_backside_check_circle = gdspy.Round([-x, y], backside_check_circle_radius, **self.Backside_Check)
            self.MainBackside.add(mirrored_y_backside_check_circle)

        return

    def add_4_KID_outers_and_arms_new(
        self,
        x,
        y,
        rel_kid_positions,
        KID_Nos,
        f0s,
        IDC_and_CC_function,
        Main_config_file_dict,
        IDC_and_frame_materials=None,
        meander_materials=None,
        trim_lengths=None,
        add_grnd_cutout=True,
        add_SiN_dep_dielectric_around=True,
        add_SiN_dep_dielectric_cutout=True,
        add_SiO_cutout=True,
        add_SiN_membrane_cutout=True,
        add_backside_check=True,
    ):
        """
        Adds the 4 KIDs geometries to the Main cell. The four kids are placed
        around the center x, y point in order top right, top left, bot right,
        bot left (as should be the order of the rel_kid_positions list).
        The KIDs geometries are defined by the dimensions within the
        Main_config_file_dict. By default it will, but optionally can choose not
        to, add all the neccessay cutouts for the structure.

        Parameters
        ----------
        x,y : float, int
            The x,y coordinates about which to center the antenna structure.

        rel_kid_positions : list
            list of [x,y] lists defining the points at which to make each KID
            relative to the center of the antenna. This KIDs are made from the very
            bottom center point of the inductive meander section. This list should
            be the positions, in order, of the top left, top right, bot left,
            bot right.

        KID_Nos : list
            list of ints which are the numbers each KID should have drawn next to
            it. This list of numbers for each KID should be the same order as the
            rel_kid_positions, TL, TR, BL, BR.

        f0s : float, int
            The resonant frequencies of the resonators. Should be in the same unit
            that the IDC_and_CC_function function takes.

        IDC_and_CC_function : function
            function for getting the IDC and CC lengths from a given f0. The
            function should take a single frequency as an argument and should
            return an array-like (28 long) and a single float **in this order**,
            the array-like should contain all the lengths for each of the 28 arms
            for the IDC and the float should be the CC length.
            A simple example funtion:

            >>> def example_IDC_CC_func(f0):
                    '''
                    Example function for IDC and CC
                    f0 : float
                        Resonant frequency in Hz.
                    '''
                    if f0 < 3.1e9:
                        IDC_lengths = np.ones(28)*1900.0
                        CC_length = 600.0
                    else:
                        IDC_lengths = np.ones(28)*1500.0
                        CC_length = 300.0

                    # return IDC array, then CC seperately
                    return IDC_lengths, CC_length

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        KwArgs
        ------
        IDC_and_frame_materials=None
            The material to make each of the IDC and frame structures out of.
            By Default this is None which will make all the KIDs out of the
            IDC_Nb material. When this is passed in it should be a length 4 list
            of strings that take any of the values "IDC_Nb", "Nb", or "Al",
            which will make the frame and IDC out of IDC_Nb, Nb_Antenna, or
            Aluminium layers respectively. The order of the values passed in
            will be attributed to each KID and should be the same order as the
            rel_kid_positions, TL, TR, BL, BR.

        meander_materials=None
            The material to make each of the inductive meander structure out of.
            By Default this is None which will make all the meanders out of the
            Aluminium material. When this is passed in it should be a length 4
            list of strings that take any of the values "Al", "IDC_Nb", or "Nb",
            which will make the inductive meander out of Aluminium, IDC_Nb, or
            Nb_Antenna layers respectively. The order of the values passed in
            will be attributed to each meander and should be the same order as
            the rel_kid_positions, TL, TR, BL, BR.

        trim_lengths=None
            Whether not not to add trim boxes to the mask. When this is defined
            it should be a list of 4 ints of floats that define how long the
            trim arms should be on the mask. Trim boxes will be made to cover
            the trim arms to bring them down to the length specified.
            **More info can be found in the add_kid method.**

        add_grnd_cutout=True
            Whether or not to add a cutout in the Nb_Groundplane layer in the
            neccary place for the KID structure.

        add_SiN_dep_dielectric_around=True
            Whether or not to add an SiN depositon layer around in neccary
            place for the KID structure, this is a box the of the pixel pitch
            around the center of the antenna structure.

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
            Whether or not to add a backside check cover in the
            neccary place for the KID structure.

        """
        config = Main_config_file_dict["resonator"]

        # Making the KID number text
        text_size = config["text_size"]  # 90
        text_x_offset = config["text_x_offset"]  # 800
        text_y_offset = config["text_y_offset"]  # 900
        text_underline_height = config["text_underline_height"]  # 3

        rot_angles = [pi / 2, -pi / 2, pi / 2, -pi / 2]
        to_mirror = [True, False, False, True]

        if not (type(trim_lengths) == list and len(trim_lengths) == 4):  # checks if trim lengths is a list of 4 ints of floats
            trim_lengths = [None, None, None, None]
        elif not all(isinstance(x, (int, float)) for x in trim_lengths):
            trim_lengths = [None, None, None, None]

        if not (type(IDC_and_frame_materials) == list and len(IDC_and_frame_materials) == 4):  # checks if materials is a list of 4 strings
            IDC_and_frame_materials = ["IDC_Nb", "IDC_Nb", "IDC_Nb", "IDC_Nb"]
        elif not all(isinstance(x, (str)) for x in IDC_and_frame_materials):
            IDC_and_frame_materials = ["IDC_Nb", "IDC_Nb", "IDC_Nb", "IDC_Nb"]

        if not (type(meander_materials) == list and len(meander_materials) == 4):  # checks if materials is a list of 4 strings
            meander_materials = ["Al", "Al", "Al", "Al"]
        elif not all(isinstance(x, (str)) for x in meander_materials):
            meander_materials = ["Al", "Al", "Al", "Al"]

        for k in range(4):
            kid_x = x + rel_kid_positions[k][0]
            kid_y = y + rel_kid_positions[k][1]

            self.add_kid(
                kid_x,
                kid_y,
                rot_angles[k],
                f0s[k],
                IDC_and_CC_function,
                Main_config_file_dict,
                mirror=to_mirror[k],
                IDC_and_frame_material=IDC_and_frame_materials[k],
                meander_material=meander_materials[k],
                trim_length=trim_lengths[k],
                add_grnd_cutout=add_grnd_cutout,
                add_SiN_dep_dielectric_cutout=add_SiN_dep_dielectric_cutout,
                add_SiO_cutout=add_SiO_cutout,
                add_SiN_membrane_cutout=add_SiN_membrane_cutout,
                add_backside_check=add_backside_check,
            )

            # Adding the KID number Text
            undln_len = (text_size / 9) * (5 * len(str(KID_Nos[k])) + 3 * (len(str(KID_Nos[k])) - 1))
            x_sign = 1 if k in [0, 2] else -1
            y_sign = 1 if k in [0, 1] else -1
            KID_text_label = gdspy.Text(
                str(KID_Nos[k]), text_size, (kid_x + x_sign * text_x_offset, kid_y + y_sign * text_y_offset), **self.Aluminium
            )
            underline = gdspy.Rectangle(
                [kid_x + x_sign * text_x_offset, kid_y + y_sign * text_y_offset],
                [kid_x + x_sign * text_x_offset + undln_len, kid_y + y_sign * text_y_offset + text_underline_height],
                **self.Aluminium,
            )
            self.Main.add(KID_text_label)
            self.ground_plane_cutouts.add([KID_text_label])
            self.silicon_nitride_cutouts.add([KID_text_label])
            self.Main.add(underline)
            self.ground_plane_cutouts.add([underline])
            self.silicon_nitride_cutouts.add([underline])

        # Making the SiN dep dielectric to go around resonators
        if add_SiN_dep_dielectric_around:
            dielectric_around_width = Main_config_file_dict["general"]["horizontal_pitch"]
            dielectric_around_height = Main_config_file_dict["general"]["vertical_pitch"]

            dielectric_around_box = gdspy.Rectangle(
                [x - dielectric_around_width / 2, y - dielectric_around_height / 2],
                [x + dielectric_around_width / 2, y + dielectric_around_height / 2],
                **self.SiN_dep,
            )
            self.silicon_nitride_positives.add(dielectric_around_box)

        return

    def add_kid(
        self,
        x,
        y,
        rot_angle,
        f0,
        IDC_and_CC_function,
        Main_config_file_dict,
        mirror=False,
        IDC_and_frame_material="IDC_Nb",
        meander_material="Al",
        trim_length=None,
        add_grnd_cutout=True,
        add_SiN_dep_dielectric_cutout=True,
        add_SiO_cutout=True,
        add_SiN_membrane_cutout=True,
        add_backside_check=True,
    ):
        """
        Adds the KID geometry to the Main cell athe the x,y cooardinate given.
        The KID is placed where the base middle of the inductive meander is at
        this x,y. The KID geometry is defined by the dimensions within the
        Main_config_file_dict. By default it will, but optionally can choose not
        to, add all the neccessay cutouts for the structure.

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
            function should take a single frequency as an argument and should
            return an array-like (28 long) and a single float **in this order**,
            the array-like should contain all the lengths for each of the 28 arms
            for the IDC and the float should be the CC length.
            A simple example funtion:

            >>> def example_IDC_CC_func(f0):
                    '''
                    Example function for IDC and CC
                    f0 : float
                        Resonant frequency in Hz.
                    '''
                    if f0 < 3.1e9:
                        IDC_lengths = np.ones(28)*1900.0
                        CC_length = 600.0
                    else:
                        IDC_lengths = np.ones(28)*1500.0
                        CC_length = 300.0

                    # return IDC array, then CC seperately
                    return IDC_lengths, CC_length

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        KwArgs
        ------
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
            Whether or not to add a backside check cover in the
            neccary place for the KID structure.

        """

        IDC_and_frame_material_lookup = {"IDC_Nb": self.IDC_Nb, "Nb": self.Nb_Antenna, "Al": self.Aluminium}
        material_idc_and_frame = IDC_and_frame_material_lookup[IDC_and_frame_material]

        meander_material_lookup = {"Al": self.Aluminium, "IDC_Nb": self.IDC_Nb, "Nb": self.Nb_Antenna}
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

        meander_path_points = [
            [
                (meander_bot_width / 2) - meander_right_width_1 + meander_right_width_2,
                meander_lw + meander_right_height_1 + meander_right_height_2 + meander_right_height_3,
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
                meander_lw + meander_left_height_1 + meander_left_height_2 + meander_left_height_3,
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

        meander_cover_box_right_points = [
            [meander_path_points[0][0] - (frame_meander_cover_box_width / 2), meander_path_points[0][1] + frame_bot_lw],
            [
                meander_path_points[0][0] - (frame_meander_cover_box_width / 2),
                meander_path_points[0][1] + frame_bot_lw - frame_meander_cover_box_height,
            ],
            [
                meander_path_points[0][0] + (frame_meander_cover_box_width / 2),
                meander_path_points[0][1] + frame_bot_lw - frame_meander_cover_box_height,
            ],
            [meander_path_points[0][0] + (frame_meander_cover_box_width / 2), meander_path_points[0][1] + frame_bot_lw],
        ]

        meander_cover_box_left_points = [
            [meander_path_points[-1][0] - (frame_meander_cover_box_width / 2), meander_path_points[-1][1] + frame_bot_lw],
            [
                meander_path_points[-1][0] - (frame_meander_cover_box_width / 2),
                meander_path_points[-1][1] + frame_bot_lw - frame_meander_cover_box_height,
            ],
            [
                meander_path_points[-1][0] + (frame_meander_cover_box_width / 2),
                meander_path_points[-1][1] + frame_bot_lw - frame_meander_cover_box_height,
            ],
            [meander_path_points[-1][0] + (frame_meander_cover_box_width / 2), meander_path_points[-1][1] + frame_bot_lw],
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
        cutout_height = (
            coupler_frame_start_y + coupler_frame_left_height + coupler_gap + coupler_lw + cutout_top_offset
        ) - cutout_start_height

        grnd_plane_cutout_width = cutout_width + 30
        grnd_plane_cutout_height = cutout_height + 30
        grnd_plane_cutout_start_height = cutout_start_height - 15

        grnd_plane_meander_cutout_poly_points = [
            [-grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height],
            [-grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height + grnd_plane_cutout_height],
            [grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height + grnd_plane_cutout_height],
            [grnd_plane_cutout_width / 2, grnd_plane_cutout_start_height],
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
            new_meander_path_points = self.mirror_points_around_yaxis(meander_path_points)
            new_meander_path_points = self.rotate_and_move_points_list(new_meander_path_points, rot_angle, x, y)
        else:
            new_meander_path_points = self.rotate_and_move_points_list(meander_path_points, rot_angle, x, y)
        meander_path = gdspy.FlexPath(
            new_meander_path_points, meander_lw, corners="circular bend", bend_radius=meander_corner_bend_radius, **material_meander
        )
        meander_path_polygons = self.get_polys_from_flexpath(meander_path)
        for i in range(len(meander_path_polygons)):
            self.Main.add(gdspy.Polygon(meander_path_polygons[i], **material_meander))

        # Adding the meander ant overlap box
        new_ant_pad_box_points = self.rotate_and_move_points_list(ant_pad_box_points, rot_angle, x, y)
        ant_pad_box = gdspy.Polygon(new_ant_pad_box_points, **material_meander)
        self.Main.add(ant_pad_box)

        # Adding the meander frame overlap boxes
        if mirror:
            new_meander_cover_box_right_points = self.mirror_points_around_yaxis(meander_cover_box_right_points)
            new_meander_cover_box_right_points = self.rotate_and_move_points_list(new_meander_cover_box_right_points, rot_angle, x, y)
        else:
            new_meander_cover_box_right_points = self.rotate_and_move_points_list(meander_cover_box_right_points, rot_angle, x, y)

        meander_cover_box_right = gdspy.Polygon(new_meander_cover_box_right_points, **material_idc_and_frame)
        self.Main.add(meander_cover_box_right)

        if mirror:
            new_meander_cover_box_left_points = self.mirror_points_around_yaxis(meander_cover_box_left_points)
            new_meander_cover_box_left_points = self.rotate_and_move_points_list(new_meander_cover_box_left_points, rot_angle, x, y)
        else:
            new_meander_cover_box_left_points = self.rotate_and_move_points_list(meander_cover_box_left_points, rot_angle, x, y)

        meander_cover_box_left = gdspy.Polygon(new_meander_cover_box_left_points, **material_idc_and_frame)
        self.Main.add(meander_cover_box_left)

        # Adding the frame left and frame right
        if mirror:
            new_frame_left_points = self.mirror_points_around_yaxis(frame_left_points)
            new_frame_left_points = self.rotate_and_move_points_list(new_frame_left_points, rot_angle, x, y)
        else:
            new_frame_left_points = self.rotate_and_move_points_list(frame_left_points, rot_angle, x, y)

        frame_left_poly = gdspy.Polygon(new_frame_left_points, **material_idc_and_frame)
        self.Main.add(frame_left_poly)

        if mirror:
            new_frame_right_points = self.mirror_points_around_yaxis(frame_right_points)
            new_frame_right_points = self.rotate_and_move_points_list(new_frame_right_points, rot_angle, x, y)
        else:
            new_frame_right_points = self.rotate_and_move_points_list(frame_right_points, rot_angle, x, y)

        frame_right_poly = gdspy.Polygon(new_frame_right_points, **material_idc_and_frame)
        self.Main.add(frame_right_poly)

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
            new_coupler_frame_points = self.mirror_points_around_yaxis(coupler_frame_points)
            new_coupler_frame_points = self.rotate_and_move_points_list(new_coupler_frame_points, rot_angle, x, y)
        else:
            new_coupler_frame_points = self.rotate_and_move_points_list(coupler_frame_points, rot_angle, x, y)

        coupler_frame_poly = gdspy.Polygon(new_coupler_frame_points, **material_idc_and_frame)
        self.Main.add(coupler_frame_poly)

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
            new_coupler_arm_points = self.mirror_points_around_yaxis(coupler_arm_points)
            new_coupler_arm_points = self.rotate_and_move_points_list(new_coupler_arm_points, rot_angle, x, y)
        else:
            new_coupler_arm_points = self.rotate_and_move_points_list(coupler_arm_points, rot_angle, x, y)

        coupler_arm_poly = gdspy.Polygon(new_coupler_arm_points, **self.Nb_Antenna)
        self.Main.add(coupler_arm_poly)

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
            self.Main.add(right_arm)

            left_arm = gdspy.Rectangle(
                [arm_start_x_left_side, arm_start_y_left_side - (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
                [arm_start_x_left_side + IDCLs[-(i + 2)], arm_start_y_left_side + (IDC_arm_lw / 2) + (i * (IDC_arm_gap + IDC_arm_lw))],
                **material_idc_and_frame,
            )
            left_arm.translate(x, y)
            if mirror:
                left_arm.mirror([x, y], [x, y + 10])
            left_arm.rotate(rot_angle, center=(x, y))
            self.Main.add(left_arm)

        # Adding the Trim arms
        right_trim_arm = gdspy.Rectangle(
            [arm_start_x_right_side, trim_arm_start_y_right_side],
            [arm_start_x_right_side - trim_arm_length_right_side, trim_arm_start_y_right_side + trim_arm_lw],
            **material_idc_and_frame,
        )
        right_trim_arm.translate(x, y)
        if mirror:
            right_trim_arm.mirror([x, y], [x, y + 10])
        right_trim_arm.rotate(rot_angle, center=(x, y))
        self.Main.add(right_trim_arm)

        left_trim_arm = gdspy.Rectangle(
            [arm_start_x_left_side, trim_arm_start_y_left_side],
            [arm_start_x_left_side + trim_arm_length_left_side, trim_arm_start_y_left_side + trim_arm_lw],
            **material_idc_and_frame,
        )
        left_trim_arm.translate(x, y)
        if mirror:
            left_trim_arm.mirror([x, y], [x, y + 10])
        left_trim_arm.rotate(rot_angle, center=(x, y))
        self.Main.add(left_trim_arm)

        # Adding the Trim boxes for the trim arms at Trim lengths specified if non-zero
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
                    **self.TrimLayer,
                )
                trim_box_for_right_trim_arm.translate(x, y)
                if mirror:
                    trim_box_for_right_trim_arm.mirror([x, y], [x, y + 10])
                trim_box_for_right_trim_arm.rotate(rot_angle, center=(x, y))
                self.Main.add(trim_box_for_right_trim_arm)

                trim_box_for_left_trim_arm = gdspy.Rectangle(
                    [
                        arm_start_x_left_side + trim_arm_length_left_side + trim_box_length_overhang,
                        trim_arm_start_y_left_side + (trim_arm_lw / 2) - (trim_box_width / 2),
                    ],
                    [
                        arm_start_x_right_side - (inner_width - trim_length),
                        trim_arm_start_y_left_side + (trim_arm_lw / 2) + (trim_box_width / 2),
                    ],
                    **self.TrimLayer,
                )
                trim_box_for_left_trim_arm.translate(x, y)
                if mirror:
                    trim_box_for_left_trim_arm.mirror([x, y], [x, y + 10])
                trim_box_for_left_trim_arm.rotate(rot_angle, center=(x, y))
                self.Main.add(trim_box_for_left_trim_arm)

        # Adding the cutout to the groundplane
        if add_grnd_cutout:
            if mirror:
                new_grnd_plane_meander_cutout_poly_points = self.mirror_points_around_yaxis(grnd_plane_meander_cutout_poly_points)
                new_grnd_plane_meander_cutout_poly_points = self.rotate_and_move_points_list(
                    new_grnd_plane_meander_cutout_poly_points, rot_angle, x, y
                )
            else:
                new_grnd_plane_meander_cutout_poly_points = self.rotate_and_move_points_list(
                    grnd_plane_meander_cutout_poly_points, rot_angle, x, y
                )

            grnd_plane_meander_cutout_poly = gdspy.Polygon(new_grnd_plane_meander_cutout_poly_points, **self.Nb_Groundplane)
            self.ground_plane_cutouts.add(grnd_plane_meander_cutout_poly)

        # Adding the cutout to the Silicon DiOxide membrane
        if add_SiO_cutout:
            if mirror:
                new_SiO_cutout_poly_points = self.mirror_points_around_yaxis(SiO_cutout_poly_points)
                new_SiO_cutout_poly_points = self.rotate_and_move_points_list(new_SiO_cutout_poly_points, rot_angle, x, y)
            else:
                new_SiO_cutout_poly_points = self.rotate_and_move_points_list(SiO_cutout_poly_points, rot_angle, x, y)

            SiO_cutout_poly = gdspy.Polygon(new_SiO_cutout_poly_points, **self.Nb_Groundplane)
            self.silicon_oxide_cutouts.add(SiO_cutout_poly)

        # Adding the cutout to the Silicon Nitride membrane
        if add_SiN_membrane_cutout:
            if mirror:
                new_SiN_membrane_cutout_poly_points = self.mirror_points_around_yaxis(SiN_membrane_cutout_poly_points)
                new_SiN_membrane_cutout_poly_points = self.rotate_and_move_points_list(new_SiN_membrane_cutout_poly_points, rot_angle, x, y)
            else:
                new_SiN_membrane_cutout_poly_points = self.rotate_and_move_points_list(SiN_membrane_cutout_poly_points, rot_angle, x, y)

            SiN_membrane_cutout_poly = gdspy.Polygon(new_SiN_membrane_cutout_poly_points, **self.Nb_Groundplane)
            self.silicon_nitride_membrane_cutouts.add(SiN_membrane_cutout_poly)

        # Adding the cutout to the SiN Dep layer
        if add_SiN_dep_dielectric_cutout:
            if mirror:
                new_SiN_dep_cutout_poly_points = self.mirror_points_around_yaxis(SiN_dep_cutout_poly_points)
                new_SiN_dep_cutout_poly_points = self.rotate_and_move_points_list(new_SiN_dep_cutout_poly_points, rot_angle, x, y)
            else:
                new_SiN_dep_cutout_poly_points = self.rotate_and_move_points_list(SiN_dep_cutout_poly_points, rot_angle, x, y)

            SiN_dep_cutout_poly = gdspy.Polygon(new_SiN_dep_cutout_poly_points, **self.SiN_dep)
            self.silicon_nitride_cutouts.add(SiN_dep_cutout_poly)

        # Adding the backside check covers
        if add_backside_check:
            if mirror:
                new_backside_check_cover_poly_points = self.mirror_points_around_yaxis(backside_check_cover_poly_points)
                new_backside_check_cover_poly_points = self.rotate_and_move_points_list(
                    new_backside_check_cover_poly_points, rot_angle, x, y
                )
            else:
                new_backside_check_cover_poly_points = self.rotate_and_move_points_list(backside_check_cover_poly_points, rot_angle, x, y)

            backside_check_cover_poly = gdspy.Polygon(new_backside_check_cover_poly_points, **self.Backside_Check)
            self.Main.add(backside_check_cover_poly)

        return

    def add_Lo_pass_filters(self, x, y, inner_ring_line_width, inner_ring_radius, init_angle, direction, Main_config_file_dict):
        """
        Adds the low pass filter arms around the inner ring of the filter bank.

        Parameters
        ----------
        x, y : float, int
            x, y coordinate of the center of the antenna to be places around.

        inner_ring_line_width : float, int
            Line width of the inner ring that the filters conect to.

        inner_ring_radius : float, int
            Radius of the inner ring that the filters conect to.

        init_angle : float, int
            Initial angle (*in radians*) to start at when making the filter arms.

        direction : str, **only "clockwise" or "anti-clockwise".**
            The direction around the inner ring in which to make the filters point.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "Lo_pass_filters".

        """

        config = Main_config_file_dict["Lo_pass_filters"]
        # LoPass arm details
        Lo_pass_arm1_linewidth = config["Lo_pass_arm1_linewidth"]  # 9.5
        Lo_pass_arm1_height = config["Lo_pass_arm1_height"]  # 20
        Lo_pass_arm1_length = config["Lo_pass_arm1_length"]  # 245
        Lo_pass_arm1_offset = config["Lo_pass_arm1_offset"]  # 152
        Lo_pass_arm1_offset_arc_angle = Lo_pass_arm1_offset / inner_ring_radius

        Lo_pass_arm2_linewidth = config["Lo_pass_arm2_linewidth"]  # 20.5
        Lo_pass_arm2_height = config["Lo_pass_arm2_height"]  # 41.5
        Lo_pass_arm2_length = config["Lo_pass_arm2_length"]  # 218
        Lo_pass_arm2_offset = config["Lo_pass_arm2_offset"]  # 278
        Lo_pass_arm2_offset_arc_angle = Lo_pass_arm2_offset / inner_ring_radius

        Lo_pass_arm3_linewidth = config["Lo_pass_arm3_linewidth"]  # 20
        Lo_pass_arm3_height = config["Lo_pass_arm3_height"]  # 40
        Lo_pass_arm3_length = config["Lo_pass_arm3_length"]  # 219
        Lo_pass_arm3_offset = config["Lo_pass_arm3_offset"]  # 278
        Lo_pass_arm3_offset_arc_angle = Lo_pass_arm3_offset / inner_ring_radius

        Lo_pass_arm4_linewidth = config["Lo_pass_arm4_linewidth"]  # 20.5
        Lo_pass_arm4_height = config["Lo_pass_arm4_height"]  # 41.5
        Lo_pass_arm4_length = config["Lo_pass_arm4_length"]  # 218
        Lo_pass_arm4_offset = config["Lo_pass_arm4_offset"]  # 278
        Lo_pass_arm4_offset_arc_angle = Lo_pass_arm4_offset / inner_ring_radius

        Lo_pass_arm5_linewidth = config["Lo_pass_arm5_linewidth"]  # 9.5
        Lo_pass_arm5_height = config["Lo_pass_arm5_height"]  # 120.5
        Lo_pass_arm5_length = config["Lo_pass_arm5_length"]  # 145
        Lo_pass_arm5_offset = config["Lo_pass_arm5_offset"]  # 278
        Lo_pass_arm5_offset_arc_angle = Lo_pass_arm5_offset / inner_ring_radius

        via_overlap_length = config["via_overlap_length"]  # 5
        via_extend_length = config["via_extend_length"]  # 5

        if direction == "clockwise":
            sign = -1
        elif direction == "anti-clockwise":
            sign = +1

        arm_linewidths = [
            Lo_pass_arm1_linewidth,
            Lo_pass_arm2_linewidth,
            Lo_pass_arm3_linewidth,
            Lo_pass_arm4_linewidth,
            Lo_pass_arm5_linewidth,
        ]
        arm_heights = [Lo_pass_arm1_height, Lo_pass_arm2_height, Lo_pass_arm3_height, Lo_pass_arm4_height, Lo_pass_arm5_height]
        arm_lengths = [Lo_pass_arm1_length, Lo_pass_arm2_length, Lo_pass_arm3_length, Lo_pass_arm4_length, Lo_pass_arm5_length]
        offset_arc_angles = [
            Lo_pass_arm1_offset_arc_angle,
            Lo_pass_arm2_offset_arc_angle,
            Lo_pass_arm3_offset_arc_angle,
            Lo_pass_arm4_offset_arc_angle,
            Lo_pass_arm5_offset_arc_angle,
        ]
        offset_arc_cumsum_angles = np.cumsum(offset_arc_angles)

        for i in range(5):
            angle = init_angle + sign * offset_arc_cumsum_angles[i]
            x_pos = x + (inner_ring_radius * cos(angle))
            y_pos = y + (inner_ring_radius * sin(angle))

            arm_geom_points = [
                [(sign) * arm_linewidths[i] / 2, 0],
                [(sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
                [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2)],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i],
                    arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
                ],
                [(-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i], arm_heights[i] + (inner_ring_line_width / 2)],
                [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
                [(-sign) * arm_linewidths[i] / 2, 0],
            ]

            arm_geom = self.rotate_and_move_points_list(arm_geom_points, angle - pi / 2, x_pos, y_pos)

            arm = gdspy.Polygon(arm_geom, **self.Nb_Antenna)
            self.Main.add(arm)

            via_points = [
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                    arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
                ],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                    arm_heights[i] + (inner_ring_line_width / 2),
                ],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                    arm_heights[i] + (inner_ring_line_width / 2),
                ],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                    arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
                ],
            ]

            via_geom = self.rotate_and_move_points_list(via_points, angle - pi / 2, x_pos, y_pos)

            via_box = gdspy.Polygon(via_geom, **self.SiN_dep)
            self.silicon_nitride_cutouts.add(via_box)

        return

    def add_Hi_pass_filters(self, x, y, inner_ring_line_width, inner_ring_radius, init_angle, direction, Main_config_file_dict):
        """
        Adds the high pass filter arms around the inner ring of the filter bank.

        Parameters
        ----------
        x, y : float, int
            x, y coordinate of the center of the antenna to be places around.

        inner_ring_line_width : float, int
            Line width of the inner ring that the filters conect to.

        inner_ring_radius : float, int
            Radius of the inner ring that the filters conect to.

        init_angle : float, int
            Initial angle (*in radians*) to start at when making the filter arms.

        direction : str, **only "clockwise" or "anti-clockwise".**
            The direction around the inner ring in which to make the filters point.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "Lo_pass_filters".

        """

        config = Main_config_file_dict["Hi_pass_filters"]
        # HiPass arm details
        Hi_pass_arm1_linewidth = config["Hi_pass_arm1_linewidth"]  # 22
        Hi_pass_arm1_height = config["Hi_pass_arm1_height"]  # 23
        Hi_pass_arm1_length = config["Hi_pass_arm1_length"]  # 120.5
        Hi_pass_arm1_offset = config["Hi_pass_arm1_offset"]  # 185
        Hi_pass_arm1_offset_arc_angle = Hi_pass_arm1_offset / inner_ring_radius

        Hi_pass_arm2_linewidth = config["Hi_pass_arm2_linewidth"]  # 41.5
        Hi_pass_arm2_height = config["Hi_pass_arm2_height"]  # 34.5
        Hi_pass_arm2_length = config["Hi_pass_arm2_length"]  # 97
        Hi_pass_arm2_offset = config["Hi_pass_arm2_offset"]  # 160
        Hi_pass_arm2_offset_arc_angle = Hi_pass_arm2_offset / inner_ring_radius

        Hi_pass_arm3_linewidth = config["Hi_pass_arm3_linewidth"]  # 40.5
        Hi_pass_arm3_height = config["Hi_pass_arm3_height"]  # 34.5
        Hi_pass_arm3_length = config["Hi_pass_arm3_length"]  # 99.5
        Hi_pass_arm3_offset = config["Hi_pass_arm3_offset"]  # 160
        Hi_pass_arm3_offset_arc_angle = Hi_pass_arm3_offset / inner_ring_radius

        Hi_pass_arm4_linewidth = config["Hi_pass_arm4_linewidth"]  # 41.5
        Hi_pass_arm4_height = config["Hi_pass_arm4_height"]  # 34.5
        Hi_pass_arm4_length = config["Hi_pass_arm4_length"]  # 97
        Hi_pass_arm4_offset = config["Hi_pass_arm4_offset"]  # 160
        Hi_pass_arm4_offset_arc_angle = Hi_pass_arm4_offset / inner_ring_radius

        Hi_pass_arm5_linewidth = config["Hi_pass_arm5_linewidth"]  # 22
        Hi_pass_arm5_height = config["Hi_pass_arm5_height"]  # 23
        Hi_pass_arm5_length = config["Hi_pass_arm5_length"]  # 120.5
        Hi_pass_arm5_offset = config["Hi_pass_arm5_offset"]  # 160
        Hi_pass_arm5_offset_arc_angle = Hi_pass_arm5_offset / inner_ring_radius

        via_overlap_length = config["via_overlap_length"]  # 5
        via_extend_length = config["via_extend_length"]  # 5

        if direction == "clockwise":
            sign = -1
        elif direction == "anti-clockwise":
            sign = +1

        arm_linewidths = [
            Hi_pass_arm1_linewidth,
            Hi_pass_arm2_linewidth,
            Hi_pass_arm3_linewidth,
            Hi_pass_arm4_linewidth,
            Hi_pass_arm5_linewidth,
        ]
        arm_heights = [Hi_pass_arm1_height, Hi_pass_arm2_height, Hi_pass_arm3_height, Hi_pass_arm4_height, Hi_pass_arm5_height]
        arm_lengths = [Hi_pass_arm1_length, Hi_pass_arm2_length, Hi_pass_arm3_length, Hi_pass_arm4_length, Hi_pass_arm5_length]
        offset_arc_angles = [
            Hi_pass_arm1_offset_arc_angle,
            Hi_pass_arm2_offset_arc_angle,
            Hi_pass_arm3_offset_arc_angle,
            Hi_pass_arm4_offset_arc_angle,
            Hi_pass_arm5_offset_arc_angle,
        ]
        offset_arc_cumsum_angles = np.cumsum(offset_arc_angles)

        for i in range(5):
            angle = init_angle + sign * offset_arc_cumsum_angles[i]
            x_pos = x + (inner_ring_radius * cos(angle))
            y_pos = y + (inner_ring_radius * sin(angle))

            arm_geom_points = [
                [(sign) * arm_linewidths[i] / 2, 0],
                [(sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
                [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2)],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i],
                    arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
                ],
                [(-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i], arm_heights[i] + (inner_ring_line_width / 2)],
                [(-sign) * arm_linewidths[i] / 2, arm_heights[i] + (inner_ring_line_width / 2)],
                [(-sign) * arm_linewidths[i] / 2, 0],
            ]

            arm_geom = self.rotate_and_move_points_list(arm_geom_points, angle - pi / 2, x_pos, y_pos)

            arm = gdspy.Polygon(arm_geom, **self.Nb_Antenna)
            self.Main.add(arm)

            via_points = [
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                    arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
                ],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (sign) * via_overlap_length,
                    arm_heights[i] + (inner_ring_line_width / 2),
                ],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                    arm_heights[i] + (inner_ring_line_width / 2),
                ],
                [
                    (-sign) * arm_linewidths[i] / 2 + (-sign) * arm_lengths[i] + (-sign) * via_extend_length,
                    arm_heights[i] + arm_linewidths[i] + (inner_ring_line_width / 2),
                ],
            ]

            via_geom = self.rotate_and_move_points_list(via_points, angle - pi / 2, x_pos, y_pos)

            via_box = gdspy.Polygon(via_geom, **self.SiN_dep)
            self.silicon_nitride_cutouts.add(via_box)

        return

    def add_combiner_section_and_get_conect_point(
        self, x, y, rot, outer_ring_conection_gap, outer_ring_linewidth, Main_config_file_dict, combiner_type="90GHZ", mirror_combiner=False
    ):
        """
        Adds the phase combiner to the outer ring of the filter bank. This is by
        default the 90GHz or optionally the 150GHz. The specific geometries
        dimensions is determined by the parameters in the config.
        This function will return the coordinate of the conection point

        Parameters
        ----------
        x, y : float, int
            x, y coordinate for the base of the combiner structure, this is the
            very bottom middle where it connects to the outer ring of the filter
            bank. i.e. the middle of the gap in the outer ring.

        rot : float, int
            The angle (**in radians**) of the rotation of the whole combiner
            geometry.

        outer_ring_conection_gap : float, int
            The gap distance in the outer ring so the combiner can connect to
            either side.

        outer_ring_linewidth : float, int
            The line width of the outer ring of the filter bank this combiner
            attaches to.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "combiner_section_90GHZ" OR "combiner_section_150GHZ"
            depending upon the "combiner_type" kwarg.

        KwArgs
        ------
        combiner_type = "90GHZ"
            This defines what type of phase combiner to add. This can take string
            values of either "90GHZ" or "150GHZ". The Main_config_file_dict should
            contain this relevant key to access the config parameters.

        mirror_combiner = False,
            This defines if the combiner should be mirrored around the line
            perpedicular to the outer ring of the filter bank. The default is
            False or can be True.

        Returns
        -------

        points_to_conect_kids_to : list
            list of [x,y] coordinate for the conection point which connects to a
            KID.

        """

        if combiner_type == "90GHZ":
            config = Main_config_file_dict["combiner_section_90GHZ"]
        elif combiner_type == "150GHZ":
            config = Main_config_file_dict["combiner_section_150GHZ"]
        else:
            raise ValueError("'combiner_type' argument should be either '90GHZ' or '150GHZ'")

        # Making First section
        first_linewidth = config["first_linewidth"]  # 5
        distance_from_outer_ring = config["distance_from_outer_ring"]  # 50
        distance_from_outer_ring += (first_linewidth / 2) + (outer_ring_linewidth / 2)

        first_section_height = config["first_section_height"]  # 560
        first_section_width = config["first_section_width"]  # 155
        first_section_back_height = config["first_section_back_height"]  # 105

        first_section_top_points = [
            (-outer_ring_linewidth / 2, outer_ring_conection_gap / 2),
            (distance_from_outer_ring, outer_ring_conection_gap / 2),
            (distance_from_outer_ring, first_section_height / 2),
            (distance_from_outer_ring + first_section_width, first_section_height / 2),
            (distance_from_outer_ring + first_section_width, first_section_height / 2 - first_section_back_height),
        ]

        first_section_bot_points = [
            (-outer_ring_linewidth / 2, -outer_ring_conection_gap / 2),
            (distance_from_outer_ring, -outer_ring_conection_gap / 2),
            (distance_from_outer_ring, -first_section_height / 2),
            (distance_from_outer_ring + first_section_width, -first_section_height / 2),
            (distance_from_outer_ring + first_section_width, -first_section_height / 2 + first_section_back_height),
        ]

        if mirror_combiner:
            mirrored_first_section_top_points = self.mirror_points_around_xaxis(first_section_top_points)
            mirrored_first_section_bot_points = self.mirror_points_around_xaxis(first_section_bot_points)

            first_section_top = gdspy.FlexPath(
                mirrored_first_section_top_points, first_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )
            first_section_bot = gdspy.FlexPath(
                mirrored_first_section_bot_points, first_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )
        else:
            first_section_top = gdspy.FlexPath(first_section_top_points, first_linewidth, corners=self.create_miter_join, **self.Nb_Antenna)
            first_section_bot = gdspy.FlexPath(first_section_bot_points, first_linewidth, corners=self.create_miter_join, **self.Nb_Antenna)

        first_section_top.rotate(rot, (0, 0))
        first_section_top.translate(x, y)
        first_section_bot.rotate(rot, (0, 0))
        first_section_bot.translate(x, y)
        self.Main.add(first_section_top)
        self.Main.add(first_section_bot)

        # Making second section
        second_section_linewidth = config["second_section_linewidth"]  # 4.5

        second_section_top_width1 = config["second_section_top_width1"]  # 25.25
        second_section_top_height1 = config["second_section_top_height1"]  # 16.5
        second_section_top_width2 = config["second_section_top_width2"]  # 20.5
        second_section_top_height2 = config["second_section_top_height2"]  # 16.5
        second_section_top_width3 = config["second_section_top_width3"]  # 29
        second_section_top_height3 = config["second_section_top_height3"]  # 16.5
        second_section_top_width4 = config["second_section_top_width4"]  # 20.5
        second_section_top_height4 = config["second_section_top_height4"]  # 11.5
        second_section_top_width5 = config["second_section_top_width5"]  # 26

        second_section_top_first_point = [first_section_top_points[-1][0], first_section_top_points[-1][1] + second_section_linewidth / 2]

        second_section_top_points = [
            second_section_top_first_point,
            [second_section_top_first_point[0] + second_section_top_width1, second_section_top_first_point[1]],
            [second_section_top_first_point[0] + second_section_top_width1, second_section_top_first_point[1] - second_section_top_height1],
            [
                second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2,
                second_section_top_first_point[1] - second_section_top_height1,
            ],
            [
                second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2,
                second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2,
            ],
            [
                second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2 + second_section_top_width3,
                second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2,
            ],
            [
                second_section_top_first_point[0] + second_section_top_width1 + second_section_top_width2 + second_section_top_width3,
                second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2 - second_section_top_height3,
            ],
            [
                second_section_top_first_point[0]
                + second_section_top_width1
                + second_section_top_width2
                + second_section_top_width3
                + second_section_top_width4,
                second_section_top_first_point[1] - second_section_top_height1 + second_section_top_height2 - second_section_top_height3,
            ],
            [
                second_section_top_first_point[0]
                + second_section_top_width1
                + second_section_top_width2
                + second_section_top_width3
                + second_section_top_width4,
                second_section_top_first_point[1]
                - second_section_top_height1
                + second_section_top_height2
                - second_section_top_height3
                + second_section_top_height4,
            ],
            [
                second_section_top_first_point[0]
                + second_section_top_width1
                + second_section_top_width2
                + second_section_top_width3
                + second_section_top_width4
                + second_section_top_width5,
                second_section_top_first_point[1]
                - second_section_top_height1
                + second_section_top_height2
                - second_section_top_height3
                + second_section_top_height4,
            ],
        ]

        if mirror_combiner:
            mirrored_second_section_top_points = self.mirror_points_around_xaxis(second_section_top_points)
            second_section_top = gdspy.FlexPath(
                mirrored_second_section_top_points, second_section_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )
        else:
            second_section_top = gdspy.FlexPath(
                second_section_top_points, second_section_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )

        second_section_top.rotate(rot, (0, 0))
        second_section_top.translate(x, y)
        self.Main.add(second_section_top)

        second_section_bot_width1 = config["second_section_bot_width1"]  # 24.25
        second_section_bot_height1 = config["second_section_bot_height1"]  # 15.5
        second_section_bot_width2 = config["second_section_bot_width2"]  # 20.5
        second_section_bot_height2 = config["second_section_bot_height2"]  # 15.5
        second_section_bot_width3 = config["second_section_bot_width3"]  # 24
        second_section_bot_height3 = config["second_section_bot_height3"]  # 15.5
        second_section_bot_width4 = config["second_section_bot_width4"]  # 20.5
        second_section_bot_height4 = config["second_section_bot_height4"]  # 15.5
        second_section_bot_width5 = config["second_section_bot_width5"]  # 32

        second_section_bot_first_point = [first_section_bot_points[-1][0], first_section_bot_points[-1][1] - second_section_linewidth / 2]

        second_section_bot_points = [
            second_section_bot_first_point,
            [second_section_bot_first_point[0] + second_section_bot_width1, second_section_bot_first_point[1]],
            [second_section_bot_first_point[0] + second_section_bot_width1, second_section_bot_first_point[1] + second_section_bot_height1],
            [
                second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2,
                second_section_bot_first_point[1] + second_section_bot_height1,
            ],
            [
                second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2,
                second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2,
            ],
            [
                second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2 + second_section_bot_width3,
                second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2,
            ],
            [
                second_section_bot_first_point[0] + second_section_bot_width1 + second_section_bot_width2 + second_section_bot_width3,
                second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2 + second_section_bot_height3,
            ],
            [
                second_section_bot_first_point[0]
                + second_section_bot_width1
                + second_section_bot_width2
                + second_section_bot_width3
                + second_section_bot_width4,
                second_section_bot_first_point[1] + second_section_bot_height1 - second_section_bot_height2 + second_section_bot_height3,
            ],
            [
                second_section_bot_first_point[0]
                + second_section_bot_width1
                + second_section_bot_width2
                + second_section_bot_width3
                + second_section_bot_width4,
                second_section_bot_first_point[1]
                + second_section_bot_height1
                - second_section_bot_height2
                + second_section_bot_height3
                - second_section_bot_height4,
            ],
            [
                second_section_bot_first_point[0]
                + second_section_bot_width1
                + second_section_bot_width2
                + second_section_bot_width3
                + second_section_bot_width4
                + second_section_bot_width5,
                second_section_bot_first_point[1]
                + second_section_bot_height1
                - second_section_bot_height2
                + second_section_bot_height3
                - second_section_bot_height4,
            ],
        ]

        if mirror_combiner:
            mirrored_second_section_bot_points = self.mirror_points_around_xaxis(second_section_bot_points)
            second_section_bot = gdspy.FlexPath(
                mirrored_second_section_bot_points, second_section_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )
        else:
            second_section_bot = gdspy.FlexPath(
                second_section_bot_points, second_section_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )

        second_section_bot.rotate(rot, (0, 0))
        second_section_bot.translate(x, y)
        self.Main.add(second_section_bot)

        start_of_second_section_vertical_linewidth = 3.5
        start_of_second_section_vertical_points = [
            [
                first_section_top_points[-1][0] - first_linewidth / 2 + start_of_second_section_vertical_linewidth / 2,
                first_section_top_points[-1][1],
            ],
            [
                first_section_bot_points[-1][0] - first_linewidth / 2 + start_of_second_section_vertical_linewidth / 2,
                first_section_bot_points[-1][1],
            ],
        ]
        if mirror_combiner:
            mirrored_start_of_second_section_vertical_points = self.mirror_points_around_xaxis(start_of_second_section_vertical_points)
            start_of_second_section_vertical = gdspy.FlexPath(
                mirrored_start_of_second_section_vertical_points, start_of_second_section_vertical_linewidth, **self.Nb_Antenna
            )
        else:
            start_of_second_section_vertical = gdspy.FlexPath(
                start_of_second_section_vertical_points, start_of_second_section_vertical_linewidth, **self.Nb_Antenna
            )

        start_of_second_section_vertical.rotate(rot, (0, 0))
        start_of_second_section_vertical.translate(x, y)
        self.Main.add(start_of_second_section_vertical)

        third_section_linewidth = config["third_section_linewidth"]  # 3

        end_of_second_section_vertical_linewidth = 6

        end_of_second_section_vertical_points = [
            [second_section_top_points[-1][0], second_section_top_points[-1][1] + second_section_linewidth / 2],
            [
                second_section_top_points[-1][0] + end_of_second_section_vertical_linewidth,
                second_section_top_points[-1][1] - second_section_linewidth / 2 + third_section_linewidth,
            ],
            [
                second_section_bot_points[-1][0] + end_of_second_section_vertical_linewidth,
                second_section_bot_points[-1][1] + second_section_linewidth / 2 - third_section_linewidth,
            ],
            [second_section_bot_points[-1][0], second_section_bot_points[-1][1] - second_section_linewidth / 2],
        ]
        if mirror_combiner:
            mirrored_end_of_second_section_vertical_points = self.mirror_points_around_xaxis(end_of_second_section_vertical_points)
            end_of_second_section_vertical = gdspy.Polygon(mirrored_end_of_second_section_vertical_points, **self.Nb_Antenna)
        else:
            end_of_second_section_vertical = gdspy.Polygon(end_of_second_section_vertical_points, **self.Nb_Antenna)

        end_of_second_section_vertical.rotate(rot, (0, 0))
        end_of_second_section_vertical.translate(x, y)
        self.Main.add(end_of_second_section_vertical)

        # Making third section
        third_section_width1 = config["third_section_width1"]  # 21.5
        third_section_height1 = config["third_section_height1"]  # 56
        third_section_width2 = config["third_section_width2"]  # 23.5
        third_section_height2 = config["third_section_height2"]  # 56
        third_section_width3 = config["third_section_width3"]  # 23.5

        third_section_start_points = [
            end_of_second_section_vertical_points[1][0],
            end_of_second_section_vertical_points[1][1] - third_section_linewidth / 2,
        ]
        third_section_end_points = [
            end_of_second_section_vertical_points[2][0],
            end_of_second_section_vertical_points[2][1] + third_section_linewidth / 2,
        ]

        third_section_points = [
            third_section_start_points,
            [third_section_start_points[0] + third_section_width1, third_section_start_points[1]],
            [third_section_start_points[0] + third_section_width1, third_section_start_points[1] - third_section_height1],
            [
                third_section_start_points[0] + third_section_width1 + third_section_width2,
                third_section_start_points[1] - third_section_height1,
            ],
            [
                third_section_start_points[0] + third_section_width1 + third_section_width2,
                third_section_start_points[1] - third_section_height1 + third_section_height2,
            ],
            [
                third_section_start_points[0] + third_section_width1 + third_section_width2 + third_section_width3,
                third_section_start_points[1] - third_section_height1 + third_section_height2,
            ],
            [
                third_section_end_points[0] + third_section_width1 + third_section_width2 + third_section_width3,
                third_section_end_points[1] + third_section_height1 - third_section_height2,
            ],
            [
                third_section_end_points[0] + third_section_width1 + third_section_width2,
                third_section_end_points[1] + third_section_height1 - third_section_height2,
            ],
            [
                third_section_end_points[0] + third_section_width1 + third_section_width2,
                third_section_end_points[1] + third_section_height1,
            ],
            [third_section_end_points[0] + third_section_width1, third_section_end_points[1] + third_section_height1],
            [third_section_end_points[0] + third_section_width1, third_section_end_points[1]],
            third_section_end_points,
        ]
        if mirror_combiner:
            mirrored_third_section_points = self.mirror_points_around_xaxis(third_section_points)
            third_section = gdspy.FlexPath(
                mirrored_third_section_points, third_section_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )
        else:
            third_section = gdspy.FlexPath(third_section_points, third_section_linewidth, corners=self.create_miter_join, **self.Nb_Antenna)

        third_section.rotate(rot, (0, 0))
        third_section.translate(x, y)
        self.Main.add(third_section)

        conection_to_kid_path_linewidth = config["conection_to_kid_path_linewidth"]  # 5
        conection_to_kid_path_start_piece_length = config["conection_to_kid_path_start_piece_length"]  # 10
        conection_to_kid_path_start_path_points = [
            [
                third_section_points[5][0] - third_section_linewidth / 2,
                third_section_points[5][1] + third_section_linewidth / 2 - conection_to_kid_path_linewidth / 2,
            ],
            [
                third_section_points[5][0] - third_section_linewidth / 2 + conection_to_kid_path_start_piece_length,
                third_section_points[5][1] + third_section_linewidth / 2 - conection_to_kid_path_linewidth / 2,
            ],
        ]

        conection_to_kid_path_start_path_points = [
            [
                third_section_points[5][0] + third_section_linewidth / 2 - conection_to_kid_path_linewidth / 2,
                third_section_points[5][1] - third_section_linewidth / 2,
            ],
            [
                third_section_points[5][0] + third_section_linewidth / 2 - conection_to_kid_path_linewidth / 2,
                third_section_points[5][1] - third_section_linewidth / 2 + conection_to_kid_path_start_piece_length,
            ],
        ]

        if mirror_combiner:
            mirrored_conection_to_kid_path_start_path_points = self.mirror_points_around_xaxis(conection_to_kid_path_start_path_points)
            conection_to_kid_path_start = gdspy.FlexPath(
                mirrored_conection_to_kid_path_start_path_points, conection_to_kid_path_linewidth, **self.Nb_Antenna
            )
        else:
            conection_to_kid_path_start = gdspy.FlexPath(
                conection_to_kid_path_start_path_points, conection_to_kid_path_linewidth, **self.Nb_Antenna
            )

        conection_to_kid_path_start.rotate(rot, (0, 0))
        conection_to_kid_path_start.translate(x, y)
        self.Main.add(conection_to_kid_path_start)

        if mirror_combiner:
            points_to_conect_kids_to = self.rotate_and_move_single_point(mirrored_conection_to_kid_path_start_path_points[-1], rot, x, y)
        else:
            points_to_conect_kids_to = self.rotate_and_move_single_point(conection_to_kid_path_start_path_points[-1], rot, x, y)

        # Making the meander section and connect to frame

        meander_height = config["meander_height"]  # 525
        meander_linewidth = config["meander_linewidth"]  # 5

        meander_initial_gap = config["meander_initial_gap"]  # 10

        meander_final_gap = config["meander_final_gap"]
        meander_gap_spacing = config["meander_gap_spacing"]

        meander_conect_bot_width = config["meander_conect_bot_width"]  # 30
        meander_conect_right_height = config["meander_conect_right_height"]  # 265
        meander_conect_final_width = config["meander_conect_final_width"]  # 7.5

        meander_conect_linewidth = config["meander_conect_linewidth"]  # 5

        no_of_full_straights = int(config["no_of_full_straights"])  # 10
        direction = 1
        meander_points = []
        meander_offset_x = distance_from_outer_ring + meander_initial_gap + (meander_linewidth / 2) + (first_linewidth / 2)

        for i in range(no_of_full_straights - 1):
            direction *= -1
            meander_points.append(
                [meander_offset_x + i * meander_gap_spacing, direction * meander_height / 2 + direction * meander_linewidth / 2]
            )
            meander_points.append(
                [meander_offset_x + i * meander_gap_spacing, -direction * meander_height / 2 - direction * meander_linewidth / 2]
            )

        meander_points.append([meander_points[-1][0] + meander_final_gap, meander_points[-1][1]])
        meander_points.append([meander_points[-1][0], meander_points[-1][1] - meander_height - (meander_linewidth * 1.5)])

        if mirror_combiner:
            mirrored_meander_points = self.mirror_points_around_xaxis(meander_points)
            meander = gdspy.FlexPath(mirrored_meander_points, meander_linewidth, corners=self.create_miter_join, **self.Aluminium)
        else:
            meander = gdspy.FlexPath(meander_points, meander_linewidth, corners=self.create_miter_join, **self.Aluminium)

        meander.rotate(rot, (0, 0))
        meander.translate(x, y)

        meander_poly_points_flexpath = self.get_polys_from_flexpath(
            meander
        )  # gets the polygon points of the outer cpw flex path via function and

        for i in range(len(meander_poly_points_flexpath)):
            meander_path_polygon = gdspy.Polygon(meander_poly_points_flexpath[i], **self.Aluminium)
            self.Main.add(meander_path_polygon)

        meander_to_frame_box_size = config["meander_to_frame_box_size"]

        meander_to_frame_box = gdspy.Rectangle(
            [
                meander_points[-1][0] - (meander_to_frame_box_size / 2),
                meander_points[-1][1] - (meander_to_frame_box_size / 2) + (meander_linewidth / 2),
            ],
            [
                meander_points[-1][0] + (meander_to_frame_box_size / 2),
                meander_points[-1][1] + (meander_to_frame_box_size / 2) + (meander_linewidth / 2),
            ],
            **self.Nb_Antenna,
        )
        if mirror_combiner:
            meander_to_frame_box.mirror([0, 0], [1, 0])

        meander_to_frame_box.rotate(rot, (0, 0))
        meander_to_frame_box.translate(x, y)
        self.Main.add(meander_to_frame_box)

        connect_meander_to_frame_points = [
            [meander_points[-1][0], meander_points[-1][1]],
            [meander_points[-1][0], -meander_height / 2 - meander_linewidth / 2],
            [meander_points[-1][0] + meander_conect_bot_width, -meander_height / 2 - meander_linewidth / 2],
            [meander_points[-1][0] + meander_conect_bot_width, -meander_height / 2 - meander_linewidth / 2 + meander_conect_right_height],
            [
                meander_points[-1][0] + meander_conect_bot_width + meander_conect_final_width,
                -meander_height / 2 - meander_linewidth / 2 + meander_conect_right_height,
            ],
        ]
        if mirror_combiner:
            mirrored_connect_meander_to_frame_points = self.mirror_points_around_xaxis(connect_meander_to_frame_points)
            connect_meander_to_frame = gdspy.FlexPath(
                mirrored_connect_meander_to_frame_points, meander_conect_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )
        else:
            connect_meander_to_frame = gdspy.FlexPath(
                connect_meander_to_frame_points, meander_conect_linewidth, corners=self.create_miter_join, **self.Nb_Antenna
            )

        connect_meander_to_frame.rotate(rot, (0, 0))
        connect_meander_to_frame.translate(x, y)
        self.Main.add(connect_meander_to_frame)

        # meander_fork
        meander_last_fork_wdith = config["meander_last_fork_wdith"]  # 5
        meander_last_fork_height = config["meander_last_fork_height"]  # 172
        meander_last_fork_linewdith = config["meander_last_fork_linewdith"]  # 2.5
        meander_last_fork_start_height = config["meander_last_fork_start_height"]  # 95

        meander_fork_start_xy = [
            connect_meander_to_frame_points[2][0] - meander_linewidth / 2,
            connect_meander_to_frame_points[2][1]
            + meander_last_fork_start_height
            + meander_linewidth / 2
            + meander_last_fork_linewdith / 2,
        ]

        meander_last_fork_points = [
            meander_fork_start_xy,
            [meander_fork_start_xy[0] - meander_last_fork_wdith - meander_last_fork_linewdith / 2, meander_fork_start_xy[1]],
            [
                meander_fork_start_xy[0] - meander_last_fork_wdith - meander_last_fork_linewdith / 2,
                meander_fork_start_xy[1] + meander_last_fork_height + meander_last_fork_linewdith / 2,
            ],
        ]
        if mirror_combiner:
            mirrored_meander_last_fork_points = self.mirror_points_around_xaxis(meander_last_fork_points)
            meander_last_fork = gdspy.FlexPath(
                mirrored_meander_last_fork_points, meander_last_fork_linewdith, corners=self.create_miter_join, **self.Nb_Antenna
            )
        else:
            meander_last_fork = gdspy.FlexPath(
                meander_last_fork_points, meander_last_fork_linewdith, corners=self.create_miter_join, **self.Nb_Antenna
            )

        meander_last_fork.rotate(rot, (0, 0))
        meander_last_fork.translate(x, y)
        self.Main.add(meander_last_fork)

        meander_last_fork_top_box_size = config["meander_last_fork_top_box_size"]

        meander_last_fork_top_box = gdspy.Rectangle(
            [meander_last_fork_points[-1][0] - meander_last_fork_top_box_size / 2, meander_last_fork_points[-1][1]],
            [
                meander_last_fork_points[-1][0] + meander_last_fork_top_box_size / 2,
                meander_last_fork_points[-1][1] + meander_last_fork_top_box_size,
            ],
            **self.Nb_Antenna,
        )
        if mirror_combiner:
            meander_last_fork_top_box.mirror([0, 0], [1, 0])

        meander_last_fork_top_box.rotate(rot, (0, 0))
        meander_last_fork_top_box.translate(x, y)
        self.Main.add(meander_last_fork_top_box)

        meander_last_fork_via_box_size_x = config["meander_last_fork_via_box_size_x"]
        meander_last_fork_via_box_size_y = config["meander_last_fork_via_box_size_y"]
        meander_last_fork_via_box_y_offset = config["meander_last_fork_via_box_y_offset"]

        meander_last_fork_via_box = gdspy.Rectangle(
            [
                meander_last_fork_points[-1][0] - meander_last_fork_via_box_size_x / 2,
                meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset,
            ],
            [
                meander_last_fork_points[-1][0] + meander_last_fork_via_box_size_x / 2,
                meander_last_fork_points[-1][1] + meander_last_fork_via_box_y_offset + meander_last_fork_via_box_size_y,
            ],
            **self.SiN_dep,
        )
        if mirror_combiner:
            meander_last_fork_via_box.mirror([0, 0], [1, 0])

        meander_last_fork_via_box.rotate(rot, (0, 0))
        meander_last_fork_via_box.translate(x, y)
        self.silicon_nitride_cutouts.add(meander_last_fork_via_box)

        return points_to_conect_kids_to

    def add_filter_bank_ring_overlap_and_get_conections(self, x, y, ant_center_x, ant_center_y, overlap_no, rot, Main_config_file_dict):
        """
        Adds a ring overlap conection to bridge the inner and outer rings over one
        another. This function will return the conection points where the inner and
        outer rings should connect to.


        Parameters
        ----------
        x, y : float, int
            The x, y coordinate of the center of the ring overlap.

        ant_center_x, ant_center_y : float, int
            The x, y coordinate of the center of the antenna structure,
            i.e. the center of the horn.

        overlap_no : int
            The number of the ring overlap. This determines where to draw it around
            the antenna. Starting at 0 for left middle and +1 for each subsequent
            overlap going anti-clockwise. Should not be more than 3, values more
            than this wrap back to 0 (left middle placement) because the overlap_no
            operates like it is modulo 4.

        rot : float
            The angle (**in radians**) which the overlap geometry should be rotated
            at. This rot angle defined as the anti-clockwise angle made with the
            positive x-axis.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "filter_bank_ring_overlap".

        Returns
        -------
        conections_dict : dict
            This dictionary contains keys that map to an [x,y] list which are the
            coordinates defining the conection points where the inner and outer
            ring should connect to this overlap structure.

            This dict has keys: **'inner_conect_0', 'inner_conect_1',
            'outer_conect_0', 'outer_conect_1'**.

        """

        config = Main_config_file_dict["filter_bank_ring_overlap"]

        outer_box_width = config["outer_box_width"]  # 100
        outer_box_height = config["outer_box_height"]  # 100
        outer_box_inner_cutout_width = config["outer_box_inner_cutout_width"]  # 16
        outer_box_inner_cutout_height = config["outer_box_inner_cutout_height"]  # 6
        linewidth = config["linewidth"]  # 5

        half_conect_taper_extra_height = config["half_conect_taper_extra_height"]  # 5
        half_conect_taper_stright_width = config["half_conect_taper_stright_width"]  # 5
        half_conect_taper_diag_width = config["half_conect_taper_diag_width"]  # 5.5
        half_conect_taper_inner_height = config["half_conect_taper_inner_height"]  # 4
        half_conect_distacne_from_center = config["half_conect_distacne_from_center"]  # 3

        half_conect_bridge_rect_width = config["half_conect_bridge_rect_width"]  # 14
        half_conect_bridge_rect_height = config["half_conect_bridge_rect_height"]  # 4

        half_conect_bridge_pad_width = config["half_conect_bridge_pad_width"]  # 2
        half_conect_bridge_pad_height = config["half_conect_bridge_pad_height"]  # 2
        half_conect_bridge_pad_offset_from_center = config["half_conect_bridge_pad_offset_from_center"]  # 5

        full_conect_taper_extra_width = config["full_conect_taper_extra_width"]  # 5
        full_conect_taper_stright_height = config["full_conect_taper_stright_height"]  # 5
        full_conect_taper_diag_height = config["full_conect_taper_diag_height"]  # 5.5
        full_conect_taper_start_from_center = config["full_conect_taper_start_from_center"]  # 3
        full_conect_center_width = config["full_conect_center_width"]  # 4

        conect_half_conect_left = self.rotate_and_move_single_point([-outer_box_width / 2, 0], rot, x, y)
        conect_half_conect_right = self.rotate_and_move_single_point([outer_box_width / 2, 0], rot, x, y)
        conect_full_conect_bot = self.rotate_and_move_single_point([0, -outer_box_height / 2], rot, x, y)
        conect_full_conect_top = self.rotate_and_move_single_point([0, outer_box_height / 2], rot, x, y)

        conections_rotated_to_quadrant_1 = [
            self.rotate(ant_center_x, ant_center_y, conect_half_conect_left[0], conect_half_conect_left[1], (overlap_no * -pi / 2)),
            self.rotate(ant_center_x, ant_center_y, conect_half_conect_right[0], conect_half_conect_right[1], (overlap_no * -pi / 2)),
            self.rotate(ant_center_x, ant_center_y, conect_full_conect_bot[0], conect_full_conect_bot[1], (overlap_no * -pi / 2)),
            self.rotate(ant_center_x, ant_center_y, conect_full_conect_top[0], conect_full_conect_top[1], (overlap_no * -pi / 2)),
        ]

        conections_sorted_list = sorted(conections_rotated_to_quadrant_1.copy(), key=lambda k: [k[1], k[0]])

        conections_dict = {}
        conections_dict["inner_conect_0"] = self.rotate(
            ant_center_x, ant_center_y, conections_sorted_list[0][0], conections_sorted_list[0][1], (overlap_no * pi / 2)
        )
        conections_dict["outer_conect_0"] = self.rotate(
            ant_center_x, ant_center_y, conections_sorted_list[1][0], conections_sorted_list[1][1], (overlap_no * pi / 2)
        )
        conections_dict["inner_conect_1"] = self.rotate(
            ant_center_x, ant_center_y, conections_sorted_list[2][0], conections_sorted_list[2][1], (overlap_no * pi / 2)
        )
        conections_dict["outer_conect_1"] = self.rotate(
            ant_center_x, ant_center_y, conections_sorted_list[3][0], conections_sorted_list[3][1], (overlap_no * pi / 2)
        )

        outer_box_poly_points = [
            [-outer_box_width / 2, -outer_box_height / 2],
            [-outer_box_width / 2, outer_box_height / 2],
            [outer_box_width / 2, outer_box_height / 2],
            [outer_box_width / 2, outer_box_inner_cutout_height / 2],
            [-outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
            [-outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
            [outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
            [outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
            [outer_box_width / 2, outer_box_inner_cutout_height / 2],
            [outer_box_width / 2, -outer_box_height / 2],
        ]

        outer_box_poly_points = self.rotate_and_move_points_list(outer_box_poly_points, rot, x, y)
        outer_box = gdspy.Polygon(outer_box_poly_points, **self.Aluminium)
        # self.Main.add(outer_box)

        outer_box_inner_cutout_points = [
            [-outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
            [-outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
            [outer_box_inner_cutout_width / 2, outer_box_inner_cutout_height / 2],
            [outer_box_inner_cutout_width / 2, -outer_box_inner_cutout_height / 2],
        ]

        outer_box_inner_cutout_points = self.rotate_and_move_points_list(outer_box_inner_cutout_points, rot, x, y)

        outer_box_inner_cutout = gdspy.Polygon(outer_box_inner_cutout_points, **self.Nb_Groundplane)
        self.ground_plane_cutouts.add(outer_box_inner_cutout)

        half_conect_poly_points_right = [
            [half_conect_distacne_from_center, half_conect_taper_inner_height / 2],
            [half_conect_distacne_from_center + half_conect_taper_diag_width, linewidth / 2 + half_conect_taper_extra_height],
            [
                half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width,
                linewidth / 2 + half_conect_taper_extra_height,
            ],
            [half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width, linewidth / 2],
            [outer_box_width / 2, linewidth / 2],
            [outer_box_width / 2, -linewidth / 2],
            [half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width, -linewidth / 2],
            [
                half_conect_distacne_from_center + half_conect_taper_diag_width + half_conect_taper_stright_width,
                -linewidth / 2 - half_conect_taper_extra_height,
            ],
            [half_conect_distacne_from_center + half_conect_taper_diag_width, -linewidth / 2 - half_conect_taper_extra_height],
            [half_conect_distacne_from_center, -half_conect_taper_inner_height / 2],
        ]

        half_conect_poly_points_right = self.rotate_and_move_points_list(half_conect_poly_points_right, rot, x, y)
        half_conect_right = gdspy.Polygon(half_conect_poly_points_right, **self.Nb_Antenna)
        self.Main.add(half_conect_right)

        half_conect_poly_points_left = [
            [-half_conect_distacne_from_center, half_conect_taper_inner_height / 2],
            [-half_conect_distacne_from_center - half_conect_taper_diag_width, linewidth / 2 + half_conect_taper_extra_height],
            [
                -half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width,
                linewidth / 2 + half_conect_taper_extra_height,
            ],
            [-half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width, linewidth / 2],
            [-outer_box_width / 2, linewidth / 2],
            [-outer_box_width / 2, -linewidth / 2],
            [-half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width, -linewidth / 2],
            [
                -half_conect_distacne_from_center - half_conect_taper_diag_width - half_conect_taper_stright_width,
                -linewidth / 2 - half_conect_taper_extra_height,
            ],
            [-half_conect_distacne_from_center - half_conect_taper_diag_width, -linewidth / 2 - half_conect_taper_extra_height],
            [-half_conect_distacne_from_center, -half_conect_taper_inner_height / 2],
        ]

        half_conect_poly_points_left = self.rotate_and_move_points_list(half_conect_poly_points_left, rot, x, y)
        half_conect_left = gdspy.Polygon(half_conect_poly_points_left, **self.Nb_Antenna)
        self.Main.add(half_conect_left)

        half_conect_bridge_rect = gdspy.Rectangle(
            [-half_conect_bridge_rect_width / 2, -half_conect_bridge_rect_height / 2],
            [half_conect_bridge_rect_width / 2, half_conect_bridge_rect_height / 2],
            **self.Nb_Groundplane,
        )
        half_conect_bridge_rect.rotate(rot, (0, 0))
        half_conect_bridge_rect.translate(x, y)
        self.Main.add(half_conect_bridge_rect)

        half_conect_bridge_pad_right = gdspy.Rectangle(
            [half_conect_bridge_pad_offset_from_center, -half_conect_bridge_pad_height / 2],
            [half_conect_bridge_pad_offset_from_center + half_conect_bridge_pad_width, half_conect_bridge_pad_height / 2],
            **self.SiN_dep,
        )
        half_conect_bridge_pad_right.rotate(rot, (0, 0))
        half_conect_bridge_pad_right.translate(x, y)
        # self.Main.add(half_conect_bridge_pad_right)
        self.silicon_nitride_cutouts.add(half_conect_bridge_pad_right)

        half_conect_bridge_pad_left = gdspy.Rectangle(
            [half_conect_bridge_pad_offset_from_center, -half_conect_bridge_pad_height / 2],
            [half_conect_bridge_pad_offset_from_center + half_conect_bridge_pad_width, half_conect_bridge_pad_height / 2],
            **self.SiN_dep,
        )
        half_conect_bridge_pad_left.rotate(pi, (0, 0))
        half_conect_bridge_pad_left.rotate(rot, (0, 0))
        half_conect_bridge_pad_left.translate(x, y)
        # self.Main.add(half_conect_bridge_pad_left)
        self.silicon_nitride_cutouts.add(half_conect_bridge_pad_left)

        full_conect_poly_points = [
            [full_conect_center_width / 2, full_conect_taper_start_from_center],
            [linewidth / 2 + full_conect_taper_extra_width, full_conect_taper_start_from_center + full_conect_taper_diag_height],
            [
                linewidth / 2 + full_conect_taper_extra_width,
                full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height,
            ],
            [linewidth / 2, full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height],
            [linewidth / 2, outer_box_height / 2],
            [-linewidth / 2, outer_box_height / 2],
            [-linewidth / 2, full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height],
            [
                -linewidth / 2 - full_conect_taper_extra_width,
                full_conect_taper_start_from_center + full_conect_taper_diag_height + full_conect_taper_stright_height,
            ],
            [-linewidth / 2 - full_conect_taper_extra_width, full_conect_taper_start_from_center + full_conect_taper_diag_height],
            [-full_conect_center_width / 2, full_conect_taper_start_from_center],
            [-full_conect_center_width / 2, -full_conect_taper_start_from_center],
            [-linewidth / 2 - full_conect_taper_extra_width, -full_conect_taper_start_from_center - full_conect_taper_diag_height],
            [
                -linewidth / 2 - full_conect_taper_extra_width,
                -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height,
            ],
            [-linewidth / 2, -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height],
            [-linewidth / 2, -outer_box_height / 2],
            [linewidth / 2, -outer_box_height / 2],
            [linewidth / 2, -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height],
            [
                linewidth / 2 + full_conect_taper_extra_width,
                -full_conect_taper_start_from_center - full_conect_taper_diag_height - full_conect_taper_stright_height,
            ],
            [linewidth / 2 + full_conect_taper_extra_width, -full_conect_taper_start_from_center - full_conect_taper_diag_height],
            [full_conect_center_width / 2, -full_conect_taper_start_from_center],
        ]

        full_conect_poly_points = self.rotate_and_move_points_list(full_conect_poly_points, rot, x, y)

        full_conect = gdspy.Polygon(full_conect_poly_points, **self.Nb_Antenna)
        self.Main.add(full_conect)

        return conections_dict

    def add_Filter_bank_and_get_conection_points(
        self, x, y, Main_config_file_dict, with_combiner=True, with_crossover=True, only_1_pol=False
    ):
        """
        Adds the filter bank structure to the chip centered at the x,y coordinate
        given. By default this will be drawn with phase combiners and ring overlap
        crossovers and with both polorizations filtered.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate to center filter bank structure around.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "filter_bank".

        KwArgs
        ------

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

        Returns
        -------
        points_to_conect_kids_dict : dict
            This dictionary contains keys that map to an [x,y] list which are the
            coordinates defining the conection points where the KIDs should connect
            to the filter bank.

            This dict has keys: **'TR', 'TL', 'BL', 'BR'**.



        """

        config = Main_config_file_dict["filter_bank"]

        # ring dimension and properties
        inner_ring_radius = config["inner_ring_radius"]  # 3920/2
        inner_ring_line_width = config["inner_ring_line_width"]  # 5
        outer_ring_radius = config["outer_ring_radius"]  # 4315/2
        outer_ring_line_width = config["outer_ring_line_width"]  # 5

        ring_overlap_distance_from_center = config["ring_overlap_distance_from_center"]  # inner_ring_radius + 150

        inner_ring_overlap_gap = config["inner_ring_overlap_gap"]  # 325
        outer_ring_overlap_gap = config["outer_ring_overlap_gap"]  # 325
        if with_combiner:
            outer_ring_conector_gap = config["outer_ring_conector_gap"]  # 20
        else:
            outer_ring_conector_gap = 0

        inner_ring_arc_length = inner_ring_overlap_gap
        inner_ring_arc_angle = inner_ring_arc_length / inner_ring_radius

        outer_ring_arc_length = outer_ring_overlap_gap
        outer_ring_arc_angle = outer_ring_arc_length / outer_ring_radius

        outer_combiner_arc_length = outer_ring_conector_gap
        outer_combiner_arc_angle = outer_combiner_arc_length / outer_ring_radius

        inner_ring_overlap_theta = inner_ring_arc_angle / 2
        outer_ring_overlap_theta = outer_ring_arc_angle / 2
        outer_combiner_theta = outer_combiner_arc_angle / 2

        inner_straight_extend_distance = config["inner_straight_extend_distance"]  # 20
        inner_bend_radius = config["inner_bend_radius"]  # 15
        outer_straight_extend_distance = config["outer_straight_extend_distance"]  # 10
        outer_bend_radius = config["outer_bend_radius"]  # 9

        Hugh_filter_thin_length1 = config["Hugh_filter_thin_length1"]  # 51
        Hugh_filter_thin_length2 = config["Hugh_filter_thin_length2"]  # 106

        Hugh_filter_thick_length1 = config["Hugh_filter_thick_length1"]  # 48
        Hugh_filter_thick_length2 = config["Hugh_filter_thick_length2"]  # 58.5

        Hugh_filter_thin_lengths = np.array([Hugh_filter_thin_length1, Hugh_filter_thin_length2])
        Hugh_filter_thin_arc_angles = Hugh_filter_thin_lengths / outer_ring_radius
        Hugh_filter_thick_lengths = np.array([Hugh_filter_thick_length1, Hugh_filter_thick_length2])
        Hugh_filter_thick_arc_angles = Hugh_filter_thick_lengths / outer_ring_radius
        Hugh_filter_thin_width = config["Hugh_filter_thin_width"]  # 2
        Hugh_filter_thick_width = config["Hugh_filter_thick_width"]  # 20
        half_Hugh_filter_total_arc_angle = (
            np.sum(Hugh_filter_thin_arc_angles) + Hugh_filter_thick_arc_angles[0] + (Hugh_filter_thick_arc_angles[1] / 2)
        )

        ring_overlap_0_rot = config["ring_overlap_0_rot"]
        ring_overlap_1_rot = config["ring_overlap_1_rot"]
        ring_overlap_2_rot = config["ring_overlap_2_rot"]
        ring_overlap_3_rot = config["ring_overlap_3_rot"]

        ring_overlap_rotations = [
            self.deg_to_rad(ring_overlap_0_rot),
            self.deg_to_rad(ring_overlap_1_rot),
            self.deg_to_rad(ring_overlap_2_rot),
            self.deg_to_rad(ring_overlap_3_rot),
        ]

        points_to_conect_kids_dict = {}
        conect_kids_dict_key_names = ["TR", "TL", "BL", "BR"]

        for i in range(4):
            if not only_1_pol or i in [0, 2]:
                inner_arc_section = gdspy.Round(
                    [x, y],
                    inner_ring_radius + (inner_ring_line_width / 2),
                    inner_ring_radius - (inner_ring_line_width / 2),
                    initial_angle=inner_ring_overlap_theta,
                    final_angle=inner_ring_overlap_theta + (pi / 2) - inner_ring_arc_angle,
                    **self.Nb_Antenna,
                )
                self.Main.add(inner_arc_section)

            if i in [0, 1] and with_crossover:
                outer_arc_section1 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (outer_ring_line_width / 2),
                    outer_ring_radius - (outer_ring_line_width / 2),
                    initial_angle=outer_ring_overlap_theta,
                    final_angle=(i * pi / 2) + (pi / 4) - outer_combiner_theta,
                    **self.Nb_Antenna,
                )
                outer_arc_section2 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (outer_ring_line_width / 2),
                    outer_ring_radius - (outer_ring_line_width / 2),
                    initial_angle=(i * pi / 2) + (pi / 4) + outer_combiner_theta,
                    final_angle=outer_ring_overlap_theta + (pi / 2) - outer_ring_arc_angle,
                    **self.Nb_Antenna,
                )
                self.Main.add(outer_arc_section1)
                self.Main.add(outer_arc_section2)

            if i in [2, 3]:
                if (not only_1_pol) or (with_crossover):
                    outer_arc_section11 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (outer_ring_line_width / 2),
                        outer_ring_radius - (outer_ring_line_width / 2),
                        initial_angle=outer_ring_overlap_theta,
                        final_angle=(i * pi / 2) + (pi / 8) - (half_Hugh_filter_total_arc_angle),
                        **self.Nb_Antenna,
                    )
                    outer_arc_section12 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (outer_ring_line_width / 2),
                        outer_ring_radius - (outer_ring_line_width / 2),
                        initial_angle=(i * pi / 2) + (pi / 8) + (half_Hugh_filter_total_arc_angle),
                        final_angle=(i * pi / 2) + (pi / 4) - outer_combiner_theta,
                        **self.Nb_Antenna,
                    )
                    self.Main.add(outer_arc_section11)
                    self.Main.add(outer_arc_section12)

                outer_arc_section21 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (outer_ring_line_width / 2),
                    outer_ring_radius - (outer_ring_line_width / 2),
                    initial_angle=(i * pi / 2) + (pi / 4) + outer_combiner_theta,
                    final_angle=(i * pi / 2) + (pi / 8) + (pi / 4) - (half_Hugh_filter_total_arc_angle),
                    **self.Nb_Antenna,
                )
                outer_arc_section22 = gdspy.Round(
                    [x, y],
                    outer_ring_radius + (outer_ring_line_width / 2),
                    outer_ring_radius - (outer_ring_line_width / 2),
                    initial_angle=(i * pi / 2) + (pi / 8) + (pi / 4) + (half_Hugh_filter_total_arc_angle),
                    final_angle=outer_ring_overlap_theta + (pi / 2) - outer_ring_arc_angle,
                    **self.Nb_Antenna,
                )
                self.Main.add(outer_arc_section21)
                self.Main.add(outer_arc_section22)

                for j in range(2):
                    mid_Hugh_angle = pi / 8 if j == 0 else pi / 8 + pi / 4
                    angle_inc = 0

                    if (only_1_pol and not with_crossover) and j == 0:
                        if i != 3:
                            continue

                        # add arc to hugh from lopass bot left
                        inner_conect_ang = inner_ring_overlap_theta - inner_ring_arc_angle
                        outer_conect_ang = outer_ring_overlap_theta - outer_ring_arc_angle
                        extra_length_connect_arc = 200
                        conection_low_pass_to_hugh_points = [
                            [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                            [
                                x + inner_ring_radius * cos(inner_conect_ang) + extra_length_connect_arc * cos(inner_conect_ang + pi / 2),
                                y + inner_ring_radius * sin(inner_conect_ang) + extra_length_connect_arc * sin(inner_conect_ang + pi / 2),
                            ],
                            [
                                x + outer_ring_radius * cos(inner_conect_ang) + extra_length_connect_arc * cos(outer_conect_ang + pi / 2),
                                y + outer_ring_radius * sin(inner_conect_ang) + extra_length_connect_arc * sin(outer_conect_ang + pi / 2),
                            ],
                            [x + outer_ring_radius * cos(outer_conect_ang), y + outer_ring_radius * sin(outer_conect_ang)],
                        ]
                        conection_low_pass_to_hugh = gdspy.FlexPath(
                            conection_low_pass_to_hugh_points,
                            inner_ring_line_width,
                            corners="circular bend",
                            bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                            **self.Nb_Antenna,
                        )
                        self.Main.add(conection_low_pass_to_hugh)
                        BL_KID_connect_1_pol = [x + outer_ring_radius * cos(-3 * pi / 4), y + outer_ring_radius * sin(-3 * pi / 4)]

                        # add from hugh bot right
                        inner_conect_ang = -(pi / 4)
                        outer_conect_ang = -(pi / 4)
                        extra_length_connect_arc = 200
                        outer_arc_conection_radius = outer_ring_radius + (outer_ring_radius - inner_ring_radius)
                        conection_from_hugh_points = [
                            [x + outer_ring_radius * cos(inner_conect_ang), y + outer_ring_radius * sin(inner_conect_ang)],
                            [
                                x + outer_ring_radius * cos(inner_conect_ang) + extra_length_connect_arc * cos(inner_conect_ang - pi / 2),
                                y + outer_ring_radius * sin(inner_conect_ang) + extra_length_connect_arc * sin(inner_conect_ang - pi / 2),
                            ],
                            [
                                x
                                + outer_arc_conection_radius * cos(inner_conect_ang)
                                + extra_length_connect_arc * cos(outer_conect_ang - pi / 2),
                                y
                                + outer_arc_conection_radius * sin(inner_conect_ang)
                                + extra_length_connect_arc * sin(outer_conect_ang - pi / 2),
                            ],
                            [
                                x + outer_arc_conection_radius * cos(outer_conect_ang),
                                y + outer_arc_conection_radius * sin(outer_conect_ang),
                            ],
                        ]
                        conection_from_hugh = gdspy.FlexPath(
                            conection_from_hugh_points,
                            inner_ring_line_width,
                            corners="circular bend",
                            bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                            **self.Nb_Antenna,
                        )
                        self.Main.add(conection_from_hugh)
                        BR_KID_connect_1_pol = [
                            x + outer_arc_conection_radius * cos(-pi / 4),
                            y + outer_arc_conection_radius * sin(-pi / 4),
                        ]

                        # add path conecting left mid inner to outer
                        inner_conect_ang = (pi / 2) + inner_ring_overlap_theta
                        outer_conect_ang = (pi / 2) + outer_ring_overlap_theta - outer_ring_arc_angle
                        extra_length_connect_path = 100
                        inner_to_outer_connect_path_points = [
                            [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                            [
                                x + inner_ring_radius * cos(inner_conect_ang) + extra_length_connect_path * cos(inner_conect_ang - pi / 2),
                                y + inner_ring_radius * sin(inner_conect_ang) + extra_length_connect_path * sin(inner_conect_ang - pi / 2),
                            ],
                            [
                                x + outer_ring_radius * cos(outer_conect_ang) + extra_length_connect_path * cos(outer_conect_ang + pi / 2),
                                y + outer_ring_radius * sin(outer_conect_ang) + extra_length_connect_path * sin(outer_conect_ang + pi / 2),
                            ],
                            [x + outer_ring_radius * cos(outer_conect_ang), y + outer_ring_radius * sin(outer_conect_ang)],
                        ]
                        inner_to_outer_connect_path = gdspy.FlexPath(
                            inner_to_outer_connect_path_points,
                            inner_ring_line_width,
                            corners="circular bend",
                            bend_radius=(outer_ring_radius - inner_ring_radius),
                            **self.Nb_Antenna,
                        )
                        self.Main.add(inner_to_outer_connect_path)

                        # add arc from hipass top
                        inner_conect_ang = pi + inner_ring_overlap_theta - inner_ring_arc_angle
                        outer_conect_ang = pi + outer_ring_overlap_theta - outer_ring_arc_angle
                        extra_inner_length_connect_arc = 400
                        extra_outer_length_connect_arc_after = 800
                        conection_high_pass_points = [
                            [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                            [
                                x
                                + inner_ring_radius * cos(inner_conect_ang)
                                + extra_inner_length_connect_arc * cos(inner_conect_ang + pi / 2),
                                y
                                + inner_ring_radius * sin(inner_conect_ang)
                                + extra_inner_length_connect_arc * sin(inner_conect_ang + pi / 2),
                            ],
                            [x + outer_ring_radius * cos(outer_conect_ang), y + outer_ring_radius * sin(outer_conect_ang)],
                            [
                                x + outer_ring_radius * cos(outer_conect_ang) + extra_outer_length_connect_arc_after,
                                y + outer_ring_radius * sin(outer_conect_ang),
                            ],
                        ]
                        conection_high_pass = gdspy.FlexPath(
                            conection_high_pass_points,
                            inner_ring_line_width,
                            corners="circular bend",
                            bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                            **self.Nb_Antenna,
                        )
                        self.Main.add(conection_high_pass)
                        TR_KID_connect_1_pol = conection_high_pass_points[-1]

                        # add small section on the left mid toward top left
                        inner_conect_ang = pi + inner_ring_arc_angle / 2
                        extra_length_straight = 100

                        small_connection_left_mid_points = [
                            [x + inner_ring_radius * cos(inner_conect_ang), y + inner_ring_radius * sin(inner_conect_ang)],
                            [
                                x + inner_ring_radius * cos(inner_conect_ang) + extra_length_straight * cos(inner_conect_ang - pi / 2),
                                y + inner_ring_radius * sin(inner_conect_ang) + extra_length_straight * sin(inner_conect_ang - pi / 2),
                            ],
                            [
                                x + inner_ring_radius * cos(inner_conect_ang) + extra_length_straight * cos(inner_conect_ang - pi / 2),
                                y
                                + inner_ring_radius * sin(inner_conect_ang)
                                + extra_length_straight * sin(inner_conect_ang - pi / 2)
                                + extra_length_straight,
                            ],
                        ]

                        small_connection_left_mid = gdspy.FlexPath(
                            small_connection_left_mid_points,
                            inner_ring_line_width,
                            corners="circular bend",
                            bend_radius=(outer_ring_radius - inner_ring_radius) / 2,
                            **self.Nb_Antenna,
                        )
                        self.Main.add(small_connection_left_mid)
                        TL_KID_connect_1_pol = small_connection_left_mid_points[-1]

                        self.Main.add(gdspy.Round(TL_KID_connect_1_pol, 10, layer=-1))
                        self.Main.add(gdspy.Round(TR_KID_connect_1_pol, 10, layer=-1))
                        self.Main.add(gdspy.Round(BL_KID_connect_1_pol, 10, layer=-1))
                        self.Main.add(gdspy.Round(BR_KID_connect_1_pol, 10, layer=-1))
                        one_pol_connect_dict = {
                            "TL": TL_KID_connect_1_pol,
                            "TR": TR_KID_connect_1_pol,
                            "BL": BL_KID_connect_1_pol,
                            "BR": BR_KID_connect_1_pol,
                        }

                        continue

                    middle_Hugh_thick_arc = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thick_width / 2),
                        outer_ring_radius - (Hugh_filter_thick_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle - (Hugh_filter_thick_arc_angles[1] / 2),
                        final_angle=(i * pi / 2) + mid_Hugh_angle + (Hugh_filter_thick_arc_angles[1] / 2),
                        **self.Nb_Antenna,
                    )
                    self.Main.add(middle_Hugh_thick_arc)
                    angle_inc += Hugh_filter_thick_arc_angles[1] / 2

                    middle_Hugh_thin_arc1 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thin_width / 2),
                        outer_ring_radius - (Hugh_filter_thin_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[1],
                        final_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc,
                        **self.Nb_Antenna,
                    )
                    middle_Hugh_thin_arc2 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thin_width / 2),
                        outer_ring_radius - (Hugh_filter_thin_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc,
                        final_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc + Hugh_filter_thin_arc_angles[1],
                        **self.Nb_Antenna,
                    )
                    self.Main.add(middle_Hugh_thin_arc1)
                    self.Main.add(middle_Hugh_thin_arc2)
                    angle_inc += Hugh_filter_thin_arc_angles[1]

                    outer_Hugh_thick_arc1 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thick_width / 2),
                        outer_ring_radius - (Hugh_filter_thick_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc - Hugh_filter_thick_arc_angles[0],
                        final_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc,
                        **self.Nb_Antenna,
                    )
                    outer_Hugh_thick_arc2 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thick_width / 2),
                        outer_ring_radius - (Hugh_filter_thick_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc,
                        final_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc + Hugh_filter_thick_arc_angles[0],
                        **self.Nb_Antenna,
                    )
                    self.Main.add(outer_Hugh_thick_arc1)
                    self.Main.add(outer_Hugh_thick_arc2)
                    angle_inc += Hugh_filter_thick_arc_angles[0]

                    outer_Hugh_thin_arc1 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thin_width / 2),
                        outer_ring_radius - (Hugh_filter_thin_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc - Hugh_filter_thin_arc_angles[0],
                        final_angle=(i * pi / 2) + mid_Hugh_angle - angle_inc,
                        **self.Nb_Antenna,
                    )
                    outer_Hugh_thin_arc2 = gdspy.Round(
                        [x, y],
                        outer_ring_radius + (Hugh_filter_thin_width / 2),
                        outer_ring_radius - (Hugh_filter_thin_width / 2),
                        initial_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc,
                        final_angle=(i * pi / 2) + mid_Hugh_angle + angle_inc + Hugh_filter_thin_arc_angles[0],
                        **self.Nb_Antenna,
                    )
                    self.Main.add(outer_Hugh_thin_arc1)
                    self.Main.add(outer_Hugh_thin_arc2)
                    angle_inc += Hugh_filter_thin_arc_angles[0]

            if i == 0:
                self.add_Lo_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "clockwise", Main_config_file_dict
                )
                self.add_Hi_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "anti-clockwise", Main_config_file_dict
                )
            if i == 1 and not only_1_pol:
                self.add_Hi_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "clockwise", Main_config_file_dict
                )
                self.add_Lo_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "anti-clockwise", Main_config_file_dict
                )
            if i == 2:
                self.add_Hi_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "clockwise", Main_config_file_dict
                )
                self.add_Lo_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "anti-clockwise", Main_config_file_dict
                )
            if i == 3 and not only_1_pol:
                self.add_Lo_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "clockwise", Main_config_file_dict
                )
                self.add_Hi_pass_filters(
                    x, y, inner_ring_line_width, inner_ring_radius, (i * pi / 2) + pi / 4, "anti-clockwise", Main_config_file_dict
                )

            combiner_xpos = x + (outer_ring_radius * cos(pi / 4 + i * (pi / 2)))
            combiner_ypos = y + (outer_ring_radius * sin(pi / 4 + i * (pi / 2)))

            if with_combiner:
                if i in [1, 3]:
                    mirror_combiner = False
                else:
                    mirror_combiner = True

                if i in [0, 1]:
                    conection_point = self.add_combiner_section_and_get_conect_point(
                        combiner_xpos,
                        combiner_ypos,
                        (pi / 4 + i * pi / 2),
                        outer_ring_conector_gap,
                        outer_ring_line_width,
                        Main_config_file_dict,
                        combiner_type="150GHZ",
                        mirror_combiner=mirror_combiner,
                    )
                else:
                    conection_point = self.add_combiner_section_and_get_conect_point(
                        combiner_xpos,
                        combiner_ypos,
                        (pi / 4 + i * pi / 2),
                        outer_ring_conector_gap,
                        outer_ring_line_width,
                        Main_config_file_dict,
                        combiner_type="90GHZ",
                        mirror_combiner=mirror_combiner,
                    )
            else:
                conection_point = [
                    x + outer_ring_radius * cos((i * pi / 2) + (pi / 4) - outer_combiner_theta),
                    y + outer_ring_radius * sin((i * pi / 2) + (pi / 4) - outer_combiner_theta),
                ]

            points_to_conect_kids_dict[conect_kids_dict_key_names[i]] = conection_point

            if with_crossover:
                ring_overlap_xpos = x + (ring_overlap_distance_from_center * cos(i * (pi / 2)))
                ring_overlap_ypos = y + (ring_overlap_distance_from_center * sin(i * (pi / 2)))
                ring_overlap_number = i

                conections_dict = self.add_filter_bank_ring_overlap_and_get_conections(
                    ring_overlap_xpos, ring_overlap_ypos, x, y, ring_overlap_number, ring_overlap_rotations[i], Main_config_file_dict
                )

                # inner conections to ring overlap
                inner_0_conection_points = [
                    conections_dict["inner_conect_0"],
                    [
                        conections_dict["inner_conect_0"][0] + inner_straight_extend_distance * cos(-3 * pi / 4 + i * pi / 2),
                        conections_dict["inner_conect_0"][1] + inner_straight_extend_distance * sin(-3 * pi / 4 + i * pi / 2),
                    ],
                    [
                        x
                        + (
                            inner_ring_radius
                            * cos((i * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)
                        ),
                        y
                        + (
                            inner_ring_radius
                            * sin((i * pi / 2) - inner_ring_arc_angle / 2 + inner_straight_extend_distance / inner_ring_radius)
                        ),
                    ],
                    [
                        x + (inner_ring_radius * cos((i * pi / 2) - inner_ring_arc_angle / 2)),
                        y + (inner_ring_radius * sin((i * pi / 2) - inner_ring_arc_angle / 2)),
                    ],
                ]

                inner_0_path = gdspy.FlexPath(
                    inner_0_conection_points,
                    inner_ring_line_width,
                    corners="circular bend",
                    bend_radius=inner_bend_radius,
                    **self.Nb_Antenna,
                )
                self.Main.add(inner_0_path)

                inner_1_conection_points = [
                    conections_dict["inner_conect_1"],
                    [
                        conections_dict["inner_conect_1"][0] + inner_straight_extend_distance * cos(3 * pi / 4 + i * pi / 2),
                        conections_dict["inner_conect_1"][1] + inner_straight_extend_distance * sin(3 * pi / 4 + i * pi / 2),
                    ],
                    [
                        x
                        + (
                            inner_ring_radius
                            * cos((i * pi / 2) + inner_ring_arc_angle / 2 - inner_straight_extend_distance / inner_ring_radius)
                        ),
                        y
                        + (
                            inner_ring_radius
                            * sin((i * pi / 2) + inner_ring_arc_angle / 2 - inner_straight_extend_distance / inner_ring_radius)
                        ),
                    ],
                    [
                        x + (inner_ring_radius * cos((i * pi / 2) + inner_ring_arc_angle / 2)),
                        y + (inner_ring_radius * sin((i * pi / 2) + inner_ring_arc_angle / 2)),
                    ],
                ]

                inner_1_path = gdspy.FlexPath(
                    inner_1_conection_points,
                    inner_ring_line_width,
                    corners="circular bend",
                    bend_radius=inner_bend_radius,
                    **self.Nb_Antenna,
                )
                self.Main.add(inner_1_path)

                # outer conections to ring overlap
                outer_0_conection_points = [
                    conections_dict["outer_conect_0"],
                    [
                        conections_dict["outer_conect_0"][0] + outer_straight_extend_distance * cos(-pi / 4 + i * pi / 2),
                        conections_dict["outer_conect_0"][1] + outer_straight_extend_distance * sin(-pi / 4 + i * pi / 2),
                    ],
                    [
                        x
                        + (
                            outer_ring_radius
                            * cos((i * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)
                        ),
                        y
                        + (
                            outer_ring_radius
                            * sin((i * pi / 2) - outer_ring_arc_angle / 2 + outer_straight_extend_distance / outer_ring_radius)
                        ),
                    ],
                    [
                        x + (outer_ring_radius * cos((i * pi / 2) - outer_ring_arc_angle / 2)),
                        y + (outer_ring_radius * sin((i * pi / 2) - outer_ring_arc_angle / 2)),
                    ],
                ]

                outer_0_path = gdspy.FlexPath(
                    outer_0_conection_points,
                    inner_ring_line_width,
                    corners="circular bend",
                    bend_radius=outer_bend_radius,
                    **self.Nb_Antenna,
                )
                self.Main.add(outer_0_path)

                outer_1_conection_points = [
                    conections_dict["outer_conect_1"],
                    [
                        conections_dict["outer_conect_1"][0] + outer_straight_extend_distance * cos(pi / 4 + i * pi / 2),
                        conections_dict["outer_conect_1"][1] + outer_straight_extend_distance * sin(pi / 4 + i * pi / 2),
                    ],
                    [
                        x
                        + (
                            outer_ring_radius
                            * cos((i * pi / 2) + outer_ring_arc_angle / 2 - outer_straight_extend_distance / outer_ring_radius)
                        ),
                        y
                        + (
                            outer_ring_radius
                            * sin((i * pi / 2) + outer_ring_arc_angle / 2 - outer_straight_extend_distance / outer_ring_radius)
                        ),
                    ],
                    [
                        x + (outer_ring_radius * cos((i * pi / 2) + outer_ring_arc_angle / 2)),
                        y + (outer_ring_radius * sin((i * pi / 2) + outer_ring_arc_angle / 2)),
                    ],
                ]

                right_path = gdspy.FlexPath(
                    outer_1_conection_points,
                    inner_ring_line_width,
                    corners="circular bend",
                    bend_radius=outer_bend_radius,
                    **self.Nb_Antenna,
                )
                self.Main.add(right_path)

            inner_ring_overlap_theta += pi / 2
            outer_ring_overlap_theta += pi / 2

        if only_1_pol and not with_crossover:
            return one_pol_connect_dict
        else:
            return points_to_conect_kids_dict

    def connect_filter_bank_to_KIDs(self, x, y, absolute_filter_bank_points, relative_kid_positions, only_1_pol_no_comb=False):
        """
        Adds a microstrip feedline connecting the filter bank connection points
        (normally the end of the connection to the phase combiner) to the KIDs.
        There is a bend coming out of both the combiner on the filter bank and
        the KID and in the middle of the microstrip feedline is a taper if the
        two needed linewidths differ. The parameters for this microstrip line
        are defined at the top of this function.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate of the center of the antenna structure,
            i.e. the center of the horn.

        absolute_filter_bank_points : dict
            This dict should have keys: **'TR', 'TL', 'BL', 'BR'** that map to
            [x,y] lists which are the conection point coordinates for the filter
            bank.

        relative_kid_positions : list
            list of [x,y] lists defining the connection point coordinates for each
            of the KIDs. This is the very bottom center point of the inductive
            meander section. This list should be the coordinates, in order, of the
            TopLeft, TopRight, BotLeft, BotRight KIDs.

        KwArgs
        ------
        only_1_pol_no_comb = False
            Default is False, when True this will alter the angles of the
            conections to the filter bank side of the conection to the KIDs.
            This is should be used when the function to add the filter bank,
            add_Filter_bank_and_get_conection_points, has the arguments
            with_combiner=False, only_1_pol=True, with_crossover=False.

        """

        KID_conection_linewidth = 3
        filter_conection_linewidth = 5
        path_line_width = filter_conection_linewidth
        KID_extend_out_distance = 100.0  # the distance over which the feedline is straight coming out of the KID and filter bank conection
        filter_extend_out_distance = 50.0
        flex_path_bend_radius = 50
        kid_side_bend_radius = 50
        filt_side_bend_radius = 50

        extra_width = 400

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

            taper = gdspy.Polygon(taper_points, **self.Nb_Antenna)  # makes the taper as a polygon

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
                **self.Nb_Antenna,
            )
            taper_to_FILT = gdspy.FlexPath(
                points_taper_to_FILT,
                filter_conection_linewidth,
                corners="circular bend",
                bend_radius=filt_side_bend_radius,
                gdsii_path=True,
                **self.Nb_Antenna,
            )

            # adds the taper and the two flex paths for each side to the main cell
            self.Main.add(taper)
            self.Main.add(KID_to_taper)
            self.Main.add(taper_to_FILT)

        return

    def connect_ants_to_Filter_Bank(
        self,
        x,
        y,
        relative_antena_conect_positions,
        antena_rot,
        Main_config_file_dict,
        terminate_ants=[],
        add_dielectric_under_conections=True,
    ):
        """
        Adds a co-planar waveguide feedline that transitons to a microstrip that
        connects the antenna pads to the inner ring of the filter bank. The
        co-planar waveguide is on the antenna pad side and the microstrip is on the
        filter bank side. The parameters of this geometry is within the config.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate of the center of the antennas,
            i.e. the center of the horn.

        relative_antena_conect_positions : list
            This list contains [x,y] lists which are the coordinates of the
            conection points of the tip of the antenna pads. These coords should be
            relative to the center of the antennas. The order of the antenna
            connections should be [bot, top, left, right], where each element is
            an [x, y].

        antena_rot : float, int
            The rotation angle (**in degrees**) of the antennas about thier center.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "antenna_cpw_microstrip_trans" AND "filter_bank" AND "antenna".

        KwArgs
        ------
        terminate_ants : list of strings,
            Default empty list. This is a list of all the antennas that should
            be terminated by an aluminium meander. This should take a list of
            any of the following strings "L", "R", "T", "B". as in left
            e.g. terminate_ants=["L", "R"].

        add_dielectric_under_conections = True
            Default is True which adds a dielectric strip under the connection
            from the antenna to the filter bank.
        """

        config = Main_config_file_dict["antenna_cpw_microstrip_trans"]

        inner_ring_radius = Main_config_file_dict["filter_bank"]["inner_ring_radius"]

        ant_line_width = Main_config_file_dict["antenna"]["top_conect_width"]
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
                ANT[0] + extend_out_distance * cos(self.deg_to_rad(antena_rot + (i * 90))),
                ANT[1] + extend_out_distance * sin(self.deg_to_rad(antena_rot + (i * 90))),
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
                **self.Nb_Antenna,
            )
            self.Main.add(ANT_to_TRANS_path)

            cutout_around_ANT_to_TANS_polys = self.get_polys_from_flexpath(
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
                **self.Nb_Antenna,
            )
            self.Main.add(TRANS_to_FILT_path)

            TRANS_cutouts = gdspy.CellReference(
                transition, tenth_pnt_straight, rotation=(self.rad_to_deg(middle_straight_angle))
            )  # adds the transition geometry between CPW and microstrip antenna a fith of the way down the middle straight
            self.ground_plane_cutouts.add(TRANS_cutouts)

            points = self.rotate_and_move_points_list(
                cpw_transition_poly_points, middle_straight_angle, tenth_pnt_straight[0], tenth_pnt_straight[1]
            )
            cpw_transition_center = gdspy.Polygon(points, **self.Nb_Antenna)
            self.Main.add(cpw_transition_center)

            if add_dielectric_under_conections:
                dielectric_path_points = [ANT, extended_ANT, extended_FILT, FILT]
                dielectric_under_connection_width = 5 * CPW_width
                dielectric_under_connection_path = gdspy.FlexPath(
                    dielectric_path_points,
                    dielectric_under_connection_width,
                    corners="circular bend",
                    bend_radius=flex_path_bend_radius,
                    **self.SiN_dep,
                )
                # self.silicon_nitride_positives.add(dielectric_under_connection_path)
                self.Main.add(dielectric_under_connection_path)

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
                **self.Aluminium,
            )

            arc_bend = gdspy.Round(
                bend_center_xy,
                radius=((meander_bend_height / 2) + (meander_width / 2)),
                inner_radius=((meander_bend_height / 2) - (meander_width / 2)),
                initial_angle=(start_ang - arc_theta_close - pi),
                final_angle=(start_ang - arc_theta_close),
                **self.Aluminium,
            )

            arc_far = gdspy.Round(
                [x, y],
                radius=(inner_ring_radius + meander_bend_height + (meander_width / 2)),
                inner_radius=(inner_ring_radius + meander_bend_height - (meander_width / 2)),
                initial_angle=(start_ang - arc_theta_close + arc_theta_far),
                final_angle=(start_ang - arc_theta_close),
                **self.Aluminium,
            )

            overlap_ang = (overlap_pad_width / 2) / inner_ring_radius
            conect_overlap = gdspy.Round(
                [x, y],
                radius=(inner_ring_radius + (overlap_pad_height / 2)),
                inner_radius=(inner_ring_radius - (overlap_pad_height / 2)),
                initial_angle=(start_ang - overlap_ang),
                final_angle=(start_ang + overlap_ang),
                **self.Aluminium,
            )

            self.Main.add(arc_close)
            self.Main.add(arc_bend)
            self.Main.add(arc_far)
            self.Main.add(conect_overlap)

        return

    def connect_ants_to_KIDs(self, x, y, relative_antena_conect_positions, relative_kid_positions, antena_rot, Main_config_file_dict):
        """
        Adds a co-planar waveguide feedline that transitons to a microstrip that
        connects the antenna pads to the meanders of the KIDs. The co-planar
        waveguide is on the antenna pad side and the microstrip is on the KID
        side. The parameters of this geometry is within the config.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate of the center of the antennas,
            i.e. the center of the horn.

        relative_antena_conect_positions : list
            This list contains [x,y] lists which are the coordinates of the
            conection points of the tip of the antenna pads. These coords should be
            relative to the center of the antennas. The order of the antenna
            connections should be [bot, top, left, right], where each element is
            an [x, y].

        relative_kid_positions : list
            list of [x,y] lists defining the connection point coordinates for each
            of the KIDs. This is the very bottom center point of the inductive
            meander section. This list should be the coordinates, in order, of the
            TopLeft, TopRight, BotLeft, BotRight KIDs.

        antena_rot : float, int
            The rotation angle (**in degrees**) of the antennas about thier center.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "antenna_cpw_microstrip_trans" AND "filter_bank" AND "antenna".

        """

        config = Main_config_file_dict["antenna_cpw_microstrip_trans"]

        inner_ring_radius = Main_config_file_dict["filter_bank"]["inner_ring_radius"]

        ant_line_width = Main_config_file_dict["antenna"]["top_conect_width"]
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

        extend_out_distance = config[
            "extend_out_distance"
        ]  # the distance ove which the feedline is straight coming out of the antenna and Filter Bank
        flex_path_bend_radius = config["flex_path_bend_radius"]

        for i in range(4):
            ANT = ANTS[i]
            KID = KIDS[i]

            extended_ANT = [
                ANT[0] + extend_out_distance * cos(self.deg_to_rad(antena_rot + (i * 90))),
                ANT[1] + extend_out_distance * sin(self.deg_to_rad(antena_rot + (i * 90))),
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
                **self.Nb_Antenna,
            )
            self.Main.add(ANT_to_TRANS_path)

            cutout_around_ANT_to_TANS_polys = self.get_polys_from_flexpath(
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
                [trans_end, extended_KID, KID], microstrip_lw, corners="circular bend", bend_radius=flex_path_bend_radius, **self.Nb_Antenna
            )
            self.Main.add(TRANS_to_KID_path)

            TRANS_cutouts = gdspy.CellReference(
                transition, tenth_pnt_straight, rotation=(self.rad_to_deg(middle_straight_angle))
            )  # adds the transition geometry between CPW and microstrip antenna a fith of the way down the middle straight
            self.ground_plane_cutouts.add(TRANS_cutouts)

            points = self.rotate_and_move_points_list(
                cpw_transition_poly_points, middle_straight_angle, tenth_pnt_straight[0], tenth_pnt_straight[1]
            )
            cpw_transition_center = gdspy.Polygon(points, **self.Nb_Antenna)
            self.Main.add(cpw_transition_center)

        return

    def get_rounded_path_from_passthough_points(self, feedline_pass_through_points, bend_radius, bend_points_gap=10):
        """
        Generates a list of [x,y] points of a path with rounded corners from an
        input list of [x,y] points of a path. Rounds corners with given bend radius.

        The bend is put before the corner point. e.g.
        bending the path [[0,1], [1,1], [1,0]] with bend radius of 1 will form a
        circular arc from [0,1] to [1,0] like the first quadrent of a circle with
        radius 1.

        Parameters
        ----------
        feedline_pass_through_points : list
            list of [x,y] lists which are the coordinates that the feedline passes
            through. **Note**, these points will not be neccessarily be included in
            the final returned rounded feedline.

        bend_radius : float, int
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
        rounded_feedline_points : list
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
        self, feed_points, Main_config_file_dict, center_material="Nb", add_bridges=True, points_to_avoid_bridging=None
    ):
        """
        Adds a Co-Planar Waveguide feedline through the points given along with the
        dielectric layer and adds brdiges across every so many units by default.
        The parameters determining the feeline geometry are the config. The CPW
        consists of the center line, cut out in the ground plane around this center
        line and a dielectric layer covering that.

        Parameters
        ----------
        feed_points : list
            list of [x,y] lists which are the coordinates for the feedline to pass
            through.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "cpw_feedline".

        KwArgs
        ----------
        center_material = "Nb"
            What material to make the center line of the CPW out of. Only takes
            string values of "Nb", "Al" or "Both", where the later will put both
            Niobium and Aluminium.

        add_bridges = True
            Adds bridges across the CPW by deafult, if set to False it will just
            add the feedline.

        points_to_avoid_bridging = None
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

        config = Main_config_file_dict["cpw_feedline"]

        bend_radius = config["bend_radius"]
        feedline_width = config["feedline_width"]
        cutout_around_feedline_width = config["cutout_around_feedline_width"]
        dielectric_under_feedline_width = config["dielectric_under_feedline_width"]

        # makes and adds the center of the cpw flex path in the material specified
        if center_material == "Nb":
            center_path = gdspy.FlexPath(
                feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Nb_Antenna
            )
            self.Main.add(center_path)

        if center_material == "Al":
            center_path = gdspy.FlexPath(
                feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Aluminium
            )
            self.Main.add(center_path)

        if center_material == "Both":
            center_path1 = gdspy.FlexPath(
                feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Aluminium
            )
            center_path2 = gdspy.FlexPath(
                feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Nb_Antenna
            )
            self.Main.add(center_path1)
            self.Main.add(center_path2)

        ground_plane_cuttout_outer_path = gdspy.FlexPath(
            feed_points, cutout_around_feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True
        )  # makes the outer of the cpw flex path
        poly_points_flexpath = self.get_polys_from_flexpath(
            ground_plane_cuttout_outer_path
        )  # gets the polygon points of the outer cpw flex path via function and

        for i in range(len(poly_points_flexpath)):
            port_to_fdln_outer_path_polygon = gdspy.Polygon(poly_points_flexpath[i], **self.Nb_Groundplane)
            self.ground_plane_cutouts.add(
                [port_to_fdln_outer_path_polygon]
            )  # makes a polygon with that so it can be cutout of the ground plane (flex path cannot be added directly as boolean operation cant cut that out)

        dielectric_under_feedline_path = gdspy.FlexPath(
            feed_points, dielectric_under_feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.SiN_dep
        )
        poly_points_flexpath = self.get_polys_from_flexpath(dielectric_under_feedline_path)

        for i in range(len(poly_points_flexpath)):
            port_to_fdln_outer_path_polygon = gdspy.Polygon(poly_points_flexpath[i], **self.SiN_dep)
            self.silicon_nitride_positives.add([port_to_fdln_outer_path_polygon])

        if not add_bridges:
            return

        if points_to_avoid_bridging is not None:
            all_avoid_points = np.array(points_to_avoid_bridging)

        rounded_feedline_points = self.get_rounded_path_from_passthough_points(feed_points, bend_radius)

        self.Main.add(gdspy.FlexPath(rounded_feedline_points, 20, layer=600))

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
                    **self.Nb_Groundplane,
                )
                bridge.rotate(angles_of_path_sections[arg], [x_bridge_pos, y_bridge_pos])
                # self.ground_plane_positives.add(gdspy.Polygon(bridge.polygons[0], **self.Nb_Groundplane)) # TODO
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
                    **self.Nb_Groundplane,
                )
                bridge.rotate(angles_of_path_sections[arg], [x_bridge_pos, y_bridge_pos])
                # self.ground_plane_positives.add(gdspy.Polygon(bridge.polygons[0], **self.Nb_Groundplane)) # TODO
                self.Main.add(bridge)

        return

    def add_arbitrary_CPW(
        self,
        feed_points,
        bend_radius,
        feedline_width,
        cutout_around_feedline_width,
        dielectric_under_feedline_width=0,
        center_material="Nb",
        add_bridges=False,
        cutout_dilectric_around_end_distance=0,
        cutout_center_around_end_distance=0,
        return_outer_poly_points=False,
    ):
        """
        Adds a Co-Planar Waveguide (CPW) feedline through the points given
        along with the dielectric layer. The parameters determining the feeline
        geometry are passed as agurments. The CPW consists of the center line,
        cut out in the ground plane around this center line and a dielectric
        layer covering that.

        Parameters
        ----------
        feed_points : list
            list of [x,y] lists which are the coordinates for the feedline to pass
            through.

        bend_radius : float, int
            The radius for the bends in the CPW to be made.

        feedline_width : float, int
            The width of the center line of the CPW.

        cutout_around_feedline_width : float, int
            The width of the ground plane cutout centered around the
            center line.

        KwArgs
        ------
        dielectric_under_feedline_width = 0:
            This is a float or int and is the width of the dielectric that
            covers the center line and the cutout around that. This is only
            added when the value is non zero.

        center_material = "Nb"
            What material to make the center line of the CPW out of. Only takes
            string values of "Nb", "Al", "Both" or "Nb_Grnd", where the "Both"
            option will put both Niobium and Aluminium.

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

        return_outer_poly_points = False
            Will return the outer polygon points if set to true.


        """

        # makes and adds the center of the cpw flex path in the material specified

        if cutout_center_around_end_distance == 0:
            if center_material == "Nb":
                center_path = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Nb_Antenna
                )
                self.Main.add(center_path)

            if center_material == "Al":
                center_path = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Aluminium
                )
                self.Main.add(center_path)

            if center_material == "Both":
                center_path1 = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Aluminium
                )
                center_path2 = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Nb_Antenna
                )
                self.Main.add(center_path1)
                self.Main.add(center_path2)

            if center_material == "Nb_Grnd":
                center_path = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, **self.Nb_Groundplane
                )
                self.Main.add(center_path)
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
                path_points, feedline_width * 1.5, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, layer=5000
            )
            poly_points_flexpath = self.get_polys_from_flexpath(center_cutout_at_end_of_path)

            for i in range(len(poly_points_flexpath)):
                path_polygon = gdspy.Polygon(poly_points_flexpath[i], layer=5000)
                feedline_center_cell_cutouts.add([path_polygon])

            if center_material == "Nb":
                center_path = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, layer=5000
                )

                new_center_line = gdspy.boolean(center_path, feedline_center_cell_cutouts.get_polygons([5000, 0]), "not", **self.Nb_Antenna)
                self.Main.add(new_center_line)

            if center_material == "Al":
                center_path = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, layer=5000
                )
                new_center_line = gdspy.boolean(center_path, feedline_center_cell_cutouts.get_polygons([5000, 0]), "not", **self.Aluminium)
                self.Main.add(new_center_line)

            if center_material == "Both":
                center_path1 = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, layer=5000
                )
                center_path2 = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, layer=5000
                )
                new_center_line1 = gdspy.boolean(
                    center_path1, feedline_center_cell_cutouts.get_polygons([5000, 0]), "not", **self.Aluminium
                )
                self.Main.add(new_center_line1)
                new_center_line2 = gdspy.boolean(
                    center_path1, feedline_center_cell_cutouts.get_polygons([5000, 0]), "not", **self.Nb_Antenna
                )
                self.Main.add(new_center_line2)

            if center_material == "Nb_Grnd":
                center_path = gdspy.FlexPath(
                    feed_points, feedline_width, corners="circular bend", bend_radius=bend_radius, gdsii_path=True, layer=5000
                )

                new_center_line = gdspy.boolean(
                    center_path, feedline_center_cell_cutouts.get_polygons([5000, 0]), "not", **self.Nb_Groundplane
                )
                self.Main.add(new_center_line)

        angle_of_last_section = np.arctan2((feed_points[-2][1] - feed_points[-1][1]), (feed_points[-2][0] - feed_points[-1][0]))
        ground_plane_cuttout_outer_path_points = feed_points + [
            [
                feed_points[-1][0] + 0.5 * cutout_around_feedline_width * cos(angle_of_last_section + pi),
                feed_points[-1][1] + 0.5 * cutout_around_feedline_width * sin(angle_of_last_section + pi),
            ]
        ]
        ground_plane_cuttout_outer_path = gdspy.FlexPath(
            ground_plane_cuttout_outer_path_points,
            cutout_around_feedline_width,
            corners="circular bend",
            bend_radius=bend_radius,
            gdsii_path=True,
        )  # makes the outer of the cpw flex path
        poly_points_flexpath = self.get_polys_from_flexpath(
            ground_plane_cuttout_outer_path
        )  # gets the polygon points of the outer cpw flex path via function and

        for i in range(len(poly_points_flexpath)):
            port_to_fdln_outer_path_polygon = gdspy.Polygon(poly_points_flexpath[i], **self.Nb_Groundplane)
            self.ground_plane_cutouts.add(
                [port_to_fdln_outer_path_polygon]
            )  # makes a polygon with that so it can be cutout of the ground plane (flex path cannot be added directly as boolean operation cant cut that out)

        exclusion_around_feedline_extra_width = 3000
        exclusion_around_feedline_path = gdspy.FlexPath(
            feed_points,
            cutout_around_feedline_width + exclusion_around_feedline_extra_width,
            corners="circular bend",
            bend_radius=bend_radius,
            gdsii_path=True,
            **self.SiN_dep,
        )
        outer_poly_points = self.get_polys_from_flexpath(exclusion_around_feedline_path)

        # Adding the dielectric under the feedlne.
        if dielectric_under_feedline_width != 0:
            dielectric_under_feedline_path = gdspy.FlexPath(
                feed_points,
                dielectric_under_feedline_width,
                corners="circular bend",
                bend_radius=bend_radius,
                gdsii_path=True,
                **self.SiN_dep,
            )
            poly_points_flexpath = self.get_polys_from_flexpath(dielectric_under_feedline_path)

            for i in range(len(poly_points_flexpath)):
                port_to_fdln_outer_path_polygon = gdspy.Polygon(poly_points_flexpath[i], **self.SiN_dep)
                self.silicon_nitride_positives.add([port_to_fdln_outer_path_polygon])

            exclusion_around_feedline_path = gdspy.FlexPath(
                feed_points,
                dielectric_under_feedline_width + exclusion_around_feedline_extra_width,
                corners="circular bend",
                bend_radius=bend_radius,
                gdsii_path=True,
                **self.SiN_dep,
            )
            outer_poly_points = self.get_polys_from_flexpath(exclusion_around_feedline_path)

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
                    corners="circular bend",
                    bend_radius=bend_radius,
                    gdsii_path=True,
                    **self.SiN_dep,
                )
                poly_points_flexpath = self.get_polys_from_flexpath(dielectric_cutout_at_end_of_path)

                for i in range(len(poly_points_flexpath)):
                    port_to_fdln_outer_path_polygon = gdspy.Polygon(poly_points_flexpath[i], **self.SiN_dep)
                    self.silicon_nitride_cutouts.add([port_to_fdln_outer_path_polygon])

        if not add_bridges:
            if not return_outer_poly_points:
                return
            else:
                return outer_poly_points

        rounded_feedline_points = self.get_rounded_path_from_passthough_points(feed_points, bend_radius)

        self.Main.add(gdspy.FlexPath(rounded_feedline_points, 20, layer=600))

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

        bridge_gap = 5300 / 2
        bridge_width = 6
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
                **self.Nb_Groundplane,
            )
            bridge.rotate(angles_of_path_sections[arg], [x_bridge_pos, y_bridge_pos])
            # self.ground_plane_positives.add(gdspy.Polygon(bridge.polygons[0], **self.Nb_Groundplane)) # TODO
            self.Main.add(bridge)

        if return_outer_poly_points:
            return outer_poly_points
        else:
            return

    def add_port_and_get_connection_point(
        self, x, y, rotation, Main_config_file_dict, center_material="Nb", add_extra_squares=True, dielectric_cutout_in_port=True
    ):
        """
        Adds a port to the mask consisting of a tapered section and a back
        straight. The port will be added where the middle base of the back
        straight section sits at the x, y given.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate for the middle back of the port to be placed.

        rotation : float, int
            The rotation angle (**in radians**) of the port.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "port" AND "cpw_feedline".

        KwArgs
        ------
        center_material="Nb"
            The material to make the center line out of. This can take any of
            the string values {"Nb", "Al"}.

        add_extra_squares=True
            By default adds a series of 3 squares either side of the port in Nb.

        dielectric_cutout_in_port=True
            By default adds a cutout in the dielectric over the center line of
            the port.
        """
        if center_material not in ["Nb", "Al"]:
            raise ValueError("center_material arg should be a str that takes the value of 'Nb' or 'Al'")

        cpw_feedline_config = Main_config_file_dict["cpw_feedline"]
        port_config = Main_config_file_dict["port"]

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

        if center_material == "Nb":
            lay = self.Nb_Antenna["layer"]
            dtype = self.Nb_Antenna["datatype"]

        if center_material == "Al":
            lay = self.Aluminium["layer"]
            dtype = self.Aluminium["datatype"]

        port_center_points = [
            [x, y + (outer_feedline_width / 2)],
            [x + outer_back_length, y + (outer_feedline_width / 2)],
            [x + outer_back_length + taper_length, y + (inner_feedline_width / 2)],
            [x + outer_back_length + taper_length, y - (inner_feedline_width / 2)],
            [x + outer_back_length, y - (outer_feedline_width / 2)],
            [x, y - (outer_feedline_width / 2)],
        ]

        port_center = gdspy.Polygon(port_center_points, layer=lay, datatype=dtype)
        port_center.rotate(rotation, [x, y])
        self.Main.add(port_center)

        port_cutout_around_points = [
            [x, y + (outer_cutout_around_feedline_width / 2)],
            [x + outer_back_length, y + (outer_cutout_around_feedline_width / 2)],
            [x + outer_back_length + taper_length, y + (inner_cutout_around_feedline_width / 2)],
            [x + outer_back_length + taper_length, y - (inner_cutout_around_feedline_width / 2)],
            [x + outer_back_length, y - (outer_cutout_around_feedline_width / 2)],
            [x, y - (outer_cutout_around_feedline_width / 2)],
        ]

        port_cutout_around = gdspy.Polygon(port_cutout_around_points, **self.Nb_Groundplane)
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

        port_dielectric = gdspy.Polygon(port_dielectric_points, **self.SiN_dep)
        port_dielectric.rotate(rotation, [x, y])
        self.silicon_nitride_positives.add(port_dielectric)

        if dielectric_cutout_in_port:
            port_dielectric_cutout_points = [
                [x, y + (dielectric_cutout_in_port_width / 2)],
                [x + dielectric_cutout_in_port_length, y + (dielectric_cutout_in_port_width / 2)],
                [x + dielectric_cutout_in_port_length, y - (dielectric_cutout_in_port_width / 2)],
                [x, y - (dielectric_cutout_in_port_width / 2)],
            ]
            port_dielectric_cutout = gdspy.Polygon(port_dielectric_cutout_points, **self.SiN_dep)
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
                    **self.Nb_Antenna,
                )
                square_top.rotate(rotation, [x, y])
                self.Main.add(square_top)

                square_bot = gdspy.Rectangle(
                    [x + square_offset_from_edge, y - offset_from_center],
                    [x + square_offset_from_edge + square_side_len, y - offset_from_center - square_side_len],
                    **self.Nb_Antenna,
                )
                square_bot.rotate(rotation, [x, y])
                self.Main.add(square_bot)

        feedline_connection_point = self.rotate(x, y, x + outer_back_length + taper_length, y, rotation)

        return feedline_connection_point

    def get_feedline_pass_through_points(self, x, y, relative_kid_positions, Main_config_file_dict):
        """
        This will get a list of coordinates where the feedline should connect to
        the KIDs. This is the end of the coupler.

        Parameters
        ----------
        x, y : float, int
            The x, y coordinate of the center of the antenna structure,
            i.e. the center of the horn.

        relative_kid_positions : list
            list of [x,y] lists defining the connection point coordinates for each
            of the KIDs. This is the very bottom center point of the inductive
            meander section. This list should be the coordinates, in order, of the
            TopLeft, TopRight, BotLeft, BotRight KIDs.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        Returns
        -------
        feed_pass_points : list
            list of [x,y] lists defining the coordinates where to feedline should
            connect to the KIDs. This list has the coonection coordinates for the
            KIDs in order top left, top right, bot left, bot right.

        """
        feed_pass_points = []
        res_config = Main_config_file_dict["resonator"]

        vert_feedline_offset_from_kid_position = (
            (res_config["meander_bot_width"] / 2)
            - (res_config["meander_lw"] / 2)
            + (res_config["meander_left_width_1"] - res_config["meander_left_width_2"])
            + res_config["frame_bot_left_width"]
            - (res_config["frame_left_lw"] / 2)
            + (res_config["coupler_frame_left_lw"] / 2)
            + res_config["left_coupler_frame_to_feed_distance"]
        )

        horizontal_feedline_offset_from_kid_position = (
            res_config["meander_lw"]
            + res_config["meander_left_height_1"]
            + res_config["meander_left_height_2"]
            + res_config["meander_left_height_3"]
            + res_config["frame_left_height"]
            + res_config["coupler_frame_left_height"]
            + res_config["coupler_gap"]
            + (res_config["coupler_lw"] / 2)
        )

        feed_pass_points.append(
            [
                x + relative_kid_positions[0][0] - horizontal_feedline_offset_from_kid_position,
                (y + relative_kid_positions[0][1] + vert_feedline_offset_from_kid_position + 36 / 2),
            ]
        )  # TL

        feed_pass_points.append(
            [
                x + relative_kid_positions[1][0] + horizontal_feedline_offset_from_kid_position,
                (y + relative_kid_positions[1][1] + vert_feedline_offset_from_kid_position + 36 / 2),
            ]
        )  # TR

        feed_pass_points.append(
            [
                x + relative_kid_positions[2][0] - horizontal_feedline_offset_from_kid_position,
                (y + relative_kid_positions[2][1] - vert_feedline_offset_from_kid_position - 36 / 2),
            ]
        )  # BL

        feed_pass_points.append(
            [
                x + relative_kid_positions[3][0] + horizontal_feedline_offset_from_kid_position,
                (y + relative_kid_positions[3][1] - vert_feedline_offset_from_kid_position - 36 / 2),
            ]
        )  # BR

        return feed_pass_points

    def get_total_height_of_resonator(self, Main_config_file_dict):
        """
        This will get the total height of the resonator from the base of the
        inductive meander to the end of the ground plane cutout at the top of the
        structure.

        Parameters
        ----------
        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        Returns
        -------
        KID_height : float
            The total height of the KID calculated from the config file.

        """

        res_config = Main_config_file_dict["resonator"]

        KID_height = (
            res_config["meander_lw"]
            + res_config["meander_left_height_1"]
            + res_config["meander_left_height_2"]
            + res_config["meander_left_height_3"]
            + res_config["frame_left_height"]
            + res_config["coupler_frame_left_height"]
            + res_config["coupler_gap"]
            + res_config["coupler_lw"]
            + res_config["cutout_top_offset"]
        )
        return KID_height

    def get_width_height_of_resonator_IDC_section(self, Main_config_file_dict):
        """
        This will get the total width and height of ground plane cutout around the
        IDC section of the KID.

        Parameters
        ----------
        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        Returns
        -------
        [KID_width, KID_height] : list
            This is a list containing, in order, the KID_width and the KID_height
            of the KID's IDC section calculated from the config file.

        """

        res_config = Main_config_file_dict["resonator"]

        KID_height = (
            res_config["cutout_bot_offset"]
            + res_config["frame_left_height"]
            + res_config["coupler_frame_left_height"]
            + res_config["coupler_gap"]
            + res_config["coupler_lw"]
            + res_config["cutout_top_offset"]
        )

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

        return [KID_width, KID_height]

    def get_total_width_of_resonator(self, Main_config_file_dict):
        """
        This will get the total width of the resonator from the far left of the
        ground plane cutout to the far right of the ground plane cutout of the
        structure. This will always be the widest part of the KID.

        Parameters
        ----------
        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator".

        Returns
        -------
        KID_width : float
            The total width of the KID calculated from the config file.

        """

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

    def get_feedline_center_to_meander_base_distance(self, Main_config_file_dict):
        """
        This will calculate the the vertical distance from the center of the
        feeline to where the base of the KIDs inductive meander will sit assuming
        the coupler arm buts up against the edge of the CPW feedline center strip.
        This distance is calculated based on the dimensions within the config file.

        Parameters
        ----------
        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator" AND "cpw_feedline".

        Returns
        -------
        feedline_to_meander_base_distance : float
            The vertical distance from the CPW feedline center to the middle of the
            base of the KIDs inductive meander calculated from the config file.
            This distance is similar to, but not the same as, the value
            "right_side_width" from within the "get_total_width_of_resonator"
            function but here this also calculates the distance to feedline
            dimension.

        """

        res_config = Main_config_file_dict["resonator"]

        feedline_to_meander_base_distance = (
            (res_config["meander_bot_width"] / 2)
            - (res_config["meander_lw"] / 2)
            + (res_config["meander_left_width_1"] - res_config["meander_left_width_2"])
            + res_config["frame_bot_left_width"]
            - (res_config["frame_left_lw"] / 2)
            + (res_config["coupler_frame_left_lw"] / 2)
            + res_config["left_coupler_frame_to_feed_distance"]
            + (Main_config_file_dict["cpw_feedline"]["feedline_width"] / 2)
        )

        return feedline_to_meander_base_distance

    def get_relative_kid_positions(self, Main_config_file_dict):
        """
        Gets the positions of the base of the KID meanders relative to the center
        of the antennas in the middle. The relative KID positions is returned as a
        list is the order, [top_left, top_right, bot_left, bot_right], where each
        element is an [x, y] list.

        Parameters
        ----------
        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "resonator" AND "general".

        Returns
        -------
        rel_kid_positions : list
            list of [x,y] lists that are the coordinates of the base of the KIDs
            meanders. These KID positons are in order top_left, top_right,
            bot_left, bot_right.

        """

        tot_kid_height = self.get_total_height_of_resonator(Main_config_file_dict)

        ground_gap_between_resonators = Main_config_file_dict["resonator"]["grndpln_gap_between_adjacent_resonators"]

        horizontal_offset_from_center_of_antenna = (
            (Main_config_file_dict["general"]["horizontal_pitch"] / 2) - tot_kid_height - (ground_gap_between_resonators / 2)
        )

        feed_to_meander_base_distance = self.get_feedline_center_to_meander_base_distance(Main_config_file_dict)

        vertical_offset_from_center_of_antenna = (Main_config_file_dict["general"]["vertical_pitch"] / 2) - feed_to_meander_base_distance

        rel_kid_positions = (
            [-horizontal_offset_from_center_of_antenna, +vertical_offset_from_center_of_antenna],
            [+horizontal_offset_from_center_of_antenna, +vertical_offset_from_center_of_antenna],
            [-horizontal_offset_from_center_of_antenna, -vertical_offset_from_center_of_antenna],
            [+horizontal_offset_from_center_of_antenna, -vertical_offset_from_center_of_antenna],
        )

        return rel_kid_positions

    def get_relative_antena_conect_positions(self, rotation, Main_config_file_dict):
        r"""
        Gets the connect positions of the end of the antennas relative to the
        center of the antennas. The relative antenna conect positions is returned
        as a list is the order, [bot, top, left, right], where each element is an
        [x, y] list. If the rotation setting in config is non-zero  this should be
        included in the rotation arg. The relative positions are only rotated by
        this rotation amount.

        Parameters
        ----------
        rotation : float
            The angle (**in radians**) of rotation of the antennas about the center
            of the antennas. A negative rotation angle is clockwise, positive is
            anti-clockwise.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "antenna".

        Return
        ------
        rel_ant_conect_positions : list
            list of [x,y] lists that describe the coordinates for the end of the
            antennas wher they will feed into a KID or filer bank. These antenna
            connection positons are in order bot, top, left, right.
            **Note** These refer to the locations of the end of the antennas before
            any rotation has been done. i.e. if rotated $\pi$ radians the top
            refers to the antenna at the bottom now, left will refer to right
            antenna and so on.

        """

        ant_config = Main_config_file_dict["antenna"]

        ant_connect_offset_from_center = ant_config["distance_from_center"] + ant_config["straight_height"] + ant_config["taper_height"]

        rot = self.deg_to_rad(rotation)

        angles = [(-pi / 2 + rot), (pi / 2 + rot), (pi + rot), (0 + rot)]

        rel_ant_conect_positions = []

        for angle in angles:
            rel_ant_conect_positions.append([ant_connect_offset_from_center * cos(angle), ant_connect_offset_from_center * sin(angle)])

        return rel_ant_conect_positions

    def add_sma_connector_and_launcher_and_get_connection(self, x, y, rot, Main_config_file_dict, bend="none"):
        """
        Creates an SMA connector where the center pin is located at the x,y given
        and returns the point where a feedline should connect.

        Parameters
        ----------
        x, y : float
            The x, y coordinate for the center of the SMA conector pin.

        rot : float
            Angle (**in radians**), the rotation of the whole assembly around the
            center x, y given. positive is anti-clockwise, negative is clockwise.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "sma_connector" AND "cpw_feedline".

        KwArgs
        ------
        bend : str, {"none", "left", "right"}
            If not "none" it will bend the cpw coming out of the sma connector
            around to either the left or right 90 degrees and then calculate the
            taper to the edge of the connector. Anything other than "left" will be regarded as "none"

        Returns
        -------
        conection : list
            The [x, y] connection point at the end of the taper for the SMA
            conector where a feedline should connect.
        """

        sma_config = Main_config_file_dict["sma_connector"]
        feed_config = Main_config_file_dict["cpw_feedline"]

        sma_square_offset_left = sma_config["sma_square_offset_left"]  # 0
        sma_square_offset_right = 7000  # sma_config["sma_square_offset_right"]#5000
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

        # TODO Make this a propper implimentation
        # making the dielectric cutout in the wide section if the launcher
        dielectirc_cutout_length = 4000
        dielectirc_cutout_height = 2500
        cutout_box = gdspy.Rectangle(
            [x + (sma_square_width / 2) + sma_square_offset_right - taper_length, y - (dielectirc_cutout_height / 2)],
            [x + (sma_square_width / 2) + sma_square_offset_right, y + (dielectirc_cutout_height / 2)],
            **self.SiN_dep,
        )
        cutout_box.rotate(rot, [x, y])
        self.silicon_nitride_cutouts.add([cutout_box])
        # self.Main.add(cutout_box)

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
            central_section_poly = gdspy.Polygon(central_section_poly_points, **self.Nb_Antenna)
            central_section_poly.rotate(rot, [x, y])
            self.Main.add(central_section_poly)

            # adding the ground plane cutout around the central path
            grnd_cutout_poly = gdspy.Polygon(grnd_cutout_poly_points, **self.Nb_Groundplane)
            grnd_cutout_poly.rotate(rot, [x, y])
            self.ground_plane_cutouts.add(grnd_cutout_poly)

            # adding the dielectric around the ground plane cutout
            dielectric_poly = gdspy.Polygon(dielectric_poly_points, **self.SiN_dep)
            dielectric_poly.rotate(rot, [x, y])
            self.silicon_nitride_positives.add(dielectric_poly)

            # getting the conection point
            connection_point_xy = self.rotate(x, y, x + (sma_square_width / 2) + sma_square_offset_right, y, rot)

            return connection_point_xy

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
            [x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end, y + (central_linewidth / 2) + gap],
            [x + (sma_square_width / 2) + sma_square_offset_right - taper_length - bend_offset_from_end, y - (central_linewidth / 2) - gap],
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

        central_section_poly = gdspy.Polygon(central_section_poly_points, **self.Nb_Antenna)
        central_section_poly.rotate(rot, [x, y])
        self.Main.add(central_section_poly)

        # adding the ground plane cutout around the central path
        grnd_cutout_poly = gdspy.Polygon(grnd_cutout_poly_points, **self.Nb_Groundplane)
        grnd_cutout_poly.rotate(rot, [x, y])
        self.ground_plane_cutouts.add(grnd_cutout_poly)

        # adding the dielectric around the ground plane cutout
        dielectric_poly = gdspy.Polygon(dielectric_poly_points, **self.SiN_dep)
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
            **self.SiN_dep,
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
            **self.Nb_Groundplane,
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
            **self.Nb_Antenna,
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
        central_taper_poly = gdspy.Polygon(central_taper_poly_points, **self.Nb_Antenna)
        central_taper_poly.rotate(rot, [x, y])
        self.Main.add(central_taper_poly)

        # adding the ground plane cutout taper the central path
        grnd_cutout_taper_poly = gdspy.Polygon(grnd_cutout_taper_poly_points, **self.Nb_Groundplane)
        grnd_cutout_taper_poly.rotate(rot, [x, y])
        self.ground_plane_cutouts.add(grnd_cutout_taper_poly)

        # adding the dielectric taper around the ground plane cutout
        dielectric_taper_poly = gdspy.Polygon(dielectric_taper_poly_points, **self.SiN_dep)
        dielectric_taper_poly.rotate(rot, [x, y])
        self.Main.add(dielectric_taper_poly)

        # getting the conection point
        connection_point_xy = self.rotate(x, y, end_of_round_x, end_of_round_y + sign * taper_length, rot)

        return connection_point_xy

    def get_sma_connector_and_laucher_bounding_box(self, x, y, rot, Main_config_file_dict):
        """
        Gets the outer bounding box around an SMA connector where the center pin is
        located at the x,y given.

        Parameters
        ----------
        x, y : float
            The x, y coordinate for the center of the SMA conector pin.

        rot : float
            Angle (**in radians**), the rotation of the whole assembly around the
            center x, y given. positive is anti-clockwise, negative is clockwise.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "sma_connector".

        Returns
        -------
        bounding_box_points : list
            Returns the bounding box points of the SMA connector. This is a list
            of [x,y] lists which are coordinates for the top left, top right,
            bot right, bot left points in this order.

        """

        sma_config = Main_config_file_dict["sma_connector"]

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

        bounding_box_points = self.rotate_and_move_points_list(bounding_box_points, rot, 0, 0, ox=x, oy=y)

        return bounding_box_points

    def remove_exclusions_from_hex_pack(self, hex_pack, exclusions_list, Main_config_file_dict):
        """
        Removes hex pack points where any part of the KID block would intersect
        with another feature in the exclusions list.

        Parameters
        ----------
        hex_pack : list
            list containing individual [x,y] lists defining the center of all the
            hex pack grid points.

        exclusions_list : list
            list containing one or many lists which each are a list of [x, y] lists
            defining the boundary points for the exclusion polygon shape.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "general".

        Returns
        -------
        new_hex_pack : list
            Returns a list containing [x, y] lists defining the center of the hex
            points (similar in form to the hex_pack argument) that no longer
            clashes with any features from the exclusion list.

        """

        kid_block_horizontal_size = Main_config_file_dict["general"]["horizontal_pitch"]
        kid_block_vertical_size = Main_config_file_dict["general"]["vertical_pitch"]

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

    def remove_exclusions_from_split_hex_pack(self, split_hex_pack, exclusions_list, Main_config_file_dict):
        """
        Removes hex pack points where any part of the KID block would intersect
        with another feature in the exclusions list.

        Parameters
        ----------
        split_hex_pack : list
            list containing two seperate lists of the hex pack grid. Each of those
            lists should be a list of [x, y] lists defining the center of the hex
            point.

        exclusions_list : list
            list containing one or many lists which each are a list of [x, y] lists
            defining the boundary points for the exclusion polygon shape.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "general".

        Returns
        -------
        [new_bot_hex_pack, new_top_hex_pack] : list
            Returns a list (similar in form to the split_hex_pack argument) of two
            lists containing [x, y] lists defining the center of the hex points
            that no longer clashes with any features from the exclusion list.

        """

        kid_block_horizontal_size = Main_config_file_dict["general"]["horizontal_pitch"]
        kid_block_vertical_size = Main_config_file_dict["general"]["vertical_pitch"]

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
        """
        combines a split hex pack grid into one long list instead of two seperate
        lists.

        Parameters
        ----------
        split_hex_pack_grid : list
            list containing two seperate lists of the hex pack grid. Each of those
            lists should be a list of [x, y] lists defining the center of the hex
            point.


        Returns
        -------

        combined_hex_pack_grid : list
            list containing individual [x,y] lists defining the center of all the
            hex pack grid points.

        """
        combined_hex_pack_grid = []

        for xy in split_hex_pack_grid[0]:
            combined_hex_pack_grid.append(xy)

        for xy in split_hex_pack_grid[1]:
            combined_hex_pack_grid.append(xy)

        return combined_hex_pack_grid

    def get_feedline_running_list(self, feed_pass_points, Main_config_file_dict, init_direction="right"):
        """
        This will take in a list of all the points that a feedline should pass
        through and then organizes them into a line running left to right on one
        row then right to left on the next row above and so on untill all points
        have been mapped like this. This can start right to left or left to right
        depending upon the init_direction. This is used to then be able to draw a
        feedline that meanders nicely through the mask passing over all the KIDs.

        Parameters
        ----------
        feed_pass_points : list
            list of [x,y] lists that define the coordinates of all the points the
            feedline should pass through.

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "cpw_feedline".

        KwArgs
        ------
        init_direction : str, {"right", "left"}
            The initial direction to order the fitst row of points. Only takes
            strings "right" or "left". Default is "right" which will take the
            lowest row and start from the left most point and go right.

        Returns
        -------
        running_list : list
            list of [x,y] lists for all the points the feedline should pass through
            ordered in rows running right to left on one row and then left to right
            on the row above and so on untill all the points in the original
            feed_pass_points list have been ordered.

        """

        sorted_list = sorted(feed_pass_points.copy(), key=lambda k: [k[1], k[0]])
        xs_sorted = np.array([a[0] for a in sorted_list])
        ys_sorted = [a[1] for a in sorted_list]
        indexes = [index for index, _ in enumerate(ys_sorted) if ys_sorted[index] != ys_sorted[index - 1]]
        indexes.append(len(ys_sorted))

        running_list = []

        bend_rad = Main_config_file_dict["cpw_feedline"]["bend_radius"]
        extra_straight_length = Main_config_file_dict["cpw_feedline"]["extra_straight_length"]

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
            little_rect = gdspy.Rectangle([xy[0] - 5, xy[1] - 5], [xy[0] + 5, xy[1] + 5], layer=200)
            self.Main.add(little_rect)

        return running_list

    def add_center_pin_and_get_bounding_points(self, center_pin_xy, center_pin_radius):
        """
        Adds a center pin hole cutout to all layers and adds a pin hole positive
        circle to the pin hole positve later. This also returns a list of [x, y]
        points defining the boundary points for the pin.

        Adds a center pin cutout to all layers and adds the pin positive to a pin
        positives layer. Also returns a list of [x, y] points defining the boundary
        points for the pin.

        Parameters
        ----------
        center_pin_xy : list
            list containing the [x,y] coordinate for the center of the pin.

        center_pin_radius : float, int
            Radius of the center pin.

        Returns
        -------
        bounding_points : list
            Returns a list containing [x, y] lists defining the boundary points of
            the center pin.

        """

        ground_plane_cutout_for_center_pin = gdspy.Round(center_pin_xy, center_pin_radius, **self.Nb_Groundplane)
        self.ground_plane_cutouts.add(ground_plane_cutout_for_center_pin)

        silicon_nitride_cutout_for_center_pin = gdspy.Round(center_pin_xy, center_pin_radius, **self.SiN_dep)
        self.silicon_nitride_cutouts.add(silicon_nitride_cutout_for_center_pin)

        silicon_oxide_cutout_for_center_pin = gdspy.Round(center_pin_xy, center_pin_radius, **self.SiO)
        self.silicon_oxide_cutouts.add(silicon_oxide_cutout_for_center_pin)

        silicon_nitride_membrane_cutout_for_center_pin = gdspy.Round(center_pin_xy, center_pin_radius, **self.SiN_Membrane)
        self.silicon_nitride_membrane_cutouts.add(silicon_nitride_membrane_cutout_for_center_pin)

        center_pin_positive = gdspy.Round(center_pin_xy, center_pin_radius, **self.Pin_hole_positives)
        self.Main.add(center_pin_positive)

        bounding_points = []
        for i in range(200):
            bounding_points.append(
                [center_pin_xy[0] + center_pin_radius * cos(i * center_pin_radius), center_pin_xy[1]]
                + center_pin_radius * sin(i * center_pin_radius)
            )

        return bounding_points

    def add_slotted_pin_and_get_bounding_points(self, slotted_pin_center, slotted_pin_length, slotted_pin_radius, center_pin_xy):
        """
        Adds a slotted pin cutout to all layers and adds the pin positive to a pin
        positives layer. Also returns a list of [x, y] points defining the boundary
        points for the pin. The slotted pin will always point towards the center of
        the chip where the center pin hole would be.

        Parameters
        ----------
        slotted_pin_center: list
            list containing the [x, y] coordinate for the center of the slotted
            pin hole.

        slotted_pin_length : float, int
            The maximal length of the slot length.
                i.e. (2*slotted_pin_radius + central straight length).

        slotted_pin_radius : float, int
            Radius of the corners of the slotted pin.

        center_pin_xy : list
            list containing the [x, y] coordinate for the center pin on the chip.
            This is used to ensure the slotted pin is pointing in the direction of
            the center pin hole.

        Returns
        -------
        bounding_points : list
            Returns a list containing [x, y] lists defining the boundary points of
            the slotted pin.

        """

        slotted_pin_angle = np.arctan2((slotted_pin_center[1] - center_pin_xy[1]), (slotted_pin_center[0] - center_pin_xy[0]))

        self.rad_to_deg(slotted_pin_angle)

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

        ang_diff = pi / 100
        loop_angle = slotted_pin_angle + pi / 2
        for i in range(100):
            bounding_points.append(
                [end_points_slot[0][0] + slotted_pin_radius * cos(loop_angle), end_points_slot[0][1] + slotted_pin_radius * sin(loop_angle)]
            )
            loop_angle -= ang_diff

        for i in range(100, 200):
            bounding_points.append(
                [end_points_slot[1][0] + slotted_pin_radius * cos(loop_angle), end_points_slot[1][1] + slotted_pin_radius * sin(loop_angle)]
            )
            loop_angle -= ang_diff

        ground_plane_cutout_for_slotted_pin = gdspy.Polygon(bounding_points, **self.Nb_Groundplane)
        self.ground_plane_cutouts.add(ground_plane_cutout_for_slotted_pin)

        silicon_nitride_cutout_for_slotted_pin = gdspy.Polygon(bounding_points, **self.SiN_dep)
        self.silicon_nitride_cutouts.add(silicon_nitride_cutout_for_slotted_pin)

        silicon_oxide_cutout_for_slotted_pin = gdspy.Polygon(bounding_points, **self.SiO)
        self.silicon_oxide_cutouts.add(silicon_oxide_cutout_for_slotted_pin)

        silicon_nitride_membrane_cutout_for_slotted_pin = gdspy.Polygon(bounding_points, **self.SiN_Membrane)
        self.silicon_nitride_membrane_cutouts.add(silicon_nitride_membrane_cutout_for_slotted_pin)

        slotted_pin_positive = gdspy.Polygon(bounding_points, **self.Pin_hole_positives)
        self.Main.add(slotted_pin_positive)

        return bounding_points

    def add_top_choke_features(self, x, y, Main_config_file_dict):
        """
        Adds the top choke anulus and waveguide hole to the mask at the xy given.
        These are added on seperate layers. The parameters controling the dimesions
        of these are defined in the config.

        Parameters
        ----------
        x, y : float, int
            the x, y position to place the top choke feature

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "top_choke".

        """
        # Getting config
        top_choke_config = Main_config_file_dict["top_choke"]
        waveguide_hole_radius = top_choke_config["waveguide_hole_radius"]  # 1200
        anulus_width = top_choke_config["anulus_width"]  # 600

        aunulus_outer_radius = waveguide_hole_radius + anulus_width
        aunulus_inner_radius = waveguide_hole_radius

        # Making the top choke waveguide hole
        top_choke_waveguide_hole = gdspy.Round([x, y], waveguide_hole_radius, **self.Top_choke_waveguide_hole)
        self.Main.add(top_choke_waveguide_hole)

        # Making the top choke anulus
        top_choke_anulus = gdspy.Round([x, y], aunulus_outer_radius, inner_radius=aunulus_inner_radius, **self.Top_choke_anulus)
        self.Main.add(top_choke_anulus)

        return

    def add_bottom_choke_features(self, x, y, rel_kid_positions, Main_config_file_dict):
        """
        At the xy given, adds the bottom choke waveguide hole and creates IDC
        cutout holes. Also adds pads to support the wafer. These are added on
        seperate layers. The parameters controling the dimesions of these are
        defined in the config.

        Parameters
        ----------
        x, y : float, int
            the x, y position to place the top choke feature

        rel_kid_positions : list
            list of [x, y] lists that describe where the base of each KIDs meander
            is placed relative to the center of the antenna.
            Expected order is top_left, top_right, bot_left, bot_right

        Main_config_file_dict : dict
            dictionary containing individual dictionarys of config settings.
            Requires "bottom_choke".

        """
        # Getting the config
        bottom_choke_config = Main_config_file_dict["bottom_choke"]
        wave_guide_hole_radius = bottom_choke_config["wave_guide_hole_radius"]  # 1200

        pad_radius = bottom_choke_config["pad_radius"]  # 50
        pad_x_offset_from_center = bottom_choke_config["pad_x_offset_from_center"]  # 1200
        pad_y_offset_from_center = bottom_choke_config["pad_y_offset_from_center"]  # 1950

        IDC_cutout_offset_top = bottom_choke_config["IDC_cutout_offset_top"]  # 100
        IDC_cutout_offset_bot = bottom_choke_config["IDC_cutout_offset_bot"]  # 100
        IDC_cutout_offset_left = bottom_choke_config["IDC_cutout_offset_left"]  # 100
        IDC_cutout_offset_right = bottom_choke_config["IDC_cutout_offset_right"]  # 100

        IDC_cutout_width_height = self.get_width_height_of_resonator_IDC_section(Main_config_file_dict)
        tot_KID_height = self.get_total_height_of_resonator(Main_config_file_dict)

        horizontal_offset_from_rel_kid_position_to_center_of_IDC_cutout = tot_KID_height - (IDC_cutout_width_height[1] / 2)

        # Making the bottom choke waveguide hole
        bottom_choke_wave_guide_hole = gdspy.Round([x, y], wave_guide_hole_radius, **self.Bottom_choke_waveguide_hole)
        self.Main.add(bottom_choke_wave_guide_hole)

        # Making the bottom choke pads
        x_sign = [+1, -1, -1, +1]
        y_sign = [+1, +1, -1, -1]
        for i in range(4):
            pad_xy = [x + x_sign[i] * pad_x_offset_from_center, y + y_sign[i] * pad_y_offset_from_center]

            pad_round = gdspy.Round(pad_xy, pad_radius, **self.Bottom_choke_pads)
            self.Main.add(pad_round)

        # Making the bottom choke IDC cutouts
        x_sign = [-1, +1, -1, +1]
        for i in range(4):
            center_of_IDC_cutout_xy = [
                x + rel_kid_positions[i][0] + x_sign[i] * horizontal_offset_from_rel_kid_position_to_center_of_IDC_cutout,
                y + rel_kid_positions[i][1],
            ]

            IDC_cutout = gdspy.Rectangle(
                [
                    center_of_IDC_cutout_xy[0] - (IDC_cutout_width_height[1] / 2) - IDC_cutout_offset_left,
                    center_of_IDC_cutout_xy[1] - (IDC_cutout_width_height[0] / 2) - IDC_cutout_offset_bot,
                ],
                [
                    center_of_IDC_cutout_xy[0] + (IDC_cutout_width_height[1] / 2) + IDC_cutout_offset_right,
                    center_of_IDC_cutout_xy[1] + (IDC_cutout_width_height[0] / 2) + IDC_cutout_offset_top,
                ],
                **self.Bottom_choke_IDC_hole,
            )
            self.Main.add(IDC_cutout)

        return

    def add_hexagonal_tabbed_dicing_line(self, hexagon_points, hex_rad, chip_center_xy, all_sma_conectors_xy_rot):
        """
        Adds a tabbed dicing line around the SOUK boundary

        Parameters
        ----------
        hexagon_points : list
            list of [x,y] lists that define the corner coordinates of the
            hexagonal boundary to create dicing line for.

        hex_rad : float, int
            The radius of the hexagon.

        chip_center_xy : list
            list of [x,y] coordinates of the center of the chip.

        all_sma_conectors_xy_rot : list
            list of [x,y,rot] lists for all the sma connectors on the mask.


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

            dice_line = gdspy.FlexPath(main_dice_line_start_end, diceline_linewidth, offset=diceline_linewidth / 2, **self.Tab_dicing_line)
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

        hex_line_test = gdspy.FlexPath(
            new_outline_points, diceline_linewidth, offset=-diceline_linewidth / 2, ends="flush", **self.Tab_dicing_line
        )
        self.Main.add(hex_line_test)

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
            **self.Tab_dicing_line,
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
                    **self.Tab_dicing_line,
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
                    **self.Tab_dicing_line,
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
                    **self.Tab_dicing_line,
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
                    **self.Tab_dicing_line,
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
                **self.Tab_dicing_line,
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
                **self.Tab_dicing_line,
            )
            tab_left.rotate(side_angle, center_point)
            self.Main.add(tab_left)

        return

    def add_hexagonal_bottom_choke_tabbed_dicing_line(
        self, hexagon_points, hex_rad, chip_center_xy, all_sma_conectors_xy_rot, bend_rad_corner
    ):
        """
        Adds a Bottom choke tabbed dicing line around the SOUK boundary

        Parameters
        ----------
        hexagon_points : list
            list of [x,y] lists that define the corner coordinates of the
            hexagonal boundary to create dicing line for.

        hex_rad : float, int
            The radius of the hexagon.

        chip_center_xy : list
            list of [x,y] coordinates of the center of the chip.

        all_sma_conectors_xy_rot : list
            list of [x,y,rot] lists for all the sma connectors on the mask.

        bend_rad_corner : float, int
            The bend radius for the outer corners

        """

        diceline_linewidth = 300
        diceline_end_extend_length = 500

        diceline_tab_length = 4000
        diceline_tab_gap = 700

        # diceline_tab_gap+=diceline_linewidth

        dice_x_offset_from_center = 8800
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
                main_dice_line_start_end, diceline_linewidth, offset=diceline_linewidth / 2, **self.Bottom_choke_Tab_dicing_line
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

        hex_line_part1_points = new_outline_points[1:11]
        hex_line_part1 = gdspy.FlexPath(
            hex_line_part1_points,
            diceline_linewidth,
            offset=-diceline_linewidth / 2,
            ends="flush",
            corners="natural",
            **self.Bottom_choke_Tab_dicing_line,
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
            **self.Bottom_choke_Tab_dicing_line,
        )
        self.Main.add(hex_line_part2)

        hex_line_part3_points = new_outline_points[12:15]
        hex_line_part3 = gdspy.FlexPath(
            hex_line_part3_points,
            diceline_linewidth,
            offset=-diceline_linewidth / 2,
            ends="flush",
            corners="natural",
            **self.Bottom_choke_Tab_dicing_line,
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
            **self.Bottom_choke_Tab_dicing_line,
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
                    **self.Bottom_choke_Tab_dicing_line,
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
                    **self.Bottom_choke_Tab_dicing_line,
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
                **self.Bottom_choke_Tab_dicing_line,
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
                **self.Bottom_choke_Tab_dicing_line,
            )
            tab_left.rotate(side_angle, center_point)
            self.Main.add(tab_left)

        return

    def add_all_positive_layers_to_mask(self):
        """
        This adds all the positive cells to the mask under new layer names.
        This is the ground_plane, silicon_nitride, silicon_oxide and the
        silicon_nitride_membrane.
        """

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

        silicon_nitride_positive_polygons = self.silicon_nitride_positives.get_polygons([self.SiN_dep["layer"], self.SiN_dep["datatype"]])
        if len(silicon_nitride_positive_polygons) != 0:
            for i in range(len(silicon_nitride_positive_polygons)):
                poly = gdspy.Polygon(silicon_nitride_positive_polygons[i], **positive_layers_dict["silicon_nitride_positives"])
                self.Main.add(poly)
        else:
            print("No silicon_nitride_positives to add to mask.")

        silicon_oxide_positive_polygons = self.silicon_oxide_positives.get_polygons([self.SiO["layer"], self.SiO["datatype"]])
        if len(silicon_oxide_positive_polygons) != 0:
            for i in range(len(silicon_oxide_positive_polygons)):
                poly = gdspy.Polygon(silicon_oxide_positive_polygons[i], **positive_layers_dict["silicon_oxide_positives"])
                self.Main.add(poly)
        else:
            print("No silicon_oxide_positives to add to mask.")

        silicon_nitride_membrane_positive_polygons = self.silicon_nitride_membrane_positives.get_polygons(
            [self.SiN_Membrane["layer"], self.SiN_Membrane["datatype"]]
        )
        if len(silicon_nitride_membrane_positive_polygons) != 0:
            for i in range(len(silicon_nitride_membrane_positive_polygons)):
                poly = gdspy.Polygon(
                    silicon_nitride_membrane_positive_polygons[i], **positive_layers_dict["silicon_nitride_membrane_positives"]
                )
                self.Main.add(poly)
        else:
            print("No silicon_nitride_membrane_positives to add to mask.")

        ground_plane_positive_polygons = self.ground_plane_positives.get_polygons(
            [self.Nb_Groundplane["layer"], self.Nb_Groundplane["datatype"]]
        )
        if len(ground_plane_positive_polygons) != 0:
            for i in range(len(ground_plane_positive_polygons)):
                poly = gdspy.Polygon(ground_plane_positive_polygons[i], **positive_layers_dict["ground_plane_positives"])
                self.Main.add(poly)
        else:
            print("No ground_plane_positives to add to mask.")

        return

    def add_all_cutout_layers_to_mask(self):
        """
        This adds all the cutout cells to the mask under new layer names.
        This is the ground_plane, silicon_nitride, silicon_oxide and the
        silicon_nitride_membrane.
        """

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

        silicon_nitride_cutout_polygons = self.silicon_nitride_cutouts.get_polygons([self.SiN_dep["layer"], self.SiN_dep["datatype"]])
        if len(silicon_nitride_cutout_polygons) != 0:
            for i in range(len(silicon_nitride_cutout_polygons)):
                poly = gdspy.Polygon(silicon_nitride_cutout_polygons[i], **cutout_layers_dict["silicon_nitride_cutouts"])
                self.Main.add(poly)
        else:
            print("No silicon_nitride_cutouts to add to mask.")

        silicon_oxide_cutout_polygons = self.silicon_oxide_cutouts.get_polygons([self.SiO["layer"], self.SiO["datatype"]])
        if len(silicon_oxide_cutout_polygons) != 0:
            for i in range(len(silicon_oxide_cutout_polygons)):
                poly = gdspy.Polygon(silicon_oxide_cutout_polygons[i], **cutout_layers_dict["silicon_oxide_cutouts"])
                self.Main.add(poly)
        else:
            print("No silicon_oxide_cutouts to add to mask.")

        silicon_nitride_membrane_cutout_polygons = self.silicon_nitride_membrane_cutouts.get_polygons(
            [self.SiN_Membrane["layer"], self.SiN_Membrane["datatype"]]
        )
        if len(silicon_nitride_membrane_cutout_polygons) != 0:
            for i in range(len(silicon_nitride_membrane_cutout_polygons)):
                poly = gdspy.Polygon(silicon_nitride_membrane_cutout_polygons[i], **cutout_layers_dict["silicon_nitride_membrane_cutouts"])
                self.Main.add(poly)
        else:
            print("No silicon_nitride_membrane_cutouts to add to mask.")

        ground_plane_cutout_polygons = self.ground_plane_cutouts.get_polygons(
            [self.Nb_Groundplane["layer"], self.Nb_Groundplane["datatype"]]
        )
        if len(ground_plane_cutout_polygons) != 0:
            for i in range(len(ground_plane_cutout_polygons)):
                poly = gdspy.Polygon(ground_plane_cutout_polygons[i], **cutout_layers_dict["ground_plane_cutouts"])
                self.Main.add(poly)
        else:
            print("No ground_plane_cutouts to add to mask.")

        return

    def add_Aluminium_Bulk_around_all_Aluminium(self, region_to_avoid=[]):
        """
        Adds a bulk aluminium block bounding box around all the aluminium in
        the Main cell of the mask. This is such that the aluminium can be etched
        out of the bulk aluminium layer.

        KwArgs
        -------
        regoin_to_avoid : list
            This by default is an empty list. When defined this should be a list
            of [x,y] list points defining a polygon (similar to defining a gdspy
            polygon) that should be enclose a regoin where no bulk aluminium
            should be added. e.g. [[x1,y1], [x2,y2], [x3,y3], [x4,y4]].
            This can be multiple different regions as seperate lists. e.g.
            [ [[x1,y1], [x2,y2] .. [xN,yN]], [[x1,y1], [x2,y2] .. [xN,yN]] ].

        """
        # TODO add a region to only copy and not expand the bulk aluminium added.

        Al_polys = self.Main.get_polygons([self.Aluminium["layer"], self.Aluminium["datatype"]])

        bulk_aluminium_positves = gdspy.Cell("BULK_ALUMINIUM_POSITVES")
        bulk_aluminium_negatives = gdspy.Cell("BULK_ALUMINIUM_NEGATIVES")

        for poly in Al_polys:
            new_poly = gdspy.Polygon(poly, **self.Aluminium_Bulk)
            [BB_xmin, BB_ymin], [BB_xmax, BB_ymax] = new_poly.get_bounding_box()

            y_offset_top = 100
            y_offset_bot = 100
            x_offset_left = 10
            x_offset_right = 10

            oversized_rectangle = gdspy.Rectangle(
                [BB_xmin - x_offset_left, BB_ymin - y_offset_bot], [BB_xmax + x_offset_right, BB_ymax + y_offset_top], **self.Aluminium_Bulk
            )
            bulk_aluminium_positves.add(oversized_rectangle)

        if len(region_to_avoid) != 0:
            if isinstance(region_to_avoid[0][0], list):  # if this element is list then there are multiple regions to avoid.
                for region in region_to_avoid:
                    negative_poly = gdspy.Polygon(region, **self.Aluminium_Bulk)
                    bulk_aluminium_negatives.add(negative_poly)
            elif isinstance(region_to_avoid[0][0], (int, float)):
                negative_poly = gdspy.Polygon(region_to_avoid, **self.Aluminium_Bulk)
                bulk_aluminium_negatives.add(negative_poly)
            else:
                raise TypeError("Incorrect format for region_to_avoid")

        bulk_aluminium = gdspy.boolean(
            bulk_aluminium_positves.get_polygons([self.Aluminium_Bulk["layer"], self.Aluminium_Bulk["datatype"]]),
            bulk_aluminium_negatives,
            "not",
            precision=0.01,
            **self.Aluminium_Bulk,
        )
        self.Main.add(bulk_aluminium)

        return

    def expected_time_progress_bar(self, expected_time, progress_bar_title, steps=100):
        """
        Shows a progress bar based on the expected_time taken to complete an
        operation. This is intended to be used in a daemon thread and get
        alongside the other long running operation, e.g. the
        do_boolean_operations() method used this for each of the longer running
        boolean operation being completed.

        Parameters
        ----------
        expected_time : int
            The expected amount of time (**in seconds**) the operation will
            last and hence how long this progress bar will run for.

        progress_bar_title : str
            The title to put in the progress bar.

        KwArgs
        ------
        steps : int
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

    def boolean_ground_plane(self):
        start_time = time.time()
        groundplane = gdspy.boolean(
            self.ground_plane_positives.get_polygons([self.Nb_Groundplane["layer"], self.Nb_Groundplane["datatype"]]),
            self.ground_plane_cutouts,
            "not",
            precision=0.01,
            **self.Nb_Groundplane,
        )
        self.Main.add(groundplane)  # then adds this groundplane geometry to the main cell

        print("    Ground Plane boolean operation DONE\n")
        end_time = time.time()
        print("    time taken to do ground = " + str(end_time - start_time) + "\n\n\n")
        return

    def boolean_SiN(self):
        start_time = time.time()
        dielectric_layer = gdspy.boolean(
            self.silicon_nitride_positives.get_polygons([self.SiN_dep["layer"], self.SiN_dep["datatype"]]),
            self.silicon_nitride_cutouts,
            "not",
            precision=0.01,
            **self.SiN_dep,
        )  # note the positives need to be a polygon set so you get polygons with tupple (3,0) which is layer 3 datatype 0 i.e. the SiN layer

        if self.silicon_nitride_positives.get_polygons([3, 0]) != []:
            self.Main.add(dielectric_layer)  # then adds this dielectric geometry to the main cell
        print("    dielectirc layer boolean operation DONE\n")
        end_time = time.time()
        print("    time taken to do sin_dep = " + str(end_time - start_time) + "\n\n\n")

        return

    def boolean_SiO(self):
        start_time = time.time()
        SiO_layer = gdspy.boolean(
            self.silicon_oxide_positives.get_polygons([self.SiO["layer"], self.SiO["datatype"]]),
            self.silicon_oxide_cutouts,
            "not",
            precision=0.01,
            **self.SiO,
        )  # note the positives need to be a polygon set so you get polygons with tupple (3,0) which is layer 3 datatype 0 i.e. the SiN layer

        self.Main.add(SiO_layer)  # then adds this dielectric geometry to the main cell
        print("    Silicon DiOxide layer boolean operation DONE\n")
        end_time = time.time()
        print("    time taken to do SiO = " + str(end_time - start_time) + "\n\n\n")
        return

    def boolean_SiN_mem(self):
        start_time = time.time()
        SiN_membrane_layer = gdspy.boolean(
            self.silicon_nitride_membrane_positives.get_polygons([self.SiN_Membrane["layer"], self.SiN_Membrane["datatype"]]),
            self.silicon_nitride_membrane_cutouts,
            "not",
            precision=0.01,
            **self.SiN_Membrane,
        )  # note the positives need to be a polygon set so you get polygons with tupple (3,0) which is layer 3 datatype 0 i.e. the SiN layer

        self.Main.add(SiN_membrane_layer)  # then adds this dielectric geometry to the main cell
        print("    Silicon Nitride membrane layer boolean operation DONE\n")
        end_time = time.time()
        print("    time taken to do SiN mem = " + str(end_time - start_time) + "\n\n\n")

    def do_boolean_operations(self):
        """
        This does all the boolean operations on the the positive cells and
        cutout cells for the mask and adds these to the Main cell ready for
        writing the final gds file. This is the ground_plane, silicon_nitride,
        silicon_oxide, silicon_nitride_membrane.
        """

        # TODO Make an alive loop for each of the steps in this process.
        # time is roughly proportioanl to number of polygons in the positives and cutouts combined.

        print("Doing the boolean operations to the mask.\n")
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

        longest_expected_time = max([ground_plane_expected_time, SiN_mem_expected_time, SiO_expected_time, SiN_mem_expected_time])

        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            bool_ground_done = executor.submit(
                self.boolean_ground_plane
            )  # , desc="Ground Plane Boolean Operation", total=ground_plane_expected_time
            bool_SiN_done = executor.submit(self.boolean_SiN)  # , desc="Silicon Nitride Boolean Operation", total=SiN_expected_time
            bool_SiO_done = executor.submit(self.boolean_SiO)  # , desc="Silicon DiOxide Boolean Operation", total=SiO_expected_time
            Boolean_SiN_mem_done = executor.submit(
                self.boolean_SiN_mem
            )  # , desc="Silicon Nitride membrane Boolean Operation", total=SiN_mem_expected_time

            # executor.submit(self.expected_time_progress_bar, longest_expected_time, "Boolean Operations")
            # executor.submit(self.expected_time_progress_bar, ground_plane_expected_time, "Ground Plane Boolean Operation")
            # executor.submit(self.expected_time_progress_bar, SiN_expected_time, "Silicon Nitride Boolean Operation")
            # executor.submit(self.expected_time_progress_bar, SiO_expected_time, "Silicon DiOxide Boolean Operation")
            # executor.submit(self.expected_time_progress_bar, SiN_mem_expected_time, "Silicon Nitride membrane Boolean Operation")

        # with concurrent.futures.ProcessPoolExecutor() as executor:
        #     bool_ground_done = executor.submit(self.boolean_ground_plane) #, desc="Ground Plane Boolean Operation", total=ground_plane_expected_time
        #     bool_SiN_done = executor.submit(self.boolean_SiN) #, desc="Silicon Nitride Boolean Operation", total=SiN_expected_time
        #     bool_SiO_done = executor.submit(self.boolean_SiO) #, desc="Silicon DiOxide Boolean Operation", total=SiO_expected_time
        #     Boolean_SiN_mem_done = executor.submit(self.boolean_SiN_mem) #, desc="Silicon Nitride membrane Boolean Operation", total=SiN_mem_expected_time
        #
        #     executor.submit(self.expected_time_progress_bar, ground_plane_expected_time, "Ground Plane Boolean Operation")
        #     executor.submit(self.expected_time_progress_bar, SiN_expected_time, "Silicon Nitride Boolean Operation")
        #     executor.submit(self.expected_time_progress_bar, SiO_expected_time, "Silicon DiOxide Boolean Operation")
        #     executor.submit(self.expected_time_progress_bar, SiN_mem_expected_time, "Silicon Nitride membrane Boolean Operation")

        # with tqdm(total=max([ground_plane_expected_time, SiN_mem_expected_time, SiO_expected_time, SiN_mem_expected_time])) as pbar:
        #     bool_done = {}
        #     with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        #         bool_done["ground"] = executor.submit(self.boolean_ground_plane) #, desc="Ground Plane Boolean Operation", total=ground_plane_expected_time
        #         bool_done["SiN"] = executor.submit(self.boolean_SiN) #, desc="Silicon Nitride Boolean Operation", total=SiN_expected_time
        #         bool_done["SiO"] = executor.submit(self.boolean_SiO) #, desc="Silicon DiOxide Boolean Operation", total=SiO_expected_time
        #         bool_done["SiN_mem"] = executor.submit(self.boolean_SiN_mem) #, desc="Silicon Nitride membrane Boolean Operation", total=SiN_mem_expected_time
        #         # Maybe this does not work assigning to dict????
        #         for key in bool_done:
        #             if not bool_done[key]:
        #                 pbar.update(very_first_time - time.time())
        #     pbar.close()
        # start_time = time.time()
        # groundplane = gdspy.boolean(self.ground_plane_positives.get_polygons([self.Nb_Groundplane["layer"],self.Nb_Groundplane["datatype"]]),
        #                             self.ground_plane_cutouts,
        #                             "not", precision=0.01, **self.Nb_Groundplane)
        # self.Main.add(groundplane) # then adds this groundplane geometry to the main cell

        # if ground_progress_thread.is_alive():
        #     ground_progress_thread.stop()
        # print("    Ground Plane boolean operation DONE\n")
        # end_time = time.time()
        # print("    time taken to do ground = " + str(end_time - start_time) + "\n\n\n")

        # SiNDep
        # combined_number_of_SiN_dep_positives_cutouts = 0
        # for li in self.silicon_nitride_positives.get_polygons():
        #     combined_number_of_SiN_dep_positives_cutouts += len(li)
        #
        # for li in self.silicon_nitride_cutouts.get_polygons():
        #     combined_number_of_SiN_dep_positives_cutouts += len(li)
        #
        # expected_time = combined_number_of_SiN_dep_positives_cutouts*2.955534e-09
        # SiN_progress_thread = threading.Thread(target=self.expected_time_progress_bar, args=(int(expected_time), "dielectric Boolean Operation"), daemon=True)
        # SiN_progress_thread.run()

        # start_time = time.time()
        # dielectric_layer = gdspy.boolean(self.silicon_nitride_positives.get_polygons([self.SiN_dep["layer"],self.SiN_dep["datatype"]]),
        #                                  self.silicon_nitride_cutouts,
        #                                  "not", precision=0.01, **self.SiN_dep) # note the positives need to be a polygon set so you get polygons with tupple (3,0) which is layer 3 datatype 0 i.e. the SiN layer
        #
        # if self.silicon_nitride_positives.get_polygons([3,0]) != []:
        #     self.Main.add(dielectric_layer) # then adds this dielectric geometry to the main cell
        #
        # if SiN_progress_thread.is_alive():
        #     SiN_progress_thread.stop()
        # print("    dielectirc layer boolean operation DONE\n")
        # end_time = time.time()
        # print("    time taken to do sin_dep = " + str(end_time - start_time) + "\n\n\n")

        # # SiOLayer
        #     combined_number_of_SiO_positives_cutouts = 0
        #     for li in self.silicon_oxide_positives.get_polygons():
        #         combined_number_of_SiO_positives_cutouts += len(li)
        #
        #     for li in self.silicon_oxide_cutouts.get_polygons():
        #         combined_number_of_SiO_positives_cutouts += len(li)
        #
        #     expected_time = combined_number_of_SiO_positives_cutouts*2.73806e-05
        #     SiO_progress_thread = threading.Thread(target=self.expected_time_progress_bar, args=(int(expected_time), "Silicon DiOxide Boolean Operation"), daemon=True)
        #     SiO_progress_thread.run()
        #
        # start_time = time.time()
        # SiO_layer = gdspy.boolean(self.silicon_oxide_positives.get_polygons([self.SiO["layer"],self.SiO["datatype"]]),
        #                           self.silicon_oxide_cutouts,
        #                           "not", precision=0.01, **self.SiO) # note the positives need to be a polygon set so you get polygons with tupple (3,0) which is layer 3 datatype 0 i.e. the SiN layer
        #
        # self.Main.add(SiO_layer) # then adds this dielectric geometry to the main cell
        # if SiO_progress_thread.is_alive():
        #     SiO_progress_thread.stop()
        # print("    Silicon DiOxide layer boolean operation DONE\n")
        # end_time = time.time()
        # print("    time taken to do SiO = " + str(end_time - start_time) + "\n\n\n")

        # # SiNmembrane
        #     combined_number_of_SiN_mem_positives_cutouts = 0
        #     for li in self.silicon_nitride_positives.get_polygons():
        #         combined_number_of_SiN_mem_positives_cutouts += len(li)
        #
        #     for li in self.silicon_oxide_cutouts.get_polygons():
        #         combined_number_of_SiN_mem_positives_cutouts += len(li)
        #
        #     expected_time = combined_number_of_SiN_mem_positives_cutouts*1.92823e-05
        #     SiN_mem_progress_thread = threading.Thread(target=self.expected_time_progress_bar, args=(int(expected_time), "SiN membrane Boolean Operation"), daemon=True)
        #     SiN_mem_progress_thread.run()
        #
        # start_time = time.time()
        # SiN_membrane_layer = gdspy.boolean(self.silicon_nitride_membrane_positives.get_polygons([self.SiN_Membrane["layer"],self.SiN_Membrane["datatype"]]),
        #                                    self.silicon_nitride_membrane_cutouts,
        #                                    "not", precision=0.01, **self.SiN_Membrane) # note the positives need to be a polygon set so you get polygons with tupple (3,0) which is layer 3 datatype 0 i.e. the SiN layer
        #
        # self.Main.add(SiN_membrane_layer) # then adds this dielectric geometry to the main cell
        # if SiN_mem_progress_thread.is_alive():
        #     SiN_mem_progress_thread.stop()
        # print("    Silicon Nitride membrane layer boolean operation DONE\n")
        # end_time = time.time()
        # print("    time taken to do SiN mem = " + str(end_time - start_time) + "\n\n\n")

        very_last_time = time.time()
        print("    The whole time taken overall = " + str(very_last_time - very_first_time) + "\n\n\n")

        return

    def append_layer_to_xml(self, layer_prop_dict):
        """
        Appends a layer property element to the xml layer properties object.

        Parameters
        ----------
        layer_prop_dict : dict
            This is a dictionary with layer property keys and values. These are
            then directly added as xml tags and values respectively.
        """
        prop = ET.SubElement(self.xml_layer_prop, "properties")

        for key, val in layer_prop_dict.items():
            layer_element = ET.SubElement(prop, key)
            layer_element.text = val

        return

    def save_layer_props(self, filename):
        """
        Save a layer properties file for the mask under the name given.
        This should be done after do_boolean_operations() has been run to
        ensure that all the layers that should exist are included within the
        generated layer properies file.

        Parameters
        ----------
        filename : str
            This is the name for the layer_properties file. If this does not
            include a ".lyp" file extention it will be automatically added.

        Output
        ------
        xml_layer_prop .lyp file.
        """

        from ._default_layer_props_info import default_layer_colors_dict
        from .colors import all_colors

        local_default_layer_colors_dict = copy.deepcopy(default_layer_colors_dict)

        self.xml_layer_prop = ET.Element("layer-properties")  # make the root of the layer props file

        layer_numbers_in_mask = list(self.Main.get_layers())
        for i in range(len(layer_numbers_in_mask)):
            layer_numbers_in_mask[i] = layer_numbers_in_mask[i] % (2**16)  # mod 65536 so layer nums are not negative.

        layer_numbers_in_lookup = list(self.all_layers_name_lookup_from_number.keys())

        for layer_number in layer_numbers_in_mask:
            if layer_number in layer_numbers_in_lookup:
                layer_name = self.all_layers_name_lookup_from_number[layer_number]

                layer_prop_dict = local_default_layer_colors_dict[layer_name]

                if layer_prop_dict["frame-color"][0] != "#":
                    layer_prop_dict["frame-color"] = all_colors[layer_prop_dict["frame-color"]]

                if layer_prop_dict["fill-color"][0] != "#":
                    layer_prop_dict["fill-color"] = all_colors[layer_prop_dict["fill-color"]]

                base_layer_name = layer_prop_dict["name"]
                layer_prop_dict["name"] = f"[{layer_number}] - {base_layer_name}"
                layer_prop_dict["source"] = f"{layer_number}/{0}@1"  # TODO 0 should be the layer_datatype

                self.append_layer_to_xml(layer_prop_dict)

            else:
                # set to some arb value.
                layer_prop_dict = local_default_layer_colors_dict["NoName"]

                color_for_no_name_layer = all_colors[list(all_colors.keys())[(layer_number * 10) % len(all_colors)]]
                layer_prop_dict["frame-color"] = color_for_no_name_layer
                layer_prop_dict["fill-color"] = color_for_no_name_layer

                layer_prop_dict["name"] = f"[{layer_number}] - NoName"
                layer_prop_dict["source"] = f"{layer_number}/{0}@1"  # TODO 0 should be the layer_datatype

                self.append_layer_to_xml(layer_prop_dict)

        tree = ET.ElementTree(self.xml_layer_prop)
        if filename[-4:] != ".lyp":
            filename += ".lyp"
        tree.write(filename, encoding="utf-8", xml_declaration=True)

        print(f'Succesfully written layerprops file "{filename}"\n')

        return

    def save_gds_file(
        self,
        filename,
        try_to_make_backside=False,
        make_mirrored_along_x=False,
        mirrored_x_filename="",
        make_mirrored_along_y=False,
        mirrored_y_filename="",
    ):
        """
        Save and write the mask. This saves everything in the Main cell to a
        gds file under the given filename. Additionally can save a mirrored
        version. If Anything exists in the MainBackside cell. this will also be
        saved with "_BACKSIDE" appended to the end.

        Parameters
        ----------
        filename : str
            This is the name for the gds file that will be written. If this
            filename does not include a ".gds" file extention it will be
            automatically added.

        KwArgs
        ------
        try_to_make_backside : Boolean
            Default False, This will try to make the backside of the mask file
            if anything exists in the MainBackside cell. If nothing exists in
            this cell, Nothing will get made.

        make_mirrored_along_x, make_mirrored_along_y : Boolean
            Whether or not a mirrored version of the mask (mirrored along the x
            or y axis respectively) should also be saved.
            See mirrored_filename for naming this.

        mirrored_x_filename, mirrored_x_filename : str
            This is the name for the mirrored (in the x or y axis respectively)
            version of the mask. Like the filename arg, if any of these do not
            include a ".gds" file extention it will be automatically added. If
            this argument is not provided the mirrored version of the mask will
            take same name as the filename provided with "_MIRRORED_ACROSS_{X/Y}"
            ,depending upon the mirror axis, appended at the end.
        """
        from phidl import Device

        if filename[-4:] != ".gds":
            filename += ".gds"

        print(f'Writing to GDS "{filename}"')
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

            print(f'Writing Mask mirrored along x-axis to GDS "{mirrored_x_filename}"')
            Main_mirrored_across_x_axis = self.Main.copy("MainCell_xax_mirror", deep_copy=True, x_reflection=True)

            MainDevice_mirrored_x = Device("MainDevice_mirrored_x")
            MainDevice_mirrored_x.add(Main_mirrored_across_x_axis)
            MainDevice_mirrored_x.write_gds(mirrored_x_filename)

        if make_mirrored_along_y:
            if mirrored_y_filename == "":
                mirrored_y_filename = filename[:-4] + "_MIRRORED_ACROSS_Y.gds"
            elif mirrored_y_filename[-4:] != ".gds":
                mirrored_y_filename += ".gds"

            print(f'Writing Mask mirrored along y-axis to GDS "{mirrored_y_filename}"')
            Main_mirrored_across_x_axis_rot90 = self.Main.copy(
                "MainCell_xax_mirror_rot90", deep_copy=True, rotation=(pi / 2), x_reflection=True
            )
            Main_mirrored_across_y_axis = Main_mirrored_across_x_axis_rot90.copy("MainCell_yax_mirror", deep_copy=True, rotation=(pi / 2))

            MainDevice_mirrored_y = Device("MainDevice_mirrored_y")
            MainDevice_mirrored_y.add(Main_mirrored_across_y_axis)
            MainDevice_mirrored_y.write_gds(mirrored_y_filename)

        print("Finished writing files to disk.")
        return
