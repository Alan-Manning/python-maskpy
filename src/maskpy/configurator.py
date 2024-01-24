import os
import sys

import gdspy
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from numpy import pi
from PyQt6.QtWidgets import QApplication, QVBoxLayout, QWidget

from maskpy._default_layer_props_info import default_layer_colors_dict
from maskpy.colors import all_colors

# from . import _default_layer_props_info, colors
from maskpy.SOUK_FUNCS_V2 import SOUK_Functions
from maskpy.V8Mux import V8_get_IDCL_and_CCL_from_f0


class MyWindow(QWidget):
    # def __init__(self, parent: typing.Optional['QWidget'] = ..., flags: QtCore.Qt.WindowType = ...) -> None:
    #     super().__init__(parent, flags)
    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)

        canvas = FigureCanvasQTAgg(fig)
        layout.addWidget(canvas)


"""
###############################################################################
Setup
###############################################################################
"""
souk = SOUK_Functions()
V8_IDCL_CCL_mux_func = V8_get_IDCL_and_CCL_from_f0
chip_center_xy = [0, 0]
"""
###############################################################################
Config
###############################################################################
"""
config_filename = r"C:\Users\Alan_\Documents\Year4\Summer Work\Code\Making Masks\6 inch reduced test mask\6Inch_reduced_test_mask_config_excel_V12_base.xlsx"
# config_filepath = r"C:\\Users\\Alan_\\Documents\\Year4\\Summer Work\\Code\\Making Masks\\6 inch reduced test mask"
config_filepath = ""
config_sheets = [
    "general",
    "test_chip_quad",
    "port",
    "top_choke",
    "bottom_choke",
    "sma_connector",
    "cpw_feedline",
    "antenna_cpw_microstrip_trans",
    "antenna",
    "resonator",
    "filter_bank",
    "filter_bank_ring_overlap",
    "Lo_pass_filters",
    "Hi_pass_filters",
    "combiner_section_90GHZ",
    "combiner_section_150GHZ",
]
# Creates the congig file as a dictionary
Main_config_file = {}
for sheet_name in config_sheets:
    Main_config_file[sheet_name] = souk.get_config_excel(config_filename, sheet_name, path=config_filepath)
"""
###############################################################################
making the wafer
###############################################################################
"""
wafer_square_width = 3000
wafer_square_height = 3000

wafer_grnd_square = gdspy.Rectangle(
    [chip_center_xy[0] - (wafer_square_width / 2), chip_center_xy[1] - (wafer_square_height / 2)],
    [chip_center_xy[0] + (wafer_square_width / 2), chip_center_xy[1] + (wafer_square_height / 2)],
    layer=souk.Nb_Groundplane["layer"],
)
souk.ground_plane_positives.add(wafer_grnd_square)

wafer_SiO_square = gdspy.Rectangle(
    [chip_center_xy[0] - (wafer_square_width / 2), chip_center_xy[1] - (wafer_square_height / 2)],
    [chip_center_xy[0] + (wafer_square_width / 2), chip_center_xy[1] + (wafer_square_height / 2)],
    layer=souk.SiO["layer"],
)
souk.silicon_oxide_positives.add(wafer_grnd_square)
wafer_SiN_mem_square = gdspy.Rectangle(
    [chip_center_xy[0] - (wafer_square_width / 2), chip_center_xy[1] - (wafer_square_height / 2)],
    [chip_center_xy[0] + (wafer_square_width / 2), chip_center_xy[1] + (wafer_square_height / 2)],
    layer=souk.SiN_Membrane["layer"],
)
souk.silicon_nitride_membrane_positives.add(wafer_grnd_square)

test_cutout = gdspy.Rectangle(
    [chip_center_xy[0] - (100 / 2), chip_center_xy[1] - (100 / 2)],
    [chip_center_xy[0] + (100 / 2), chip_center_xy[1] + (100 / 2)],
    layer=souk.Nb_Groundplane,
)
souk.ground_plane_cutouts.add(test_cutout)

souk.add_kid(0, 0, 0, 3.0e9, V8_IDCL_CCL_mux_func, Main_config_file)

souk.do_boolean_operations()
"""
###############################################################################
Plotting the result
###############################################################################
"""
fig = Figure()
ax = fig.add_subplot()

# default_layer_colors_dict = _default_layer_props_info.default_layer_colors_dict
# all_colors = colors.all_colors

layer_numbers_in_mask = list(souk.Main.get_layers())
for i, lay_no in enumerate(layer_numbers_in_mask):
    # mod 65536 so layer nums are not negative.
    layer_numbers_in_mask[i] = lay_no % (2**16)

layer_numbers_in_lookup = list(souk.all_layers_name_lookup_from_number.keys())

for current_layer_number in layer_numbers_in_mask:
    current_layer_name = souk.all_layers_name_lookup_from_number[current_layer_number]
    polys_in_layer = souk.Main.get_polygons([current_layer_number, 0])
    frame_color_for_layer = default_layer_colors_dict[current_layer_name]["frame-color"]
    fill_color_for_layer = default_layer_colors_dict[current_layer_name]["fill-color"]
    fill_alpha_for_layer = 0.2
    name_for_layer = default_layer_colors_dict[current_layer_name]["name"]

    ax.fill([], [], color=fill_color_for_layer, edgecolor=frame_color_for_layer, alpha=fill_alpha_for_layer, label=name_for_layer)

    if frame_color_for_layer[0] != "#":
        frame_color_for_layer = all_colors[frame_color_for_layer]

    if fill_color_for_layer[0] != "#":
        fill_color_for_layer = all_colors[fill_color_for_layer]

    for poly_points in polys_in_layer:
        coords = list(poly_points)
        coords.append(list(poly_points[0]))
        xs, ys = zip(*coords)
        ax.fill(
            xs, ys, facecolor=fill_color_for_layer, edgecolor=frame_color_for_layer, alpha=fill_alpha_for_layer, zorder=current_layer_number
        )


app = QApplication(sys.argv)
window = MyWindow()
window.show()
sys.exit(app.exec())
