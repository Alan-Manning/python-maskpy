import os
from typing import Literal

import gdspy
from numpy import pi as pi

from ..layers import Layer
from ..logging import TextColor, styled_text


def _validate_layer_map(layer_map: None | dict[tuple[int, int] | list[int], Layer]) -> None:
    """Validate the layer map passed in is in the correct form and has correct
    types."""

    if not isinstance(layer_map, dict):
        msg = f"layer_map should be a dictionary with keys that should be of type tuple[int, int] | list[int] (len(2)), and values that should be of type Layer. Recieved type {type(layer_map)}."
        raise TypeError(msg)

    for key, val in layer_map.items():
        if not isinstance(key, (tuple, list)):
            msg = f"dict key should be tuple or list. Got {key} of type {type(key)}."
            raise TypeError(msg)
        if not len(key) == 2:
            msg = f"tuple key should be of len 2. Got {key}."
            raise TypeError(msg)
        if not all(isinstance(x, int) for x in key):
            msg = f"dict key should be of type tuple[int, int]. Got tuple[{type(key[0])}, {type(key[1])}]."
            raise TypeError(msg)
        if not isinstance(val, Layer):
            msg = f"layer should be of type Layer. Got {val} of type {type(val)}."
            raise TypeError(msg)
    return


def add_gds_file_to_mask(
    mask_builder,
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
    mask_builder: SoukMaskBuilder
        mask_builder object.

    file_name: str,
        The filename for the gds mask file. If this does not include the
        `.gds` file extention then it will be automatically inserted.

    xy : list[float | int] | tuple[float | int, float | int]
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
    _validate_layer_map(layer_map)

    # Append the .gds file ext if it doesn't already exist.
    if not file_name.endswith(".gds"):
        file_name += ".gds"

    if file_path is None:
        file_path = os.getcwd()

    full_file = os.path.join(file_path, file_name)

    # Check file exists
    if not os.path.isfile(full_file):
        raise FileNotFoundError(styled_text(f"Warning: Could Not find the file: '{full_file}'.", color=TextColor.ERROR))

    library = gdspy.GdsLibrary(infile=full_file)
    top_level_cell: gdspy.Cell = library.top_level()[0]

    polygon_dict = top_level_cell.get_polygons(by_spec=True)

    if placement == "center":
        dx = xy[0]
        dy = xy[1]

    else:
        bbox = top_level_cell.get_bounding_box()
        if bbox is None:
            print(
                f"Could not determine the bounding box for {file_name} so could not perform `{placement}` placement. Defaulting to `center` placement."
            )
            dx = xy[0]
            dy = xy[1]
        else:
            ((bb_x_min, bb_y_min), (bb_x_max, bb_y_max)) = bbox
            width = bb_x_max - bb_x_min
            height = bb_y_max - bb_y_min
            match placement:
                case "bot_left":
                    dx = xy[0] + (width / 2)
                    dy = xy[1] + (height / 2)
                case "bot_center":
                    dx = xy[0]
                    dy = xy[1] + (height / 2)
                case "bot_right":
                    dx = xy[0] - (width / 2)
                    dy = xy[1] + (height / 2)
                case "center_left":
                    dx = xy[0] + (width / 2)
                    dy = xy[1]
                case "center_right":
                    dx = xy[0] - (width / 2)
                    dy = xy[1]
                case "top_left":
                    dx = xy[0] + (width / 2)
                    dy = xy[1] - (height / 2)
                case "top_center":
                    dx = xy[0]
                    dy = xy[1] - (height / 2)
                case "top_right":
                    dx = xy[0] - (width / 2)
                    dy = xy[1] - (height / 2)
                case _:
                    dx = xy[0]
                    dy = xy[1]

    for from_layer, to_layer in layer_map.items():
        if not isinstance(from_layer, tuple):
            raise TypeError(f"layer_map keys should be of type tuple[int, int]. Got key {from_layer}.")
        elif not isinstance(from_layer[0], int):
            raise TypeError(f"layer_map keys should be of type tuple[int, int]. Got key {from_layer}.")
        elif not isinstance(from_layer[1], int):
            raise TypeError(f"layer_map keys should be of type tuple[int, int]. Got key {from_layer}.")
        elif not len(from_layer) == 2:
            raise TypeError(f"layer_map keys should be of type tuple[int, int]. Got key {from_layer}.")

        polygons_in_layer = polygon_dict.get(from_layer)
        if polygons_in_layer is None:
            raise LookupError(f"Could not find layer {from_layer} in mask {file_name}.")

        for polygon in polygons_in_layer:
            poly = gdspy.Polygon(
                polygon,
                layer=to_layer.number,
                datatype=to_layer.datatype,
            )
            poly.translate(dx, dy)
            mask_builder.Main.add(poly)

    print(f"Succesfully added file `{file_name}` to mask at {xy}.")
    return
