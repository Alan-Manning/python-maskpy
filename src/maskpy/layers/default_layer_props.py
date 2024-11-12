import os
from dataclasses import asdict, dataclass
from enum import StrEnum
from typing import TypedDict
from xml.etree.ElementTree import Element, ElementTree, SubElement

import yaml


class LayerPropType(StrEnum):
    """Type of layer properties

    Enum Members
    ------------
    ANL
    SOUK

    """

    ANL = "anl"
    SOUK = "souk"


@dataclass
class SingleLayerProperty:
    """Single layer properties data class.

    Parameters
    ----------
    frame_color: str,
        string of Either color name from colors.all_colors or hex "#rrggbb".

    fill_color: str
        string of Either color name from colors.all_colors or hex "#rrggbb".

    frame_brightness: str
    fill_brightness: str
    dither_pattern: str
    line_style: str
    valid: str
    visible: str
    transparent: str
    width: str
    marked: str
    xfill: str
    animation: str
    name: str
    source: str
    """

    frame_color: str
    fill_color: str
    frame_brightness: str
    fill_brightness: str
    dither_pattern: str
    line_style: str
    valid: str
    visible: str
    transparent: str
    width: str
    marked: str
    xfill: str
    animation: str
    name: str
    source: str

    def add_to_xml_layer_properties(self, layer_properties_xml: Element) -> Element:
        """Convert the SingleLayerProperty object to an xml layer property for
        use in klayout .lyp file.

        Returns
        -------
        xml_layer_property: ElementTree.Element
        """
        single_xml_layer_property = SubElement(layer_properties_xml, "properties")

        for key, val in asdict(self).items():
            layer_element = SubElement(single_xml_layer_property, key.replace("_", "-"))
            layer_element.text = val

        return layer_properties_xml


def get_default_layer_properties(layer_prop_type: LayerPropType) -> dict[str, SingleLayerProperty]:
    """Get the default layer properties for the specified layer_prop_type.

    Parameters
    ----------
    layer_prop_type: LayerPropType
        This is the type to get the default for. Takes any member of LayerPropType.
    """
    if not isinstance(layer_prop_type, LayerPropType):
        raise TypeError(f"layer_prop_type should be of type LayerPropType, not {type(layer_prop_type)}.")

    match layer_prop_type:
        case LayerPropType.ANL:
            default_layer_properties_filename = "default_ANL_layer_properties"
        case LayerPropType.SOUK:
            default_layer_properties_filename = "default_SOUK_layer_properties"
        case _:
            raise ValueError(f"Unable to locate defuatl layer properties for {layer_prop_type}.")

    try:
        file_path = os.path.dirname(os.path.realpath(__file__))
        full_file = os.path.join(file_path, f"default_props/{default_layer_properties_filename}.yaml")
        with open(full_file, "r") as f:
            default_layer_props_dict: dict[str, dict[str, str]] = yaml.safe_load(f)
    except Exception as e:
        raise e

    default_layer_props: dict[str, SingleLayerProperty] = {}

    for key, val in default_layer_props_dict.items():
        single_layer_prop = SingleLayerProperty(
            frame_color=val["frame-color"],
            fill_color=val["fill-color"],
            frame_brightness=val["frame-brightness"],
            fill_brightness=val["fill-brightness"],
            dither_pattern=val["dither-pattern"],
            line_style=val["line-style"],
            valid=val["valid"],
            visible=val["visible"],
            transparent=val["transparent"],
            width=val["width"],
            marked=val["marked"],
            xfill=val["xfill"],
            animation=val["animation"],
            name=val["name"],
            source=val["source"],
        )
        default_layer_props.update({key: single_layer_prop})

    return default_layer_props
