# from .default_layer_props import LayerPropType, get_default_layer_properties
# from .layer_set import Layer, SoukMaskLayerSet
import xml.etree.ElementTree as ET

from ..colors import all_colors
from . import ANLMaskLayerSet, LayerPropType, SoukMaskLayerSet, get_default_layer_properties


def generate_layer_properties(
    filename: str,
    layer_numbers_in_mask: list[int],
    layer_prop_type: LayerPropType,
    layer_set: SoukMaskLayerSet | ANLMaskLayerSet,
    add_layer_number_to_name: bool = True,
) -> None:
    """Generate layer properties file.

    Parameters
    ----------
    filename : str
        This is the name for the layer_properties file. If this does not
        include a ".lyp" file extention it will be automatically added.

    layer_numbers_in_mask: list[int]
        This is a list of integers that represent all the layer numbers that
        exist in the mask to generate a layer properties file for.

    layer_prop_type: LayerPropType
        The LayerPropType.

    layer_set: SoukMaskLayerSet | ANLMaskLayerSet
        The SoukMaskLayerSet.

    KwArgs
    ------
    add_layer_number_to_name: bool = True
        This by default will add the layer number to the name of the layer in
        the layer properties file. This takes the form, `[num] - layer_name`.
        When False the layer name will just be `layer_name`.

    Outputs
    -------
    xml_layer_props: `filename.lyp`
    """

    default_layer_props = get_default_layer_properties(layer_prop_type)

    layer_properties = ET.Element("layer-properties")  # make the root of the layer props file

    layer_numbers_in_mask.sort()

    for layer_number in layer_numbers_in_mask:
        layer_number = layer_number % (2**16)  # mod 65536 so layer nums are not negative.

        layer_name = layer_set.get_layer_name_from_number(layer_number)

        # print(f"default_layer_props: {default_layer_props}")
        # print(f"layer_set: {layer_set}")
        # print(f"layer_name: {layer_name}")
        # print(f"layer_number: {layer_number}")
        #
        # print("")
        # test = default_layer_props.get(layer_name)
        # print(f"test: {test}")
        #

        # if layer_name is None:
        #     single_layer_property = default_layer_props["NoName"]
        # else:
        #     single_layer_property = default_layer_props.get(layer_name)
        #     if single_layer_property is None:
        #         single_layer_property = default_layer_props["NoName"]
        #

        # print(f"single_layer_property: {single_layer_property}")

        single_layer_property = (
            default_layer_props["NoName"] if layer_name is None else default_layer_props.get(layer_name, default_layer_props["NoName"])
        )

        if layer_name is None:
            frame_color = fill_color = all_colors[list(all_colors.keys())[(layer_number * 10) % len(all_colors)]]
            name = "NoName"
        else:
            frame_color = single_layer_property.frame_color
            if not frame_color.startswith("#"):
                frame_color = all_colors[frame_color]

            fill_color = single_layer_property.fill_color
            if not fill_color.startswith("#"):
                fill_color = all_colors[fill_color]

            name = single_layer_property.name

        if add_layer_number_to_name:
            name = f"[{layer_number}] - {name}"

        source = f"{layer_number}/{0}@1"  # TODO 0 should be the layer_datatype

        single_layer_property.frame_color = frame_color
        single_layer_property.fill_color = fill_color
        single_layer_property.name = name
        single_layer_property.source = source

        layer_properties = single_layer_property.add_to_xml_layer_properties(layer_properties)

    # print(f"filename before: {filename}")
    if filename.split(".")[-1] == "lyp":
        filename = str(filename.split(".")[0])
    filename += ".lyp"
    # print(f"filename after: {filename}")

    layer_properties_root = ET.ElementTree(layer_properties)
    layer_properties_root.write(filename, encoding="utf-8", xml_declaration=True)

    print(f'Succesfully written layerprops file "{filename}"\n')

    return


# def append_layer_to_xml(
#     xml_layer_prop: ET.Element,
#     single_xml_layer_props: ET.Element,
# ) -> ET.Element:
#     """Appends a layer property element to the xml layer properties object.
#
#     Parameters
#     ----------
#     single_xml_layer_props : ET.Element
#     """
#     prop = ET.SubElement(xml_layer_prop, "properties")
#
#     for key, val in single_layer_props.items():
#         layer_element = ET.SubElement(prop, key)
#         layer_element.text = val
#
#     return xml_layer_prop
