import os
from collections.abc import Sequence
from typing import Literal

import gdspy


def _get_logo_gds_files_path() -> str:
    """Get the path for where all the gds logo files are stored within maskpy.

    Returns
    -------
    file_path: str
        The file path for logo gds files.
    """
    file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "logo_gds_files")
    return file_path


def add_cardiff_logo(
    mask_builder,
    logo_xy: Sequence[float | int],
    size: float | int = 100,
) -> None:
    """Adds the cardiff logo as a cutout in the ng groundplane layer.
    Requires the logo gds file. Logo will be placed centered on the logo_xy
    given.

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
    """
    x = logo_xy[0]
    y = logo_xy[1]

    default_size = 100
    scale_factor = size / default_size

    logo_file_name = "Cardiff_Uni_Logo.gds"
    logo_file_path = _get_logo_gds_files_path()

    full_file = os.path.join(logo_file_path, logo_file_name)

    try:
        lib = gdspy.GdsLibrary(infile=full_file)
    except FileNotFoundError:
        print(f"Logo gds file '{logo_file_name}' not found.")
        print(f"full file and path: {full_file}")
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
        poly = gdspy.Polygon(
            text_polygons[i],
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.Main.add(poly)

    for i in range(len(box_polygons)):
        poly = gdspy.Polygon(
            box_polygons[i],
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.ground_plane_cutouts.add(poly)

    return


def add_so_logo(
    mask_builder,
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
    x = logo_xy[0]
    y = logo_xy[1]

    logo_file_name = "SO_Logo.gds"
    logo_file_path = _get_logo_gds_files_path()

    full_file = os.path.join(logo_file_path, logo_file_name)

    try:
        lib = gdspy.GdsLibrary(infile=full_file)
    except FileNotFoundError:
        print(f"Logo gds file '{logo_file_name}' not found.")
        print(f"full file and path: {full_file}")
        return
    except Exception:
        print("Issue with file " + logo_file_name + ".")
        return

    main_cell = lib.top_level()[0]
    pol_dict = main_cell.get_polygons(by_spec=True)

    main_text = pol_dict[(1, 0)]
    mountains = pol_dict[(2, 0)]
    telescope = pol_dict[(3, 0)]
    outer_ring = pol_dict[(4, 0)]
    light_beam = pol_dict[(5, 0)]
    stars = pol_dict[(6, 0)]
    outer_text = pol_dict[(7, 0)]
    mountains_without_telescope = pol_dict[(8, 0)]
    all_one_layer = pol_dict[(9, 0)]
    all_one_layer_without_light_beam = pol_dict[(10, 0)]
    full_cover_circle = pol_dict[(11, 0)]

    if draw_all_in_one_layer:
        default_size = 1800
        scale_factor = size / default_size

        for i in range(len(full_cover_circle)):
            poly = gdspy.Polygon(
                full_cover_circle[i],
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.ground_plane_cutouts.add(poly)

        for i in range(len(all_one_layer)):
            poly = gdspy.Polygon(
                all_one_layer[i],
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.Main.add(poly)

    return


def add_souk_logo(
    mask_builder,
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
    x = logo_xy[0]
    y = logo_xy[1]

    logo_file_name = "SOUK_Logo.gds"
    logo_file_path = _get_logo_gds_files_path()

    full_file = os.path.join(logo_file_path, logo_file_name)

    try:
        lib = gdspy.GdsLibrary(infile=full_file)
    except FileNotFoundError:
        print(f"Logo gds file '{logo_file_name}' not found.")
        print(f"full file and path: {full_file}")
        return
    except Exception:
        print("Issue with file " + logo_file_name + ".")
        return

    main_cell = lib.top_level()[0]
    pol_dict = main_cell.get_polygons(by_spec=True)

    all_one_layer = pol_dict[(0, 0)]
    all_one_layer_no_outer_text = pol_dict[(1, 0)]
    outer_cover_circle = pol_dict[(2, 0)]
    inner_cover_circle = pol_dict[(3, 0)]
    inner_text = pol_dict[(4, 0)]
    outer_text = pol_dict[(5, 0)]
    ground_and_circle_with_telescope_cutout = pol_dict[(6, 0)]
    ground_and_circle = pol_dict[(7, 0)]
    telescope = pol_dict[(8, 0)]
    light_beams = pol_dict[(9, 0)]
    stars = pol_dict[(10, 0)]

    if include_outer_ring_text:
        default_size = 231
    else:
        default_size = 201.9525

    scale_factor = size / default_size

    if draw_all_in_one_layer and include_outer_ring_text:
        for poly_points in outer_cover_circle:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.ground_plane_cutouts.add(poly)

        for poly_points in all_one_layer:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.Main.add(poly)

        return

    if draw_all_in_one_layer and (not include_outer_ring_text):
        for poly_points in inner_cover_circle:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.ground_plane_cutouts.add(poly)

        for poly_points in all_one_layer_no_outer_text:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.Main.add(poly)

        return

    if include_outer_ring_text:
        for poly_points in outer_cover_circle:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.ground_plane_cutouts.add(poly)

        for poly_points in outer_text:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.Main.add(poly)
    else:
        for poly_points in inner_cover_circle:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.ground_plane_cutouts.add(poly)

    for poly_points in inner_text:
        poly = gdspy.Polygon(
            poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.Main.add(poly)

    for poly_points in ground_and_circle_with_telescope_cutout:
        poly = gdspy.Polygon(
            poly_points,
            layer=mask_builder.layers.Aluminium.number,
            datatype=mask_builder.layers.Aluminium.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.Main.add(poly)

    for poly_points in telescope:
        poly = gdspy.Polygon(
            poly_points,
            layer=mask_builder.layers.Aluminium.number,
            datatype=mask_builder.layers.Aluminium.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.Main.add(poly)

    for poly_points in light_beams:
        poly = gdspy.Polygon(
            poly_points,
            layer=mask_builder.layers.Aluminium.number,
            datatype=mask_builder.layers.Aluminium.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.Main.add(poly)

    for poly_points in stars:
        poly = gdspy.Polygon(
            poly_points,
            layer=mask_builder.layers.Nb_Groundplane.number,
            datatype=mask_builder.layers.Nb_Groundplane.datatype,
        )
        poly.scale(scale_factor, scale_factor)
        poly.translate(x, y)
        mask_builder.Main.add(poly)

    return


def add_AM_signature(
    mask_builder,
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
    x = signature_xy[0]
    y = signature_xy[1]

    signature_file_name = "AM_signature.gds"
    signature_file_path = _get_logo_gds_files_path()

    full_file = os.path.join(signature_file_path, signature_file_name)

    try:
        lib = gdspy.GdsLibrary(infile=full_file)
    except FileNotFoundError:
        print(f"Logo gds file '{signature_file_name}' not found.")
        print(f"full file and path: {full_file}")
        return
    except Exception:
        print("Issue with file " + signature_file_name + ".")
        return

    main_cell = lib.top_level()[0]
    pol_dict = main_cell.get_polygons(by_spec=True)

    full_sig = pol_dict[(0, 0)]
    stacked_sig = pol_dict[(1, 0)]
    forename = pol_dict[(2, 0)]
    surname = pol_dict[(3, 0)]

    default_size = 889
    polygons_for_variant = full_sig

    match variant:
        case 1:
            default_size = 493
            polygons_for_variant = stacked_sig
        case 2:
            default_size = 492
            polygons_for_variant = forename
        case 3:
            default_size = 404
            polygons_for_variant = surname
        # case 0 or anything else is default.

    scale_factor = size / default_size

    if draw_as_cutout:
        for poly_points in polygons_for_variant:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Nb_Groundplane.number,
                datatype=mask_builder.layers.Nb_Groundplane.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.ground_plane_cutouts.add(poly)
    else:
        for poly_points in polygons_for_variant:
            poly = gdspy.Polygon(
                poly_points,
                layer=mask_builder.layers.Aluminium.number,
                datatype=mask_builder.layers.Aluminium.datatype,
            )
            poly.scale(scale_factor, scale_factor)
            poly.translate(x, y)
            mask_builder.Main.add(poly)

    return
