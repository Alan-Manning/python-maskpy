from numpy import cos, pi, sin

from ...souk_mask_builder import SoukMaskBuilder
from ...souk_mask_configs import SoukMaskConfig, get_mask_default_config


def add_hex_pack_base(
    mask_builder: SoukMaskBuilder,
    add_sma_connectors: bool = True,
) -> tuple[
    list[list[float]],
    list[list[float]],
    list[float],
    list[float],
    list[float],
    list[float],
]:
    """This adds a preconfigured hex pack which includes many things. This
    generates the exact split hex pack used in the horn block V1 and adds the
    tabbed dice line to fit with that layer sandwich. This also adds all the
    SMA connectors and mask logos.

    KwArgs
    ------
    add_sma_connectors: bool = True
        This will add the cpw structure under the sma connectors if True, Default.


    Returns
    -------
    This returns a lot of things in one tuple:
    (
        bot_hex_grid: list[list[float]],
        top_hex_grid: list[list[float]],
        sma_bot_right_connect_xyr: list[float],
        sma_mid_left_connect_xyr: list[float],
        sma_mid_right_connect_xyr: list[float],
        sma_top_right_connect_xyr: list[float],
    )
    """

    # The specific values used in milling horn block V1
    BIG_HEXAGON_RADIUS: float = 67587.17967697246
    CHIP_CENTER_XY: list[float] = [0.0, 0.0]
    HORIZONTAL_PITCH: float = 5306.0
    VERTICAL_PITCH: float = 4590.0
    MIDDLE_GAP_DISTANCE: int = 600

    bot_hex_grid, top_hex_grid = mask_builder.generate_half_split_hex_pack(
        BIG_HEXAGON_RADIUS,
        CHIP_CENTER_XY,
        HORIZONTAL_PITCH,
        VERTICAL_PITCH,
        middle_gap_distance=MIDDLE_GAP_DISTANCE,
        adjust_center_hex_grid=True,
        return_two_grids=True,
        return_split_point=False,
    )
    ###########################################################################
    # define a big hexagon for getting sma connector center positions
    big_hexagon_points: list[list[float]] = []
    for i in range(6):
        big_hexagon_points.append(
            [
                CHIP_CENTER_XY[0] + BIG_HEXAGON_RADIUS * cos(i * (pi / 3)),
                CHIP_CENTER_XY[1] + BIG_HEXAGON_RADIUS * sin(i * (pi / 3)),
            ]
        )
    ###########################################################################
    # coords_for_sma_centers
    sma_bot_right_xyr = [
        big_hexagon_points[5][0] - 2500,
        big_hexagon_points[5][1] + 2500,
        pi,
    ]
    sma_mid_left_xyr = [
        big_hexagon_points[3][0] + 6750,
        big_hexagon_points[3][1] + 4590,
        0,
    ]
    sma_mid_right_xyr = [
        big_hexagon_points[0][0] - 4300,
        big_hexagon_points[0][1],
        pi,
    ]
    sma_top_right_xyr = [
        big_hexagon_points[1][0] - 2500,
        big_hexagon_points[1][1] - 2500,
        pi,
    ]

    # connection positions for the feedline to the sma connectors
    if add_sma_connectors:
        sma_bot_right_connect_xyr = mask_builder.add_sma_connector_and_launcher_and_get_connection(
            sma_bot_right_xyr[0],
            sma_bot_right_xyr[1],
            sma_bot_right_xyr[2],
        )
        sma_bot_right_connect_xyr.append(sma_bot_right_xyr[2])

        sma_mid_left_connect_xyr = mask_builder.add_sma_connector_and_launcher_and_get_connection(
            sma_mid_left_xyr[0],
            sma_mid_left_xyr[1],
            sma_mid_left_xyr[2],
        )
        sma_mid_left_connect_xyr.append(sma_mid_left_xyr[2])

        sma_mid_right_connect_xyr = mask_builder.add_sma_connector_and_launcher_and_get_connection(
            sma_mid_right_xyr[0],
            sma_mid_right_xyr[1],
            sma_mid_right_xyr[2],
        )
        sma_mid_right_connect_xyr.append(sma_mid_right_xyr[2])

        sma_top_right_connect_xyr = mask_builder.add_sma_connector_and_launcher_and_get_connection(
            sma_top_right_xyr[0],
            sma_top_right_xyr[1],
            sma_top_right_xyr[2],
        )
        sma_top_right_connect_xyr.append(sma_top_right_xyr[2])
    else:
        sma_config = get_mask_default_config(SoukMaskConfig.SMA_CONNECTOR, config_override=None)

        offset_from_sma_center = (sma_config["sma_square_width"] / 2) + sma_config["sma_square_offset_right"]

        sma_bot_right_connect_xyr = [
            sma_bot_right_xyr[0] + offset_from_sma_center * cos(sma_bot_right_xyr[2]),
            sma_bot_right_xyr[1] + offset_from_sma_center * sin(sma_bot_right_xyr[2]),
            sma_bot_right_xyr[2],
        ]
        sma_mid_left_connect_xyr = [
            sma_mid_left_xyr[0] + offset_from_sma_center * cos(sma_mid_left_xyr[2]),
            sma_mid_left_xyr[1] + offset_from_sma_center * sin(sma_mid_left_xyr[2]),
            sma_mid_left_xyr[2],
        ]
        sma_mid_right_connect_xyr = [
            sma_mid_right_xyr[0] + offset_from_sma_center * cos(sma_mid_right_xyr[2]),
            sma_mid_right_xyr[1] + offset_from_sma_center * sin(sma_mid_right_xyr[2]),
            sma_mid_right_xyr[2],
        ]
        sma_top_right_connect_xyr = [
            sma_top_right_xyr[0] + offset_from_sma_center * cos(sma_top_right_xyr[2]),
            sma_top_right_xyr[1] + offset_from_sma_center * sin(sma_top_right_xyr[2]),
            sma_top_right_xyr[2],
        ]

    # Adding the bounding box of the sma connectors to the hex pack exclusions list
    sma_bot_right_bounding_box_points = mask_builder.get_sma_connector_and_laucher_bounding_box(
        sma_bot_right_xyr[0],
        sma_bot_right_xyr[1],
        sma_bot_right_xyr[2],
    )
    sma_mid_left_bounding_box_points = mask_builder.get_sma_connector_and_laucher_bounding_box(
        sma_mid_left_xyr[0],
        sma_mid_left_xyr[1],
        sma_mid_left_xyr[2],
    )
    sma_mid_right_bounding_box_points = mask_builder.get_sma_connector_and_laucher_bounding_box(
        sma_mid_right_xyr[0],
        sma_mid_right_xyr[1],
        sma_mid_right_xyr[2],
    )
    sma_top_right_bounding_box_points = mask_builder.get_sma_connector_and_laucher_bounding_box(
        sma_top_right_xyr[0],
        sma_top_right_xyr[1],
        sma_top_right_xyr[2],
    )

    ###########################################################################
    mask_builder.add_hexagonal_tabbed_dicing_line(
        big_hexagon_points,
        BIG_HEXAGON_RADIUS,
        CHIP_CENTER_XY,
        [sma_bot_right_xyr, sma_mid_left_xyr, sma_mid_right_xyr, sma_top_right_xyr],
    )

    # Adding a Bottom choke Tabbed dice line around the chip
    TABBED_DICT_LINE_CORNER_BEND_RAD: float = 5000.0

    mask_builder.add_hexagonal_bottom_choke_tabbed_dicing_line(
        big_hexagon_points,
        BIG_HEXAGON_RADIUS,
        CHIP_CENTER_XY,
        [sma_bot_right_xyr, sma_mid_left_xyr, sma_mid_right_xyr, sma_top_right_xyr],
        TABBED_DICT_LINE_CORNER_BEND_RAD,
    )

    ###########################################################################
    # Defining the center and slotted pin hole dimensions and addind them to
    # chip and adding their bounding box to the hex exclusions list.

    CENTER_PIN_RADIUS: float = 2020.0 / 2
    center_pin_points = mask_builder.add_center_pin_and_get_bounding_points(CHIP_CENTER_XY, CENTER_PIN_RADIUS)

    slotted_pin_center = top_hex_grid[-25]
    SLOTTED_PIN_LENGTH: float = 1000.0
    SLOTTED_PIN_RADIUS: float = 2020.0 / 2

    slotted_pin_points = mask_builder.add_slotted_pin_and_get_bounding_points(
        slotted_pin_center, SLOTTED_PIN_LENGTH, SLOTTED_PIN_RADIUS, CHIP_CENTER_XY
    )

    ###########################################################################
    # Making the hex exclusions list and adding all things to be exclusded

    hex_exclusions = []

    hex_exclusions.append(sma_bot_right_bounding_box_points)
    hex_exclusions.append(sma_mid_left_bounding_box_points)
    hex_exclusions.append(sma_mid_right_bounding_box_points)
    hex_exclusions.append(sma_top_right_bounding_box_points)

    hex_exclusions.append(center_pin_points)
    hex_exclusions.append(slotted_pin_points)

    extra_exclusion_top_left = [
        [26000, 49000],
        [26000, 54000],
        [35000, 54000],
        [35000, 49000],
    ]
    hex_exclusions.append(extra_exclusion_top_left)
    extra_exclusion_mid_right = [
        [57000, -5000],
        [57000, 5000],
        [61000, 5000],
        [61000, -5000],
    ]
    hex_exclusions.append(extra_exclusion_mid_right)
    extra_exclusion_mid_left = [
        [-62000, -1000],
        [-62000, 9500],
        [-54000, 9500],
        [-54000, -1000],
    ]
    hex_exclusions.append(extra_exclusion_mid_left)
    extra_exclusion_bot_right = [
        [26000, -57000],
        [26000, -50000],
        [36000, -50000],
        [36000, -57000],
    ]
    hex_exclusions.append(extra_exclusion_bot_right)

    bot_hex_grid, top_hex_grid = mask_builder.remove_exclusions_from_split_hex_pack([bot_hex_grid, top_hex_grid], hex_exclusions)

    return (
        bot_hex_grid,
        top_hex_grid,
        sma_bot_right_connect_xyr,
        sma_mid_left_connect_xyr,
        sma_mid_right_connect_xyr,
        sma_top_right_connect_xyr,
    )
