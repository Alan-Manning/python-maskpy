from .css4_colors import CSS4_COLORS
from .material_colors import MATERIAL_COLORS
from .matplotlib_colors import MATPLOTLIB_TABLEAU_COLORS

all_colors: dict[str, str] = {}

all_colors.update(MATPLOTLIB_TABLEAU_COLORS)
all_colors.update(MATERIAL_COLORS)
all_colors.update(CSS4_COLORS)

all_colors_in_palettes: dict[str, dict[str, str]] = {
    "MATPLOTLIB_TABLEAU_COLORS": MATPLOTLIB_TABLEAU_COLORS,
    "MATERIAL_COLORS": MATERIAL_COLORS,
    "CSS4_COLORS": CSS4_COLORS,
}


def hex_to_dec(hex_str: str) -> int:
    """Convert a hexadecimal string value to an integer.

    Parameters
    ----------
    hex_str : str
        Hexidecimal string, can include leading "#" or not.

    Returns
    -------
    integer : int
        The int representation of the hex string input.
    """
    hex_str = hex_str[1:] if hex_str[0] == "#" else hex_str
    integer = int(sum(int(x, 16) * pow(16, len(hex_str) - i - 1) for i, x in enumerate(hex_str)))
    return integer


def hex_to_rbg(hex_col: str) -> tuple[int, int, int]:
    """Convert a hexadecimal color string value to an rgb color.

    Parameters
    ----------
    hex_col : str
        Hexidecimal color string, can include leading "#" or not.
        eg, "#fefefe" or "fefefe".

    Returns
    -------
    rbg_col : tuple[int, int, int]
        List of [red, green, blue], each being int values 0 to 255 representing
        the hex string input.
    """
    hex_col = hex_col[1:] if hex_col[0] == "#" else hex_col
    hex_r = hex_col[:2]
    hex_g = hex_col[2:4]
    hex_b = hex_col[4:6]
    rbg_col = (hex_to_dec(hex_r), hex_to_dec(hex_g), hex_to_dec(hex_b))
    return rbg_col


def get_color_palettes() -> list[str]:
    """Get all the color palettes defined.

    Returns
    -------
    all_palettes : list[str]
        List of strings whcih are the name of all the color palettes defined.
    """
    return list(all_colors_in_palettes.keys())


def get_colors(palette: str) -> dict[str, str]:
    """Get a color pallete of hex codes from defined palletes.

    Parameters
    ----------
    palette : str
        String of a palette name, to see all availible palettes use
        get_color_palettes().

    Returns
    -------
    palette : dict
        A dictionary containing keys with color names with values of hex codes.
        Hex code is a string of format "#xxyyzz" including the "#".
    """
    match palette:
        case "MATPLOTLIB_TABLEAU_COLORS":
            return MATPLOTLIB_TABLEAU_COLORS
        case "MATERIAL_COLORS":
            return MATERIAL_COLORS
        case "CSS4_COLORS":
            return CSS4_COLORS
    raise (ValueError(f"Unable to find pallete `{palette}`"))


def get_closest_color(hex_color: str) -> tuple[str, str, str]:
    """Gets the closest color in all the color palettes defined. Take in a hex
    code and returns the closest color name, hex, and the pallete its in.

    Parameters
    ----------
    hex_color : str
        String containing a hex color including the leading hash.
        eg. "#fefefe".

    Returns
    -------
    (closest_color_name, closest_color_hex, in_palette): tuple[str, str, str]
        closest_color_name: str
            The closest color's name
        closest_color_hex: str
            The closest color's Hex code in format "#xxyyzz", including the #.
        in_palette: str
            The pallete who the closest color belongs to.
    """
    rgb_color = hex_to_rbg(hex_color)

    for palette in all_colors_in_palettes.keys():
        for i, key in enumerate(all_colors_in_palettes[palette]):
            hex_col = all_colors_in_palettes[palette][key][1:]
            rgb = hex_to_rbg(hex_col)
            diff = ((rgb_color[0] - rgb[0]) ** 2 + (rgb_color[1] - rgb[1]) ** 2 + (rgb_color[2] - rgb[2]) ** 2) ** 0.5

            if i == 0:
                prev_lowest_diff = diff + 1

            if diff < prev_lowest_diff:
                prev_lowest_diff = diff
                closest_color_name = key
                in_palette = palette
                closest_color_hex = all_colors[key]

                if diff == 0:
                    return closest_color_name, closest_color_hex, in_palette

    return closest_color_name, closest_color_hex, in_palette
