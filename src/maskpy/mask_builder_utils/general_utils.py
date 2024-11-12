from datetime import datetime
from typing import NoReturn

from numpy import pi

from maskpy.logging.pretty_print import TextColor, styled_text, styled_type_error


def deg_to_rad(deg: int | float) -> float:
    """Converts an angle given in degrees to radians."""
    return (deg / 180.0) * pi


def rad_to_deg(rad: int | float) -> float:
    """Converts an angle given in radians to degrees."""
    return rad * (180.0 / pi)


def get_time(
    time: datetime,
    time_format: str = "%H#:%M# #-# %D# %B# %Y",
    usetex: bool = False,
) -> str:
    """Get the datetime as a string formatted with the character directives
    and.

    codes given. i.e. convert a datetime object to something like
    `14:05 - 1st October 2024` with `time_format="%H#:%M# #-# %D# %B# %Y"`.
    This can generate some parts in latex super sub script when usetex is True.

    Parameters
    ----------
    time: datetime
        A datetime object to get the human date format for.

    KwArgs
    ------
    time_format: str = "%H#:%M# #-# %D# %B# %Y"
        This will by default print `14:05 - 1st October 2024`.
        However this can be formatted any way by providing your own directive
        and code pairs. These are all detailed below in the `Notes` but to give
        an example using the same datetime above.
        >>> get_time(
        >>>     datetime.now(),
        >>>     time_format="%Y#-%m#-%d# %H#:%M#:%S",
        >>> )
        >>> # "2024-10-01 14:05:00"

    usetex: bool = False,
        When True This will use latex to put specific parts in supersciprt.
        This is namely the day suffix i.e. st in 1st, nd in 2nd...

    Returns
    -------
    text: str
        The formatted time string.

    Notes
    -----
    Below are all the directives and codes
    +-----+-------------------------------------------------------------------+
    | D C | Directive and Code - meaning - Example                            |
    +=====+===================================================================+
    | #?  | This is a literal character. Anything followed by a a # will be   |
    |     | literally interpreted as that character written after the #.      |
    | %a  | Weekday as locale’s abbreviated name.                             |
    |     |     Example - Sun, Mon, …, Sat                                    |
    | %A  | Weekday as locale’s full name.                                    |
    |     |     Example - Sunday, Monday, …, Saturday                         |
    | %w  | Weekday as a decimal number, where 0 is Sunday and 6 is Saturday. |
    |     |     Example - 0, 1, …, 6                                          |
    | %d  | Day of the month as a zero-padded decimal number.                 |
    |     |     Example - 01, 02, …, 31                                       |
    | %D  | Day of the month as a decimal number with day suffix.             |
    |     |     Example - 1st, 2nd, 3rd, 4th …, 31st                          |
    | %b  | Month as locale’s abbreviated name.                               |
    |     |     Example - Jan, Feb, …, Dec                                    |
    | %B  | Month as locale’s full name.                                      |
    |     |     Example - January, February, …, December                      |
    | %m  | Month as a zero-padded decimal number.                            |
    |     |     Example - 01, 02, …, 12                                       |
    | %y  | Year without century as a zero-padded decimal number.             |
    |     |     Example - 00, 01, …, 99                                       |
    | %Y  | Year with century as a decimal number.                            |
    |     |     Example - 0001, 0002, …, 2013, 2014, …, 9998, 9999            |
    | %H  | Hour (24-hour clock) as a zero-padded decimal number.             |
    |     |     Example - 00, 01, …, 23                                       |
    | %I  | Hour (12-hour clock) as a zero-padded decimal number.             |
    |     |     Example - 01, 02, …, 12                                       |
    | %p  | Locale’s equivalent of either AM or PM.                           |
    |     |     Example - AM, PM                                              |
    | %M  | Minute as a zero-padded decimal number.                           |
    |     |     Example - 00, 01, …, 59                                       |
    | %S  | Second as a zero-padded decimal number.                           |
    |     |     Example - 00, 01, …, 59                                       |
    | %f  | Microsecond as a decimal number, zero-padded to 6 digits.         |
    |     |     Example - 000000, 000001, …, 999999                           |
    | %z  | UTC offset in the form; ±HHMM[SS[.ffffff]];                       |
    |     |     Example - (empty), +0000, -0400, +1030, +063415, -030712.345  |
    | %Z  | Time zone name (empty string if the object is naive).             |
    |     |     Example - (empty), UTC, GMT                                   |
    | %j  | Day of the year as a zero-padded decimal number.                  |
    |     |     Example - 001, 002, …, 366                                    |
    | %U  | Week number of the year as a zero-padded decimal number           |
    |     |     Example - 00, 01, …, 53                                       |
    | %W  | Week number of the year as a zero-padded decimal number           |
    |     |     Example - 00, 01, …, 53                                       |
    | %c  | Locale’s appropriate date and time representation.                |
    |     |     Example - Tue Aug 16 21:30:00 1988                            |
    | %x  | Locale’s appropriate date representation.                         |
    |     |     Example - 08/16/88 (None);08/16/1988                          |
    | %X  | Locale’s appropriate time representation.                         |
    |     |     Example - 21:30:00                                            |
    +-----+-------------------------------------------------------------------+
    """
    if not isinstance(time_format, str):
        styled_type_error(time_format, "time_format", str)

    if len(time_format) % 2 == 1:
        raise ValueError("Malformed time_format. Every part should be either `%` or `#` followed a single character.")

    def time_format_error(
        time_format: str,
        error_char: str,
        directive_or_code: str,
        valid_options: set,
        pos: int,
    ) -> NoReturn:
        """Raise a time format error for an incorrectly formatted time code
        string.

        Parameters
        ----------
        error_char: str,
            The character that is the error.
        directive_or_code: str,
            Whether this character is part of the directive or the code.
        valid_options: set,
            The valid options this could be
        pos: int,
            The position in the time_format string
        """
        sub = ""
        if directive_or_code == "directive":
            sub += styled_text(time_format[: (2 * pos)], color=TextColor.GREEN)
            sub += styled_text(time_format[(2 * pos) : (2 * pos) + 1], color=TextColor.ERROR)
            sub += time_format[(2 * pos) + 1 :]
            sub += '"'
            sub += "\n"
            sub += "_"
            sub += styled_text(("__" * pos), color=TextColor.GREEN)
            sub += styled_text("^", color=TextColor.RED)
            sub += "_"
        if directive_or_code == "code":
            sub += styled_text(time_format[: (2 * pos) + 1], color=TextColor.GREEN)
            sub += styled_text(time_format[(2 * pos) + 1 : (2 * pos) + 2], color=TextColor.ERROR)
            sub += time_format[(2 * pos) + 2 :]
            sub += '"'
            sub += "\n"
            sub += "_"
            sub += styled_text(("__" * pos), color=TextColor.GREEN)
            sub += "_"
            sub += styled_text("^", color=TextColor.RED)
        raise ValueError(
            f"Encountered invaid {directive_or_code} `"
            + styled_text(f"{error_char}", color=TextColor.ERROR)
            + "` in `time_format`. Expected one of "
            + styled_text(f"{valid_options}", color=TextColor.GREEN)
            + ".\n"
            + '"'
            + sub
            + "__" * ((int(len(time_format) / 2)) - pos - 1)
            + "_"
        )

    valid_time_directives = set("%#")
    valid_time_codes = set("aAwdDbBmyYHIpMSfzZjUWcxX")
    directives = time_format[::2]
    codes = time_format[1::2]

    time_text = ""
    for pos, (directive, code) in enumerate(zip(directives, codes, strict=False)):
        if directive not in valid_time_directives:
            time_format_error(
                time_format,
                directive,
                "directive",
                valid_time_directives,
                pos,
            )

        if directive == "#":
            time_text += code
            continue

        if code not in valid_time_codes:
            time_format_error(
                time_format,
                code,
                "code",
                valid_time_codes,
                pos,
            )

        # Handle %D code as its not in native strftime.
        if code != "D":
            time_text += time.strftime(f"{directive}{code}")
        else:
            directive_text = ""
            day = time.strftime("%d")
            if day[0] == "0":
                day = day[1]
            directive_text += day

            if usetex:
                directive_text += r"$^{"
            match day[-1]:
                case "1":
                    directive_text += "st"
                case "2":
                    directive_text += "nd"
                case "3":
                    directive_text += "rd"
                case _:
                    directive_text += "th"
            if usetex:
                directive_text += r"}$"

            time_text += directive_text

    return time_text
