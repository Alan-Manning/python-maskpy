from enum import StrEnum


class AmberResonatorType(StrEnum):
    """Type of a Resonator.

    Enum Members
    ------------
    ORIGINAL_Q50K
    """

    ORIGINAL_Q50K = "original_q50k"

    @classmethod
    def _missing_(cls, value):
        value = value.upper()
        for member in cls:
            if member == value:
                return member
        return None
