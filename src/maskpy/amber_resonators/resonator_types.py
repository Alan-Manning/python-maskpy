from enum import StrEnum


class AmberResonatorType(StrEnum):
    """Type of a Resonator.

    Enum Members
    ------------
    ORIGINAL_NB_IND2
    ORIGINAL_NB_IND4
    ORIGINAL_AL_IND2
    ORIGINAL_AL_IND4
    """

    # ORIGINAL_Q50K = "original_q50k"
    ORIGINAL_NB_IND2 = "original_nb_ind2"
    ORIGINAL_NB_IND4 = "original_nb_ind4"
    ORIGINAL_AL_IND2 = "original_al_ind2"
    ORIGINAL_AL_IND4 = "original_al_ind4"

    @classmethod
    def _missing_(cls, value):
        value = value.upper()
        for member in cls:
            if member == value:
                return member
        return None
