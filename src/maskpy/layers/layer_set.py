from __future__ import annotations

from dataclasses import dataclass, fields


@dataclass(init=False, repr=False, eq=False, frozen=True)
class Layer:
    """Mask layer. This is an immutable layer for use in mask_builder.

    Parameters
    ----------
    layer_name: str
        The name for the layer.

    number: int
        The number for the layer. This should be within 0->2^16 inclusive.

    datatype: int
        The datatype for the layer. This should be within 0->255 inclusive.

    """

    name: str
    number: int
    datatype: int

    def __init__(self, layer_name: str, layer_number: int, datatype: int):
        if not isinstance(layer_name, str):
            raise TypeError(f"Layer's name should be of type str, not {type(layer_name)}.")

        if not isinstance(layer_number, int):
            raise TypeError(f"Layer's number should be of type int, not {type(layer_number)}.")
        self._validate_number(layer_number)

        if not isinstance(datatype, int):
            raise TypeError(f"Layer's datatype should be of type int, not {type(datatype)}.")
        self._validate_datatype(datatype)

        object.__setattr__(self, "name", layer_name)
        object.__setattr__(self, "number", layer_number)
        object.__setattr__(self, "datatype", datatype)

    def __repr__(self) -> str:
        rep = f"Layer:\n   name   = {self.name}\n   number = {self.number}\n   datatype  = {self.datatype}"
        return rep

    def __str__(self) -> str:
        return str(f"{self.name} [{self.number}/{self.datatype}]")

    def __eq__(self, other) -> bool:
        if not isinstance(other, Layer):
            return False

        return (self.number == other.number) and (self.datatype == other.datatype)

    def _validate_number(self, layer_number: int) -> None:
        min = 0
        max = 32767
        # max = 65535

        if min <= layer_number <= max:
            return
        else:
            raise ValueError(f"Layer number should be within inclusive range ({min} -> {max}). Current layer number is {layer_number}")

    def _validate_datatype(self, datatype: int) -> None:
        min = 0
        max = 32767
        # max = 65535

        if min <= datatype <= max:
            return
        else:
            raise ValueError(f"Layer datatype should be within inclusive range ({min} -> {max}). Current datatype is {datatype}")


DEFAULT_SOUK_LAYERS = {
    "Aluminium": (1, 0),
    "Nb_Antenna": (2, 0),
    "SiN_dep": (3, 0),
    "IDC_Nb": (4, 0),
    "SiO": (5, 0),
    "SiN_Membrane": (6, 0),
    "Aluminium_Bulk": (7, 0),
    "Aluminium_Direct": (8, 0),
    "Backside_Check": (9, 0),
    "Backside_Flipped": (29, 0),
    "Nb_Groundplane": (30, 0),
    "Pin_hole_positives": (31, 0),
    "TrimLayer": (90, 0),
    "Aluminium_Patch": (100, 0),
    "Aluminium_Etch": (101, 0),
    "Nb_Patch": (102, 0),
    "Nb_Etch": (103, 0),
    "Top_choke_waveguide_hole": (150, 0),
    "Top_choke_anulus": (151, 0),
    "Bottom_choke_waveguide_hole": (155, 0),
    "Bottom_choke_IDC_hole": (156, 0),
    "Bottom_choke_pads": (157, 0),
    "Bottom_choke_backshort": (158, 0),
    "Tab_dicing_line": (160, 0),
    "Bottom_choke_Tab_dicing_line": (161, 0),
    "Chip_holder": (175, 0),
    "General_labeling": (32767, 0),
}
# "General_labeling": (65535, 0),

DEFAULT_ANL_LAYERS: dict[str, tuple[int, int]] = {
    "nb_wiring": (103, 0),
    "aluminium": (104, 0),
    "dielectric": (105, 0),
    "groundplane": (106, 0),
    "oxide": (107, 0),
    "membrane": (108, 0),
    "al_patch": (109, 0),
    "al_etch": (110, 0),
    "nb_patch": (111, 0),
    "nb_etch": (112, 0),
}


@dataclass(init=False, repr=True, eq=True, frozen=True)
class SoukMaskLayerSet:
    """SoukMaskLayerSet. This is an immutable dataclass filled with Layer
    atributes for use in a mask.

    KwArgs
    ------
    Any of the Atributes can be overriden with an altername Layer instance.
    These are all by default None which will set the layer to the default
    layer number and data type. Some example constructions would be:
    >>> layer_set = SoukMaskLayerSet()
    >>>
    >>> al_lay_name = "Ali"
    >>> al_lay_no = 200
    >>> al_lay_dtype = 1
    >>> layer_set2 = SoukMaskLayerSet(
    >>>     Aluminium=Layer(al_lay_name, al_lay_no, al_lay_dtype),
    >>> )
    >>>
    >>> al_layer = Layer("Ali", 200, 0)
    >>> nb_layer = Layer("Niob", 300, 0)
    >>> layer_set3 = SoukMaskLayerSet(
    >>>     Aluminium=al_layer,
    >>>     Nb_Antenna=nb_layer,
    >>> )

    Atributes
    ---------
    Aluminium: Layer
    Nb_Antenna: Layer
    SiN_dep: Layer
    IDC_Nb: Layer
    SiO: Layer
    SiN_Membrane: Layer
    Aluminium_Bulk: Layer
    Aluminium_Direct: Layer
    Backside_Check: Layer
    Backside_Flipped: Layer
    Nb_Groundplane: Layer
    Pin_hole_positives: Layer
    TrimLayer: Layer
    Aluminium_Patch: Layer
    Aluminium_Etch: Layer
    Nb_Patch: Layer
    Nb_Etch: Layer
    Top_choke_waveguide_hole: Layer
    Top_choke_anulus: Layer
    Bottom_choke_waveguide_hole: Layer
    Bottom_choke_IDC_hole: Layer
    Bottom_choke_pads: Layer
    Bottom_choke_backshort: Layer
    Tab_dicing_line: Layer
    Bottom_choke_Tab_dicing_line: Layer
    Chip_holder: Layer
    General_labeling: Layer
    """

    Aluminium: Layer
    Nb_Antenna: Layer
    SiN_dep: Layer
    IDC_Nb: Layer
    SiO: Layer
    SiN_Membrane: Layer
    Aluminium_Bulk: Layer
    Aluminium_Direct: Layer
    Backside_Check: Layer
    Backside_Flipped: Layer
    Nb_Groundplane: Layer
    Pin_hole_positives: Layer
    TrimLayer: Layer
    Aluminium_Patch: Layer
    Aluminium_Etch: Layer
    Nb_Patch: Layer
    Nb_Etch: Layer
    Top_choke_waveguide_hole: Layer
    Top_choke_anulus: Layer
    Bottom_choke_waveguide_hole: Layer
    Bottom_choke_IDC_hole: Layer
    Bottom_choke_pads: Layer
    Bottom_choke_backshort: Layer
    Tab_dicing_line: Layer
    Bottom_choke_Tab_dicing_line: Layer
    Chip_holder: Layer
    General_labeling: Layer

    def __init__(
        self,
        Aluminium: Layer | None = None,
        Nb_Antenna: Layer | None = None,
        SiN_dep: Layer | None = None,
        IDC_Nb: Layer | None = None,
        SiO: Layer | None = None,
        SiN_Membrane: Layer | None = None,
        Aluminium_Bulk: Layer | None = None,
        Aluminium_Direct: Layer | None = None,
        Backside_Check: Layer | None = None,
        Backside_Flipped: Layer | None = None,
        Nb_Groundplane: Layer | None = None,
        Pin_hole_positives: Layer | None = None,
        TrimLayer: Layer | None = None,
        Aluminium_Patch: Layer | None = None,
        Aluminium_Etch: Layer | None = None,
        Nb_Patch: Layer | None = None,
        Nb_Etch: Layer | None = None,
        Top_choke_waveguide_hole: Layer | None = None,
        Top_choke_anulus: Layer | None = None,
        Bottom_choke_waveguide_hole: Layer | None = None,
        Bottom_choke_IDC_hole: Layer | None = None,
        Bottom_choke_pads: Layer | None = None,
        Bottom_choke_backshort: Layer | None = None,
        Tab_dicing_line: Layer | None = None,
        Bottom_choke_Tab_dicing_line: Layer | None = None,
        Chip_holder: Layer | None = None,
        General_labeling: Layer | None = None,
    ):

        temp_init_loop_vars = [
            ("Aluminium", Aluminium),
            ("Nb_Antenna", Nb_Antenna),
            ("SiN_dep", SiN_dep),
            ("IDC_Nb", IDC_Nb),
            ("SiO", SiO),
            ("SiN_Membrane", SiN_Membrane),
            ("Aluminium_Bulk", Aluminium_Bulk),
            ("Aluminium_Direct", Aluminium_Direct),
            ("Backside_Check", Backside_Check),
            ("Backside_Flipped", Backside_Flipped),
            ("Nb_Groundplane", Nb_Groundplane),
            ("Pin_hole_positives", Pin_hole_positives),
            ("TrimLayer", TrimLayer),
            ("Aluminium_Patch", Aluminium_Patch),
            ("Aluminium_Etch", Aluminium_Etch),
            ("Nb_Patch", Nb_Patch),
            ("Nb_Etch", Nb_Etch),
            ("Top_choke_waveguide_hole", Top_choke_waveguide_hole),
            ("Top_choke_anulus", Top_choke_anulus),
            ("Bottom_choke_waveguide_hole", Bottom_choke_waveguide_hole),
            ("Bottom_choke_IDC_hole", Bottom_choke_IDC_hole),
            ("Bottom_choke_pads", Bottom_choke_pads),
            ("Bottom_choke_backshort", Bottom_choke_backshort),
            ("Tab_dicing_line", Tab_dicing_line),
            ("Bottom_choke_Tab_dicing_line", Bottom_choke_Tab_dicing_line),
            ("Chip_holder", Chip_holder),
            ("General_labeling", General_labeling),
        ]

        lay_num_dt_to_name_map: dict[tuple[int, int], str] = {}
        lay_name_to_layer_map: dict[str, Layer] = {}

        for name, lay in temp_init_loop_vars:
            if lay is None:
                def_layer_num, def_layer_dtype = DEFAULT_SOUK_LAYERS[name]
                object.__setattr__(
                    self,
                    name,
                    Layer(name, def_layer_num, def_layer_dtype),
                )
                lay_num_dt_to_name_map.update({(def_layer_num, def_layer_dtype): name})
                lay_name_to_layer_map.update({name: Layer(name, def_layer_num, def_layer_dtype)})
            elif isinstance(lay, Layer):
                lay_num_dt_to_name_map = self._check_layer_clash(lay, lay_num_dt_to_name_map)
                object.__setattr__(self, name, lay)
            else:
                raise TypeError(
                    f"Failed to create SoukMaskLayerSet with {name} set to type {type(lay)}. Should be of type `Layer` or None for defualt values."
                )

        object.__setattr__(self, "_layer_num_dt_to_name_map", lay_num_dt_to_name_map)
        object.__setattr__(self, "_layer_name_to_layer_map", lay_name_to_layer_map)

    # def _check_layer_clash(self, layer: Layer, lays_already_in: set[tuple[int, int]]):
    def _check_layer_clash(self, layer: Layer, lay_num_dt_to_name_map: dict[tuple[int, int], str]):
        if (layer.number, layer.datatype) in lay_num_dt_to_name_map.keys():
            raise ValueError(f"Cannot create SoukMaskLayerSet with layer: {layer} as it collides with another layer.")

        return lay_num_dt_to_name_map.update({(layer.number, layer.datatype): layer.name})

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        string = "SoukMaskLayerSet:\n"
        for field in fields(self):
            string += f"    {field.name}:\n"
            string += f"    ├─  number : {getattr(self, field.name).number}\n"
            string += f"    └─  dtype  : {getattr(self, field.name).datatype}\n"

        return string

    def get_layer_name_from_number(self, layer_number: int, datatype: int = 0, not_found_return_val: str | None = None) -> str | None:
        """Get the layer name for a given layer_number. If the layer number
        does not exist None will be returned.

        Parameters
        ----------
        layer_number: int
            The layer number who's name to get.

        KwArgs
        ------
        datatype: int = 0
            The datatype for the layer number. This is Default 0.

        not_found_return_val: str | None = None
            This is a value to return if that layer name cannot be found from
            the layer number provided. This defaults to None but can optionally
            be any other string.

        Returns
        -------
        layer_name: str | None
            The layer name if the layer number exists else None.
        """

        if not isinstance(layer_number, int):
            raise TypeError(f"layer_number should be of type int not {type(layer_number)}")
        if not isinstance(datatype, int):
            raise TypeError(f"datatype should be of type int not {type(datatype)}")

        layer_name: str | None = self._layer_num_dt_to_name_map.get((layer_number, datatype))
        if layer_name is None:
            return not_found_return_val

        return layer_name

    def get_layer_from_layer_name(self, layer_name: str) -> Layer:
        """Get the layer object from a layer name str.

        Parameters
        ----------
        layer_name: str
            The name of the layer to get the Layer object for.

        Returns
        -------
        layer: Layer
            The Layer object for the layer name given.
        """
        if not isinstance(layer_name, str):
            raise TypeError(f"layer_name should be of type str not {type(layer_name)}")

        layer = self._layer_name_to_layer_map.get(layer_name, None)

        if layer is None:
            raise LookupError(f"Could not find Layer object with layer name `{layer_name}` in SoukMaskLayerSet.")

        return layer


@dataclass(init=False, repr=True, eq=True, frozen=True)
class ANLMaskLayerSet:
    """ANLMaskLayerSet. This is an immutable dataclass filled with Layer
    atributes for use in a mask.

    KwArgs
    ------
    Any of the Atributes can be overriden with an altername Layer instance.
    These are all by default None which will set the layer to the default
    layer number and data type. Some example constructions would be:
    >>> layer_set = ANLMaskLayerSet()
    >>>
    >>> al_lay_name = "Ali"
    >>> al_lay_no = 200
    >>> al_lay_dtype = 1
    >>> layer_set2 = ANLMaskLayerSet(
    >>>     aluminium=Layer(al_lay_name, al_lay_no, al_lay_dtype),
    >>> )
    >>>
    >>> al_layer = Layer("Ali", 200, 0)
    >>> nb_layer = Layer("Niob", 300, 0)
    >>> layer_set3 = ANLMaskLayerSet(
    >>>     aluminium=al_layer,
    >>>     nb_wiring=nb_layer,
    >>> )

    Atributes
    ---------
    nb_wiring: Layer
    aluminium: Layer
    dielectric: Layer
    groundplane: Layer
    oxide: Layer
    membrane: Layer
    al_patch: Layer
    al_etch: Layer
    nb_patch: Layer
    nb_etch: Layer
    """

    nb_wiring: Layer
    aluminium: Layer
    dielectric: Layer
    groundplane: Layer
    oxide: Layer
    membrane: Layer
    al_patch: Layer
    al_etch: Layer
    nb_patch: Layer
    nb_etch: Layer

    def __init__(
        self,
        nb_wiring: Layer | None = None,
        aluminium: Layer | None = None,
        dielectric: Layer | None = None,
        groundplane: Layer | None = None,
        oxide: Layer | None = None,
        membrane: Layer | None = None,
        al_patch: Layer | None = None,
        al_etch: Layer | None = None,
        nb_patch: Layer | None = None,
        nb_etch: Layer | None = None,
    ):

        temp_init_loop_vars = [
            ("nb_wiring", nb_wiring),
            ("aluminium", aluminium),
            ("dielectric", dielectric),
            ("groundplane", groundplane),
            ("oxide", oxide),
            ("membrane", membrane),
            ("al_patch", al_patch),
            ("al_etch", al_etch),
            ("nb_patch", nb_patch),
            ("nb_etch", nb_etch),
        ]

        lay_num_dt_to_name_map: dict[tuple[int, int], str] = {}

        for name, lay in temp_init_loop_vars:
            if lay is None:
                def_layer_num, def_layer_dtype = DEFAULT_ANL_LAYERS[name]
                object.__setattr__(
                    self,
                    name,
                    Layer(name, def_layer_num, def_layer_dtype),
                )
                lay_num_dt_to_name_map.update({(def_layer_num, def_layer_dtype): name})
            elif isinstance(lay, Layer):
                lay_num_dt_to_name_map = self._check_layer_clash(lay, lay_num_dt_to_name_map)
                object.__setattr__(self, name, lay)
            else:
                raise TypeError(
                    f"Failed to create ANLMaskLayerSet with {name} set to type {type(lay)}. Should be of type `Layer` or None for defualt values."
                )

        object.__setattr__(self, "_layer_num_dt_to_name_map", lay_num_dt_to_name_map)

    # def _check_layer_clash(self, layer: Layer, lays_already_in: set[tuple[int, int]]):
    def _check_layer_clash(self, layer: Layer, lay_num_dt_to_name_map: dict[tuple[int, int], str]):
        if (layer.number, layer.datatype) in lay_num_dt_to_name_map.keys():
            raise ValueError(f"Cannot create ANLMaskLayerSet with layer: {layer} as it collides with another layer.")

        return lay_num_dt_to_name_map.update({(layer.number, layer.datatype): layer.name})

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        string = "ANLMaskLayerSet:\n"
        for field in fields(self):
            string += f"    {field.name}:\n"
            string += f"    ├─  number = {getattr(self, field.name).number}\n"
            string += f"    └─  dtype  = {getattr(self, field.name).datatype}\n"

        return string

    def get_layer_name_from_number(self, layer_number: int, datatype: int = 0, not_found_return_val: str | None = None) -> str | None:
        """Get the layer name for a given layer_number. If the layer number
        does not exist None will be returned.

        Parameters
        ----------
        layer_number: int
            The layer number who's name to get.

        KwArgs
        ------
        datatype: int = 0
            The datatype for the layer number. This is Default 0.

        not_found_return_val: str | None = None
            This is a value to return if that layer name cannot be found from
            the layer number provided. This defaults to None but can optionally
            be any other string.

        Returns
        -------
        layer_name: str | None
            The layer name if the layer number exists else None.
        """

        if not isinstance(layer_number, int):
            raise TypeError(f"layer_number should be of type int not {type(layer_number)}")
        if not isinstance(datatype, int):
            raise TypeError(f"datatype should be of type int not {type(datatype)}")

        layer_name: str | None = self._layer_num_dt_to_name_map.get((layer_number, datatype))
        if layer_name is None:
            return not_found_return_val

        return layer_name
