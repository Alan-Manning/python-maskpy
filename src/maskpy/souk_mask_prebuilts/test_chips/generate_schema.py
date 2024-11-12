from __future__ import annotations

import json
from enum import Enum
from pydoc import locate
from typing import Any, Callable, Literal, NotRequired, Required, TypedDict, get_type_hints

from maskpy.souk_mask_configs import SoukMaskConfig, get_mask_default_config
from maskpy.souk_mask_prebuilts.test_chips.types import FilterBankStructureSettings, QuadSettings, SoukResonatorType, TestChipSettings


def _get_types_in_tuple(types_str: str) -> list[str]:
    """."""
    # strip "tuple[" and "]" from begining and end
    types_str = types_str[6:-1]

    count = 0
    indent_count = 0
    start = 0
    types_list: list[str] = []
    for i, char in enumerate(types_str):
        if char == "[":
            indent_count += 1
        if char == "]":
            indent_count -= 1
        if (char == ",") and (indent_count == 0):
            count += 1
            types_list.append(types_str[start:i])
            start = i + 2

    types_list.append(types_str[start:])
    count += 1

    return types_list


def _get_config_properties_requireds_extra_defs(definition_name: str) -> tuple[dict[str, dict[str, Any]], list[str], dict[str, str]]:
    """."""
    props, extra_defs = _get_confg_properties(definition_name)
    req = []
    return props, req, extra_defs


def _get_confg_properties(definition_name: str) -> tuple[dict[str, dict[str, Any]], set[str]]:
    """."""
    config_name = definition_name.split("_config_override")[0]

    souk_mask_config_values = [member.value for member in SoukMaskConfig]

    properties_dict: dict[str, dict[str, Any]] = {}
    extra_defs: dict[str, str] = {}

    basic_type_to_schema_type_dict: dict[str, str] = {
        "<class 'int'>": "number",
        "<class 'float'>": "number",
        "<class 'str'>": "string",
        "<class 'bool'>": "boolean",
    }

    if config_name not in souk_mask_config_values:
        raise LookupError(f"Couldnt find default config for {config_name}, from required_definition {definition_name}")

    default_config = get_mask_default_config(
        SoukMaskConfig(config_name),
    )

    for key, val in default_config.items():
        schema_type = basic_type_to_schema_type_dict.get(str(type(val)))
        properties_dict[key] = {
            "title": key,
            "type": schema_type,
        }
    return properties_dict, extra_defs


def _get_properties_requireds_extra_defs(
    obj: type | object,
) -> tuple[dict[str, dict[str, Any]], list[str], dict[str, str]]:
    """."""

    print(f"obj: {obj}")
    print(f"type(obj): {type(obj)}")

    stringed_typed_of_obj = str(type(obj))
    print(f"stringed_typed_of_obj: {stringed_typed_of_obj}")

    match stringed_typed_of_obj:
        case "<class 'typing._TypedDictMeta'>":
            props, extra_defs = _get_properties_for_typeddict(obj)
            req = _get_required_for_typeddict(obj)
        case "<class 'enum.EnumType'>":
            props, extra_defs = _get_properties_for_enum(obj)
            req = _get_required_for_enum(obj)
        case _:
            raise LookupError(f"unable to get properties and extra_defs for {obj} of type {stringed_typed_of_obj}")

    req = _get_required_for_typeddict(obj)

    return props, req, extra_defs


def _get_required_for_enum(obj: Enum) -> list[str]:
    return []


def _get_required_for_typeddict(type_dict: type) -> list[str]:
    """."""
    type_dict_name = str(type_dict).split("'")[1].split(".")[-1]

    match type_dict_name:
        case "TestChipSettings":
            required_keys = [
                "chip_id",
                "top_right_text",
                "top_left_text",
                "bottom_left_text",
                "add_ports",
                "quad_0_settings",
                "quad_1_settings",
                "quad_2_settings",
                "quad_3_settings",
            ]
        case "QuadSettings":
            required_keys = [
                "resonator_type",
                "add_antenna",
                "antenna_rotation",
                "text_under_quad",
                "couple_KID_to_ANT",
                "add_filter_bank",
                "add_top_choke_features",
                "add_bot_choke_features",
                "meander_materials",
                "IDC_and_frame_materials",
            ]
        case "FilterBankStructureSettings":
            required_keys = [
                "with_combiner",
                "with_crossover",
                "only_1_pol",
            ]
        case _:
            required_keys = []

    return required_keys


def _get_properties_for_enum(obj: Enum) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    """."""

    properties_dict: dict[str, dict[str, Any]] = {}
    extra_defs: dict[str, str] = {}

    basic_type_to_schema_type_dict: dict[str, str] = {
        "int": "integer",
        "float": "number",
        "str": "string",
        "bool": "boolean",
    }

    enum_members = [member.value for member in SoukMaskConfig]
    enum_members_types = set([type(el) for el in enum_members])
    if len(enum_members_types) == 1:
        schema_enum_type = enum_members_types.pop()
    if enum_members_types == set({int, float}):
        schema_enum_type = "number"
    else:
        schema_enum_type = "str"

    properties_dict[obj.__class__.__name__] = {
        "type": schema_enum_type,
        "enum": enum_members,
    }

    return properties_dict, extra_defs


def _get_properties_for_typeddict(type_dict: type) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    """."""
    basic_type_to_schema_type_dict: dict[str, str] = {
        "int": "integer",
        "float": "number",
        "str": "string",
        "bool": "boolean",
    }

    properties_dict: dict[str, dict[str, Any]] = {}
    # extra_defs: set[str] = set()
    extra_defs: dict[str, str] = {}

    for prop_name, raw_prop_type in get_type_hints(type_dict).items():

        properties_dict[prop_name] = {
            "title": prop_name,
        }
        type_string = str(raw_prop_type)

        schema_type = None
        schema_ref = None
        schema_items = None
        schema_max_items = None
        schema_min_items = None

        match type_string:
            case s if s.startswith("<class"):
                type_string = type_string.split("'")[1]
                schema_type = basic_type_to_schema_type_dict.get(type_string)

                # this means it is another TypedDict
                if schema_type is None:
                    type_dict_name = s.split("'")[1].split(".")[-1]
                    schema_ref = f"#/definitions/{type_dict_name}"
                    # extra_defs.add(type_dict_name)
                    extra_defs[type_dict_name] = type_string

            case s if s.startswith("<enum"):
                type_string = type_string.split("'")[1]
                schema_ref = f"#/definitions/{type_string}"
                # extra_defs.add(type_string)
                extra_defs[type_string] = type_string

            case s if s.startswith("typing.Union"):
                schema_ref = f"#/definitions/{prop_name}"
                # extra_defs.add(prop_name)
                extra_defs[prop_name] = s

            case s if s.startswith("dict["):
                schema_ref = f"#/definitions/{prop_name}"
                # extra_defs.add(prop_name)
                extra_defs[prop_name] = s

            case s if s.startswith("tuple["):
                tuple_types = _get_types_in_tuple(type_string)
                schema_type = "array"
                schema_items = {"$ref": f"#/definitions/{prop_name}"}
                # extra_defs.add(prop_name)
                extra_defs[prop_name] = s
                schema_min_items = len(tuple_types)
                schema_max_items = len(tuple_types)

            case _:
                print(f"{prop_name}: {raw_prop_type}")
                print(f"{prop_name}: {type_string}")
                raise RuntimeError("Could not parse TypedDict to schema")

        if schema_type:
            properties_dict[prop_name]["type"] = schema_type
        if schema_ref:
            properties_dict[prop_name]["$ref"] = schema_ref
        if schema_items:
            properties_dict[prop_name]["items"] = schema_items
        if schema_min_items:
            properties_dict[prop_name]["minItems"] = schema_min_items
        if schema_max_items:
            properties_dict[prop_name]["maxItems"] = schema_max_items

    return properties_dict, extra_defs


def _test_print_schema(schema_dict, indent_level=1, first_call=True):
    """."""
    tab = 4 * " "
    if first_call:
        print("{")

    for item_no, (key, val) in enumerate(schema_dict.items()):
        if isinstance(val, dict):
            print(f'{tab*indent_level}"{key}": ' + "{")
            _test_print_schema(val, indent_level=indent_level + 1, first_call=False)
        else:
            if isinstance(val, list):
                print(f'{tab*indent_level}"{key}": ' + "[")
                for el_no, el in enumerate(val):
                    if isinstance(val, dict):
                        _test_print_schema(val, indent_level=indent_level + 1, first_call=False)
                    else:
                        print(f'{tab*(indent_level+1)}"{el}",')

                print(f"{tab*indent_level}],")
            else:
                print(f'{tab*indent_level}"{key}": "{val}",')

    print(f"{tab*( indent_level - 1 )}" + "},")


def generate_schema() -> None:
    """Generates a schema for help in defining TestChipSettings in a
    yaml/json file."""

    main_schema: dict[str, Any] = {}

    schema_type = "http://json-schema.org/draft-07/schema#"
    main_ref = "#/definitions/TestChipsSettings"

    main_schema["$schema"] = schema_type
    main_schema["$ref"] = main_ref

    definitions: dict[str, dict[str, Any]] = {}

    test_chips_settings: dict[str, Any] = {
        "type": "object",
        "additionalProperties": True,
        "patternProperties": {
            "^chip_\\d+$": {
                "$ref": "#/definitions/TestChipSettings",
            },
        },
    }
    definitions["TestChipsSettings"] = test_chips_settings

    (properties, required_keys, extra_definitions_required) = _get_properties_requireds_extra_defs(TestChipSettings)

    test_chip_settings: dict[str, Any] = {
        "type": "object",
        "additionalProperties": False,
        "properties": properties,
        "required": required_keys,
    }
    definitions["TestChipSettings"] = test_chip_settings

    print("extra_definitions_required:")
    print(extra_definitions_required)

    # more_extra_definitions_required: dict[str, str] = {}

    # if no extra_definitions_required were done
    if extra_definitions_required:
        finished_with_extra_definitions = False
    else:
        finished_with_extra_definitions = True

    print(f"finished_with_extra_definitions: {finished_with_extra_definitions}")

    while not finished_with_extra_definitions:
        print("in while")

        # clear any defs from past run
        more_extra_definitions_required: dict[str, str] = {}

        for extra_req_name, extra_req_type in extra_definitions_required.items():
            print(f"extra_req_name: {extra_req_name}")
            print(f"extra_req_type: {extra_req_type}")

            if extra_req_name.__contains__("config"):
                if extra_req_name == "resonator_config_overrides":
                    properties: dict[str, dict[str, Any]] = {}
                    required_keys: list[str] = []
                    inner_extra_definitions_required: dict[str, str] = {}
                else:
                    (
                        properties,
                        required_keys,
                        inner_extra_definitions_required,
                    ) = _get_config_properties_requireds_extra_defs(extra_req_name)

            elif extra_req_name.__contains__("SoukResonatorType"):
                (
                    properties,
                    required_keys,
                    inner_extra_definitions_required,
                ) = _get_properties_requireds_extra_defs(SoukResonatorType)

            elif extra_req_name.__contains__("meander_materials"):
                properties: dict[str, dict[str, Any]] = {
                    "meander_materials": {
                        "title": "meander_materials",
                        "type": "string",
                        "enum": ["Al", "Nb"],
                    }
                }
                required_keys: list[str] = []
                inner_extra_definitions_required: dict[str, str] = {}

            elif extra_req_name.__contains__("IDC_and_frame_materials"):
                properties: dict[str, dict[str, Any]] = {
                    "IDC_and_frame_materials": {
                        "title": "IDC_and_frame_materials",
                        "type": "string",
                        "enum": ["Al", "Nb"],
                    }
                }
                required_keys: list[str] = []
                inner_extra_definitions_required: dict[str, str] = {}

            elif extra_req_name.__contains__("mux_func_overrides"):
                properties: dict[str, dict[str, Any]] = {
                    "mux_func_overrides": {
                        "title": "mux_func_overrides",
                        "type": "array",
                        "items": {"type": "string"},
                        "minItems": 4,
                        "maxItems": 4,
                    }
                }
                required_keys: list[str] = []
                inner_extra_definitions_required: dict[str, str] = {}

            elif extra_req_name.__contains__("trim_lengths"):
                properties: dict[str, dict[str, Any]] = {
                    "trim_lengths": {
                        "title": "trim_lengths",
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 4,
                        "maxItems": 4,
                    }
                }
                required_keys: list[str] = []
                inner_extra_definitions_required: dict[str, str] = {}

            elif extra_req_name.__contains__("add_filter_bank"):
                (
                    properties_FilterBankStructureSettings,
                    required_keys_FilterBankStructureSettings,
                    inner_extra_definitions_required,
                ) = _get_properties_requireds_extra_defs(FilterBankStructureSettings)

                properties: dict[str, dict[str, Any]] = {
                    "add_filter_bank": {
                        "oneOf": [
                            {
                                "type": "object",
                                "additionalProperties": False,
                                "properties": properties_FilterBankStructureSettings,
                                "required": required_keys_FilterBankStructureSettings,
                            },
                            {
                                "type": "boolean",
                                "enum": [False],
                            },
                        ],
                    }
                }

            elif extra_req_type.__contains__("tuple") and not extra_req_name.__contains__("config"):
                properties: dict[str, dict[str, Any]] = {}
                required_keys: list[str] = []
                inner_extra_definitions_required: dict[str, str] = {}

            else:
                prop_type = locate(extra_req_type)
                (
                    properties,
                    required_keys,
                    inner_extra_definitions_required,
                ) = _get_properties_requireds_extra_defs(prop_type)

            for k, v in inner_extra_definitions_required.items():
                more_extra_definitions_required[k] = v

            extra_definitions: dict[str, Any] = {
                "type": "object",
                "additionalProperties": False,
                "properties": properties,
                "required": required_keys,
            }

            definitions[extra_req_name] = extra_definitions

        print(f"more_extra_definitions_required: {more_extra_definitions_required}")

        # if no extra defs required then just break out and finish up
        if not more_extra_definitions_required:
            finished_with_extra_definitions = True
            break

        # set any defs from current run for next
        extra_definitions_required: dict[str, str] = {}
        for k, v in more_extra_definitions_required.items():
            extra_definitions_required[k] = v

    main_schema["definitions"] = definitions

    _test_print_schema(main_schema)

    # result = json.dumps(main_schema)

    with open("schema.json", "w", encoding="utf-8") as f:
        json.dump(main_schema, f, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    generate_schema()
