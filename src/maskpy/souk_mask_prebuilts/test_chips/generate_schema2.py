from __future__ import annotations

import json
from enum import Enum, StrEnum
from pydoc import locate
from typing import Any, Callable, Literal, NotRequired, Required, TypedDict, get_type_hints

from maskpy.souk_mask_configs import SoukMaskConfig, get_mask_default_config
from maskpy.souk_mask_prebuilts.test_chips.types import FilterBankStructureSettings, QuadSettings, SoukResonatorType, TestChipSettings


class NumberValidator(Enum):
    """Validate types for any number"""

    fle = "float value less than or equal to"
    flt = "float value less than"
    fge = "float value greater than or equal to"
    fgt = "float value greater than"
    _float = "any float"
    ile = "integer value less than or equal to"
    ilt = "integer value less than"
    ige = "integer value greater than or equal to"
    igt = "integer value greater than"
    _int = "any integer"


class SchemaObj:
    """Schema object using draft-07.
    http://json-schema.org/draft-07/schema#

    References
    ----------
    https://datatracker.ietf.org/doc/html/draft-wright-json-schema-validation-00#section-5
    """

    def __init__(
        self,
        name: str,
        ref: Any | None = None,
        additionalItems: bool | SchemaObj | None = None,
        additionalProperties: bool | SchemaObj | None = None,
        allOf: set[SchemaObj] | None = None,
        anyOf: set[SchemaObj] | None = None,
        const: Any | None = None,
        contains: SchemaObj | None = None,
        default: Any | None = None,
        dependencies: SchemaObj | None = None,
        desctiption: str | None = None,
        _else: SchemaObj | None = None,
        enum: set[Any] | None = None,
        exclusiveMaximum: bool | None = None,
        exclusiveMinimum: bool | None = None,
        format: str | None = None,
        _if: SchemaObj | None = None,
        items: SchemaObj | None = None,
        maxItems: int | None = None,
        maxLength: int | None = None,
        maxProperties: int | None = None,
        maximum: float | int | None = None,
        minItems: int | None = None,
        minLength: int | None = None,
        minProperties: int | None = None,
        minimum: float | int | None = None,
        multipleOf: int | None = None,
        _not: SchemaObj | None = None,
        oneOf: set[SchemaObj] | None = None,
        pattern: str | None = None,
        patternProperties: SchemaObj | None = None,
        properties: SchemaObj | None = None,
        propertyNames: SchemaObj | None = None,
        required: set[str] | None = None,
        then: SchemaObj | None = None,
        title: str | None = None,
        _type: SchemaType | set[str] | None = None,
        uniqueItems: bool | None = None,
    ):
        self.name: str = name
        self.ref = ref
        self.additionalItems = additionalItems
        self.additionalProperties = additionalProperties
        self.allOf = allOf
        self.anyOf = anyOf
        self.const = const
        self.contains = contains
        self.default = default
        self.dependencies = dependencies
        self.desctiption = desctiption
        self._else = _else
        self.enum = enum
        self.exclusiveMaximum = exclusiveMaximum
        self.exclusiveMinimum = exclusiveMinimum
        self.format = format
        self._if = _if
        self.items = items
        self.maxItems = maxItems
        self.maxLength = maxLength
        self.maxProperties = maxProperties
        self.maximum = maximum
        self.minItems = minItems
        self.minLength = minLength
        self.minProperties = minProperties
        self.minimum = minimum
        self.multipleOf = multipleOf
        self._not = _not
        self.oneOf = oneOf
        self.pattern = pattern
        self.patternProperties = patternProperties
        self.properties = properties
        self.propertyNames = propertyNames
        self.required = required
        self.then = then
        self.title = title
        self._type = _type
        self.uniqueItems = uniqueItems

    def _type_error(
        self,
        prop_name: str,
        correct_type: type | list[type] | str | list[str],
        current_type: type,
    ) -> None:
        """raise type error for property"""
        raise (TypeError(f"{prop_name} is not of correct type. Should be {correct_type}. Not {current_type}"))

    def validate_number(
        self,
        prop_name: str,
        prop_value: Any,
        validator: NumberValidator,
        limit: float | int | None = None,
    ) -> bool:
        """Validate a number fits some validator constraints."""
        if not isinstance(prop_value, (float | int)):
            self._type_error(prop_name, [float, int], type(prop_value))

        if limit is None:
            if not (validator == NumberValidator._float or validator == NumberValidator._int):
                raise TypeError(
                    "Either limit should be a float or an int, or validator should be anything but NumberValidator._float or NumberValidator._int."
                )
            else:
                return True

        valid: bool = False
        match validator:
            case NumberValidator.fle:
                valid = (prop_value <= limit) and isinstance(prop_value, float)
            case NumberValidator.flt:
                valid = (prop_value < limit) and isinstance(prop_value, float)
            case NumberValidator.fge:
                valid = (prop_value >= limit) and isinstance(prop_value, float)
            case NumberValidator.fgt:
                valid = (prop_value > limit) and isinstance(prop_value, float)
            case NumberValidator.ile:
                valid = (prop_value <= limit) and isinstance(prop_value, int)
            case NumberValidator.ilt:
                valid = (prop_value < limit) and isinstance(prop_value, int)
            case NumberValidator.ige:
                valid = (prop_value >= limit) and isinstance(prop_value, int)
            case NumberValidator.igt:
                valid = (prop_value > limit) and isinstance(prop_value, int)

        if valid:
            return True
        else:
            raise ValueError(f"{prop_name} should be {validator.value} {limit if limit else ''}. Current value is {prop_value}")

    def as_dict(self) -> dict[str, Any]:
        """."""
        schema_object_dict: dict[str, Any] = {}

        if self.ref is not None:
            schema_object_dict["$ref"] = f"#/definitions/{self.ref}"

        if self.additionalItems is not None:
            if isinstance(self.additionalItems, bool):
                schema_object_dict["additionalItems"] = self.additionalItems
            elif isinstance(self.additionalItems, SchemaObj):
                schema_object_dict["additionalItems"] = self.additionalItems.as_dict()
            else:
                self._type_error("additionalItems", [bool, SchemaObj], type(self.additionalItems))

        if self.additionalProperties:
            if isinstance(self.additionalProperties, bool):
                schema_object_dict["additionalProperties"] = self.additionalProperties
            elif isinstance(self.additionalProperties, SchemaObj):
                schema_object_dict["additionalProperties"] = self.additionalProperties.as_dict()
            else:
                self._type_error("additionalProperties", [bool, SchemaObj], type(self.additionalProperties))

        if self.allOf:
            items: list[dict[str, Any]] = []
            for item in self.allOf:
                items.append(item.as_dict())
            schema_object_dict["allOf"] = items

        if self.anyOf:
            items: list[dict[str, Any]] = []
            for item in self.anyOf:
                items.append(item.as_dict())
            schema_object_dict["anyOf"] = items

        if self.const:
            schema_object_dict["const"] = self.const

        if self.contains:
            schema_object_dict["contains"] = self.contains.as_dict()

        if self.default:
            schema_object_dict["default"] = self.default

        if self.dependencies:
            schema_object_dict["dependencies"] = self.dependencies.as_dict()

        if self.desctiption:
            schema_object_dict["desctiption"] = self.desctiption

        if self._else:
            schema_object_dict["else"] = self._else.as_dict()

        if self.enum:
            items: list[dict[str, Any]] = []
            for item in self.enum:
                items.append(item)
            schema_object_dict["enum"] = items

        if self.exclusiveMaximum is not None:
            if isinstance(self.exclusiveMaximum, bool):
                schema_object_dict["exclusiveMaximum"] = self.exclusiveMaximum
            else:
                self._type_error("exclusiveMaximum", bool, type(self.exclusiveMaximum))

        if self.exclusiveMinimum is not None:
            if isinstance(self.exclusiveMinimum, bool):
                schema_object_dict["exclusiveMinimum"] = self.exclusiveMinimum
            else:
                self._type_error("exclusiveMinimum", bool, type(self.exclusiveMinimum))

        if self.format:
            schema_object_dict["format"] = self.format

        if self._if:
            schema_object_dict["if"] = self._if.as_dict()

        if self.items:
            if isinstance(self.items, SchemaObj):
                schema_object_dict["items"] = self.items.as_dict()
            else:
                self._type_error("items", SchemaObj, type(self.items))

        if self.maxItems:
            if self.validate_number(
                "maxItems",
                self.maxItems,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["maxItems"] = self.maxItems

        if self.maxLength:
            if self.validate_number(
                "maxLength",
                self.maxLength,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["maxLength"] = self.maxLength

        if self.maxProperties:
            if self.validate_number(
                "maxProperties",
                self.maxProperties,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["maxProperties"] = self.maxProperties

        if self.maximum:
            if self.validate_number(
                "maxProperties",
                self.maxProperties,
                NumberValidator._float,
            ):
                schema_object_dict["maximum"] = self.maximum
            else:
                self._type_error("maximum", [float, int], type(self.maximum))

        if self.minItems:
            if self.validate_number(
                "minItems",
                self.minItems,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["minItems"] = self.minItems

        if self.minLength:
            if self.validate_number(
                "minLength",
                self.minLength,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["minLength"] = self.minLength

        if self.minProperties:
            if self.validate_number(
                "minProperties",
                self.minProperties,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["minProperties"] = self.minProperties

        if self.minimum:
            if self.validate_number(
                "minimum",
                self.minimum,
                NumberValidator._float,
            ):
                schema_object_dict["minimum"] = self.minimum

        if self.multipleOf:
            if self.validate_number(
                "multipleOf",
                self.multipleOf,
                NumberValidator.ige,
                limit=0,
            ):
                schema_object_dict["minimum"] = self.minimum

        if self._not:
            if isinstance(self._not, SchemaObj):
                schema_object_dict["not"] = self._not.as_dict()
            else:
                self._type_error("_not", SchemaObj, type(self._not))

        if self.oneOf:
            if isinstance(self.oneOf, set):
                items: list[dict[str, Any]] = []
                for item in self.oneOf:
                    if isinstance(item, SchemaObj):
                        items.append(item.as_dict())
                    else:
                        self._type_error(f"oneOf[{item}]", SchemaObj, type(item))
                schema_object_dict["not"] = items
            else:
                self._type_error("oneOf", set[SchemaObj], type(self._not))

        if self.pattern:
            if isinstance(self.pattern, str):
                schema_object_dict["pattern"] = self.pattern
            else:
                self._type_error("pattern", str, type(self.pattern))

        if self.patternProperties:
            if isinstance(self.patternProperties, SchemaObj):
                schema_object_dict["patternProperties"] = self.patternProperties.as_dict()
            else:
                self._type_error("patternProperties", SchemaObj, type(self.patternProperties))

        if self.properties:
            if isinstance(self.properties, SchemaObj):
                schema_object_dict["properties"] = self.properties.as_dict()
            else:
                self._type_error("properties", SchemaObj, type(self.properties))

        if self.propertyNames:
            if isinstance(self.propertyNames, SchemaObj):
                schema_object_dict["propertyNames"] = self.propertyNames.as_dict()
            else:
                self._type_error("propertyNames", SchemaObj, type(self.propertyNames))

        if self.required:
            if isinstance(self.required, set):
                items: list[str] = []
                for item in self.required:
                    if isinstance(item, str):
                        items.append(item)
                    else:
                        self._type_error(f"required[{item}]", str, type(item))
                schema_object_dict["required"] = items
            else:
                self._type_error("required", set[SchemaObj], type(self.required))

        if self.then:
            if isinstance(self.then, SchemaObj):
                schema_object_dict["then"] = self.then.as_dict()
            else:
                self._type_error("then", SchemaObj, type(self.then))

        if self.title:
            if isinstance(self.title, str):
                schema_object_dict["title"] = self.title
            else:
                self._type_error("title", str, type(self.title))

        if self._type:
            if isinstance(self._type, SchemaType):
                schema_object_dict["type"] = self._type.value
            else:
                self._type_error("_type", SchemaType, type(self._type))

        if self.uniqueItems is not None:

            if isinstance(self.uniqueItems, bool):
                schema_object_dict["uniqueItems"] = self.uniqueItems
            else:
                self._type_error("uniqueItems", bool, type(self.uniqueItems))

        return schema_object_dict


class SchemaVersion(StrEnum):
    Draft07 = "http://json-schema.org/draft-07/schema#"


class SchemaType(StrEnum):
    type_string = "string"
    type_number = "number"
    type_integer = "integer"
    type_object = "object"
    type_array = "array"
    type_boolean = "boolean"
    type_null = "null"


class Schema:
    """."""

    def __init__(
        self,
        ref: str,
        definitions: SchemaObj | set(SchemaObj) | None = None,
        schema_version: SchemaVersion = SchemaVersion.Draft07,
        schema_objects: list[SchemaObj] | None = None,
    ):
        """."""
        self.schema_version = schema_version
        self.objects: list[SchemaObj] = []

        self.root_schema: dict[str, Any] = {
            "$schema": self.schema_version.value,
            "$ref": f"#/definitions/{ref}",
            "definitions": {},
        }

    def _check_is_schema_obj(self, schema_object):
        """."""
        if not isinstance(schema_object, SchemaObj):
            raise TypeError(f"schema_object should be of type SchemaObj. Not of type {type(schema_object)}")

    def add_schema_object(self, schema_object: SchemaObj) -> None:
        """Add an object to the schema.

        Parameters
        ----------
        schema_object: SchemaObj
            The SchemaObj to add to the Schema.
        """
        self._check_is_schema_obj(schema_object)
        self.objects.append(schema_object)

    def add_definition(self, definition: SchemaObj):
        """Add a definition to the schema.

        Parameters
        ----------
        definition: SchemaObj
            The definition SchemaObj to add to the Schema.
        """
        self._check_is_schema_obj(definition)
        self.root_schema["definition"][definition.name] = definition.as_dict()


def main():
    return


if __name__ == "__main__":
    main()
