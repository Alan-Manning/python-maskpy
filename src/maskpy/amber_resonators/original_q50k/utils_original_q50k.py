from .config_original_q50k import get_resonator_config


def _get_config_checking_override(config_override: dict[str, float | int] | None) -> dict[str, float | int]:

    if config_override is not None:
        if isinstance(config_override, dict):
            return config_override
        else:
            print(f"\033[93mWarning: config_override not of correct format. Defaulting to using base config\033[0m")
            return get_resonator_config()
    else:
        return get_resonator_config()
