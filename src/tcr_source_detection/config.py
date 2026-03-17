import json
import os
from typing import Dict, Any, Optional

CONFIG_FILE = "tcr_configs.json"

def load_configs() -> Dict[str, Any]:
    """Loads all saved configurations from the JSON file."""
    if not os.path.exists(CONFIG_FILE):
        return {}
    try:
        with open(CONFIG_FILE, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading configs: {e}")
        return {}

def save_config(name: str, config_data: Dict[str, Any]):
    """Saves or updates a named configuration."""
    configs = load_configs()
    configs[name] = config_data
    try:
        with open(CONFIG_FILE, "w") as f:
            json.dump(configs, f, indent=4)
        print(f"Configuration '{name}' saved to {CONFIG_FILE}")
    except Exception as e:
        print(f"Error saving config: {e}")

def get_config(name: str) -> Optional[Dict[str, Any]]:
    """Retrieves a specific configuration by name."""
    configs = load_configs()
    return configs.get(name)

def list_configs():
    """Prints all available configuration names."""
    configs = load_configs()
    if not configs:
        print("No saved configurations found.")
    else:
        print("Available configurations:")
        for name in configs.keys():
            print(f"  - {name}")
