#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Config file reader blatantly copied from ASE

@author: roncofaber
"""

import os
import configparser
from platformdirs import user_config_dir

def load_config():
    # Determine the path to the user's configuration file
    MDINT_CONFIG_DIR = user_config_dir("mdinterface")
    user_config_file = MDINT_CONFIG_DIR + "/config.ini"
    
    os.environ["MDINT_CONFIG_DIR"] = MDINT_CONFIG_DIR
    
    # Check if the user's configuration file exists
    if os.path.exists(user_config_file):
        config_file = user_config_file
    else:
        # Fall back to the default configuration file in the package directory
        print(f"No config.ini file found:\n{user_config_file}")
        config_file = os.path.join(os.path.dirname(__file__), 'config.ini')
        print(f"Defaulting on reading file:\n{config_file}")
        print("Make sure it makes sense.")
    
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)
    
    if 'settings' in config:
        for key in config['settings']:
            os.environ[key] = config['settings'][key]

# Load the configuration when the module is imported
if __name__ == "__main__":
    load_config()
