#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Config file reader blatantly copied from ASE

@author: roncofaber
"""

import configparser
import os
from platformdirs import user_config_dir

def load_config():
    config_dir = user_config_dir("mdinterface")
    user_config_file = os.path.join(config_dir, "config.ini")
    os.environ.setdefault("MDINT_CONFIG_DIR", config_dir)

    if os.path.exists(user_config_file):
        config_file = user_config_file
    else:
        config_file = os.path.join(os.path.dirname(__file__), "config.ini")

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)

    if "settings" in config:
        for key, value in config["settings"].items():
            os.environ.setdefault(key, value)

    return config

if __name__ == "__main__":
    load_config()
