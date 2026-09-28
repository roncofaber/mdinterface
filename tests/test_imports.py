import os
import subprocess
import sys
import importlib

import pytest

from mdinterface.config import load_config


@pytest.mark.parametrize("module, name", [
    ("mdinterface", "SimulationBox"),
    ("mdinterface", "BoxBuilder"),
    ("mdinterface.build", "BoxBuilder"),
    ("mdinterface.build.builder", "BoxBuilder"),
])
def test_removed_builders_are_not_exported(module, name):
    assert not hasattr(importlib.import_module(module), name)


def test_legacy_builder_module_removed():
    with pytest.raises(ModuleNotFoundError, match="mdinterface.simulationbox"):
        importlib.import_module("mdinterface.simulationbox")


def test_package_import_is_silent_and_does_not_load_optional_modules(tmp_path):
    env = os.environ.copy()
    env["XDG_CONFIG_HOME"] = str(tmp_path)
    command = [
        sys.executable,
        "-c",
        (
            "import sys; import mdinterface; "
            "assert 'mdinterface.externals.aimd' not in sys.modules; "
            "assert 'matplotlib.pyplot' not in sys.modules"
        ),
    ]

    result = subprocess.run(command, env=env, capture_output=True, text=True, check=True)

    assert result.stdout == ""
    assert result.stderr == ""


def test_root_load_config_remains_available():
    import mdinterface

    assert callable(mdinterface.load_config)


def test_load_config_reads_user_file(tmp_path, monkeypatch):
    config_dir = tmp_path / "mdinterface"
    config_dir.mkdir()
    (config_dir / "config.ini").write_text("[settings]\nBOSSdir = configured-boss\n")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.delenv("BOSSdir", raising=False)
    monkeypatch.delenv("MDINT_CONFIG_DIR", raising=False)

    config = load_config()

    assert config["settings"]["BOSSdir"] == "configured-boss"
    assert os.environ["BOSSdir"] == "configured-boss"
    assert os.environ["MDINT_CONFIG_DIR"] == str(config_dir)


def test_load_config_preserves_environment_values(tmp_path, monkeypatch):
    config_dir = tmp_path / "mdinterface"
    config_dir.mkdir()
    (config_dir / "config.ini").write_text("[settings]\nBOSSdir = configured-boss\n")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setenv("BOSSdir", "environment-boss")

    load_config()

    assert os.environ["BOSSdir"] == "environment-boss"
