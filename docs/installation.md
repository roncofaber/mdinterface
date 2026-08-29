# Installation

## Requirements

- **Python** 3.10-3.14

Core dependencies are handled automatically by `pip`, including the upstream PACKMOL package and executable. See [requirements.txt](https://github.com/roncofaber/mdinterface/blob/main/requirements.txt).

## Installing mdinterface

=== "Stable release"

    ```bash
    pip install mdinterface
    ```

=== "Development version"

    ```bash
    git clone https://github.com/roncofaber/mdinterface.git
    cd mdinterface
    pip install -e .
    ```

## Optional extras

```bash
pip install mdinterface[resp]   # RESP charge analysis (PySCF / gpu4pyscf)
pip install mdinterface[aimd]   # FAIRChem AIMD
pip install mdinterface[all]    # everything
```

### LigParGen (automatic OPLS-AA parameters)

Follow the instructions on the [LigParGen GitHub](https://github.com/Isra3l/ligpargen) (or try [this fork](https://github.com/roncofaber/ligpargen) if you hit installation issues).

LigParGen requires [BOSS](http://zarbi.chem.yale.edu/software.html), a 32-bit binary. Point `mdinterface` to it via `BOSSdir` in the config file. Three modes are supported depending on how BOSS is available:

The configuration file is read when LigParGen is invoked. An existing `BOSSdir` environment variable takes precedence over the file value.

=== "Apptainer / Singularity (HPC)"

    Build the container with [boss-container](https://github.com/roncofaber/boss-container), then:

    ```ini
    # ~/.config/mdinterface/config.ini
    [settings]
    BOSSdir = /path/to/boss-container.sif
    ```

=== "Docker (local)"

    Build the container with [boss-container](https://github.com/roncofaber/boss-container), then:

    ```ini
    # ~/.config/mdinterface/config.ini
    [settings]
    BOSSdir = boss-container:latest
    ```

=== "Native BOSS"

    Requires `csh` installed on the host and a working 32-bit BOSS binary:

    ```ini
    # ~/.config/mdinterface/config.ini
    [settings]
    BOSSdir = /path/to/boss
    ```

### RESP charges with PySCF

Install [PySCF](https://github.com/pyscf/pyscf) and [PyMBXAS](https://gitlab.com/roncofaber/pymbxas). RESP fitting currently requires [gpu4pyscf](https://github.com/pyscf/gpu4pyscf).

### AIMD with FAIRChem

```bash
pip install fairchem-core
```
