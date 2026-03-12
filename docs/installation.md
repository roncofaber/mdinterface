# Installation

## Requirements

- **Python** 3.8+
- **PACKMOL** (see below)

Core Python dependencies are handled automatically by `pip`, see [requirements.txt](https://github.com/roncofaber/mdinterface/blob/main/requirements.txt).

## Installing PACKMOL

PACKMOL must be installed separately and available on your `PATH`:

```bash
conda install -c conda-forge packmol
```

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

Point `mdinterface` to your BOSS directory via a config file:

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
