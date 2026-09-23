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

Install the [mdinterface-compatible LigParGen fork](https://github.com/roncofaber/ligpargen) in the same environment, then verify that its command is available:

```bash
python -m pip install "git+https://github.com/roncofaber/ligpargen.git"
conda install -c conda-forge openbabel
ligpargen -h
obabel -V
```

Open Babel is required because `mdinterface` passes molecules to LigParGen as XYZ files.

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

The [boss-container](https://github.com/roncofaber/boss-container) repository contains only build recipes. Each user must supply a licensed BOSS installation when building an image. A built Docker image, exported archive, or Apptainer file contains BOSS and must remain private unless the license explicitly permits sharing it.

The container supplies the 32-bit userspace libraries that are commonly missing on modern x86-64 Linux systems. Because Docker and Apptainer share the host kernel, they do not by themselves solve an unsupported CPU architecture, a kernel without the 32-bit x86 ABI, or a security policy that blocks the required system calls.

### Full parameterization development environment

Use the checked-in environment definition when developing LigParGen or the planned OpenFF integration:

```bash
mamba env create -f environment-full.yml
mamba activate mdinterface-full
```

This installs Python 3.12, the local mdinterface checkout, the mdinterface-compatible LigParGen fork, Open Babel, AmberTools, OpenFF Toolkit and Interchange, CPU-only NAGL and PyTorch, PACKMOL, the test suite, and documentation tooling. The full environment uses Python 3.12 and NumPy 1.x to satisfy the current AmberTools dependency stack, while the core mdinterface package continues to support Python 3.10-3.14. BOSS is not included and must still be configured through `BOSSdir`.

If LigParGen is cloned next to mdinterface, replace the installed Git version with the local checkout for development:

```bash
python -m pip install -e ../ligpargen
```

OpenFF 2.2.1 parameterization has been validated through AmberTools, and OpenFF 2.3 has been validated through CPU-only NAGL on Python 3.14. OpenFF-to-LAMMPS export and mdinterface import have also been validated. OpenFF is not yet a supported `Specie` backend because mdinterface must first preserve Fourier torsion styles, mixing rules, switching behavior, and 1-4 scaling metadata across assembly and output.

### RESP charges with PySCF

Install [PySCF](https://github.com/pyscf/pyscf) and [PyMBXAS](https://gitlab.com/roncofaber/pymbxas). RESP fitting currently requires [gpu4pyscf](https://github.com/pyscf/gpu4pyscf).

### AIMD with FAIRChem

```bash
pip install fairchem-core
```
