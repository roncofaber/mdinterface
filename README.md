# mdinterface: build interface systems for Molecular Dynamics simulations

`mdinterface` is a Python package for building systems to use in Molecular Dynamics (MD) simulations. The package was primarly developed to help construct electrolyte/electrode interfaces, but is well suited to generate suitable MD boxes of liquids, electrolyte systems and polymer networks.

## Features

- Create and configure molecular systems with solvents, solutes, and interfaces.
- Generate simulation boxes for MD simulations.
- Populate boxes with ions and solvents using PACKMOL.
- Automatically write LAMMPS data files and coefficients.
- Database of common liquids/gases with pre-defined classical force fields parameters.
- Generate polymer chains from given monomers.

## Requirements

Check the file [requirements.txt](requirements.txt) to see which packages are needed. Installing the package using `pip` should already take care of all dependencies.

## Installation

### Install using `pip`

You can simply install the latest release of the package and all dependencies using:

```bash
pip install mdinterface
```

### Install directly the source code

Alternatively you can obtain `mdinterface` directly from the repository by following those steps:

Clone the repository in the desired location:

```bash
git clone git@gitlab.com:roncofaber/mdinterface.git
```

Install the required packages:

```bash
cd mdinterface
conda install -c conda-forge --file requirements.txt
```

Install the package with pip:

```bash
pip install .
```

### Install a development environment

If you plan of making changes, clone the package and install the requirements but then add it to your development environment with:

```bash
pip install --no-build-isolation --no-deps -e .
```
