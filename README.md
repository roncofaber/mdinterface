# mdinterface: Build Interface Systems for Molecular Dynamics Simulations

`mdinterface` is a Python package designed to build systems for Molecular Dynamics (MD) simulations. Initially developed to construct electrolyte/electrode solid-liquid interfaces, it is also well-suited for generating MD boxes of liquids, electrolyte systems, and polymer networks.

## Features

Using `mdinterface` you can:

- Create and configure molecular systems with solvents, solutes, and interfaces.
- Populate your system with ions and solvents using PACKMOL.
- Provide the starting concentration profile of the species to get the MD where you want it to be, but faster.
- Automatically write LAMMPS data files and coefficients.
- Import common molecules/metals/polymers with pre-defined classical force fields parameters in the database.
- Generate polymer chains of any length from a starting monomer.

## Requirements

Check the file [requirements.txt](requirements.txt) to see which packages are needed. Installing the package using `pip` should already take care of all dependencies.

Additionaly, you need the `packmol` utility installed to generate MD boxes. You can follow the instructions at [https://m3g.github.io/packmol/](https://m3g.github.io/packmol/) to install it. Alternatively, you can install it using `conda`:

```bash
conda install -c conda-forge packmol
```

## Installation

### Install using `pip`

You can simply install the latest release of the package and all dependencies using:

```bash
pip install mdinterface
```

### Install directly the source code

Alternatively you can obtain `mdinterface` directly from the repository by following these steps:

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

## Usage

Check the files in [examples](mdinterface/examples/) to learn how to use the package.