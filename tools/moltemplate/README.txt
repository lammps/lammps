
Moltemplate
===========

##  Description

This folder used to contain a distribution of Moltemplate, a general purpose,
cross-platform, text-based molecule and topology builder for LAMMPS.
Moltemplate was originally conceived for building custom coarse-grained
molecular models, but it has since been generalized for all-atom simulations
as well.  It currently supports the OPLS, COMPASS, AMBER(GAFF,GAFF2),
MARTINI, SPICA(SDK), LOPLS(2015), and TraPPE(1998) force fields
(new force fields and examples are added continually through user
contributions).

Moltemplate is now distributed as a Python package through PyPI and
can also be downloaded from https://github.com/jewettaij/moltemplate/releases
The most up-to-date version is usually available through GitHub.

## Typical usage

    moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd] system.lt

## Web page

Documentation, examples, and supporting code can be found at:

https://moltemplate.org

## Tutorial files

The folder tutorial-files contains the files used in the "Moltemplate Tutorial"
in the LAMMPS manual (https://docs.lammps.org/Howto_moltemplate.html).

## Requirements

Moltemplate requires a Bourne-compatible shell (e.g. bash) and Python 3.
It runs on Linux and macOS, and on Windows within the Windows Subsystem for
Linux (WSL, see https://docs.lammps.org/Howto_wsl.html).

## Installation

The recommended way is to install moltemplate with pip into a Python
virtual environment:

    python3 -m venv $HOME/moltemplate-env
    source $HOME/moltemplate-env/bin/activate
    pip install moltemplate

The virtual environment has to be activated (with the "source" command
above) every time before using moltemplate.  To install a version
downloaded from GitHub instead, unpack the archive and run "pip install ."
in its top-level folder (with the virtual environment active).

Later, you can uninstall moltemplate using:

    pip uninstall moltemplate

Please see the moltemplate documentation for alternative ways of installing it.
