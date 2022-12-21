This directory contains the lepton library from the OpenMM software
which allows to efficiently evaluate mathematical expressions from
strings.  This library is used with the LEPTON package that support
force styles within LAMMPS that make use of this library.

You can type "make lib-lepton" from the src directory to see help on how
to build this library via make commands, or you can do the same thing
by typing "python Install.py" from within this directory, or you can
do it manually by following the instructions below.

---------------------

Lepton (short for “lightweight expression parser”) is a C++ library for
parsing, evaluating, differentiating, and analyzing mathematical
expressions. It takes expressions in the form of text strings, then
converts them to an internal representation suitable for evaluation or
analysis. Here are some of its major features:

- Support for a large number of standard mathematical functions and operations.
- Support for user defined custom functions.
- A variety of optimizations for automatically simplifying expressions.
- Computing analytic derivatives.
- Representing parsed expressions in two different forms (tree or program) suitable for
  further analysis or processing.

Lepton was originally created for use in the [OpenMM project](https://openmm.org)
ch5md is developed by Pierre de Buyl and is released under the 3-clause BSD
license that can be found in the file LICENSE.

To use the h5md dump style in lammps, execute
make -f Makefile.h5cc
in this directory then
make yes-h5md
in the src directory of LAMMPS to rebuild LAMMPS.

Note that you must have the h5cc compiler installed to use
Makefile.h5cc.  It should be part

If HDF5 is not in a standard system location, edit Makefile.lammps accordingly.

In the case of 2015 and more recent debian and ubuntu systems where concurrent
serial and mpi are possible, use the full platform depedent path, i.e.
`HDF5_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial`
