Unit test collection
====================

This part of the LAMMPS distribution contains scripts and inputs to run
unit tests for different parts of the package. The tests are organized
into folders based on which part of the distribution they are testing.


* c-library:     tests of the C-library interface
* commands:      tests for simple LAMMPS input commands
* cplusplus:     tests for using the LAMMPS library from C++
* force-styles:  tests for styles that compute or modify forces
                 i.e pair, bond, angle, kspace styles and some fixes.
* formats:       tests related to file formats: atom styles, dump styles
                 tokenizers potential file readers, dump styles
* fortran:       tests for the LAMMPS Fortran 03 module
* python:        tests for the LAMMPS Python module and python command
* testing:       convenience and support functions for test tools
* tools:         tests for tools form the tools folder
* utils:         tests for utility functions, e.g. from the utils namespace
