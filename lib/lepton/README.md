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
- Support for just-in-time compilation via asmjit library on x86 (autodetected)
  This should make evaluation about 2 times faster

Lepton was originally created for use in the [OpenMM project](https://openmm.org)
