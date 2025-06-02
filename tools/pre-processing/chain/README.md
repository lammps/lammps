
# Introduction

The chain tool is a simple fortran code that produces polymer chains with a
single monomer type per chain. The input is read through the standard input
and the output is written to the standard output. The code uses LAMMPS data
file format, `molecular` atom style and `lj` units.

# Compiling

Compile using `make` in the directory. Requires the gfortran compiler. Edit the
Makefile otherwise.

`make clean` will clean the compilation file and executable.

# Usage

You can use the code by providing the input files through the standard input:
`./chain < def.chain > data.lmp`

Two example files are provided:

* `def.chain` defines a system with a single chain type.
* `def.chain.ab` defines a system with 2 chain types each with their own
  caracteristics.

See each input file for comments.
