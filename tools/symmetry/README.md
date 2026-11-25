# find-symmetry

The `find-symmetry` program writes symmetry files for the LAMMPS
`fix symmetry` command.  It uses the [spglib](https://spglib.readthedocs.io/)
library, which is downloaded and compiled automatically.

## Building

The tool is built as part of LAMMPS when configuring with CMake using
`-D BUILD_TOOLS=on`.  It can also be compiled on its own:

```bash
cmake -S tools/symmetry -B build-symmetry
cmake --build build-symmetry
```

## Usage

```
find-symmetry [options] <data file>
find-symmetry [options] -g <space group number>
```

Given a LAMMPS data file, e.g. written by the `write_data` command, the
space group of the atoms in the simulation box is determined and its
operators are written in fractional coordinates of the box.  When the box
consists of several unit cells, the operators of one unit cell are written
together with a "supercell" entry.  If the positions deviate from exact
symmetry by more than the default tolerance of fix symmetry, a suitable
value for its `tol` keyword is printed.

With `-g` and a space group number from 1 to 230, the operators of the
space group in its standard setting are written.  `-H` selects a different
setting by its Hall number from 1 to 530.  The simulation box must then be
the conventional unit cell of that setting or, with `-s nx ny nz`, a box
of nx x ny x nz such cells.

Options:

- `-g <number>` space group number from 1 to 230
- `-H <number>` Hall number from 1 to 530
- `-s <nx> <ny> <nz>` with -g or -H: add a supercell entry
- `-t <distance>` tolerance for finding the symmetry in distance units of
  the data file (default 1.0e-5)
- `-n` do not reduce the operators to a unit cell; this can produce a very
  large number of operators
- `-o <file>` write to file instead of the screen
- `-h` print help

The atom types are used to distinguish atoms.  Data files with type labels
are supported.  Data files for general triclinic boxes are not supported.
