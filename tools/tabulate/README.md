# Python scripts to generate tabulated potential files for LAMMPS

This directory contains a Python module 'tabulate' that can be used to
generate tabulated potential files for pair style table, bond style
table, and angle style table

To create tables, you need to define your energy and - optionally -
force functions and then an instance of either the PairTabulate(),
BondTabulate(), AngleTabulate(), DihedralTabulate(), or WallTabulate()
class from the tabulate module and call its run() method to generate the
table.  Most of the settings (number of points, minimum, maximum etc.)
are provided as command line flags.  The run() method may be called
multiple times to generate multiple tables, for instance after changing
parameters of the energy/force functions.  If the force function is not
provided, the force will be determined from the energy function through
numerical differentiation.

Please see the individual tabulation scripts in this folder for examples:

| File                          | Description                                                                   |
| ------------------------------|-------------------------------------------------------------------------------|
| pair_lj_tabulate.py           | creates two Lennard-Jones pair potential tables with different parameters     |
| bond_morse_tabulate.py        | creates a table for a Morse bond potential table                              |
| angle_harmonic_tabulate.py    | creates a table for a harmonic angle potential table                          |
| dihedral_harmonic_tabulate.py | creates a table for a harmonic dihedral potential table                       |
| pair_hybrid_tabulate.py       | creates a Morse/Lennard-Jones hybrid potential table with smooth switching    |
| wall_harmonic_tabulate.py     | creates a table for fix wall/table with a simple repulsive harmonic potential |
| wall_multi_tabulate.py        | creates a table for fix wall/table with multiple tables                       |
| pair_bi_tabulate.py           | creates a table from radial distribution file using Boltzmann Inversion       |

Common command line flags:

```
options:
  -h, --help                          show this help message and exit
  --num-points NUM, -n NUM            Number of tabulated points (default: 1000)
  --filename FILENAME, -f FILENAME    Name of output file (default: -)
  --diff-num, -d                      Differentiate energy function numerically
  --inner XMIN, -i XMIN               Inner cutoff of table (required for pair)
  --outer XMAX, -o XMAX               Outer cutoff of table (required)
```
