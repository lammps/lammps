This directory shows how to use the `fix born` command
to calculate the full matrix of elastic constants
for cubic diamond at finite temperature
running the Stillinger-Weber potential.

The input script `in.elastic` can be run
directly from LAMMPS, or via a Python wrapper
script.

to run directly from LAMMPS, use:

mpirun -np 4 lmp.exe -in in.elastic

This simulates an orthorhombic box with the cubic crystal axes
aligned with x, y, and z.
The default settings replicate the 1477~K benchmark of
Kluge, Ray, and Rahman (1986) that is Ref.[15] in:
Y. Zhen, C. Chu, Computer Physics Communications 183(2012) 261-265

The script contains many adjustable parameters that can be used
to generate different crystal structures, supercell sizes,
and sampling rates.

to run via the Python wrapper, use:

mpirun -np 4 python elastic.py

This will first run the orthorhombic supercell as before,
follows by an equivalent simulation using a triclinic structure.
The script shows how the standard triclinic primitive cell for cubic diamond
can be rotated in to the LAMMPS upper triangular frame. The resultant
elastic constant matrix does not exhibit the standard symmetries of cubic crystals.
However, the matrix is then rotated back to the standard orientation
to recover the cubic symmetry form of the elastic matrix,
resulting in elastic constants that are the same for both
simulations, modulo statistical uncertainty.

