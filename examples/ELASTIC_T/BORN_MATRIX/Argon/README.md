This repository is a test case for the compute born/matrix. It provides short
scripts creating argon fcc crystal and computing the Born term using the two
methods described in the documentation.

In the __Analytical__ directory the terms are computed using the analytical
derivation of the Born term for the lj/cut pair style.

In the __Numdiff__ directory, the Born term is evaluated through small
numerical differences of the stress tensor. This method can be used with any
interaction potential.

Both script show examples on how to compute the full Cij elastic stiffness
tensor in LAMMPS.
