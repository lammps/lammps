The ``LIBLAMMPS`` Fortran Module
********************************

The ``LIBLAMMPS`` module provides an interface to call LAMMPS from a
Fortran code.  It is based on the LAMMPS C-library interface and
requires a Fortran 2003 compatible compiler to be compiled.

While C libraries have a defined binary interface (ABI) and can thus be
used from multiple compiler versions from different vendors for as long
as they are compatible with the hosting operating system, the same is
not true for Fortran codes.  Thus the LAMMPS Fortran module needs to be
compiled alongside the code using it from the source code in
``examples/COUPLE/fortran/lammps.f90``.  When linking, you also need to
:doc:`link to the LAMMPS library <Build_link>`.  A typical command line
for a simple program using the Fortran interface would be:

.. code-block:: bash

   mpifort -o testlib.x  lammps.f90 testlib.f90 -L. -llammps

Please note, that the MPI compiler wrapper is only required when the
calling the library from an MPI parallel code.

