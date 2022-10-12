Parallel algorithms
-------------------

LAMMPS is designed to enable running simulations in parallel using the
MPI parallel communication standard with distributed data via domain
decomposition.  The parallelization aims to be efficient result in good
strong scaling (= good speedup for the same system) and good weak
scaling (= the computational cost of enlarging the system is
proportional to the system size).  Additional parallelization using GPUs
or OpenMP can also be applied within the sub-domain assigned to an MPI
process.  For clarity, most of the following illustrations show the 2d
simulation case. The underlying algorithms in those cases, however,
apply to both 2d and 3d cases equally well.

.. note::

   The text and most of the figures in this chapter were adapted
   for the manual from the section on parallel algorithms in the
   :ref:`new LAMMPS paper <lammps_paper>`.

.. toctree::
   :maxdepth: 1

   Developer_par_part
   Developer_par_comm
   Developer_par_neigh
   Developer_par_long
   Developer_par_openmp
