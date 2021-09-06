Parallel algorithms
-------------------

LAMMPS is from ground up designed to be running in parallel using the
MPI standard with distributed data via domain decomposition.  The
parallelization has to be efficient to enable good strong scaling (=
good speedup for the same system) and good weak scaling (= the
computational cost of enlarging the system is proportional to the system
size).  Additional parallelization using GPUs or OpenMP can then be
applied within the sub-domain assigned to an MPI process.  For clarity,
most of the following illustrations show the 2d simulation case. The
underlying algorithms in those cases, however, apply to both 2d and 3d
cases equally well.

.. toctree::
   :maxdepth: 1

   Developer_par_part
   Developer_par_comm
   Developer_par_neigh
   Developer_par_long
   Developer_par_openmp
