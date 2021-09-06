OpenMP Parallelism
^^^^^^^^^^^^^^^^^^

The styles in the INTEL, KOKKOS, and OPENMP package offer to use OpenMP
thread parallelism to predominantly distribute loops over local data
and thus follow an orthogonal parallelization strategy to the
decomposition into spatial domains used by the :doc:`MPI partitioning
<Developer_par_part>`.  For clarity, this section discusses only the
implementation in the OPENMP package as it is the simplest. The INTEL
and KOKKOS package offer additional options and are more complex since
they support more features and different hardware like co-processors
or GPUs.

Avoiding data races
-------------------

