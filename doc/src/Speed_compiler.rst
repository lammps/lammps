Compiler optimization
=====================

The CMake presets that select a compiler (``gcc.cmake``,
``clang.cmake``, ``oneapi.cmake``) and the default settings of the
CMake build use moderate optimization (``-O2``) and generate code that
runs on any CPU of the same architecture, e.g. on any x86_64 CPU.
These settings are chosen for reliability: the LAMMPS executable runs
on every compatible CPU, and its results agree closely with the
reference data of the LAMMPS unit tests.

For users who want more speed and accept the limitations discussed
below, there are three presets with more aggressive optimization:
``gcc-opt.cmake``, ``clang-opt.cmake``, and ``oneapi-opt.cmake``.  Each
is used instead of the corresponding compiler preset, for example:

.. code-block:: bash

   cmake -S cmake -B build -C cmake/presets/gcc-opt.cmake -C cmake/presets/most.cmake

These presets select ``-O3`` and, on x86_64 CPUs, the x86-64-v3
instruction set (``-march=x86-64-v3``), which includes the AVX2 vector
instructions and fused multiply-add (FMA).  ``gcc-opt.cmake`` also adds
``-fno-semantic-interposition``.

Expected speedup
""""""""""""""""

The table below lists the average speedup with the three presets over
LAMMPS compiled with GCC and the default settings, measured on three
different CPUs, for the plain styles on one CPU core (first number) and
for the OPENMP package styles with 8 threads (second number).  The
averages are geometric means over 23 inputs (9 inputs with threads)
with pairwise, manybody, and ReaxFF potentials and with long-range or
cutoff Coulomb interactions.  The oneAPI compilers were tested only on
the Cascade Lake CPU.

.. table::

   +--------------------------------------+-----------+-----------+------------+
   | CPU                                  | GCC       | Clang     | oneAPI     |
   +======================================+===========+===========+============+
   | AMD Ryzen Threadripper PRO 7985WX    | 7% / 6%   | 5% / 5%   |            |
   | (Zen 4)                              |           |           |            |
   +--------------------------------------+-----------+-----------+------------+
   | Intel Xeon E5-2620 v4 (Broadwell)    | 5% / 7%   | 6% / 5%   |            |
   +--------------------------------------+-----------+-----------+------------+
   | Intel Xeon Gold 6248 (Cascade Lake)  | 10% / 9%  | 11% / 8%  | 10% / 8%   |
   +--------------------------------------+-----------+-----------+------------+

The compilers were GCC 15 and 16, Clang 21 and 22, and Intel oneAPI
2026.2.  For individual inputs, the speedup ranged from 3 percent
slower to 23 percent faster.  With their default settings, the three
compilers were within 2 percent of each other.  When running one copy
of the input on every CPU core at the same time, the speedup was about
the same (within 1 to 3 percent), except for the ReaxFF input, which is
then limited by the memory bandwidth and gains less.

Compiling for the specific CPU of the machine (``-march=native`` with
GCC and Clang or ``-xHost`` with the Intel compilers) instead of
x86-64-v3 added at most 1.5 percent on average on these CPUs, even on
those with AVX-512 instructions.  Most of the floating-point work of the plain
and OPENMP package styles is not vectorized (see the section on
floating-point operations on the :doc:`Benchmarks <Speed_bench>` page).
Also, an executable compiled for a specific CPU can stop with an
"illegal instruction" error on a different CPU, e.g. when LAMMPS was
compiled on the login node of a cluster with a newer CPU than its
compute nodes.  The presets therefore use x86-64-v3.

Limitations
"""""""""""

- An executable compiled with these presets runs only on CPUs that
  support the x86-64-v3 instruction set: most Intel CPUs since Haswell
  (2013) and AMD CPUs since Zen (2017).  On older CPUs it stops with an
  "illegal instruction" error.

- With FMA, a multiplication and an addition are combined into one
  operation with a single rounding step.  This changes results in the
  last digits compared to LAMMPS compiled with the default settings.
  Since MD trajectories are chaotic, such small differences grow over
  time, so trajectories from differently compiled executables
  eventually diverge, while averaged properties agree within their
  statistical uncertainty.  This is the same kind of difference as
  between different compilers or different numbers of MPI processes.
  Even with the default settings, the same executable can produce
  slightly different results on different CPUs: the math library of the
  GNU C library selects optimized versions of functions like ``exp()``,
  ``log()``, or ``pow()`` that use FMA instructions on CPUs that support
  them.

- The presets do not use options like ``-ffast-math`` that allow the
  compiler to ignore the IEEE 754 rules for floating-point math.  With
  the Intel compilers, which use such optimizations by default, this
  gained only about 2 percent and caused more than 20 unit tests to fail,
  some with NaN results.  Therefore, LAMMPS adds ``-fp-model precise``
  when compiling with the Intel compilers.

- With the oneAPI compilers, the :doc:`INTEL package <Speed_intel>`
  compiles LAMMPS for the CPU of the build machine (``-xHost``), unless
  a target CPU is already selected.  With ``oneapi-opt.cmake``, the
  INTEL package styles are thus compiled for x86-64-v3 and do not use
  AVX-512 instructions.

- Some Linux distributions change the defaults of their compilers.  For
  example, GCC in Red Hat Enterprise Linux 10 and compatible
  distributions generates code for x86-64-v3 by default, so that
  LAMMPS compiled with ``gcc.cmake`` on those distributions also uses
  AVX2 and FMA instructions.

Other factors
"""""""""""""

The gains from compiler optimization are modest compared to those from
the :doc:`accelerator packages <Speed_packages>` or from the choice of
model and settings in the input (see :doc:`Measuring performance
<Speed_measure>`).  All the measurements above were done on single
nodes.  On clusters, the speed and latency of the network between the
nodes usually affects the parallel performance far more than compiler
optimization.

The Linux kernel enables mitigations for security vulnerabilities of
the CPU, which can slow down specific instructions.  For example, the
mitigation for "Gather Data Sampling" (also known as "Downfall") on
Intel CPUs from the Skylake to the Ice Lake generations makes vector
gather instructions much slower.  On the Cascade Lake CPU from the table
above, the INTEL package styles, which use these instructions, ran on
average 55 percent faster with this mitigation disabled, while the
plain, OPENMP, and KOKKOS styles were not affected.  Whether such
mitigations may be disabled is a decision for the system
administrators.
