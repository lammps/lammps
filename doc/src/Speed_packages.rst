Accelerator packages
====================

Accelerated versions of various :doc:`pair_style <pair_style>`,
:doc:`fixes <fix>`, :doc:`computes <compute>`, and other commands have
been added to LAMMPS, which will typically run faster than the
standard non-accelerated versions.  Some require appropriate hardware
to be present on your system, e.g. GPUs or Intel Xeon Phi
co-processors.

All of these commands are in packages provided with LAMMPS.  An
overview of packages is give on the :doc:`Packages <Packages>` doc
pages.

These are the accelerator packages currently in LAMMPS, either as
standard or user packages:

+-----------------------------------------+-------------------------------------------------------+
| :doc:`GPU Package <Speed_gpu>`          | for NVIDIA GPUs as well as OpenCL support             |
+-----------------------------------------+-------------------------------------------------------+
| :doc:`USER-INTEL Package <Speed_intel>` | for Intel CPUs and Intel Xeon Phi                     |
+-----------------------------------------+-------------------------------------------------------+
| :doc:`KOKKOS Package <Speed_kokkos>`    | for Nvidia GPUs, Intel Xeon Phi, and OpenMP threading |
+-----------------------------------------+-------------------------------------------------------+
| :doc:`USER-OMP Package <Speed_omp>`     | for OpenMP threading and generic CPU optimizations    |
+-----------------------------------------+-------------------------------------------------------+
| :doc:`OPT Package <Speed_opt>`          | generic CPU optimizations                             |
+-----------------------------------------+-------------------------------------------------------+


.. toctree::
   :maxdepth: 1
   :hidden:

   Speed_gpu
   Speed_intel
   Speed_kokkos
   Speed_omp
   Speed_opt

Inverting this list, LAMMPS currently has acceleration support for
three kinds of hardware, via the listed packages:

+----------------+-----------------------------------------------------------------------------------------------------------------------------+
| Many-core CPUs | :doc:`USER-INTEL <Speed_intel>`, :doc:`KOKKOS <Speed_kokkos>`, :doc:`USER-OMP <Speed_omp>`, :doc:`OPT <Speed_opt>` packages |
+----------------+-----------------------------------------------------------------------------------------------------------------------------+
| NVIDIA GPUs    | :doc:`GPU <Speed_gpu>`, :doc:`KOKKOS <Speed_kokkos>` packages                                                               |
+----------------+-----------------------------------------------------------------------------------------------------------------------------+
| Intel Phi      | :doc:`USER-INTEL <Speed_intel>`, :doc:`KOKKOS <Speed_kokkos>` packages                                                      |
+----------------+-----------------------------------------------------------------------------------------------------------------------------+

Which package is fastest for your hardware may depend on the size
problem you are running and what commands (accelerated and
non-accelerated) are invoked by your input script.  While these doc
pages include performance guidelines, there is no substitute for
trying out the different packages appropriate to your hardware.

Any accelerated style has the same name as the corresponding standard
style, except that a suffix is appended.  Otherwise, the syntax for
the command that uses the style is identical, their functionality is
the same, and the numerical results it produces should also be the
same, except for precision and round-off effects.

For example, all of these styles are accelerated variants of the
Lennard-Jones :doc:`pair_style lj/cut <pair_lj>`:

* :doc:`pair_style lj/cut/gpu <pair_lj>`
* :doc:`pair_style lj/cut/intel <pair_lj>`
* :doc:`pair_style lj/cut/kk <pair_lj>`
* :doc:`pair_style lj/cut/omp <pair_lj>`
* :doc:`pair_style lj/cut/opt <pair_lj>`

To see what accelerate styles are currently available for a particular
style, find the style name in the `Commands\_all <lc_>`_
style pages (fix,compute,pair,etc) and see what suffixes are listed
(g,i,k,o,t) with it.  The doc pages for individual commands
(e.g. :doc:`pair lj/cut <pair_lj>` or :doc:`fix nve <fix_nve>`) also list
any accelerated variants available for that style.

To use an accelerator package in LAMMPS, and one or more of the styles
it provides, follow these general steps.  Details vary from package to
package and are explained in the individual accelerator doc pages,
listed above:

+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| build the accelerator library                                                                                                  | only for GPU package                                                 |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| install the accelerator package                                                                                                | make yes-opt, make yes-user-intel, etc                               |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| add compile/link flags to Makefile.machine in src/MAKE                                                                         | only for USER-INTEL, KOKKOS, USER-OMP, OPT packages                  |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| re-build LAMMPS                                                                                                                | make machine                                                         |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| prepare and test a regular LAMMPS simulation                                                                                   | lmp\_machine -in in.script; mpirun -np 32 lmp\_machine -in in.script |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| enable specific accelerator support via '-k on' :doc:`command-line switch <Run_options>`,                                      | only needed for KOKKOS package                                       |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| set any needed options for the package via "-pk" :doc:`command-line switch <Run_options>` or :doc:`package <package>` command, | only if defaults need to be changed                                  |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+
| use accelerated styles in your input via "-sf" :doc:`command-line switch <Run_options>` or :doc:`suffix <suffix>` command      | lmp\_machine -in in.script -sf gpu                                   |
+--------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------+

Note that the first 4 steps can be done as a single command with
suitable make command invocations. This is discussed on the
:doc:`Packages <Packages>` doc pages, and its use is illustrated in the
individual accelerator sections.  Typically these steps only need to
be done once, to create an executable that uses one or more
accelerator packages.

The last 4 steps can all be done from the command-line when LAMMPS is
launched, without changing your input script, as illustrated in the
individual accelerator sections.  Or you can add
:doc:`package <package>` and :doc:`suffix <suffix>` commands to your input
script.

.. note::

   With a few exceptions, you can build a single LAMMPS executable
   with all its accelerator packages installed.  Note however that the
   USER-INTEL and KOKKOS packages require you to choose one of their
   hardware options when building for a specific platform.  I.e. CPU or
   Phi option for the USER-INTEL package.  Or the OpenMP, Cuda, or Phi
   option for the KOKKOS package.

These are the exceptions.  You cannot build a single executable with:

* both the USER-INTEL Phi and KOKKOS Phi options
* the USER-INTEL Phi or Kokkos Phi option, and the GPU package

See the examples/accelerate/README and make.list files for sample
Make.py commands that build LAMMPS with any or all of the accelerator
packages.  As an example, here is a command that builds with all the
GPU related packages installed (GPU, KOKKOS with Cuda), including
settings to build the needed auxiliary GPU libraries for Kepler GPUs:


.. parsed-literal::

   Make.py -j 16 -p omp gpu kokkos -cc nvcc wrap=mpi   -gpu mode=double arch=35 -kokkos cuda arch=35 lib-all file mpi

The examples/accelerate directory also has input scripts that can be
used with all of the accelerator packages.  See its README file for
details.

Likewise, the bench directory has FERMI and KEPLER and PHI
sub-directories with Make.py commands and input scripts for using all
the accelerator packages on various machines.  See the README files in
those dirs.

As mentioned above, the `Benchmark page <http://lammps.sandia.gov/bench.html>`_ of the LAMMPS web site gives
performance results for the various accelerator packages for several
of the standard LAMMPS benchmark problems, as a function of problem
size and number of compute nodes, on different hardware platforms.

Here is a brief summary of what the various packages provide.  Details
are in the individual accelerator sections.

* Styles with a "gpu" suffix are part of the GPU package, and can be run
  on NVIDIA GPUs.  The speed-up on a GPU depends on a variety of
  factors, discussed in the accelerator sections.
* Styles with an "intel" suffix are part of the USER-INTEL
  package. These styles support vectorized single and mixed precision
  calculations, in addition to full double precision.  In extreme cases,
  this can provide speedups over 3.5x on CPUs.  The package also
  supports acceleration in "offload" mode to Intel(R) Xeon Phi(TM)
  co-processors.  This can result in additional speedup over 2x depending
  on the hardware configuration.
* Styles with a "kk" suffix are part of the KOKKOS package, and can be
  run using OpenMP on multicore CPUs, on an NVIDIA GPU, or on an Intel
  Xeon Phi in "native" mode.  The speed-up depends on a variety of
  factors, as discussed on the KOKKOS accelerator page.
* Styles with an "omp" suffix are part of the USER-OMP package and allow
  a pair-style to be run in multi-threaded mode using OpenMP.  This can
  be useful on nodes with high-core counts when using less MPI processes
  than cores is advantageous, e.g. when running with PPPM so that FFTs
  are run on fewer MPI processors or when the many MPI tasks would
  overload the available bandwidth for communication.
* Styles with an "opt" suffix are part of the OPT package and typically
  speed-up the pairwise calculations of your simulation by 5-25% on a
  CPU.


The individual accelerator package doc pages explain:

* what hardware and software the accelerated package requires
* how to build LAMMPS with the accelerated package
* how to run with the accelerated package either via command-line switches or modifying the input script
* speed-ups to expect
* guidelines for best performance
* restrictions
