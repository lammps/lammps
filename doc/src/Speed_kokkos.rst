KOKKOS package
==============

Kokkos is a templated C++ library that provides abstractions to allow
a single implementation of an application kernel (e.g. a pair style)
to run efficiently on different kinds of hardware, such as GPUs, Intel
Xeon Phis, or many-core CPUs. Kokkos maps the C++ kernel onto
different back end languages such as CUDA, OpenMP, or Pthreads.  The
Kokkos library also provides data abstractions to adjust (at compile
time) the memory layout of data structures like 2d and 3d arrays to
optimize performance on different hardware. For more information on
Kokkos, see `the Kokkos GitHub page <https://github.com/kokkos/kokkos>`_.

The LAMMPS KOKKOS package contains versions of pair, fix, and atom
styles that use data structures and macros provided by the Kokkos
library, which is included with LAMMPS in /lib/kokkos. The KOKKOS
package was developed primarily by Christian Trott (Sandia) and Stan
Moore (Sandia) with contributions of various styles by others,
including Sikandar Mashayak (UIUC), Ray Shan (Sandia), and Dan Ibanez
(Sandia). For more information on developing using Kokkos abstractions
see the Kokkos `Wiki <https://github.com/kokkos/kokkos/wiki>`_.

Kokkos currently provides support for 4 modes of execution (per MPI
task). These are Serial (MPI-only for CPUs and Intel Phi), OpenMP
(threading for many-core CPUs and Intel Phi), CUDA (for NVIDIA
GPUs) and HIP (for AMD GPUs). You choose the mode at build time to
produce an executable compatible with a specific hardware.

.. admonition:: C++14 support
   :class: note

   Kokkos requires using a compiler that supports the c++14 standard. For
   some compilers, it may be necessary to add a flag to enable c++14 support.
   For example, the GNU compiler uses the -std=c++14 flag. For a list of
   compilers that have been tested with the Kokkos library, see the Kokkos
   `README <https://github.com/kokkos/kokkos/blob/master/README.md>`_.

.. admonition:: NVIDIA CUDA support
   :class: note

   To build with Kokkos support for NVIDIA GPUs, the NVIDIA CUDA toolkit
   software version 9.0 or later must be installed on your system. See
   the discussion for the :doc:`GPU package <Speed_gpu>` for details of
   how to check and do this.

.. admonition:: CUDA and MPI library compatibility
   :class: note

   Kokkos with CUDA currently implicitly assumes that the MPI library is
   GPU-aware. This is not always the case, especially when using
   pre-compiled MPI libraries provided by a Linux distribution. This is
   not a problem when using only a single GPU with a single MPI
   rank. When running with multiple MPI ranks, you may see segmentation
   faults without GPU-aware MPI support. These can be avoided by adding
   the flags :doc:`-pk kokkos gpu/aware off <Run_options>` to the
   LAMMPS command line or by using the command :doc:`package kokkos
   gpu/aware off <package>` in the input file.

.. admonition:: AMD GPU support
   :class: note

   To build with Kokkos the HIPCC compiler from the AMD ROCm software
   version 3.5 or later is required.  Supporting this Kokkos mode in
   LAMMPS is still work in progress.  Please contact the LAMMPS developers
   if you run into problems.

Building LAMMPS with the KOKKOS package
"""""""""""""""""""""""""""""""""""""""

See the :ref:`Build extras <kokkos>` page for instructions.

Running LAMMPS with the KOKKOS package
""""""""""""""""""""""""""""""""""""""

All Kokkos operations occur within the context of an individual MPI task
running on a single node of the machine. The total number of MPI tasks
used by LAMMPS (one or multiple per compute node) is set in the usual
manner via the ``mpirun`` or ``mpiexec`` commands, and is independent of
Kokkos. E.g. the mpirun command in OpenMPI does this via its ``-np`` and
``-npernode`` switches. Ditto for MPICH via ``-np`` and ``-ppn``.

Running on a multicore CPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a quick overview of how to use the KOKKOS package
for CPU acceleration, assuming one or more 16-core nodes.

.. code-block:: bash

   mpirun -np 16 lmp_kokkos_mpi_only -k on -sf kk -in in.lj        # 1 node, 16 MPI tasks/node, no multi-threading
   mpirun -np 2 -ppn 1 lmp_kokkos_omp -k on t 16 -sf kk -in in.lj  # 2 nodes, 1 MPI task/node, 16 threads/task
   mpirun -np 2 lmp_kokkos_omp -k on t 8 -sf kk -in in.lj          # 1 node,  2 MPI tasks/node, 8 threads/task
   mpirun -np 32 -ppn 4 lmp_kokkos_omp -k on t 4 -sf kk -in in.lj  # 8 nodes, 4 MPI tasks/node, 4 threads/task

To run using the KOKKOS package, use the "-k on", "-sf kk" and "-pk
kokkos" :doc:`command-line switches <Run_options>` in your mpirun
command.  You must use the "-k on" :doc:`command-line switch <Run_options>` to enable the KOKKOS package. It takes
additional arguments for hardware settings appropriate to your system.
For OpenMP use:

.. parsed-literal::

   -k on t Nt

The "t Nt" option specifies how many OpenMP threads per MPI task to
use with a node. The default is Nt = 1, which is MPI-only mode.  Note
that the product of MPI tasks \* OpenMP threads/task should not exceed
the physical number of cores (on a node), otherwise performance will
suffer. If Hyper-Threading (HT) is enabled, then the product of MPI
tasks \* OpenMP threads/task should not exceed the physical number of
cores \* hardware threads.  The "-k on" switch also issues a
"package kokkos" command (with no additional arguments) which sets
various KOKKOS options to default values, as discussed on the
:doc:`package <package>` command doc page.

The "-sf kk" :doc:`command-line switch <Run_options>` will automatically
append the "/kk" suffix to styles that support it.  In this manner no
modification to the input script is needed. Alternatively, one can run
with the KOKKOS package by editing the input script as described
below.

.. note::

   When using a single OpenMP thread, the Kokkos Serial back end (i.e.
   Makefile.kokkos_mpi_only) will give better performance than the OpenMP
   back end (i.e. Makefile.kokkos_omp) because some of the overhead to make
   the code thread-safe is removed.

.. note::

   Use the "-pk kokkos" :doc:`command-line switch <Run_options>` to
   change the default :doc:`package kokkos <package>` options. See its doc
   page for details and default settings. Experimenting with its options
   can provide a speed-up for specific calculations. For example:

.. code-block:: bash

   mpirun -np 16 lmp_kokkos_mpi_only -k on -sf kk -pk kokkos newton on neigh half comm no -in in.lj       # Newton on, Half neighbor list, non-threaded comm

If the :doc:`newton <newton>` command is used in the input
script, it can also override the Newton flag defaults.

For half neighbor lists and OpenMP, the KOKKOS package uses data
duplication (i.e. thread-private arrays) by default to avoid
thread-level write conflicts in the force arrays (and other data
structures as necessary). Data duplication is typically fastest for
small numbers of threads (i.e. 8 or less) but does increase memory
footprint and is not scalable to large numbers of threads. An
alternative to data duplication is to use thread-level atomic operations
which do not require data duplication. The use of atomic operations can
be enforced by compiling LAMMPS with the "-DLMP_KOKKOS_USE_ATOMICS"
pre-processor flag. Most but not all Kokkos-enabled pair_styles support
data duplication. Alternatively, full neighbor lists avoid the need for
duplication or atomic operations but require more compute operations per
atom.  When using the Kokkos Serial back end or the OpenMP back end with
a single thread, no duplication or atomic operations are used. For CUDA
and half neighbor lists, the KOKKOS package always uses atomic operations.

CPU Cores, Sockets and Thread Affinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using multi-threading, it is important for performance to bind
both MPI tasks to physical cores, and threads to physical cores, so
they do not migrate during a simulation.

If you are not certain MPI tasks are being bound (check the defaults
for your MPI installation), binding can be forced with these flags:

.. parsed-literal::

   OpenMPI 1.8:  mpirun -np 2 --bind-to socket --map-by socket ./lmp_openmpi ...
   Mvapich2 2.0: mpiexec -np 2 --bind-to socket --map-by socket ./lmp_mvapich ...

For binding threads with KOKKOS OpenMP, use thread affinity environment
variables to force binding. With OpenMP 3.1 (gcc 4.7 or later, intel 12
or later) setting the environment variable ``OMP_PROC_BIND=true`` should
be sufficient. In general, for best performance with OpenMP 4.0 or later
set ``OMP_PROC_BIND=spread`` and ``OMP_PLACES=threads``.  For binding
threads with the KOKKOS pthreads option, compile LAMMPS with the hwloc
or libnuma support enabled as described in the :ref:`extra build options page <kokkos>`.

Running on Knight's Landing (KNL) Intel Xeon Phi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a quick overview of how to use the KOKKOS package for the
Intel Knight's Landing (KNL) Xeon Phi:

KNL Intel Phi chips have 68 physical cores. Typically 1 to 4 cores are
reserved for the OS, and only 64 or 66 cores are used. Each core has 4
Hyper-Threads,so there are effectively N = 256 (4\*64) or N = 264 (4\*66)
cores to run on. The product of MPI tasks \* OpenMP threads/task should
not exceed this limit, otherwise performance will suffer. Note that
with the KOKKOS package you do not need to specify how many KNLs there
are per node; each KNL is simply treated as running some number of MPI
tasks.

Examples of mpirun commands that follow these rules are shown below.

.. code-block:: bash

   # Running on an Intel KNL node with 68 cores (272 threads/node via 4x hardware threading):
   mpirun -np 64 lmp_kokkos_phi -k on t 4 -sf kk -in in.lj      # 1 node, 64 MPI tasks/node, 4 threads/task
   mpirun -np 66 lmp_kokkos_phi -k on t 4 -sf kk -in in.lj      # 1 node, 66 MPI tasks/node, 4 threads/task
   mpirun -np 32 lmp_kokkos_phi -k on t 8 -sf kk -in in.lj      # 1 node, 32 MPI tasks/node, 8 threads/task
   mpirun -np 512 -ppn 64 lmp_kokkos_phi -k on t 4 -sf kk -in in.lj  # 8 nodes, 64 MPI tasks/node, 4 threads/task

The -np setting of the mpirun command sets the number of MPI
tasks/node. The "-k on t Nt" command-line switch sets the number of
threads/task as Nt. The product of these two values should be N, i.e.
256 or 264.

.. note::

   The default for the :doc:`package kokkos <package>` command when
   running on KNL is to use "half" neighbor lists and set the Newton
   flag to "on" for both pairwise and bonded interactions. This will
   typically be best for many-body potentials. For simpler pairwise
   potentials, it may be faster to use a "full" neighbor list with
   Newton flag to "off".  Use the "-pk kokkos" :doc:`command-line switch
   <Run_options>` to change the default :doc:`package kokkos <package>`
   options. See its documentation page for details and default
   settings. Experimenting with its options can provide a speed-up for
   specific calculations. For example:

.. code-block:: bash

   mpirun -np 64 lmp_kokkos_phi -k on t 4 -sf kk -pk kokkos comm host -in in.reax      #  Newton on, half neighbor list, threaded comm
   mpirun -np 64 lmp_kokkos_phi -k on t 4 -sf kk -pk kokkos newton off neigh full comm no -in in.lj      # Newton off, full neighbor list, non-threaded comm

.. note::

   MPI tasks and threads should be bound to cores as described
   above for CPUs.

.. note::

   To build with Kokkos support for Intel Xeon Phi co-processors
   such as Knight's Corner (KNC), your system must be configured to use
   them in "native" mode, not "offload" mode like the INTEL package
   supports.

Running on GPUs
^^^^^^^^^^^^^^^

Use the "-k" :doc:`command-line switch <Run_options>` to specify the
number of GPUs per node. Typically the -np setting of the mpirun command
should set the number of MPI tasks/node to be equal to the number of
physical GPUs on the node. You can assign multiple MPI tasks to the same
GPU with the KOKKOS package, but this is usually only faster if some
portions of the input script have not been ported to use Kokkos. In this
case, also packing/unpacking communication buffers on the host may give
speedup (see the KOKKOS :doc:`package <package>` command). Using CUDA MPS
is recommended in this scenario.

Using a GPU-aware MPI library is highly recommended. GPU-aware MPI use can be
avoided by using :doc:`-pk kokkos gpu/aware off <package>`. As above for
multicore CPUs (and no GPU), if N is the number of physical cores/node,
then the number of MPI tasks/node should not exceed N.

.. parsed-literal::

   -k on g Ng

Here are examples of how to use the KOKKOS package for GPUs, assuming
one or more nodes, each with two GPUs:

.. code-block:: bash

   mpirun -np 2 lmp_kokkos_cuda_openmpi -k on g 2 -sf kk -in in.lj          # 1 node,   2 MPI tasks/node, 2 GPUs/node
   mpirun -np 32 -ppn 2 lmp_kokkos_cuda_openmpi -k on g 2 -sf kk -in in.lj  # 16 nodes, 2 MPI tasks/node, 2 GPUs/node (32 GPUs total)

.. note::

   The default for the :doc:`package kokkos <package>` command when
   running on GPUs is to use "full" neighbor lists and set the Newton
   flag to "off" for both pairwise and bonded interactions, along with
   threaded communication. When running on Maxwell or Kepler GPUs, this
   will typically be best. For Pascal GPUs and beyond, using "half"
   neighbor lists and setting the Newton flag to "on" may be faster. For
   many pair styles, setting the neighbor binsize equal to twice the CPU
   default value will give speedup, which is the default when running on
   GPUs. Use the "-pk kokkos" :doc:`command-line switch <Run_options>`
   to change the default :doc:`package kokkos <package>` options. See
   its documentation page for details and default
   settings. Experimenting with its options can provide a speed-up for
   specific calculations. For example:

.. code-block:: bash

   mpirun -np 2 lmp_kokkos_cuda_openmpi -k on g 2 -sf kk -pk kokkos newton on neigh half binsize 2.8 -in in.lj      # Newton on, half neighbor list, set binsize = neighbor ghost cutoff

.. note::

   When using a GPU, you will achieve the best performance if your
   input script does not use fix or compute styles which are not yet
   Kokkos-enabled. This allows data to stay on the GPU for multiple
   timesteps, without being copied back to the host CPU. Invoking a
   non-Kokkos fix or compute, or performing I/O for
   :doc:`thermo <thermo_style>` or :doc:`dump <dump>` output will cause data
   to be copied back to the CPU incurring a performance penalty.

.. note::

   To get an accurate timing breakdown between time spend in pair,
   kspace, etc., you must set the environment variable CUDA_LAUNCH_BLOCKING=1.
   However, this will reduce performance and is not recommended for production runs.

Run with the KOKKOS package by editing an input script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively the effect of the "-sf" or "-pk" switches can be
duplicated by adding the :doc:`package kokkos <package>` or :doc:`suffix kk <suffix>` commands to your input script.

The discussion above for building LAMMPS with the KOKKOS package, the
``mpirun`` or ``mpiexec`` command, and setting appropriate thread
properties are the same.

You must still use the "-k on" :doc:`command-line switch <Run_options>`
to enable the KOKKOS package, and specify its additional arguments for
hardware options appropriate to your system, as documented above.

You can use the :doc:`suffix kk <suffix>` command, or you can explicitly add a
"kk" suffix to individual styles in your input script, e.g.

.. code-block:: LAMMPS

   pair_style lj/cut/kk 2.5

You only need to use the :doc:`package kokkos <package>` command if you
wish to change any of its option defaults, as set by the "-k on"
:doc:`command-line switch <Run_options>`.

**Using OpenMP threading and CUDA together:**

With the KOKKOS package, both OpenMP multi-threading and GPUs can be
compiled and used together in a few special cases. In the makefile for
the conventional build, the KOKKOS_DEVICES variable must include both,
"Cuda" and "OpenMP", as is the case for ``/src/MAKE/OPTIONS/Makefile.kokkos_cuda_mpi``.

.. code-block:: bash

   KOKKOS_DEVICES=Cuda,OpenMP

When building with CMake you need to enable both features as it is done
in the ``kokkos-cuda.cmake`` CMake preset file.

.. code-block:: bash

   cmake ../cmake -DKokkos_ENABLE_CUDA=yes -DKokkos_ENABLE_OPENMP=yes

The suffix "/kk" is equivalent to "/kk/device", and for Kokkos CUDA,
using the "-sf kk" in the command line gives the default CUDA version
everywhere.  However, if the "/kk/host" suffix is added to a specific
style in the input script, the Kokkos OpenMP (CPU) version of that
specific style will be used instead.  Set the number of OpenMP threads
as "t Nt" and the number of GPUs as "g Ng"

.. parsed-literal::

   -k on t Nt g Ng

For example, the command to run with 1 GPU and 8 OpenMP threads is then:

.. code-block:: bash

   mpiexec -np 1 lmp_kokkos_cuda_openmpi -in in.lj -k on g 1 t 8 -sf kk

Conversely, if the "-sf kk/host" is used in the command line and then
the "/kk" or "/kk/device" suffix is added to a specific style in your
input script, then only that specific style will run on the GPU while
everything else will run on the CPU in OpenMP mode. Note that the
execution of the CPU and GPU styles will NOT overlap, except for a
special case:

A kspace style and/or molecular topology (bonds, angles, etc.) running
on the host CPU can overlap with a pair style running on the
GPU. First compile with "--default-stream per-thread" added to CCFLAGS
in the Kokkos CUDA Makefile.  Then explicitly use the "/kk/host"
suffix for kspace and bonds, angles, etc.  in the input file and the
"kk" suffix (equal to "kk/device") on the command line.  Also make
sure the environment variable CUDA_LAUNCH_BLOCKING is not set to "1"
so CPU/GPU overlap can occur.

Performance to expect
"""""""""""""""""""""

The performance of KOKKOS running in different modes is a function of
your hardware, which KOKKOS-enable styles are used, and the problem
size.

Generally speaking, the following rules of thumb apply:

* When running on CPUs only, with a single thread per MPI task,
  performance of a KOKKOS style is somewhere between the standard
  (un-accelerated) styles (MPI-only mode), and those provided by the
  OPENMP package. However the difference between all 3 is small (less
  than 20%).
* When running on CPUs only, with multiple threads per MPI task,
  performance of a KOKKOS style is a bit slower than the OPENMP
  package.
* When running large number of atoms per GPU, KOKKOS is typically faster
  than the GPU package when compiled for double precision. The benefit
  of using single or mixed precision with the GPU package depends
  significantly on the hardware in use and the simulated system and pair
  style.
* When running on Intel hardware, KOKKOS is not as fast as
  the INTEL package, which is optimized for x86 hardware (not just
  from Intel) and compilation with the Intel compilers.  The INTEL
  package also can increase the vector length of vector instructions
  by switching to single or mixed precision mode.

See the `Benchmark page <https://www.lammps.org/bench.html>`_ of the
LAMMPS website for performance of the KOKKOS package on different
hardware.

Advanced Kokkos options
"""""""""""""""""""""""

There are other allowed options when building with the KOKKOS package
that can improve performance or assist in debugging or profiling.
They are explained on the :ref:`KOKKOS section of the build extras <kokkos>` doc page,

Restrictions
""""""""""""

Currently, there are no precision options with the KOKKOS package. All
compilation and computation is performed in double precision.
