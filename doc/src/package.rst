.. index:: package

package command
===============

Syntax
""""""

.. code-block:: LAMMPS

   package style args

* style = *gpu* or *intel* or *kokkos* or *omp*
* args = arguments specific to the style

  .. parsed-literal::

       *gpu* args = Ngpu keyword value ...
         Ngpu = # of GPUs per node
         zero or more keyword/value pairs may be appended
         keywords = *neigh* or *newton* or *pair/only* or *binsize* or *split* or *gpuID* or *tpa* or *blocksize* or *omp* or *platform* or *device_type* or *ocl_args*
           *neigh* value = *yes* or *no*
             *yes* = neighbor list build on GPU (default)
             *no* = neighbor list build on CPU
           *newton* = *off* or *on*
             *off* = set Newton pairwise flag off (default and required)
             *on* = set Newton pairwise flag on (currently not allowed)
           *pair/only* = *off* or *on*
             *off* = apply "gpu" suffix to all available styles in the GPU package (default)
             *on* = apply "gpu" suffix only pair styles
           *binsize* value = size
             size = bin size for neighbor list construction (distance units)
           *split* = fraction
             fraction = fraction of atoms assigned to GPU (default = 1.0)
           *tpa* value = Nlanes
             Nlanes = # of GPU vector lanes (CUDA threads) used per atom
           *blocksize* value = size
             size = thread block size for pair force computation
           *omp* value = Nthreads
             Nthreads = number of OpenMP threads to use on CPU (default = 0)
           *platform* value = id
             id = For OpenCL, platform ID for the GPU or accelerator
           *gpuID* values = id
             id = ID of first GPU to be used on each node
           *device_type* value = *intelgpu* or *nvidiagpu* or *amdgpu* or *applegpu* or *generic* or *custom*,val1,val2,...
             val1,val2,... = custom OpenCL accelerator configuration parameters (see below for details)
           *ocl_args* value = args
             args = List of additional OpenCL compiler arguments delimited by colons
       *intel* args = NPhi keyword value ...
         Nphi = # of co-processors per node
         zero or more keyword/value pairs may be appended
         keywords = *mode* or *omp* or *lrt* or *balance* or *ghost* or *tpc* or *tptask* or *pppm_table* or *no_affinity*
           *mode* value = *single* or *mixed* or *double*
             single = perform force calculations in single precision
             mixed = perform force calculations in mixed precision
             double = perform force calculations in double precision
           *omp* value = Nthreads
             Nthreads = number of OpenMP threads to use on CPU (default = 0)
           *lrt* value = *yes* or *no*
             *yes* = use additional thread dedicated for some PPPM calculations
             *no* = do not dedicate an extra thread for some PPPM calculations
           *balance* value = split
             split = fraction of work to offload to co-processor, -1 for dynamic
           *ghost* value = *yes* or *no*
             *yes* = include ghost atoms for offload
             *no* = do not include ghost atoms for offload
           *tpc* value = Ntpc
             Ntpc = max number of co-processor threads per co-processor core (default = 4)
           *tptask* value = Ntptask
             Ntptask = max number of co-processor threads per MPI task (default = 240)
           *pppm_table* value = *yes* or *no*
             *yes* = Precompute pppm values in table (doesn't change accuracy)
             *no* = Compute pppm values on the fly
           *no_affinity* values = none
       *kokkos* args = keyword value ...
         zero or more keyword/value pairs may be appended
         keywords = *neigh* or *neigh/qeq* or *neigh/thread* or *neigh/transpose* or *newton* or *binsize* or *comm* or *comm/exchange* or *comm/forward* or *comm/pair/forward* or *comm/fix/forward* or *comm/reverse* or *comm/pair/reverse* or *sort* or *atom/map* or *gpu/aware* or *pair/only*
           *neigh* value = *full* or *half*
             full = full neighbor list
             half = half neighbor list built in thread-safe manner
           *neigh/qeq* value = *full* or *half*
             full = full neighbor list
             half = half neighbor list built in thread-safe manner
           *neigh/thread* value = *off* or *on*
             *off* = thread only over atoms
             *on* = thread over both atoms and neighbors
           *neigh/transpose* value = *off* or *on*
             *off* = use same memory layout for GPU neigh list build as pair style
             *on* = use transposed memory layout for GPU neigh list build
           *newton* = *off* or *on*
             *off* = set Newton pairwise and bonded flags off
             *on* = set Newton pairwise and bonded flags on
           *binsize* value = size
             size = bin size for neighbor list construction (distance units)
           *comm* value = *no* or *host* or *device*
             use value for comm/exchange and comm/forward and comm/pair/forward and comm/fix/forward and comm/reverse
           *comm/exchange* value = *no* or *host* or *device*
           *comm/forward* value = *no* or *host* or *device*
           *comm/pair/forward* value = *no* or *device*
           *comm/fix/forward* value = *no* or *device*
           *comm/reverse* value = *no* or *host* or *device*
             *no* = perform communication pack/unpack in non-KOKKOS mode
             *host* = perform pack/unpack on host (e.g. with OpenMP threading)
             *device* = perform pack/unpack on device (e.g. on GPU)
           *comm/pair/reverse* value = *no* or *device*
             *no* = perform communication pack/unpack in non-KOKKOS mode
             *device* = perform pack/unpack on device (e.g. on GPU)
           *sort* value = *no* or *device*
             *no* = perform atom sorting in non-KOKKOS mode
             *device* = perform atom sorting on device (e.g. on GPU)
           *atom/map* value = *no* or *device*
             *no* = build atom map in non-KOKKOS mode
             *device* = build atom map on device (e.g. on GPU)
           *gpu/aware* = *off* or *on*
             *off* = do not use GPU-aware MPI
             *on* = use GPU-aware MPI (default)
           *pair/only* = *off* or *on*
             *off* = use device acceleration (e.g. GPU) for all available styles in the KOKKOS package (default)
             *on*  = use device acceleration only for pair styles (and host acceleration for others)
       *omp* args = Nthreads keyword value ...
         Nthreads = # of OpenMP threads to associate with each MPI process
         zero or more keyword/value pairs may be appended
         keywords = *neigh*
           *neigh* value = *yes* or *no*
             *yes* = threaded neighbor list build (default)
             *no* = non-threaded neighbor list build

Examples
""""""""

.. code-block:: LAMMPS

   package gpu 0
   package gpu 1 split 0.75
   package gpu 2 split -1.0
   package gpu 0 omp 2 device_type intelgpu
   package kokkos neigh half comm device
   package omp 0 neigh no
   package omp 4
   package intel 1
   package intel 2 omp 4 mode mixed balance 0.5

Description
"""""""""""

This command invokes package-specific settings for the various
accelerator packages available in LAMMPS.  Currently the following
packages use settings from this command: GPU, INTEL, KOKKOS, and
OPENMP.

If this command is specified in an input script, it must be near the
top of the script, before the simulation box has been defined.  This
is because it specifies settings that the accelerator packages use in
their initialization, before a simulation is defined.

This command can also be specified from the command-line when
launching LAMMPS, using the "-pk" :doc:`command-line switch <Run_options>`.  The syntax is exactly the same as when used
in an input script.

Note that all of the accelerator packages require the package command
to be specified (except the OPT package), if the package is to be used
in a simulation (LAMMPS can be built with an accelerator package
without using it in a particular simulation).  However, in all cases,
a default version of the command is typically invoked by other
accelerator settings.

The KOKKOS package requires a "-k on" :doc:`command-line switch <Run_options>` respectively, which invokes a "package
kokkos" command with default settings.

For the GPU, INTEL, and OPENMP packages, if a "-sf gpu" or "-sf
intel" or "-sf omp" :doc:`command-line switch <Run_options>` is used to
auto-append accelerator suffixes to various styles in the input
script, then those switches also invoke a "package gpu", "package
intel", or "package omp" command with default settings.

.. note::

   A package command for a particular style can be invoked multiple
   times when a simulation is setup, e.g. by the :doc:`-c on, -k on, -sf, and -pk command-line switches <Run_options>`, and by using this command
   in an input script.  Each time it is used all of the style options are
   set, either to default values or to specified settings.  I.e. settings
   from previous invocations do not persist across multiple invocations.

See the :doc:`Accelerator packages <Speed_packages>` page for more details
about using the various accelerator packages for speeding up LAMMPS
simulations.

----------

The *gpu* style invokes settings associated with the use of the GPU
package.

The *Ngpu* argument sets the number of GPUs per node. If *Ngpu* is 0
and no other keywords are specified, GPU or accelerator devices are
auto-selected. In this process, all platforms are searched for
accelerator devices and GPUs are chosen if available. The device with
the highest number of compute cores is selected. The number of devices
is increased to be the number of matching accelerators with the same
number of compute cores. If there are more devices than MPI tasks,
the additional devices will be unused. The auto-selection of GPUs/
accelerator devices and platforms can be restricted by specifying
a non-zero value for *Ngpu* and / or using the *gpuID*, *platform*,
and *device_type* keywords as described below. If there are more MPI
tasks (per node) than GPUs, multiple MPI tasks will share each GPU.

Optional keyword/value pairs can also be specified.  Each has a
default value as listed below.

The *neigh* keyword specifies where neighbor lists for pair style
computation will be built.  If *neigh* is *yes*, which is the default,
neighbor list building is performed on the GPU.  If *neigh* is *no*,
neighbor list building is performed on the CPU.  GPU neighbor list
building currently cannot be used with a triclinic box.  GPU neighbor
lists are not compatible with commands that are not GPU-enabled.  When
a non-GPU enabled command requires a neighbor list, it will also be
built on the CPU.  In these cases, it will typically be more efficient
to only use CPU neighbor list builds.

The *newton* keyword sets the Newton flags for pairwise (not bonded)
interactions to *off* or *on*, the same as the :doc:`newton <newton>`
command allows.  Currently, only an *off* value is allowed, since all
the GPU package pair styles require this setting.  This means more
computation is done, but less communication.  In the future a value of
*on* may be allowed, so the *newton* keyword is included as an option
for compatibility with the package command for other accelerator
styles.  Note that the newton setting for bonded interactions is not
affected by this keyword.

The *pair/only* keyword can change how any "gpu" suffix is applied.
By default a suffix is applied to all styles for which an accelerated
variant is available.  However, that is not always the most effective
way to use an accelerator.  With *pair/only* set to *on* the suffix
will only by applied to supported pair styles, which tend to be the
most effective in using an accelerator and their operation can be
overlapped with all other computations on the CPU.

The *binsize* keyword sets the size of bins used to bin atoms in
neighbor list builds performed on the GPU, if *neigh* = *yes* is set.
If *binsize* is set to 0.0 (the default), then the binsize is set
automatically using heuristics in the GPU package.

The *split* keyword can be used for load balancing force calculations
between CPU and GPU cores in GPU-enabled pair styles. If 0 < *split* <
1.0, a fixed fraction of particles is offloaded to the GPU while force
calculation for the other particles occurs simultaneously on the CPU.
If *split* < 0.0, the optimal fraction (based on CPU and GPU timings)
is calculated every 25 timesteps, i.e. dynamic load-balancing across
the CPU and GPU is performed.  If *split* = 1.0, all force
calculations for GPU accelerated pair styles are performed on the GPU.
In this case, other :doc:`hybrid <pair_hybrid>` pair interactions,
:doc:`bond <bond_style>`, :doc:`angle <angle_style>`,
:doc:`dihedral <dihedral_style>`, :doc:`improper <improper_style>`, and
:doc:`long-range <kspace_style>` calculations can be performed on the
CPU while the GPU is performing force calculations for the GPU-enabled
pair style.  If all CPU force computations complete before the GPU
completes, LAMMPS will block until the GPU has finished before
continuing the timestep.

As an example, if you have two GPUs per node and 8 CPU cores per node,
and would like to run on 4 nodes (32 cores) with dynamic balancing of
force calculation across CPU and GPU cores, you could specify

.. code-block:: bash

   mpirun -np 32 -sf gpu -in in.script    # launch command
   package gpu 2 split -1                 # input script command

In this case, all CPU cores and GPU devices on the nodes would be
utilized.  Each GPU device would be shared by 4 CPU cores. The CPU
cores would perform force calculations for some fraction of the
particles at the same time the GPUs performed force calculation for
the other particles.

The *gpuID* keyword is used to specify the first ID for the GPU or
other accelerator that LAMMPS will use. For example, if the ID is
1 and *Ngpu* is 3, GPUs 1-3 will be used. Device IDs should be
determined from the output of nvc_get_devices, ocl_get_devices,
or hip_get_devices
as provided in the lib/gpu directory. When using OpenCL with
accelerators that have main memory NUMA, the accelerators can be
split into smaller virtual accelerators for more efficient use
with MPI.

The *tpa* keyword sets the number of GPU vector lanes per atom used to
perform force calculations.  With a default value of 1, the number of
lanes will be chosen based on the pair style, however, the value can
be set explicitly with this keyword to fine-tune performance.  For
large cutoffs or with a small number of particles per GPU, increasing
the value can improve performance. The number of lanes per atom must
be a power of 2 and currently cannot be greater than the SIMD width
for the GPU / accelerator. In the case it exceeds the SIMD width, it
will automatically be decreased to meet the restriction.

The *blocksize* keyword allows you to tweak the number of threads used
per thread block. This number should be a multiple of 32 (for GPUs)
and its maximum depends on the specific GPU hardware. Typical choices
are 64, 128, or 256. A larger block size increases occupancy of
individual GPU cores, but reduces the total number of thread blocks,
thus may lead to load imbalance. On modern hardware, the sensitivity
to the blocksize is typically low.

The *Nthreads* value for the *omp* keyword sets the number of OpenMP
threads allocated for each MPI task. This setting controls OpenMP
parallelism only for routines run on the CPUs. For more details on
setting the number of OpenMP threads, see the discussion of the
*Nthreads* setting on this page for the "package omp" command.
The meaning of *Nthreads* is exactly the same for the GPU, INTEL,
and GPU packages.

The *platform* keyword is only used with OpenCL to specify the ID for
an OpenCL platform. See the output from ocl_get_devices in the lib/gpu
directory. In LAMMPS only one platform can be active at a time and by
default (id=-1) the platform is auto-selected to find the GPU with the
most compute cores. When *Ngpu* or other keywords are specified, the
auto-selection is appropriately restricted. For example, if *Ngpu* is
3, only platforms with at least 3 accelerators are considered. Similar
restrictions can be enforced by the *gpuID* and *device_type* keywords.

The *device_type* keyword can be used for OpenCL to specify the type of
GPU to use or specify a custom configuration for an accelerator. In most
cases this selection will be automatic and there is no need to use the
keyword. The *applegpu* type is not specific to a particular GPU vendor,
but is separate due to the more restrictive Apple OpenCL implementation.
For expert users, to specify a custom configuration, the *custom* keyword
followed by the next parameters can be specified:

CONFIG_ID, SIMD_SIZE, MEM_THREADS, SHUFFLE_AVAIL, FAST_MATH,
THREADS_PER_ATOM, THREADS_PER_CHARGE, THREADS_PER_THREE, BLOCK_PAIR,
BLOCK_BIO_PAIR, BLOCK_ELLIPSE, PPPM_BLOCK_1D, BLOCK_NBOR_BUILD,
BLOCK_CELL_2D, BLOCK_CELL_ID, MAX_SHARED_TYPES, MAX_BIO_SHARED_TYPES,
PPPM_MAX_SPLINE, NBOR_PREFETCH.

CONFIG_ID can be 0. SHUFFLE_AVAIL in {0,1} indicates that inline-PTX
(NVIDIA) or OpenCL extensions (Intel) should be used for horizontal
vector operations. FAST_MATH in {0,1} indicates that OpenCL fast math
optimizations are used during the build and hardware-accelerated
transcendental functions are used when available. THREADS_PER_* give the
default *tpa* values for ellipsoidal models, styles using charge, and
any other styles. The BLOCK_* parameters specify the block sizes for
various kernel calls and the MAX_*SHARED*_ parameters are used to
determine the amount of local shared memory to use for storing model
parameters.

For OpenCL, the routines are compiled at runtime for the specified GPU
or accelerator architecture. The *ocl_args* keyword can be used to
specify additional flags for the runtime build.

----------

The *intel* style invokes settings associated with the use of the INTEL
package.  The keywords *balance*, *ghost*, *tpc*, and *tptask* are
**only** applicable if LAMMPS was built with Xeon Phi co-processor
support and are otherwise ignored.

The *Nphi* argument sets the number of co-processors per node.
This can be set to any value, including 0, if LAMMPS was not
built with co-processor support.

Optional keyword/value pairs can also be specified.  Each has a
default value as listed below.

The *Nthreads* value for the *omp* keyword sets the number of OpenMP
threads allocated for each MPI task. This setting controls OpenMP
parallelism only for routines run on the CPUs. For more details on
setting the number of OpenMP threads, see the discussion of the
*Nthreads* setting on this page for the "package omp" command.
The meaning of *Nthreads* is exactly the same for the GPU, INTEL,
and GPU packages.

The *mode* keyword determines the precision mode to use for
computing pair style forces, either on the CPU or on the co-processor,
when using a INTEL supported :doc:`pair style <pair_style>`.  It
can take a value of *single*, *mixed* which is the default, or
*double*\ .  *Single* means single precision is used for the entire
force calculation.  *Mixed* means forces between a pair of atoms are
computed in single precision, but accumulated and stored in double
precision, including storage of forces, torques, energies, and virial
quantities.  *Double* means double precision is used for the entire
force calculation.

The *lrt* keyword can be used to enable "Long Range Thread (LRT)"
mode. It can take a value of *yes* to enable and *no* to disable.
LRT mode generates an extra thread (in addition to any OpenMP threads
specified with the OMP_NUM_THREADS environment variable or the *omp*
keyword). The extra thread is dedicated for performing part of the
:doc:`PPPM solver <kspace_style>` computations and communications. This
can improve parallel performance on processors supporting
Simultaneous Multithreading (SMT) such as Hyper-Threading (HT) on Intel
processors. In this mode, one additional thread is generated per MPI
process. LAMMPS will generate a warning in the case that more threads
are used than available in SMT hardware on a node. If the PPPM solver
from the INTEL package is not used, then the LRT setting is
ignored and no extra threads are generated. Enabling LRT will replace
the :doc:`run_style <run_style>` with the *verlet/lrt/intel* style that
is identical to the default *verlet* style aside from supporting the
LRT feature. This feature requires setting the pre-processor flag
-DLMP_INTEL_USELRT in the makefile when compiling LAMMPS.

The *balance* keyword sets the fraction of :doc:`pair style <pair_style>` work
offloaded to the co-processor for split values between 0.0 and 1.0 inclusive.
While this fraction of work is running on the co-processor, other calculations
will run on the host, including neighbor and pair calculations that are not
offloaded, as well as angle, bond, dihedral, kspace, and some MPI
communications.  If *split* is set to -1, the fraction of work is dynamically
adjusted automatically throughout the run.  This typically give performance
within 5 to 10 percent of the optimal fixed fraction.

The *ghost* keyword determines whether or not ghost atoms, i.e. atoms
at the boundaries of processor subdomains, are offloaded for neighbor
and force calculations.  When the value = "no", ghost atoms are not
offloaded.  This option can reduce the amount of data transfer with
the co-processor and can also overlap MPI communication of forces with
computation on the co-processor when the :doc:`newton pair <newton>`
setting is "on".  When the value = "yes", ghost atoms are offloaded.
In some cases this can provide better performance, especially if the
*balance* fraction is high.

The *tpc* keyword sets the max # of co-processor threads *Ntpc* that
will run on each core of the co-processor.  The default value = 4,
which is the number of hardware threads per core supported by the
current generation Xeon Phi chips.

The *tptask* keyword sets the max # of co-processor threads (Ntptask*
assigned to each MPI task.  The default value = 240, which is the
total # of threads an entire current generation Xeon Phi chip can run
(240 = 60 cores \* 4 threads/core).  This means each MPI task assigned
to the Phi will enough threads for the chip to run the max allowed,
even if only 1 MPI task is assigned.  If 8 MPI tasks are assigned to
the Phi, each will run with 30 threads.  If you wish to limit the
number of threads per MPI task, set *tptask* to a smaller value.
E.g. for *tptask* = 16, if 8 MPI tasks are assigned, each will run
with 16 threads, for a total of 128.

Note that the default settings for *tpc* and *tptask* are fine for
most problems, regardless of how many MPI tasks you assign to a Phi.

.. versionadded:: 15Jun2023

The *pppm_table* keyword with the argument yes allows to use a
pre-computed table to efficiently spread the charge to the PPPM grid.
This feature is enabled by default but can be turned off using the
keyword with the argument *no*.

The *no_affinity* keyword will turn off automatic setting of core
affinity for MPI tasks and OpenMP threads on the host when using
offload to a co-processor. Affinity settings are used when possible
to prevent MPI tasks and OpenMP threads from being on separate NUMA
domains and to prevent offload threads from interfering with other
processes/threads used for LAMMPS.

----------

The *kokkos* style invokes settings associated with the use of the
KOKKOS package.

All of the settings are optional keyword/value pairs. Each has a default
value as listed below.

The *neigh* keyword determines how neighbor lists are built. A value of
*half* uses a thread-safe variant of half-neighbor lists, the same as
used by most pair styles in LAMMPS, which is the default when running on
CPUs (i.e. the Kokkos CUDA back end is not enabled).

A value of *full* uses a full neighbor lists and is the default when
running on GPUs. This performs twice as much computation as the *half*
option, however that is often a win because it is thread-safe and
does not require atomic operations in the calculation of pair forces. For
that reason, *full* is the default setting for GPUs. However, when
running on CPUs, a *half* neighbor list is the default because it are
often faster, just as it is for non-accelerated pair styles. Similarly,
the *neigh/qeq* keyword determines how neighbor lists are built for
:doc:`fix qeq/reaxff/kk <fix_qeq_reaxff>`.

If the *neigh/thread* keyword is set to *off*, then the KOKKOS package
threads only over atoms. However, for small systems, this may not expose
enough parallelism to keep a GPU busy. When this keyword is set to *on*,
the KOKKOS package threads over both atoms and neighbors of atoms. When
using *neigh/thread* *on*, the :doc:`newton pair <newton>` setting must
be "off". Using *neigh/thread* *on* may be slower for large systems, so
this this option is turned on by default only when running on one or
more GPUs and there are 16k atoms or less owned by an MPI rank. Not all
KOKKOS-enabled potentials support this keyword yet, and only thread over
atoms. Many simple pairwise potentials such as Lennard-Jones do support
threading over both atoms and neighbors.

If the *neigh/transpose* keyword is set to *off*, then the KOKKOS
package will use the same memory layout for building the neighbor list on
GPUs as used for the pair style. When this keyword is set to *on* it
will use a different (transposed) memory layout to build the neighbor
list on GPUs. This can be faster in some cases (e.g. ReaxFF HNS
benchmark) but slower in others (e.g. Lennard Jones benchmark). The
copy between different memory layouts is done out of place and
therefore doubles the memory overhead of the neighbor list, which can
be significant.

The *newton* keyword sets the Newton flags for pairwise and bonded
interactions to *off* or *on*, the same as the :doc:`newton <newton>`
command allows. The default for GPUs is *off* because this will almost
always give better performance for the KOKKOS package. This means more
computation is done, but less communication. However, when running on
CPUs a value of *on* is the default since it can often be faster, just
as it is for non-accelerated pair styles

The *binsize* keyword sets the size of bins used to bin atoms during
neighbor list builds. The same value can be set by the
:doc:`neigh_modify binsize <neigh_modify>` command. Making it an option
in the package kokkos command allows it to be set from the command line.
The default value for CPUs is 0.0, which means the LAMMPS default will be
used, which is bins = 1/2 the size of the pairwise cutoff + neighbor skin
distance. This is fine when neighbor lists are built on the CPU. For GPU
builds, a 2x larger binsize equal to the pairwise cutoff + neighbor skin
is often faster, which is the default. Note that if you use a
longer-than-usual pairwise cutoff, e.g. to allow for a smaller fraction
of KSpace work with a :doc:`long-range Coulombic solver <kspace_style>`
because the GPU is faster at performing pairwise interactions, then this
rule of thumb may give too large a binsize and the default should be
overridden with a smaller value.

The *comm* and *comm/exchange* and *comm/forward* and *comm/pair/forward*
and *comm/fix/forward* and *comm/reverse* and *comm/pair/reverse*
keywords determine whether the host or device performs the packing and
unpacking of data when communicating per-atom data between processors.
"Exchange" communication happens only on timesteps that neighbor lists
are rebuilt. The data is only for atoms that migrate to new processors.
"Forward" communication happens every timestep. "Reverse" communication
happens every timestep if the *newton* option is on. The data is for
atom coordinates and any other atom properties that needs to be updated
for ghost atoms owned by each processor. "Pair/comm" controls additional
communication in pair styles, such as pair_style EAM. "Fix/comm" controls
additional communication in fixes, such as fix SHAKE.

The *comm* keyword is simply a short-cut to set the same value for all
the comm keywords.

The value options for the keywords are *no* or *host* or *device*\ . A
value of *no* means to use the standard non-KOKKOS method of
packing/unpacking data for the communication. A value of *host* means to
use the host, typically a multicore CPU, and perform the
packing/unpacking in parallel with threads. A value of *device* means to
use the device, typically a GPU, to perform the packing/unpacking
operation.

For the *comm/pair/forward* or *comm/fix/forward* or *comm/pair/reverse*
keywords, if a value of *host* is used it will be automatically
be changed to *no* since these keywords don't support *host* mode. The
value of *no* will also always be used when running on the CPU, i.e. setting
the value to *device* will have no effect if the pair/fix style is
running on the CPU. For the *comm/fix/forward* or *comm/pair/reverse*
keywords, not all styles support *device* mode and in that case will run
in *no* mode instead.

The optimal choice for these keywords depends on the input script and
the hardware used. The *no* value is useful for verifying that the
Kokkos-based *host* and *device* values are working correctly. It is the
default when running on CPUs since it is usually the fastest.

When running on CPUs or Xeon Phi, the *host* and *device* values work
identically. When using GPUs, the *device* value is the default since it
will typically be optimal if all of your styles used in your input
script are supported by the KOKKOS package. In this case data can stay
on the GPU for many timesteps without being moved between the host and
GPU, if you use the *device* value. If your script uses styles (e.g.
fixes) which are not yet supported by the KOKKOS package, then data has
to be moved between the host and device anyway, so it is typically faster
to let the host handle communication, by using the *host* value. Using
*host* instead of *no* will enable use of multiple threads to
pack/unpack communicated data. When running small systems on a GPU,
performing the exchange pack/unpack on the host CPU can give speedup
since it reduces the number of CUDA kernel launches.

The *sort* keyword determines whether the host or device performs atom
sorting, see the :doc:`atom_modify sort <atom_modify>` command.  The value
options for the *sort* keyword are *no* or *device* similar to the *comm*
keywords above. If a value of *host* is used it will be automatically be
changed to *no* since the *sort* keyword does not support *host* mode. Not
all fix styles with extra atom data support *device* mode and in that case
a warning will be given and atom sorting will run in *no* mode instead.

.. versionadded:: 17Apr2024

The *atom/map* keyword determines whether the host or device builds the
atom_map, see the :doc:`atom_modify map <atom_modify>` command.  The
value options for the *atom/map* keyword are identical to the *sort*
keyword above.

The *gpu/aware* keyword chooses whether GPU-aware MPI will be used. When
this keyword is set to *on*, buffers in GPU memory are passed directly
through MPI send/receive calls. This reduces overhead of first copying
the data to the host CPU. However GPU-aware MPI is not supported on all
systems, which can lead to segmentation faults and would require using a
value of *off*\ . If LAMMPS can safely detect that GPU-aware MPI is not
available (currently only possible with OpenMPI v2.0.0 or later), then
the *gpu/aware* keyword is automatically set to *off* by default. When
the *gpu/aware* keyword is set to *off* while any of the *comm*
keywords are set to *device*, the value for these *comm* keywords will
be automatically changed to *no*\ . This setting has no effect if not
running on GPUs or if using only one MPI rank. GPU-aware MPI is available
for OpenMPI 1.8 (or later versions), Mvapich2 1.9 (or later) when the
"MV2_USE_CUDA" environment variable is set to "1", CrayMPI, and IBM
Spectrum MPI when the "-gpu" flag is used.

The *pair/only* keyword can change how the KOKKOS suffix "kk" is applied
when using an accelerator device.  By default device acceleration is always
used for all available styles.  With *pair/only* set to *on* the suffix
setting will choose device acceleration only for pair styles and run all
other force computations on the host CPU.  The *comm* flags, along with the
*sort* and *atom/map* keywords will also automatically be changed to *no*\ .
This can result in better performance for certain configurations and
system sizes.

----------

The *omp* style invokes settings associated with the use of the
OPENMP package.

The *Nthreads* argument sets the number of OpenMP threads allocated for
each MPI task.  For example, if your system has nodes with dual
quad-core processors, it has a total of 8 cores per node.  You could
use two MPI tasks per node (e.g. using the -ppn option of the mpirun
command in MPICH or -npernode in OpenMPI), and set *Nthreads* = 4.
This would use all 8 cores on each node.  Note that the product of MPI
tasks \* threads/task should not exceed the physical number of cores
(on a node), otherwise performance will suffer.

Setting *Nthreads* = 0 instructs LAMMPS to use whatever value is the
default for the given OpenMP environment. This is usually determined
via the *OMP_NUM_THREADS* environment variable or the compiler
runtime.  Note that in most cases the default for OpenMP capable
compilers is to use one thread for each available CPU core when
*OMP_NUM_THREADS* is not explicitly set, which can lead to poor
performance.

Here are examples of how to set the environment variable when
launching LAMMPS:

.. code-block:: bash

   env OMP_NUM_THREADS=4 lmp_machine -sf omp -in in.script
   env OMP_NUM_THREADS=2 mpirun -np 2 lmp_machine -sf omp -in in.script
   mpirun -x OMP_NUM_THREADS=2 -np 2 lmp_machine -sf omp -in in.script

or you can set it permanently in your shell's start-up script.
All three of these examples use a total of 4 CPU cores.

Note that different MPI implementations have different ways of passing
the OMP_NUM_THREADS environment variable to all MPI processes.  The
second example line above is for MPICH; the third example line with -x is
for OpenMPI.  Check your MPI documentation for additional details.

What combination of threads and MPI tasks gives the best performance
is difficult to predict and can depend on many components of your
input.  Not all features of LAMMPS support OpenMP threading via the
OPENMP package and the parallel efficiency can be very different,
too.

.. note::

   If you build LAMMPS with the GPU, INTEL, and / or OPENMP
   packages, be aware these packages all allow setting of the *Nthreads*
   value via their package commands, but there is only a single global
   *Nthreads* value used by OpenMP.  Thus if multiple package commands are
   invoked, you should ensure the values are consistent.  If they are
   not, the last one invoked will take precedence, for all packages.
   Also note that if the :doc:`-sf hybrid intel omp command-line switch <Run_options>` is used, it invokes a "package intel" command, followed by a
   "package omp" command, both with a setting of *Nthreads* = 0. Likewise
   for a hybrid suffix for gpu and omp. Note that KOKKOS also supports
   setting the number of OpenMP threads from the command line using the
   "-k on" :doc:`command-line switch <Run_options>`. The default for
   KOKKOS is 1 thread per MPI task, so any other number of threads should
   be explicitly set using the "-k on" command-line switch (and this
   setting should be consistent with settings from any other packages
   used).

Optional keyword/value pairs can also be specified.  Each has a
default value as listed below.

The *neigh* keyword specifies whether neighbor list building will be
multi-threaded in addition to force calculations.  If *neigh* is set
to *no* then neighbor list calculation is performed only by MPI tasks
with no OpenMP threading.  If *mode* is *yes* (the default), a
multi-threaded neighbor list build is used.  Using *neigh* = *yes* is
almost always faster and should produce identical neighbor lists at the
expense of using more memory.  Specifically, neighbor list pages are
allocated for all threads at the same time and each thread works
within its own pages.

----------

Restrictions
""""""""""""

This command cannot be used after the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command.

The *gpu* style of this command can only be invoked if LAMMPS was built
with the GPU package.  See the :doc:`Build package <Build_package>` doc
page for more info.

The *intel* style of this command can only be invoked if LAMMPS was
built with the INTEL package.  See the :doc:`Build package <Build_package>` page for more info.

The *kokkos* style of this command can only be invoked if LAMMPS was built
with the KOKKOS package.  See the :doc:`Build package <Build_package>`
doc page for more info.

The *omp* style of this command can only be invoked if LAMMPS was built
with the OPENMP package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`suffix <suffix>`, :doc:`-pk command-line switch <Run_options>`

Defaults
""""""""

For the GPU package, the default parameters and settings are:

.. parsed-literal::

   Ngpu = 0, neigh = yes, newton = off, binsize = 0.0, split = 1.0, gpuID = 0 to Ngpu-1, tpa = 1, omp = 0, platform=-1.

These settings are made automatically if the "-sf gpu"
:doc:`command-line switch <Run_options>` is used.  If it is not used,
you must invoke the package gpu command in your input script or via the
"-pk gpu" :doc:`command-line switch <Run_options>`.

For the INTEL package, the default parameters and settings are:

.. parsed-literal::

   Nphi = 1, omp = 0, mode = mixed, lrt = no, balance = -1, tpc = 4, tptask = 240, pppm_table = yes

The default ghost option is determined by the pair style being used.
This value is output to the screen in the offload report at the end of each
run.  Note that all of these settings, except "omp" and "mode", are ignored if
LAMMPS was not built with Xeon Phi co-processor support.  These settings are
made automatically if the "-sf intel" :doc:`command-line switch <Run_options>`
is used.  If it is not used, you must invoke the package intel command in your
input script or via the "-pk intel" :doc:`command-line switch <Run_options>`.

For the KOKKOS package when using GPUs, the option defaults are:

.. parsed-literal::

   neigh = full, neigh/qeq = full, newton = off, binsize = 2x LAMMPS default value, comm = device, sort = device, atom/map = device, neigh/transpose = off, gpu/aware = on

For GPUs, option neigh/thread = on when there are 16k atoms or less on
an MPI rank, otherwise it is "off". When LAMMPS can safely detect that
GPU-aware MPI is not available, the default value of gpu/aware becomes
"off".

For the KOKKOS package when using CPUs or Xeon Phis, the option defaults are:

.. parsed-literal::

   neigh = half, neigh/qeq = half, newton = on, binsize = 0.0, comm = no, sort = no, atom/map = no

These settings are made automatically by
the required "-k on" :doc:`command-line switch <Run_options>`.  You can
change them by using the package kokkos command in your input script or
via the :doc:`-pk kokkos command-line switch <Run_options>`.

For the OMP package, the defaults are

.. parsed-literal::

   Nthreads = 0, neigh = yes

These settings are made automatically if the "-sf omp"
:doc:`command-line switch <Run_options>` is used.  If it is not used,
you must invoke the package omp command in your input script or via the
"-pk omp" :doc:`command-line switch <Run_options>`.
