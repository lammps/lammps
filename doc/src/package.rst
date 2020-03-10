.. index:: package

package command
===============

Syntax
""""""


.. parsed-literal::

   package style args

* style = *gpu* or *intel* or *kokkos* or *omp*
* args = arguments specific to the style

  .. parsed-literal::

       *gpu* args = Ngpu keyword value ...
         Ngpu = # of GPUs per node
         zero or more keyword/value pairs may be appended
         keywords = *neigh* or *newton* or *binsize* or *split* or *gpuID* or *tpa* or *device* or *blocksize*
           *neigh* value = *yes* or *no*
             yes = neighbor list build on GPU (default)
             no = neighbor list build on CPU
           *newton* = *off* or *on*
             off = set Newton pairwise flag off (default and required)
             on = set Newton pairwise flag on (currently not allowed)
           *binsize* value = size
             size = bin size for neighbor list construction (distance units)
           *split* = fraction
             fraction = fraction of atoms assigned to GPU (default = 1.0)
           *gpuID* values = first last
             first = ID of first GPU to be used on each node
             last = ID of last GPU to be used on each node
           *tpa* value = Nthreads
             Nthreads = # of GPU threads used per atom
           *device* value = device_type or platform_id:device_type or platform_id:custom,val1,val2,val3,..,val13
             platform_id = numerical OpenCL platform id (default: -1)
             device_type = *kepler* or *fermi* or *cypress* or *intel* or *phi* or *generic* or *custom*
             val1,val2,... = custom OpenCL tune parameters (see below for details)
           *blocksize* value = size
             size = thread block size for pair force computation
       *intel* args = NPhi keyword value ...
         Nphi = # of co-processors per node
         zero or more keyword/value pairs may be appended
         keywords = *mode* or *omp* or *lrt* or *balance* or *ghost* or *tpc* or *tptask* or *no_affinity*
           *mode* value = *single* or *mixed* or *double*
             single = perform force calculations in single precision
             mixed = perform force calculations in mixed precision
             double = perform force calculations in double precision
           *omp* value = Nthreads
             Nthreads = number of OpenMP threads to use on CPU (default = 0)
           *lrt* value = *yes* or *no*
             yes = use additional thread dedicated for some PPPM calculations
             no = do not dedicate an extra thread for some PPPM calculations
           *balance* value = split
             split = fraction of work to offload to co-processor, -1 for dynamic
           *ghost* value = *yes* or *no*
             yes = include ghost atoms for offload
             no = do not include ghost atoms for offload
           *tpc* value = Ntpc
             Ntpc = max number of co-processor threads per co-processor core (default = 4)
           *tptask* value = Ntptask
             Ntptask = max number of co-processor threads per MPI task (default = 240)
           *no_affinity* values = none
       *kokkos* args = keyword value ...
         zero or more keyword/value pairs may be appended
         keywords = *neigh* or *neigh/qeq* or *neigh/thread* or *newton* or *binsize* or *comm* or *comm/exchange* or *comm/forward* or *comm/reverse* or *cuda/aware*
           *neigh* value = *full* or *half*
             full = full neighbor list
             half = half neighbor list built in thread-safe manner
           *neigh/qeq* value = *full* or *half*
             full = full neighbor list
             half = half neighbor list built in thread-safe manner
           *neigh/thread* value = *off* or *on*
             off = thread only over atoms
             on = thread over both atoms and neighbors
           *newton* = *off* or *on*
             off = set Newton pairwise and bonded flags off
             on = set Newton pairwise and bonded flags on
           *binsize* value = size
             size = bin size for neighbor list construction (distance units)
           *comm* value = *no* or *host* or *device*
             use value for comm/exchange and comm/forward and comm/reverse
           *comm/exchange* value = *no* or *host* or *device*
           *comm/forward* value = *no* or *host* or *device*
           *comm/reverse* value = *no* or *host* or *device*
             no = perform communication pack/unpack in non-KOKKOS mode
             host = perform pack/unpack on host (e.g. with OpenMP threading)
             device = perform pack/unpack on device (e.g. on GPU)
           *cuda/aware* = *off* or *on*
             off = do not use CUDA-aware MPI
             on = use CUDA-aware MPI (default)
       *omp* args = Nthreads keyword value ...
         Nthread = # of OpenMP threads to associate with each MPI process
         zero or more keyword/value pairs may be appended
         keywords = *neigh*
           *neigh* value = *yes* or *no*
             yes = threaded neighbor list build (default)
             no = non-threaded neighbor list build



Examples
""""""""


.. parsed-literal::

   package gpu 1
   package gpu 1 split 0.75
   package gpu 2 split -1.0
   package gpu 1 device kepler
   package gpu 1 device 2:generic
   package gpu 1 device custom,32,4,8,256,11,128,256,128,32,64,8,128,128
   package kokkos neigh half comm device
   package omp 0 neigh no
   package omp 4
   package intel 1
   package intel 2 omp 4 mode mixed balance 0.5

Description
"""""""""""

This command invokes package-specific settings for the various
accelerator packages available in LAMMPS.  Currently the following
packages use settings from this command: GPU, USER-INTEL, KOKKOS, and
USER-OMP.

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

For the GPU, USER-INTEL, and USER-OMP packages, if a "-sf gpu" or "-sf
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

See the :doc:`Speed packages <Speed_packages>` doc page for more details
about using the various accelerator packages for speeding up LAMMPS
simulations.


----------


The *gpu* style invokes settings associated with the use of the GPU
package.

The *Ngpu* argument sets the number of GPUs per node.  There must be
at least as many MPI tasks per node as GPUs, as set by the mpirun or
mpiexec command.  If there are more MPI tasks (per node)
than GPUs, multiple MPI tasks will share each GPU.

Optional keyword/value pairs can also be specified.  Each has a
default value as listed below.

The *neigh* keyword specifies where neighbor lists for pair style
computation will be built.  If *neigh* is *yes*\ , which is the default,
neighbor list building is performed on the GPU.  If *neigh* is *no*\ ,
neighbor list building is performed on the CPU.  GPU neighbor list
building currently cannot be used with a triclinic box.  GPU neighbor
lists are not compatible with commands that are not GPU-enabled.  When
a non-GPU enabled command requires a neighbor list, it will also be
built on the CPU.  In these cases, it will typically be more efficient
to only use CPU neighbor list builds.

The *newton* keyword sets the Newton flags for pairwise (not bonded)
interactions to *off* or *on*\ , the same as the :doc:`newton <newton>`
command allows.  Currently, only an *off* value is allowed, since all
the GPU package pair styles require this setting.  This means more
computation is done, but less communication.  In the future a value of
*on* may be allowed, so the *newton* keyword is included as an option
for compatibility with the package command for other accelerator
styles.  Note that the newton setting for bonded interactions is not
affected by this keyword.

The *binsize* keyword sets the size of bins used to bin atoms in
neighbor list builds performed on the GPU, if *neigh* = *yes* is set.
If *binsize* is set to 0.0 (the default), then bins = the size of the
pairwise cutoff + neighbor skin distance.  This is 2x larger than the
LAMMPS default used for neighbor list building on the CPU.  This will
be close to optimal for the GPU, so you do not normally need to use
this keyword.  Note that if you use a longer-than-usual pairwise
cutoff, e.g. to allow for a smaller fraction of KSpace work with a
:doc:`long-range Coulombic solver <kspace_style>` because the GPU is
faster at performing pairwise interactions, then it may be optimal to
make the *binsize* smaller than the default.  For example, with a
cutoff of 20\*sigma in LJ :doc:`units <units>` and a neighbor skin
distance of sigma, a *binsize* = 5.25\*sigma can be more efficient than
the default.

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


.. parsed-literal::

   mpirun -np 32 -sf gpu -in in.script    # launch command
   package gpu 2 split -1                 # input script command

In this case, all CPU cores and GPU devices on the nodes would be
utilized.  Each GPU device would be shared by 4 CPU cores. The CPU
cores would perform force calculations for some fraction of the
particles at the same time the GPUs performed force calculation for
the other particles.

The *gpuID* keyword allows selection of which GPUs on each node will
be used for a simulation.  The *first* and *last* values specify the
GPU IDs to use (from 0 to Ngpu-1).  By default, first = 0 and last =
Ngpu-1, so that all GPUs are used, assuming Ngpu is set to the number
of physical GPUs.  If you only wish to use a subset, set Ngpu to a
smaller number and first/last to a sub-range of the available GPUs.

The *tpa* keyword sets the number of GPU thread per atom used to
perform force calculations.  With a default value of 1, the number of
threads will be chosen based on the pair style, however, the value can
be set explicitly with this keyword to fine-tune performance.  For
large cutoffs or with a small number of particles per GPU, increasing
the value can improve performance. The number of threads per atom must
be a power of 2 and currently cannot be greater than 32.

The *device* keyword can be used to tune parameters optimized for a
specific accelerator and platform when using OpenCL. OpenCL supports
the concept of a **platform**\ , which represents one or more devices that
share the same driver (e.g. there would be a different platform for
GPUs from different vendors or for CPU based accelerator support).
In LAMMPS only one platform can be active at a time and by default
the first platform with an accelerator is selected. This is equivalent
to using a platform ID of -1. The platform ID is a number corresponding
to the output of the ocl\_get\_devices tool. The platform ID is passed
to the GPU library, by prefixing the *device* keyword with that number
separated by a colon. For CUDA, the *device* keyword is ignored.
Currently, the device tuning support is limited to NVIDIA Kepler, NVIDIA
Fermi, AMD Cypress, Intel x86\_64 CPU, Intel Xeon Phi, or a generic device.
More devices may be added later.  The default device type can be
specified when building LAMMPS with the GPU library, via setting a
variable in the lib/gpu/Makefile that is used.

In addition, a device type *custom* is available, which is followed by
13 comma separated numbers, which allows to set those tweakable parameters
from the package command. It can be combined with the (colon separated)
platform id. The individual settings are:

* MEM\_THREADS
* THREADS\_PER\_ATOM
* THREADS\_PER\_CHARGE
* BLOCK\_PAIR
* MAX\_SHARED\_TYPES
* BLOCK\_NBOR\_BUILD
* BLOCK\_BIO\_PAIR
* BLOCK\_ELLIPSE
* WARP\_SIZE
* PPPM\_BLOCK\_1D
* BLOCK\_CELL\_2D
* BLOCK\_CELL\_ID
* MAX\_BIO\_SHARED\_TYPES

The *blocksize* keyword allows you to tweak the number of threads used
per thread block. This number should be a multiple of 32 (for GPUs)
and its maximum depends on the specific GPU hardware. Typical choices
are 64, 128, or 256. A larger block size increases occupancy of
individual GPU cores, but reduces the total number of thread blocks,
thus may lead to load imbalance.


----------


The *intel* style invokes settings associated with the use of the
USER-INTEL package.  All of its settings, except the *omp* and *mode*
keywords, are ignored if LAMMPS was not built with Xeon Phi
co-processor support.  All of its settings, including the *omp* and
*mode* keyword are applicable if LAMMPS was built with co-processor
support.

The *Nphi* argument sets the number of co-processors per node.
This can be set to any value, including 0, if LAMMPS was not
built with co-processor support.

Optional keyword/value pairs can also be specified.  Each has a
default value as listed below.

The *omp* keyword determines the number of OpenMP threads allocated
for each MPI task when any portion of the interactions computed by a
USER-INTEL pair style are run on the CPU.  This can be the case even
if LAMMPS was built with co-processor support; see the *balance*
keyword discussion below.  If you are running with less MPI tasks/node
than there are CPUs, it can be advantageous to use OpenMP threading on
the CPUs.

.. note::

   The *omp* keyword has nothing to do with co-processor threads on
   the Xeon Phi; see the *tpc* and *tptask* keywords below for a
   discussion of co-processor threads.

The *Nthread* value for the *omp* keyword sets the number of OpenMP
threads allocated for each MPI task.  Setting *Nthread* = 0 (the
default) instructs LAMMPS to use whatever value is the default for the
given OpenMP environment. This is usually determined via the
*OMP\_NUM\_THREADS* environment variable or the compiler runtime, which
is usually a value of 1.

For more details, including examples of how to set the OMP\_NUM\_THREADS
environment variable, see the discussion of the *Nthreads* setting on
this doc page for the "package omp" command.  Nthreads is a required
argument for the USER-OMP package.  Its meaning is exactly the same
for the USER-INTEL package.

.. note::

   If you build LAMMPS with both the USER-INTEL and USER-OMP
   packages, be aware that both packages allow setting of the *Nthreads*
   value via their package commands, but there is only a single global
   *Nthreads* value used by OpenMP.  Thus if both package commands are
   invoked, you should insure the two values are consistent.  If they are
   not, the last one invoked will take precedence, for both packages.
   Also note that if the :doc:`-sf hybrid intel omp command-line switch <Run_options>` is used, it invokes a "package intel"
   command, followed by a "package omp" command, both with a setting of
   *Nthreads* = 0.

The *mode* keyword determines the precision mode to use for
computing pair style forces, either on the CPU or on the co-processor,
when using a USER-INTEL supported :doc:`pair style <pair_style>`.  It
can take a value of *single*\ , *mixed* which is the default, or
*double*\ .  *Single* means single precision is used for the entire
force calculation.  *Mixed* means forces between a pair of atoms are
computed in single precision, but accumulated and stored in double
precision, including storage of forces, torques, energies, and virial
quantities.  *Double* means double precision is used for the entire
force calculation.

The *lrt* keyword can be used to enable "Long Range Thread (LRT)"
mode. It can take a value of *yes* to enable and *no* to disable.
LRT mode generates an extra thread (in addition to any OpenMP threads
specified with the OMP\_NUM\_THREADS environment variable or the *omp*
keyword). The extra thread is dedicated for performing part of the
:doc:`PPPM solver <kspace_style>` computations and communications. This
can improve parallel performance on processors supporting
Simultaneous Multithreading (SMT) such as Hyper-Threading (HT) on Intel
processors. In this mode, one additional thread is generated per MPI
process. LAMMPS will generate a warning in the case that more threads
are used than available in SMT hardware on a node. If the PPPM solver
from the USER-INTEL package is not used, then the LRT setting is
ignored and no extra threads are generated. Enabling LRT will replace
the :doc:`run_style <run_style>` with the *verlet/lrt/intel* style that
is identical to the default *verlet* style aside from supporting the
LRT feature. This feature requires setting the pre-processor flag
-DLMP\_INTEL\_USELRT in the makefile when compiling LAMMPS.

The *balance* keyword sets the fraction of :doc:`pair style <pair_style>` work offloaded to the co-processor for split
values between 0.0 and 1.0 inclusive.  While this fraction of work is
running on the co-processor, other calculations will run on the host,
including neighbor and pair calculations that are not offloaded, as
well as angle, bond, dihedral, kspace, and some MPI communications.
If *split* is set to -1, the fraction of work is dynamically adjusted
automatically throughout the run.  This typically give performance
within 5 to 10 percent of the optimal fixed fraction.

The *ghost* keyword determines whether or not ghost atoms, i.e. atoms
at the boundaries of processor sub-domains, are offloaded for neighbor
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

The *no\_affinity* keyword will turn off automatic setting of core
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
doesn't require atomic operations in the calculation of pair forces. For
that reason, *full* is the default setting for GPUs. However, when
running on CPUs, a *half* neighbor list is the default because it are
often faster, just as it is for non-accelerated pair styles. Similarly,
the *neigh/qeq* keyword determines how neighbor lists are built for :doc:`fix qeq/reax/kk <fix_qeq_reax>`. If not explicitly set, the value of
*neigh/qeq* will match *neigh*\ .

If the *neigh/thread* keyword is set to *off*\ , then the KOKKOS package
threads only over atoms. However, for small systems, this may not expose
enough parallelism to keep a GPU busy. When this keyword is set to *on*\ ,
the KOKKOS package threads over both atoms and neighbors of atoms. When
using *neigh/thread* *on*\ , a full neighbor list must also be used. Using
*neigh/thread* *on* may be slower for large systems, so this this option
is turned on by default only when there are 16K atoms or less owned by
an MPI rank and when using a full neighbor list. Not all KOKKOS-enabled
potentials support this keyword yet, and only thread over atoms. Many
simple pair-wise potentials such as Lennard-Jones do support threading
over both atoms and neighbors.

The *newton* keyword sets the Newton flags for pairwise and bonded
interactions to *off* or *on*\ , the same as the :doc:`newton <newton>`
command allows. The default for GPUs is *off* because this will almost
always give better performance for the KOKKOS package. This means more
computation is done, but less communication. However, when running on
CPUs a value of *on* is the default since it can often be faster, just
as it is for non-accelerated pair styles

The *binsize* keyword sets the size of bins used to bin atoms in
neighbor list builds. The same value can be set by the :doc:`neigh_modify binsize <neigh_modify>` command. Making it an option in the package
kokkos command allows it to be set from the command line. The default
value for CPUs is 0.0, which means the LAMMPS default will be used,
which is bins = 1/2 the size of the pairwise cutoff + neighbor skin
distance. This is fine when neighbor lists are built on the CPU. For GPU
builds, a 2x larger binsize equal to the pairwise cutoff + neighbor skin
is often faster, which is the default. Note that if you use a
longer-than-usual pairwise cutoff, e.g. to allow for a smaller fraction
of KSpace work with a :doc:`long-range Coulombic solver <kspace_style>`
because the GPU is faster at performing pairwise interactions, then this
rule of thumb may give too large a binsize and the default should be
overridden with a smaller value.

The *comm* and *comm/exchange* and *comm/forward* and *comm/reverse*
keywords determine whether the host or device performs the packing and
unpacking of data when communicating per-atom data between processors.
"Exchange" communication happens only on timesteps that neighbor lists
are rebuilt. The data is only for atoms that migrate to new processors.
"Forward" communication happens every timestep. "Reverse" communication
happens every timestep if the *newton* option is on. The data is for
atom coordinates and any other atom properties that needs to be updated
for ghost atoms owned by each processor.

The *comm* keyword is simply a short-cut to set the same value for both
the *comm/exchange* and *comm/forward* and *comm/reverse* keywords.

The value options for all 3 keywords are *no* or *host* or *device*\ . A
value of *no* means to use the standard non-KOKKOS method of
packing/unpacking data for the communication. A value of *host* means to
use the host, typically a multi-core CPU, and perform the
packing/unpacking in parallel with threads. A value of *device* means to
use the device, typically a GPU, to perform the packing/unpacking
operation.

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
to be move between the host and device anyway, so it is typically faster
to let the host handle communication, by using the *host* value. Using
*host* instead of *no* will enable use of multiple threads to
pack/unpack communicated data. When running small systems on a GPU,
performing the exchange pack/unpack on the host CPU can give speedup
since it reduces the number of CUDA kernel launches.

The *cuda/aware* keyword chooses whether CUDA-aware MPI will be used. When
this keyword is set to *on*\ , buffers in GPU memory are passed directly
through MPI send/receive calls. This reduces overhead of first copying
the data to the host CPU. However CUDA-aware MPI is not supported on all
systems, which can lead to segmentation faults and would require using a
value of *off*\ . If LAMMPS can safely detect that CUDA-aware MPI is not
available (currently only possible with OpenMPI v2.0.0 or later), then
the *cuda/aware* keyword is automatically set to *off* by default. When
the *cuda/aware* keyword is set to *off* while any of the *comm*
keywords are set to *device*\ , the value for these *comm* keywords will
be automatically changed to *host*\ . This setting has no effect if not
running on GPUs or if using only one MPI rank. CUDA-aware MPI is available
for OpenMPI 1.8 (or later versions), Mvapich2 1.9 (or later) when the
"MV2\_USE\_CUDA" environment variable is set to "1", CrayMPI, and IBM
Spectrum MPI when the "-gpu" flag is used.


----------


The *omp* style invokes settings associated with the use of the
USER-OMP package.

The *Nthread* argument sets the number of OpenMP threads allocated for
each MPI task.  For example, if your system has nodes with dual
quad-core processors, it has a total of 8 cores per node.  You could
use two MPI tasks per node (e.g. using the -ppn option of the mpirun
command in MPICH or -npernode in OpenMPI), and set *Nthreads* = 4.
This would use all 8 cores on each node.  Note that the product of MPI
tasks \* threads/task should not exceed the physical number of cores
(on a node), otherwise performance will suffer.

Setting *Nthread* = 0 instructs LAMMPS to use whatever value is the
default for the given OpenMP environment. This is usually determined
via the *OMP\_NUM\_THREADS* environment variable or the compiler
runtime.  Note that in most cases the default for OpenMP capable
compilers is to use one thread for each available CPU core when
*OMP\_NUM\_THREADS* is not explicitly set, which can lead to poor
performance.

Here are examples of how to set the environment variable when
launching LAMMPS:


.. parsed-literal::

   env OMP_NUM_THREADS=4 lmp_machine -sf omp -in in.script
   env OMP_NUM_THREADS=2 mpirun -np 2 lmp_machine -sf omp -in in.script
   mpirun -x OMP_NUM_THREADS=2 -np 2 lmp_machine -sf omp -in in.script

or you can set it permanently in your shell's start-up script.
All three of these examples use a total of 4 CPU cores.

Note that different MPI implementations have different ways of passing
the OMP\_NUM\_THREADS environment variable to all MPI processes.  The
2nd example line above is for MPICH; the 3rd example line with -x is
for OpenMPI.  Check your MPI documentation for additional details.

What combination of threads and MPI tasks gives the best performance
is difficult to predict and can depend on many components of your
input.  Not all features of LAMMPS support OpenMP threading via the
USER-OMP package and the parallel efficiency can be very different,
too.

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

The gpu style of this command can only be invoked if LAMMPS was built
with the GPU package.  See the :doc:`Build package <Build_package>` doc
page for more info.

The intel style of this command can only be invoked if LAMMPS was
built with the USER-INTEL package.  See the :doc:`Build package <Build_package>` doc page for more info.

The kk style of this command can only be invoked if LAMMPS was built
with the KOKKOS package.  See the :doc:`Build package <Build_package>`
doc page for more info.

The omp style of this command can only be invoked if LAMMPS was built
with the USER-OMP package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`suffix <suffix>`, :doc:`-pk command-line switch <Run_options>`

Default
"""""""

For the GPU package, the default is Ngpu = 1 and the option defaults
are neigh = yes, newton = off, binsize = 0.0, split = 1.0, gpuID = 0
to Ngpu-1, tpa = 1, and device = not used.  These settings are made
automatically if the "-sf gpu" :doc:`command-line switch <Run_options>`
is used.  If it is not used, you must invoke the package gpu command
in your input script or via the "-pk gpu" :doc:`command-line switch <Run_options>`.

For the USER-INTEL package, the default is Nphi = 1 and the option
defaults are omp = 0, mode = mixed, lrt = no, balance = -1, tpc = 4,
tptask = 240.  The default ghost option is determined by the pair
style being used.  This value is output to the screen in the offload
report at the end of each run.  Note that all of these settings,
except "omp" and "mode", are ignored if LAMMPS was not built with Xeon
Phi co-processor support.  These settings are made automatically if the
"-sf intel" :doc:`command-line switch <Run_options>` is used.  If it is
not used, you must invoke the package intel command in your input
script or via the "-pk intel" :doc:`command-line switch <Run_options>`.

For the KOKKOS package, the option defaults for GPUs are neigh = full,
neigh/qeq = full, newton = off, binsize for GPUs = 2x LAMMPS default
value, comm = device, cuda/aware = on. When LAMMPS can safely detect
that CUDA-aware MPI is not available, the default value of cuda/aware
becomes "off". For CPUs or Xeon Phis, the option defaults are neigh =
half, neigh/qeq = half, newton = on, binsize = 0.0, and comm = no. The
option neigh/thread = on when there are 16K atoms or less on an MPI
rank, otherwise it is "off". These settings are made automatically by
the required "-k on" :doc:`command-line switch <Run_options>`. You can
change them by using the package kokkos command in your input script or
via the :doc:`-pk kokkos command-line switch <Run_options>`.

For the OMP package, the default is Nthreads = 0 and the option
defaults are neigh = yes.  These settings are made automatically if
the "-sf omp" :doc:`command-line switch <Run_options>` is used.  If it
is not used, you must invoke the package omp command in your input
script or via the "-pk omp" :doc:`command-line switch <Run_options>`.
