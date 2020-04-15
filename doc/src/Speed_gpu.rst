GPU package
===========

The GPU package was developed by Mike Brown while at SNL and ORNL
and his collaborators, particularly Trung Nguyen (now at Northwestern).
It provides GPU versions of many pair styles and for parts of the
:doc:`kspace_style pppm <kspace_style>` for long-range Coulombics.
It has the following general features:

* It is designed to exploit common GPU hardware configurations where one
  or more GPUs are coupled to many cores of one or more multi-core CPUs,
  e.g. within a node of a parallel machine.
* Atom-based data (e.g. coordinates, forces) are moved back-and-forth
  between the CPU(s) and GPU every timestep.
* Neighbor lists can be built on the CPU or on the GPU
* The charge assignment and force interpolation portions of PPPM can be
  run on the GPU.  The FFT portion, which requires MPI communication
  between processors, runs on the CPU.
* Force computations of different style (pair vs. bond/angle/dihedral/improper)
  can be performed concurrently on the GPU and CPU(s), respectively.
* It allows for GPU computations to be performed in single or double
  precision, or in mixed-mode precision, where pairwise forces are
  computed in single precision, but accumulated into double-precision
  force vectors.
* LAMMPS-specific code is in the GPU package.  It makes calls to a
  generic GPU library in the lib/gpu directory.  This library provides
  NVIDIA support as well as more general OpenCL support, so that the
  same functionality is supported on a variety of hardware.

**Required hardware/software:**

To compile and use this package in CUDA mode, you currently need
to have an NVIDIA GPU and install the corresponding NVIDIA CUDA
toolkit software on your system (this is primarily tested on Linux
and completely unsupported on Windows):

* Check if you have an NVIDIA GPU: cat /proc/driver/nvidia/gpus/\*/information
* Go to http://www.nvidia.com/object/cuda_get.html
* Install a driver and toolkit appropriate for your system (SDK is not necessary)
* Run lammps/lib/gpu/nvc_get_devices (after building the GPU library, see below) to
  list supported devices and properties

To compile and use this package in OpenCL mode, you currently need
to have the OpenCL headers and the (vendor neutral) OpenCL library installed.
In OpenCL mode, the acceleration depends on having an `OpenCL Installable Client Driver (ICD) <https://www.khronos.org/news/permalink/opencl-installable-client-driver-icd-loader>`_
installed. There can be multiple of them for the same or different hardware
(GPUs, CPUs, Accelerators) installed at the same time. OpenCL refers to those
as 'platforms'.  The GPU library will select the **first** suitable platform,
but this can be overridden using the device option of the :doc:`package <package>`
command. run lammps/lib/gpu/ocl_get_devices to get a list of available
platforms and devices with a suitable ICD available.

**Building LAMMPS with the GPU package:**

See the :ref:`Build extras <gpu>` doc page for
instructions.

**Run with the GPU package from the command line:**

The mpirun or mpiexec command sets the total number of MPI tasks used
by LAMMPS (one or multiple per compute node) and the number of MPI
tasks used per node.  E.g. the mpirun command in MPICH does this via
its -np and -ppn switches.  Ditto for OpenMPI via -np and -npernode.

When using the GPU package, you cannot assign more than one GPU to a
single MPI task.  However multiple MPI tasks can share the same GPU,
and in many cases it will be more efficient to run this way.  Likewise
it may be more efficient to use less MPI tasks/node than the available
# of CPU cores.  Assignment of multiple MPI tasks to a GPU will happen
automatically if you create more MPI tasks/node than there are
GPUs/mode.  E.g. with 8 MPI tasks/node and 2 GPUs, each GPU will be
shared by 4 MPI tasks.

Use the "-sf gpu" :doc:`command-line switch <Run_options>`, which will
automatically append "gpu" to styles that support it.  Use the "-pk
gpu Ng" :doc:`command-line switch <Run_options>` to set Ng = # of
GPUs/node to use.

.. code-block:: bash

   lmp_machine -sf gpu -pk gpu 1 -in in.script                         # 1 MPI task uses 1 GPU
   mpirun -np 12 lmp_machine -sf gpu -pk gpu 2 -in in.script           # 12 MPI tasks share 2 GPUs on a single 16-core (or whatever) node
   mpirun -np 48 -ppn 12 lmp_machine -sf gpu -pk gpu 2 -in in.script   # ditto on 4 16-core nodes

Note that if the "-sf gpu" switch is used, it also issues a default
:doc:`package gpu 1 <package>` command, which sets the number of
GPUs/node to 1.

Using the "-pk" switch explicitly allows for setting of the number of
GPUs/node to use and additional options.  Its syntax is the same as
same as the "package gpu" command.  See the :doc:`package <package>`
command doc page for details, including the default values used for
all its options if it is not specified.

Note that the default for the :doc:`package gpu <package>` command is to
set the Newton flag to "off" pairwise interactions.  It does not
affect the setting for bonded interactions (LAMMPS default is "on").
The "off" setting for pairwise interaction is currently required for
GPU package pair styles.

**Or run with the GPU package by editing an input script:**

The discussion above for the mpirun/mpiexec command, MPI tasks/node,
and use of multiple MPI tasks/GPU is the same.

Use the :doc:`suffix gpu <suffix>` command, or you can explicitly add an
"gpu" suffix to individual styles in your input script, e.g.

.. code-block:: LAMMPS

   pair_style lj/cut/gpu 2.5

You must also use the :doc:`package gpu <package>` command to enable the
GPU package, unless the "-sf gpu" or "-pk gpu" :doc:`command-line switches <Run_options>` were used.  It specifies the number of
GPUs/node to use, as well as other options.

**Speed-ups to expect:**

The performance of a GPU versus a multi-core CPU is a function of your
hardware, which pair style is used, the number of atoms/GPU, and the
precision used on the GPU (double, single, mixed). Using the GPU package
in OpenCL mode on CPUs (which uses vectorization and multithreading) is
usually resulting in inferior performance compared to using LAMMPS' native
threading and vectorization support in the USER-OMP and USER-INTEL packages.

See the `Benchmark page <https://lammps.sandia.gov/bench.html>`_ of the
LAMMPS web site for performance of the GPU package on various
hardware, including the Titan HPC platform at ORNL.

You should also experiment with how many MPI tasks per GPU to use to
give the best performance for your problem and machine.  This is also
a function of the problem size and the pair style being using.
Likewise, you should experiment with the precision setting for the GPU
library to see if single or mixed precision will give accurate
results, since they will typically be faster.

**Guidelines for best performance:**

* Using multiple MPI tasks per GPU will often give the best performance,
  as allowed my most multi-core CPU/GPU configurations.
* If the number of particles per MPI task is small (e.g. 100s of
  particles), it can be more efficient to run with fewer MPI tasks per
  GPU, even if you do not use all the cores on the compute node.
* The :doc:`package gpu <package>` command has several options for tuning
  performance.  Neighbor lists can be built on the GPU or CPU.  Force
  calculations can be dynamically balanced across the CPU cores and
  GPUs.  GPU-specific settings can be made which can be optimized
  for different hardware.  See the :doc:`package <package>` command
  doc page for details.
* As described by the :doc:`package gpu <package>` command, GPU
  accelerated pair styles can perform computations asynchronously with
  CPU computations. The "Pair" time reported by LAMMPS will be the
  maximum of the time required to complete the CPU pair style
  computations and the time required to complete the GPU pair style
  computations. Any time spent for GPU-enabled pair styles for
  computations that run simultaneously with :doc:`bond <bond_style>`,
  :doc:`angle <angle_style>`, :doc:`dihedral <dihedral_style>`,
  :doc:`improper <improper_style>`, and :doc:`long-range <kspace_style>`
  calculations will not be included in the "Pair" time.
* When the *mode* setting for the package gpu command is force/neigh,
  the time for neighbor list calculations on the GPU will be added into
  the "Pair" time, not the "Neigh" time.  An additional breakdown of the
  times required for various tasks on the GPU (data copy, neighbor
  calculations, force computations, etc) are output only with the LAMMPS
  screen output (not in the log file) at the end of each run.  These
  timings represent total time spent on the GPU for each routine,
  regardless of asynchronous CPU calculations.
* The output section "GPU Time Info (average)" reports "Max Mem / Proc".
  This is the maximum memory used at one time on the GPU for data
  storage by a single MPI process.

Restrictions
""""""""""""

None.
