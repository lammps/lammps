Packages with extra build options
=================================

When building with some packages, additional steps may be required,
in addition to

.. list-table::
   :align: center
   :header-rows: 1
   :widths: 50 50
   :width: 80%

   * - CMake build
     - Traditional make
   * - .. code-block:: bash

          cmake -D PKG_NAME=yes

     - .. code-block:: bash

          make yes-name

as described on the :doc:`Build_package <Build_package>` page.

For a CMake build there may be additional optional or required
variables to set.

.. versionchanged:: TBD

The traditional build system with GNU make no longer supports packages
that require extra steps in the ``lammps/lib`` directory.

This is the list of packages that may require additional steps.

.. this list must be kept in sync with its counterpart in Build_package.rst
.. table_from_list::
   :columns: 6

   * :ref:`ADIOS <adios>`
   * :ref:`APIP <apip>`
   * :ref:`COLVARS <colvar>`
   * :ref:`COMPRESS <compress>`
   * :ref:`ELECTRODE <electrode>`
   * :ref:`GPU <gpu>`
   * :ref:`H5MD <h5md>`
   * :ref:`INTEL <intel>`
   * :ref:`KIM <kim>`
   * :ref:`KOKKOS <kokkos>`
   * :ref:`LEPTON <lepton>`
   * :ref:`MACHDYN <machdyn>`
   * :ref:`MDI <mdi>`
   * :ref:`MISC <misc>`
   * :ref:`ML-HDNNP <ml-hdnnp>`
   * :ref:`ML-IAP <mliap>`
   * :ref:`ML-PACE <ml-pace>`
   * :ref:`ML-POD <ml-pod>`
   * :ref:`ML-QUIP <ml-quip>`
   * :ref:`MOLFILE <molfile>`
   * :ref:`NETCDF <netcdf>`
   * :ref:`OPENMP <openmp>`
   * :ref:`OPT <opt>`
   * :ref:`PLUMED <plumed>`
   * :ref:`PYTHON <python>`
   * :ref:`QMMM <qmmm>`
   * :ref:`RHEO <rheo>`
   * :ref:`SCAFACOS <scafacos>`
   * :ref:`VORONOI <voronoi>`
   * :ref:`VTK <vtk>`

----------

.. _compress:

COMPRESS package
----------------

To build with this package you must have the `zlib compression library
<https://zlib.net>`_ available on your system to build dump styles with
a ``/gz`` suffix.  There are also styles using the
`Zstandard <https://facebook.github.io/zstd/>`_ library which have a
'/zstd' suffix.  The zstd library version must be at least 1.4.  Older
versions use an incompatible API and thus LAMMPS will fail to compile.

.. tabs::

   .. tab:: CMake build

      If CMake cannot find the zlib library or include files, you can set
      these variables:

      .. code-block:: bash

         -D ZLIB_INCLUDE_DIR=path    # path to zlib.h header file
         -D ZLIB_LIBRARY=path        # path to libz.a (.so) file

      Support for Zstandard compression is auto-detected and for that
      CMake depends on the `pkg-config
      <https://www.freedesktop.org/wiki/Software/pkg-config/>`_ tool to
      identify the necessary flags to compile with this library, so the
      corresponding ``libzstandard.pc`` file must be in a folder where
      ``pkg-config`` can find it, which may require adding it to the
      ``PKG_CONFIG_PATH`` environment variable.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The COMPRESS package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _gpu:

GPU package
---------------------

To build with this package, you must choose options for precision and
which GPU hardware to build for. The GPU package currently supports
three different types of back ends: OpenCL, CUDA and HIP.

CMake build
^^^^^^^^^^^

.. code-block:: bash

   -D GPU_API=value             # value = opencl (default) or cuda or hip
   -D GPU_PREC=value            # precision setting
                                # value = double or mixed (default) or single
   -D GPU_ARCH=value            # primary GPU hardware choice for GPU_API=cuda
                                # value = sm_XX (see below, default is sm_75)
   -D GPU_DEBUG=value           # enable debug code in the GPU package library,
                                # mostly useful for developers
                                # value = yes or no (default)
   -D HIP_PATH=value            # value = path to HIP installation. Must be set if
                                # GPU_API=HIP
   -D HIP_ARCH=value            # primary GPU hardware choice for GPU_API=hip
                                # value depends on selected HIP_PLATFORM
                                # default is 'gfx906' for HIP_PLATFORM=amd and 'sm_75' for
                                # HIP_PLATFORM=nvcc
   -D HIP_USE_DEVICE_SORT=value # enables GPU sorting
                                # value = yes (default) or no
   -D CUDPP_OPT=value           # use GPU binning with CUDA (should be off for modern GPUs)
                                # enables CUDA Performance Primitives, must be "no" for
                                # CUDA_MPS_SUPPORT=yes
                                # value = yes or no (default)
   -D CUDA_MPS_SUPPORT=value    # enables some tweaks required to run with active
                                # nvidia-cuda-mps daemon
                                # value = yes or no (default)
   -D CUDA_BUILD_MULTIARCH=value  # enables building CUDA kernels for all supported GPU
                                  # architectures
                                  # value = yes (default) or no
   -D USE_STATIC_OPENCL_LOADER=value  # downloads/includes OpenCL ICD loader library,
                                      # no local OpenCL headers/libs needed
                                      # value = yes (default) or no

The GPU package supports 3 precision modes: single, double, and mixed, with
the latter being the default.  In the double precision mode, atom positions,
forces and energies are stored, computed and accumulated in double precision.
In the mixed precision mode, forces and energies are accumulated in double precision
while atom coordinates are stored and arithmetic operations are performed
in single precision. In the single precision mode, all are stored, executed
and accumulated in single precision.

To specify the precision mode (output to the screen before LAMMPS runs for
verification), set ``GPU_PREC`` to one of ``single``, ``double``, or ``mixed``.

Some accelerators or OpenCL implementations only support single precision.
This mode should be used with care and appropriate validation as the errors
can scale with system size in this implementation. This can be useful for
accelerating test runs when setting up a simulation for production runs on
another machine. In the case where only single precision is supported, either
LAMMPS must be compiled with ``-DFFT_SINGLE`` to use PPPM with GPU acceleration
or GPU acceleration should be disabled for PPPM (e.g. suffix off or ``pair/only``
as described in the LAMMPS documentation).

``GPU_ARCH`` settings for different GPU hardware is as follows:

* ``sm_30`` for Kepler (supported since CUDA 5 and until CUDA 10.x)
* ``sm_35`` or ``sm_37`` for Kepler (supported since CUDA 5 and until CUDA 11.x)
* ``sm_50`` or ``sm_52`` for Maxwell (supported since CUDA 6)
* ``sm_60`` or ``sm_61`` for Pascal (supported since CUDA 8)
* ``sm_70`` for Volta (supported since CUDA 9)
* ``sm_75`` for Turing (supported since CUDA 10)
* ``sm_80`` or ``sm_86`` for Ampere (supported since CUDA 11, ``sm_86`` since CUDA 11.1)
* ``sm_89`` for Lovelace (supported since CUDA 11.8)
* ``sm_90`` or ``sm_90a`` for Hopper (supported since CUDA 12.0)
* ``sm_100`` or ``sm_103`` for Blackwell B100/B200/B300 (supported since CUDA 12.8)
* ``sm_120`` for Blackwell B20x/B40 (supported since CUDA 12.8)
* ``sm_121`` for Blackwell (supported since CUDA 12.9)

A more detailed list can be found, for example,
at `Wikipedia's CUDA article <https://en.wikipedia.org/wiki/CUDA#GPUs_supported>`_

CMake can detect which version of the CUDA toolkit is used and thus will
try to include support for **all** major GPU architectures supported by
this toolkit.  Thus the ``GPU_ARCH`` setting is merely an optimization, to
have code for the preferred GPU architecture directly included rather
than having to wait for the JIT compiler of the CUDA driver to translate
it.  This behavior can be turned off (e.g. to speed up compilation) by
setting ``CUDA_ENABLE_MULTIARCH`` to ``no``.

When compiling for CUDA or HIP with CUDA, version 8.0 or later of the
CUDA toolkit is required and a GPU architecture of Kepler or later,
which must *also* be supported by the CUDA toolkit in use **and** the
CUDA driver in use.  When compiling for OpenCL, OpenCL version 1.2 or
later is required and the GPU must be supported by the GPU driver and
OpenCL runtime bundled with the driver.

Please note that the GPU library accesses the CUDA driver library
directly, so it needs to be linked with the CUDA driver library
(``libcuda.so``) that ships with the Nvidia driver.  If you are
compiling LAMMPS on the head node of a GPU cluster, this library may not
be installed, so you may need to copy it over from one of the compute
nodes (best into this directory).  Recent versions of the CUDA toolkit
starting from CUDA 9 provide a dummy ``libcuda.so`` library (typically
under ``$(CUDA_HOME)/lib64/stubs``), that can be used for linking.

To support the CUDA multi-process server (MPS) you can set the define
``-DCUDA_MPS_SUPPORT``.  Please note that in this case you must **not**
use the CUDA performance primitives and thus set the variable
``CUDPP_OPT`` to empty.

If you are compiling for OpenCL, the default setting is to download,
build, and link with a static OpenCL ICD loader library and standard
OpenCL headers.  This way no local OpenCL development headers or library
needs to be present and only OpenCL compatible drivers need to be
installed to use OpenCL.  If this is not desired, you can set
``USE_STATIC_OPENCL_LOADER`` to ``no``.

If ``GERYON_NUMA_FISSION`` is defined at build time (``-DGPU_DEBUG=no``),
LAMMPS will consider separate NUMA nodes on GPUs or accelerators as
separate devices.  For example, a 2-socket CPU would appear as two separate
devices for OpenCL (and LAMMPS would require two MPI processes to use both
sockets with the GPU library - each with its own device ID as output by
ocl_get_devices).  OpenCL version 1.2 or later is required.

If you are compiling with HIP, note that before running CMake you will
have to set appropriate environment variables. Some variables such as
``HCC_AMDGPU_TARGET`` (for ROCm <= 4.0) or ``CUDA_PATH`` are
necessary for ``hipcc`` and the linker to work correctly.

When compiling for HIP ROCm, GPU sorting with ``-D
HIP_USE_DEVICE_SORT=on`` requires installing the ``hipcub`` library
(https://github.com/ROCmSoftwarePlatform/hipCUB).  The HIP CUDA-backend
additionally requires cub (https://nvlabs.github.io/cub).  Setting
``-DDOWNLOAD_CUB=yes`` will download and compile CUB.

The GPU library has some multi-thread support using OpenMP.  If LAMMPS
is built with ``-D BUILD_OMP=on`` this will also be enabled.

For a debug build, set ``GPU_DEBUG`` to be ``yes``.

.. versionadded:: 3Aug2022

Using the CHIP-SPV implementation of HIP is supported. It allows one to
run HIP code on Intel GPUs via the OpenCL or Level Zero back ends. To use
CHIP-SPV, you must set ``-DHIP_USE_DEVICE_SORT=OFF`` in your CMake
command-line as CHIP-SPV does not yet support hipCUB. As of Summer 2022,
the use of HIP for Intel GPUs is experimental. You should only use this
option in preparations to run on Aurora system at Argonne.

.. code:: bash

   # AMDGPU target (ROCm <= 4.0)
   export HIP_PLATFORM=hcc
   export HIP_PATH=/path/to/HIP/install
   export HCC_AMDGPU_TARGET=gfx906
   cmake -D PKG_GPU=on -D GPU_API=HIP -D HIP_ARCH=gfx906 -D CMAKE_CXX_COMPILER=hipcc ..
   make -j 4

.. code:: bash

   # AMDGPU target (ROCm >= 4.1)
   export HIP_PLATFORM=amd
   export HIP_PATH=/path/to/HIP/install
   cmake -D PKG_GPU=on -D GPU_API=HIP -D HIP_ARCH=gfx906 -D CMAKE_CXX_COMPILER=hipcc ..
   make -j 4

.. code:: bash

   # CUDA target (not recommended, use GPU_API=cuda)
   # !!! DO NOT set CMAKE_CXX_COMPILER !!!
   export HIP_PLATFORM=nvcc
   export HIP_PATH=/path/to/HIP/install
   export CUDA_PATH=/usr/local/cuda
   cmake -D PKG_GPU=on -D GPU_API=HIP -D HIP_ARCH=sm_70 ..
   make -j 4

.. code:: bash

   # SPIR-V target (Intel GPUs)
   export HIP_PLATFORM=spirv
   export HIP_PATH=/path/to/HIP/install
   export CMAKE_CXX_COMPILER=<hipcc/clang++>
   cmake -D PKG_GPU=on -D GPU_API=HIP ..
   make -j 4

----------

.. _kim:

KIM package
---------------------

To build with this package, the KIM library with API v2 must be downloaded
and built on your system. It must include the KIM models that you want to
use with LAMMPS.

If you would like to use the :doc:`kim query <kim_commands>`
command, you also need to have libcurl installed with the matching
development headers and the curl-config tool.

If you would like to use the :doc:`kim property <kim_commands>`
command, you need to build LAMMPS with the PYTHON package installed
and linked to Python 3.6 or later. See the :ref:`PYTHON package build info <python>`
for more details on this. After successfully building LAMMPS with Python, you
also need to install the ``kim-property`` Python package, which can be easily
done using *pip* as ``pip install kim-property``, or from the *conda-forge*
channel as ``conda install kim-property`` if LAMMPS is built in Conda. More
detailed information is available at:
`kim-property installation <https://github.com/openkim/kim-property#installing-kim-property>`_.

In addition to installing the KIM API, it is also necessary to install the
library of KIM models (interatomic potentials).
See `Obtaining KIM Models <https://openkim.org/doc/usage/obtaining-models>`_ to
learn how to install a pre-build binary of the OpenKIM Repository of Models.
See the list of all KIM models here: https://openkim.org/browse/models

(Also note that when downloading and installing from source
the KIM API library with all its models, may take a long time (tens of
minutes to hours) to build.  Of course you only need to do that once.)

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_KIM=value           # download OpenKIM API v2 for build
                                         # value = no (default) or yes
         -D LMP_DEBUG_CURL=value         # set libcurl verbose mode on/off
                                         # value = off (default) or on
         -D LMP_NO_SSL_CHECK=value       # tell libcurl to not verify the peer
                                         # value = no (default) or yes
         -D KIM_EXTRA_UNITTESTS=value    # enables extra unit tests
                                         # value = no (default) or yes

      If ``DOWNLOAD_KIM`` is set to ``yes`` (or ``on``), the KIM API library
      will be downloaded and built inside the CMake build directory.  If
      the KIM library is already installed on your system (in a location
      where CMake cannot find it), you may need to set the
      ``PKG_CONFIG_PATH`` environment variable so that libkim-api can be
      found, or run the command ``source kim-api-activate``.

      Extra unit tests can only be available if they are explicitly requested
      (``KIM_EXTRA_UNITTESTS`` is set to ``yes`` (or ``on``)) and the prerequisites
      are met. See :ref:`KIM Extra unit tests <kim_extra_unittests>` for
      more details on this.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The KIM package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

Debugging OpenKIM web queries in LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``LMP_DEBUG_CURL`` is set, the libcurl verbose mode will be turned
on, and any libcurl calls within the KIM web query display a lot of
information about libcurl operations.  You hardly ever want this set in
production use, you will almost always want this when you debug or
report problems.

The libcurl library performs peer SSL certificate verification by
default.  This verification is done using a CA certificate store that
the SSL library can use to make sure the peer's server certificate is
valid.  If SSL reports an error ("certificate verify failed") during the
handshake and thus refuses further communicate with that server, you can
set ``LMP_NO_SSL_CHECK`` to override that behavior.  When LAMMPS is
compiled with ``LMP_NO_SSL_CHECK`` set, libcurl does not verify the peer
and connection attempts will succeed regardless of the names in the
certificate. This option is insecure.  As an alternative, you can
specify your own CA cert path by setting the environment variable
``CURL_CA_BUNDLE`` to the path of your choice.  A call to the KIM web
query would get this value from the environment variable.

.. _kim_extra_unittests:

KIM Extra unit tests (CMake only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During development, testing, or debugging, if
:doc:`unit testing <Build_development>` is enabled in LAMMPS, one can also
enable extra tests on :doc:`KIM commands <kim_commands>` by setting the
``KIM_EXTRA_UNITTESTS`` to ``yes`` (or ``on``).

Enabling the extra unit tests have some requirements,

* It requires to have internet access.
* It requires to have libcurl installed with the matching development headers
  and the curl-config tool.
* It requires to build LAMMPS with the PYTHON package installed and linked to
  Python 3.6 or later. See the :ref:`PYTHON package build info <python>` for
  more details on this.
* It requires to have ``kim-property`` Python package installed, which can be
  easily done using *pip* as ``pip install kim-property``, or from the
  *conda-forge* channel as ``conda install kim-property`` if LAMMPS is built in
  Conda. More detailed information is available at:
  `kim-property installation <https://github.com/openkim/kim-property#installing-kim-property>`_.
* It is also necessary to install the following KIM models:

  * ``EAM_Dynamo_MendelevAckland_2007v3_Zr__MO_004835508849_000``
  * ``EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005``
  * ``LennardJones612_UniversalShifted__MO_959249795837_003``

  See `Obtaining KIM Models <https://openkim.org/doc/usage/obtaining-models>`_
  to learn how to install a pre-built binary of the OpenKIM Repository of
  Models or see
  `Installing KIM Models <https://openkim.org/doc/usage/obtaining-models/#installing_models>`_
  to learn how to install the specific KIM models.

----------

.. _kokkos:

KOKKOS package
--------------

Using the KOKKOS package requires choosing several settings.  You have
to select whether you want to compile with parallelization on the host
and whether you want to include offloading of calculations to a device
(e.g. a GPU).  The default setting is to have no host parallelization
and no device offloading.  In addition, you can select the hardware
architecture to select the instruction set.  Since most hardware is
backward compatible, you may choose settings for an older architecture
to have an executable that will run on this and newer architectures.

.. note::

   If you run Kokkos on a different GPU architecture than what LAMMPS
   was compiled with, there will be a delay during device initialization
   while the just-in-time compiler is recompiling all GPU kernels for
   the new hardware.  This is, however, only supported for GPUs of the
   **same** major hardware version and different minor hardware versions,
   e.g. 5.0 and 5.2 but not 5.2 and 6.0.  LAMMPS will abort with an
   error message indicating a mismatch, if the major version differs.

The settings discussed below have been tested with LAMMPS and are
confirmed to work.  Kokkos is an active project with ongoing improvements
and projects working on including support for additional architectures.
More information on Kokkos can be found on the
`Kokkos GitHub project <https://github.com/kokkos>`_.

Available Architecture settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are the possible choices for the Kokkos architecture ID.
They must be specified in uppercase.

.. list-table::
   :header-rows: 0
   :widths: auto

   *  - **Arch-ID**
      - **HOST or GPU**
      - **Description**
   *  - NATIVE
      - HOST
      - Local machine
   *  - AMDAVX
      - HOST
      - AMD chip
   *  - ARMV80
      - HOST
      - ARMv8.0 Compatible CPU
   *  - ARMV81
      - HOST
      - ARMv8.1 Compatible CPU
   *  - ARMV8_THUNDERX
      - HOST
      - ARMv8 Cavium ThunderX CPU
   *  - ARMV8_THUNDERX2
      - HOST
      - ARMv8 Cavium ThunderX2 CPU
   *  - A64FX
      - HOST
      - ARMv8.2 with SVE Support
   *  - ARMV9_GRACE
      - HOST
      - ARMv9 NVIDIA Grace CPU
   *  - SNB
      - HOST
      - Intel Sandy/Ivy Bridge CPUs
   *  - HSW
      - HOST
      - Intel Haswell CPUs
   *  - BDW
      - HOST
      - Intel Broadwell Xeon E-class CPUs
   *  - ICL
      - HOST
      - Intel Ice Lake Client CPUs (AVX512)
   *  - ICX
      - HOST
      - Intel Ice Lake Xeon Server CPUs (AVX512)
   *  - SKL
      - HOST
      - Intel Skylake Client CPUs
   *  - SKX
      - HOST
      - Intel Skylake Xeon Server CPUs (AVX512)
   *  - KNC
      - HOST
      - Intel Knights Corner Xeon Phi
   *  - KNL
      - HOST
      - Intel Knights Landing Xeon Phi
   *  - SPR
      - HOST
      - Intel Sapphire Rapids Xeon Server CPUs (AVX512)
   *  - POWER8
      - HOST
      - IBM POWER8 CPUs
   *  - POWER9
      - HOST
      - IBM POWER9 CPUs
   *  - ZEN
      - HOST
      - AMD Zen architecture
   *  - ZEN2
      - HOST
      - AMD Zen2 architecture
   *  - ZEN3
      - HOST
      - AMD Zen3 architecture
   *  - ZEN4
      - HOST
      - AMD Zen4 architecture
   *  - ZEN5
      - HOST
      - AMD Zen5 architecture
   *  - RISCV_SG2042
      - HOST
      - SG2042 (RISC-V) CPUs
   *  - RISCV_RVA22V
      - HOST
      - RVA22V (RISC-V) CPUs
   *  - KEPLER30
      - GPU
      - NVIDIA Kepler generation CC 3.0
   *  - KEPLER32
      - GPU
      - NVIDIA Kepler generation CC 3.2
   *  - KEPLER35
      - GPU
      - NVIDIA Kepler generation CC 3.5
   *  - KEPLER37
      - GPU
      - NVIDIA Kepler generation CC 3.7
   *  - MAXWELL50
      - GPU
      - NVIDIA Maxwell generation CC 5.0
   *  - MAXWELL52
      - GPU
      - NVIDIA Maxwell generation CC 5.2
   *  - MAXWELL53
      - GPU
      - NVIDIA Maxwell generation CC 5.3
   *  - PASCAL60
      - GPU
      - NVIDIA Pascal generation CC 6.0
   *  - PASCAL61
      - GPU
      - NVIDIA Pascal generation CC 6.1
   *  - VOLTA70
      - GPU
      - NVIDIA Volta generation CC 7.0
   *  - VOLTA72
      - GPU
      - NVIDIA Volta generation CC 7.2
   *  - TURING75
      - GPU
      - NVIDIA Turing generation CC 7.5
   *  - AMPERE80
      - GPU
      - NVIDIA Ampere generation CC 8.0
   *  - AMPERE86
      - GPU
      - NVIDIA Ampere generation CC 8.6
   *  - ADA89
      - GPU
      - NVIDIA Ada generation CC 8.9
   *  - HOPPER90
      - GPU
      - NVIDIA Hopper generation CC 9.0
   *  - BLACKWELL100
      - GPU
      - NVIDIA Blackwell generation CC 10.0
   *  - BLACKWELL120
      - GPU
      - NVIDIA Blackwell generation CC 12.0
   *  - AMD_GFX906
      - GPU
      - AMD GPU MI50/60
   *  - AMD_GFX908
      - GPU
      - AMD GPU MI100
   *  - AMD_GFX90A
      - GPU
      - AMD GPU MI200
   *  - AMD_GFX940
      - GPU
      - AMD GPU MI300
   *  - AMD_GFX942
      - GPU
      - AMD GPU MI300
   *  - AMD_GFX942_APU
      - GPU
      - AMD APU MI300A
   *  - AMD_GFX1030
      - GPU
      - AMD GPU V620/W6800
   *  - AMD_GFX1100
      - GPU
      - AMD GPU RX7900XTX
   *  - AMD_GFX1103
      - GPU
      - AMD APU Phoenix
   *  - INTEL_GEN
      - GPU
      - SPIR64-based devices, e.g. Intel GPUs, using JIT
   *  - INTEL_DG1
      - GPU
      - Intel Iris XeMAX GPU
   *  - INTEL_GEN9
      - GPU
      - Intel GPU Gen9
   *  - INTEL_GEN11
      - GPU
      - Intel GPU Gen11
   *  - INTEL_GEN12LP
      - GPU
      - Intel GPU Gen12LP
   *  - INTEL_XEHP
      - GPU
      - Intel GPU Xe-HP
   *  - INTEL_PVC
      - GPU
      - Intel GPU Ponte Vecchio
   *  - INTEL_DG2
      - GPU
      - Intel GPU DG2

This list was last updated for version 4.6.2 of the Kokkos library.

.. tabs::

   .. tab:: Basic CMake build settings:

      For multicore CPUs using OpenMP, set these 2 variables.

      .. code-block:: bash

         -D Kokkos_ARCH_HOSTARCH=yes  # HOSTARCH = HOST from list above
         -D Kokkos_ENABLE_OPENMP=yes
         -D BUILD_OMP=yes

      Please note that enabling OpenMP for KOKKOS requires that OpenMP is
      also :ref:`enabled for the rest of LAMMPS <serial>`.

      For Intel KNLs using OpenMP, set these variables:

      .. code-block:: bash

         -D Kokkos_ARCH_KNL=yes
         -D Kokkos_ENABLE_OPENMP=yes

      For NVIDIA GPUs using CUDA, set these variables:

      .. code-block:: bash

         -D Kokkos_ARCH_HOSTARCH=yes   # HOSTARCH = HOST from list above
         -D Kokkos_ARCH_GPUARCH=yes    # GPUARCH = GPU from list above
         -D Kokkos_ENABLE_CUDA=yes
         -D Kokkos_ENABLE_OPENMP=yes

      This will also enable executing FFTs on the GPU, either via the
      internal KISSFFT library, or - by preference - with the cuFFT
      library bundled with the CUDA toolkit, depending on whether CMake
      can identify its location.

      For AMD or NVIDIA GPUs using HIP, set these variables:

      .. code-block:: bash

         -D Kokkos_ARCH_HOSTARCH=yes   # HOSTARCH = HOST from list above
         -D Kokkos_ARCH_GPUARCH=yes    # GPUARCH = GPU from list above
         -D Kokkos_ENABLE_HIP=yes
         -D Kokkos_ENABLE_OPENMP=yes

      This will enable FFTs on the GPU, either by the internal KISSFFT library
      or with the hipFFT wrapper library, which will call out to the
      platform-appropriate vendor library: rocFFT on AMD GPUs or cuFFT on
      NVIDIA GPUs.

      For Intel GPUs using SYCL, set these variables:

      .. code-block:: bash

         -D Kokkos_ARCH_HOSTARCH=yes   # HOSTARCH = HOST from list above
         -D Kokkos_ARCH_GPUARCH=yes    # GPUARCH = GPU from list above
         -D Kokkos_ENABLE_SYCL=yes
         -D Kokkos_ENABLE_OPENMP=yes
         -D FFT_KOKKOS=MKL_GPU

      This will enable FFTs on the GPU using the oneMKL library.

      To simplify compilation, seven preset files are included in the
      ``cmake/presets`` folder, ``kokkos-serial.cmake``,
      ``kokkos-openmp.cmake``, ``kokkos-cuda.cmake``,
      ``kokkos-cuda-nowrapper.cmake``, ``kokkos-hip.cmake``,
      ``kokkos-sycl-nvidia.cmake``, and ``kokkos-sycl-intel.cmake``.
      They will enable the KOKKOS package and enable some hardware
      choices.  For GPU support those preset files may need to be
      customized to match the hardware used.  For some platforms,
      e.g. CUDA, the Kokkos library will try to auto-detect a suitable
      configuration.  So to compile with CUDA device parallelization
      with some common packages enabled, you can do the following:

      .. code-block:: bash

         mkdir build-kokkos-cuda
         cd build-kokkos-cuda
         cmake -C ../cmake/presets/basic.cmake \
               -C ../cmake/presets/kokkos-cuda-nowrapper.cmake ../cmake
         cmake --build .

      The ``kokkos-openmp.cmake`` preset can be combined with any of the
      others, but it is not possible to combine multiple GPU
      acceleration settings (CUDA, HIP, SYCL) into a single executable.

   .. tab:: Basic traditional make settings:

      Choose which hardware to support in ``Makefile.machine`` via
      ``KOKKOS_DEVICES`` and ``KOKKOS_ARCH`` settings.  See the
      ``src/MAKE/OPTIONS/Makefile.kokkos*`` files for examples.

      For multicore CPUs using OpenMP:

      .. code-block:: make

         KOKKOS_DEVICES = OpenMP
         KOKKOS_ARCH = HOSTARCH          # HOSTARCH = HOST from list above

      For Intel KNLs using OpenMP:

      .. code-block:: make

         KOKKOS_DEVICES = OpenMP
         KOKKOS_ARCH = KNL

      For NVIDIA GPUs using CUDA:

      .. code-block:: make

         KOKKOS_DEVICES = Cuda
         KOKKOS_ARCH = HOSTARCH,GPUARCH  # HOSTARCH = HOST from list above that is
                                         #            hosting the GPU
                                         # GPUARCH = GPU from list above
         KOKKOS_CUDA_OPTIONS = "enable_lambda"
         FFT_INC = -DFFT_CUFFT           # enable use of cuFFT (optional)
         FFT_LIB = -lcufft               # link to cuFFT library

      For GPUs, you also need the following lines in your
      ``Makefile.machine`` before the CC line is defined.  They tell
      ``mpicxx`` to use an ``nvcc`` compiler wrapper, which will use
      ``nvcc`` for compiling CUDA files and a C++ compiler for
      non-Kokkos, non-CUDA files.

      .. code-block:: make

         # For OpenMPI
         KOKKOS_ABSOLUTE_PATH = $(shell cd $(KOKKOS_PATH); pwd)
         export OMPI_CXX = $(KOKKOS_ABSOLUTE_PATH)/config/nvcc_wrapper
         CC = mpicxx

      .. code-block:: make

         # For MPICH and derivatives
         KOKKOS_ABSOLUTE_PATH = $(shell cd $(KOKKOS_PATH); pwd)
         CC = mpicxx -cxx=$(KOKKOS_ABSOLUTE_PATH)/config/nvcc_wrapper

      For AMD or NVIDIA GPUs using HIP:

      .. code-block:: make

         KOKKOS_DEVICES = HIP
         KOKKOS_ARCH = HOSTARCH,GPUARCH  # HOSTARCH = HOST from list above that is
                                         #            hosting the GPU
                                         # GPUARCH = GPU from list above
         FFT_INC = -DFFT_HIPFFT          # enable use of hipFFT (optional)
         FFT_LIB = -lhipfft              # link to hipFFT library

      For Intel GPUs using SYCL:

      .. code-block:: make

         KOKKOS_DEVICES = SYCL
         KOKKOS_ARCH = HOSTARCH,GPUARCH  # HOSTARCH = HOST from list above that is
                                         #            hosting the GPU
                                         # GPUARCH = GPU from list above
         FFT_INC = -DFFT_KOKKOS_MKL_GPU  # enable use of oneMKL for Intel GPUs (optional)
                                         # link to oneMKL FFT library
         FFT_LIB = -lmkl_sycl_dft -lmkl_intel_ilp64 -lmkl_tbb_thread -mkl_core -ltbb

Advanced KOKKOS compilation settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are other allowed options when building with the KOKKOS package
that can improve performance or assist in debugging or profiling. Below
are some examples that may be useful in combination with LAMMPS.  For
the full list (which keeps changing as the Kokkos package itself evolves),
please consult the Kokkos library documentation.

As alternative to using multi-threading via OpenMP
(``-DKokkos_ENABLE_OPENMP=on`` or ``KOKKOS_DEVICES=OpenMP``) it is also
possible to use Posix threads directly (``-DKokkos_ENABLE_PTHREAD=on``
or ``KOKKOS_DEVICES=Pthread``).  While binding of threads to individual
or groups of CPU cores is managed in OpenMP with environment variables,
you need assistance from either the "hwloc" or "libnuma" library for the
Pthread thread parallelization option. To enable use with CMake:
``-DKokkos_ENABLE_HWLOC=on`` or ``-DKokkos_ENABLE_LIBNUMA=on``; and with
conventional make: ``KOKKOS_USE_TPLS=hwloc`` or
``KOKKOS_USE_TPLS=libnuma``.

The CMake option ``-DKokkos_ENABLE_LIBRT=on`` or the makefile setting
``KOKKOS_USE_TPLS=librt`` enables the use of a more accurate timer
mechanism on many Unix-like platforms for internal profiling.

The CMake option ``-DKokkos_ENABLE_DEBUG=on`` or the makefile setting
``KOKKOS_DEBUG=yes`` enables printing of run-time
debugging information that can be useful. It also enables runtime
bounds checking on Kokkos data structures.  As to be expected, enabling
this option will negatively impact the performance and thus is only
recommended when developing a Kokkos-enabled style in LAMMPS.

The CMake option ``-DKokkos_ENABLE_CUDA_UVM=on`` or the makefile
setting ``KOKKOS_CUDA_OPTIONS=enable_lambda,force_uvm`` enables the
use of CUDA "Unified Virtual Memory" (UVM) in Kokkos.  UVM allows to
transparently use RAM on the host to supplement the memory used on the
GPU (with some performance penalty) and thus enables running larger
problems that would otherwise not fit into the RAM on the GPU.

The CMake option ``-D KOKKOS_PREC=value`` sets the floating point
precision of the calculations, where ``value`` can be one of:
``double`` (FP64, default) or ``mixed`` (FP64 for accumulation of
forces, energy, and virial, FP32 otherwise) or ``single`` (FP32).
Similarly the makefile settings ``-DLMP_KOKKOS_DOUBLE_DOUBLE``
(default), ``-DLMP_KOKKOS_SINGLE_DOUBLE``, and
``-DLMP_KOKKOS_SINGLE_SINGLE`` set double, mixed, single precision
respectively. When using reduced precision (single or mixed), the
simulation should be carefully checked to ensure it is stable and that
energy is acceptably conserved.

The CMake option ``-D KOKKOS_LAYOUT=value`` sets the array layout of
Kokkos views (e.g. forces, velocities, etc.) on GPUs, where ``value``
can be one of: ``legacy`` (mostly LayoutRight, default) or ``default``
(mostly LayoutLeft). Similarly the makefile settings
``-DLMP_KOKKOS_LAYOUT_LEGACY`` (default) and
``-DLMP_KOKKOS_LAYOUT_DEFAULT`` set legacy or default layouts
respectively. Using the default layout (LayoutLeft) can give speedup
on GPUs for some models, but a slowdown for others. LayoutRight is
always used for positions on GPUs since it has been found to be
faster, and when compiling exclusively for CPUs.

----------

.. _lepton:

LEPTON package
--------------

To build with this package, you must build the Lepton library which is
included in the LAMMPS source distribution in the ``lib/lepton`` folder.

.. tabs::

   .. tab:: CMake build

      This is the recommended build procedure for using Lepton in
      LAMMPS. No additional settings are normally needed besides
      ``-D PKG_LEPTON=yes``.

      On x86 hardware the Lepton library will also include a just-in-time
      compiler for faster execution.  This is auto detected but can
      be explicitly disabled by setting ``-D LEPTON_ENABLE_JIT=no``
      (or enabled by setting it to yes).

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The LEPTON package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _machdyn:

MACHDYN package
-------------------------------

To build with this package, you must download the Eigen3 library.
Eigen3 is a template library, so you do not need to build it.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_EIGEN3            # download Eigen3, value = no (default) or yes
         -D EIGEN3_INCLUDE_DIR=path    # path to Eigen library (only needed if a
                                       # custom location)

      If ``DOWNLOAD_EIGEN3`` is set, the Eigen3 library will be
      downloaded and inside the CMake build directory.  If the Eigen3
      library is already on your system (in a location where CMake
      cannot find it), set ``EIGEN3_INCLUDE_DIR`` to the directory the
      ``Eigen3`` include file is in.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The MACHDYN package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _mliap:

ML-IAP package
---------------------------

Building the ML-IAP package requires including the :ref:`ML-SNAP
<PKG-ML-SNAP>` package.  There will be an error message if this requirement
is not satisfied.  Using the *mliappy* model also requires enabling
Python support, which in turn requires to include the :ref:`PYTHON
<PKG-PYTHON>` package **and** requires to have the `cython
<https://cython.org>`_ software installed and with it a working
``cythonize`` command.  This feature requires compiling LAMMPS with
Python version 3.6 or later.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D MLIAP_ENABLE_PYTHON=value   # enable mliappy model (default is autodetect)

      Without this setting, CMake will check whether it can find a
      suitable Python version and the ``cythonize`` command and choose
      the default accordingly.  During the build procedure the provided
      .pyx file(s) will be automatically translated to C++ code and compiled.
      Please do **not** run ``cythonize`` manually in the ``src/ML-IAP`` folder,
      as that can lead to compilation errors if Python support is not enabled.
      If you did it by accident, please remove the generated .cpp and .h files.

   .. tab:: Traditional make

      The build uses the ``lib/python/Makefile.mliap_python`` file in the
      compile/link process to add a rule to update the files generated by
      the ``cythonize`` command in case the corresponding .pyx file(s) were
      modified.  You may need to modify ``lib/python/Makefile.lammps``
      if the LAMMPS build fails.

      To enable building the ML-IAP package with Python support enabled,
      you need to add ``-DMLIAP_PYTHON`` to the ``LMP_INC`` variable in
      your machine makefile.  You may have to manually run the
      ``cythonize`` command on .pyx file(s) in the ``src`` folder, if
      this is not automatically done during installing the ML-IAP
      package.  Please do **not** run ``cythonize`` in the ``src/ML-IAP``
      folder, as that can lead to compilation errors if Python support
      is not enabled.  If you did this by accident, please remove the
      generated .cpp and .h files.

----------

.. _opt:

OPT package
---------------------

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_OPT=yes``

   .. tab:: Traditional make

      The compiler flag ``-restrict`` must be used to build LAMMPS with
      the OPT package when using Intel compilers.  It should be added to
      the ``CCFLAGS`` line of your ``Makefile.machine``.  See
      ``src/MAKE/OPTIONS/Makefile.opt`` for an example.

----------

.. _python:

PYTHON package
---------------------------

Building with the PYTHON package requires you have a the Python
development headers and library available on your system, which
needs to be Python version 3.6 or later.  See ``lib/python/README``
for additional details.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D Python_EXECUTABLE=path   # path to Python executable to use

      Without this setting, CMake will guess the default Python version
      on your system.  To use a different Python version, you can either
      create a virtualenv, activate it and then run cmake.  Or you can
      set the Python_EXECUTABLE variable to specify which Python
      interpreter should be used.  Note note that you will also need to
      have the development headers installed for this version,
      e.g. python3-devel.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The PYTHON package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _voronoi:

VORONOI package
-----------------------------

To build with this package, you must download and build the
`Voro++ library <https://math.lbl.gov/voro++/>`_ or install a
binary package provided by your operating system.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_VORO=value    # download Voro++ for build
                                   # value = no (default) or yes
         -D VORO_LIBRARY=path      # Voro++ library file
                                   # (only needed if at custom location)
         -D VORO_INCLUDE_DIR=path  # Voro++ include directory
                                   # (only needed if at custom location)

      If ``DOWNLOAD_VORO`` is set, the Voro++ library will be downloaded
      and built inside the CMake build directory.  If the Voro++ library
      is already on your system (in a location CMake cannot find it),
      ``VORO_LIBRARY`` is the filename (plus path) of the Voro++ library
      file, not the directory the library file is in.
      ``VORO_INCLUDE_DIR`` is the directory the Voro++ include file is
      in.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The VORONOI package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _adios:

ADIOS package
-----------------------------------

The ADIOS package requires the `ADIOS I/O library
<https://github.com/ornladios/ADIOS2>`_, version 2.3.1 or newer. Make
sure that you have ADIOS built either with or without MPI to match if
you build LAMMPS with or without MPI.  ADIOS compilation settings for
LAMMPS are automatically detected, if the PATH and LD_LIBRARY_PATH
environment variables have been updated for the local ADIOS installation
and the instructions below are followed for the respective build
systems.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D ADIOS2_DIR=path        # path is where ADIOS 2.x is installed
         -D PKG_ADIOS=yes

   .. tab:: Traditional make

      Turn on the ADIOS package before building LAMMPS. If the
      ADIOS 2.x software is installed in PATH, there is nothing else to
      do:

      .. code-block:: bash

         make yes-adios

      otherwise, set ADIOS2_DIR environment variable when turning on the package:

      .. code-block:: bash

         ADIOS2_DIR=path make yes-adios   # path is where ADIOS 2.x is installed

----------

.. _apip:

APIP package
-----------------------------

The APIP package depends on the library of the
:ref:`ML-PACE <ml-pace>` package.
The code for the library can be found
at: `https://github.com/ICAMS/lammps-user-pace/ <https://github.com/ICAMS/lammps-user-pace/>`_

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_APIP=yes``
      and ``-D PKG_ML-PACE=yes``.
      One can use a local version of the ML-PACE library instead of
      automatically downloading the library as described :ref:`here <ml-pace>`.


   .. tab:: Traditional make

      .. versionchanged:: TBD

      The APIP package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _colvar:

COLVARS package
---------------

This package enables the use of the `Colvars <https://colvars.github.io/>`_
module included in the LAMMPS source distribution.


.. tabs::

   .. tab:: CMake build

      This is the recommended build procedure for using Colvars in
      LAMMPS. No additional settings are normally needed besides ``-D
      PKG_COLVARS=yes``. The following CMake variables are available.

      .. code-block:: bash

         -D PKG_COLVARS=yes          # enable the package itself
         -D COLVARS_LEPTON=yes       # use the Lepton library for custom expression (on by defaul)
         -D COLVARS_DEBUG=no         # eneable debugging message (verbose, off by default)

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The COLVARS package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _electrode:

ELECTRODE package
-----------------

This package depends on the KSPACE package.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D PKG_ELECTRODE=yes          # enable the package itself
         -D PKG_KSPACE=yes             # the ELECTRODE package requires KSPACE
         -D USE_INTERNAL_LINALG=value  #

      Features in the ELECTRODE package are dependent on code in the
      KSPACE package so the latter one *must* be enabled.

      The ELECTRODE package also requires LAPACK (and BLAS) and CMake
      can identify their locations and pass that info to the ELECTRODE
      build script.  But on some systems this may cause problems when
      linking or the dependency is not desired.  Try enabling
      ``USE_INTERNAL_LINALG`` in those cases to use the bundled linear
      algebra library and work around the limitation.

   .. tab:: Traditional make

      The ELECTRODE package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _ml-pace:

ML-PACE package
-----------------------------

This package requires a library that can be downloaded and built
in lib/pace or somewhere else, which must be done before building
LAMMPS with this package. The code for the library can be found
at: `https://github.com/ICAMS/lammps-user-pace/ <https://github.com/ICAMS/lammps-user-pace/>`_

Instead of including the ML-PACE package directly into LAMMPS, it
is also possible to skip this step and build the ML-PACE package as
a plugin using the CMake script files in the ``examples/PACKAGE/pace/plugin``
folder and then load this plugin at runtime with the :doc:`plugin command <plugin>`.

.. tabs::

   .. tab:: CMake build

      By default the library will be downloaded from the git repository
      and built automatically when the ML-PACE package is enabled with
      ``-D PKG_ML-PACE=yes``.  The location for the sources may be
      customized by setting the variable ``PACELIB_URL`` when
      configuring with CMake (e.g. to use a local archive on machines
      without internet access).  Since CMake checks the validity of the
      archive with ``md5sum`` you may also need to set ``PACELIB_MD5``
      if you provide a different library version than what is downloaded
      automatically.


   .. tab:: Traditional make

      .. versionchanged:: TBD

      The ML-PACE package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _ml-pod:

ML-POD package
-----------------------------

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_ML-POD=yes``.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The ML-POD package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _ml-quip:

ML-QUIP package
---------------------------------

To build with this package, you must download and build the QUIP
library.  It can be obtained from GitHub.  For support of GAP
potentials, additional files with specific licensing conditions need
to be downloaded and configured.  The automatic download will from
within CMake will download the non-commercial use version.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_QUIP=value       # download QUIP library for build
                                      # value = no (default) or yes
         -D QUIP_LIBRARY=path         # path to libquip.a
                                      # (only needed if a custom location)
         -D USE_INTERNAL_LINALG=value # Use the internal linear algebra library
                                      # instead of LAPACK
                                      # value = no (default) or yes

      CMake will try to download and build the QUIP library from GitHub,
      if it is not found on the local machine. This requires to have git
      installed. It will use the same compilers and flags as used for
      compiling LAMMPS.  Currently this is only supported for the GNU
      and the Intel compilers. Set the ``QUIP_LIBRARY`` variable if you
      want to use a previously compiled and installed QUIP library and
      CMake cannot find it.

      The QUIP library requires LAPACK (and BLAS) and CMake can identify
      their locations and pass that info to the QUIP build script. But
      on some systems this triggers a (current) limitation of CMake and
      the configuration will fail. Try enabling ``USE_INTERNAL_LINALG`` in
      those cases to use the bundled linear algebra library and work around
      the limitation.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The ML-QUIP package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _plumed:

PLUMED package
-------------------------------------

.. _plumedinstall: https://www.plumed.org/doc-v2.9/user-doc/html/_installation.html

Before building LAMMPS with this package, you must first build PLUMED.
PLUMED can be built as part of the LAMMPS build or installed separately
from LAMMPS using the generic `PLUMED installation instructions <plumedinstall_>`_.
The PLUMED package has been tested to work with Plumed versions
2.4.x, to 2.9.x and will error out, when trying to run calculations
with a different version of the Plumed kernel.

PLUMED can be linked into MD codes in three different modes: static,
shared, and runtime.  With the "static" mode, all the code that PLUMED
requires is linked statically into LAMMPS. LAMMPS is then fully
independent from the PLUMED installation, but you have to rebuild/relink
it in order to update the PLUMED code inside it.  With the "shared"
linkage mode, LAMMPS is linked to a shared library that contains the
PLUMED code.  This library should preferably be installed in a globally
accessible location. When PLUMED is linked in this way the same library
can be used by multiple MD packages.  Furthermore, the PLUMED library
LAMMPS uses can be updated without the need for a recompile of LAMMPS
for as long as the shared PLUMED library is ABI-compatible.

The third linkage mode is "runtime" which allows the user to specify
which PLUMED kernel should be used at runtime by using the PLUMED_KERNEL
environment variable. This variable should point to the location of the
libplumedKernel.so dynamical shared object, which is then loaded at
runtime. This mode of linking is particularly convenient for doing
PLUMED development and comparing multiple PLUMED versions as these sorts
of comparisons can be done without recompiling the hosting MD code. All
three linkage modes are supported by LAMMPS on selected operating
systems (e.g. Linux) and using either CMake or traditional make
build. The "static" mode should be the most portable, while the
"runtime" mode support in LAMMPS makes the most assumptions about
operating system and compiler environment. If one mode does not work,
try a different one, switch to a different build system, consider a
global PLUMED installation or consider downloading PLUMED during the
LAMMPS build.

Instead of including the PLUMED package directly into LAMMPS, it
is also possible to skip this step and build the PLUMED package as
a plugin using the CMake script files in the ``examples/PACKAGE/plumed/plugin``
folder and then load this plugin at runtime with the :doc:`plugin command <plugin>`.

.. tabs::

   .. tab:: CMake build

      When the ``-D PKG_PLUMED=yes`` flag is included in the cmake
      command you must ensure that `the GNU Scientific Library (GSL)
      <https://www.gnu.org/software/gsl/>` is installed in locations
      that are accessible in your environment.  There are then two
      additional variables that control the manner in which PLUMED is
      obtained and linked into LAMMPS.

      .. code-block:: bash

         -D DOWNLOAD_PLUMED=value   # download PLUMED for build
                                    # value = no (default) or yes
         -D PLUMED_MODE=value       # Linkage mode for PLUMED
                                    # value = static (default), shared,
                                    #         or runtime

      If ``DOWNLOAD_PLUMED`` is set to ``yes``, the PLUMED library will be
      downloaded (the version of PLUMED that will be downloaded is
      hard-coded to a vetted version of PLUMED, usually a recent stable
      release version) and built inside the CMake build directory.  If
      ``DOWNLOAD_PLUMED`` is set to "no" (the default), CMake will try
      to detect and link to an installed version of PLUMED.  For this to
      work, the PLUMED library has to be installed into a location where
      the ``pkg-config`` tool can find it or the ``PKG_CONFIG_PATH``
      environment variable has to be set up accordingly.  PLUMED should
      be installed in such a location if you compile it using the
      default make; make install commands.

      The ``PLUMED_MODE`` setting determines the linkage mode for the
      PLUMED library.  The allowed values for this flag are "static"
      (default), "shared", or "runtime".  If you want to switch the
      linkage mode, just re-run CMake with a different setting. For a
      discussion of PLUMED linkage modes, please see above.  When
      ``DOWNLOAD_PLUMED`` is enabled the static linkage mode is
      recommended.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The PLUMED package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _h5md:

H5MD package
---------------------------------

To build with this package you must have the HDF5 software package
installed on your system, which should include the h5cc compiler and
the HDF5 library.

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_H5MD=yes``.

      This should auto-detect the H5MD library on your system.  Several
      advanced CMake H5MD options exist if you need to specify where it
      is installed.  Use the ccmake (terminal window) or cmake-gui
      (graphical) tools to see these options and set them interactively
      from their user interfaces.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The H5MD package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _ml-hdnnp:

ML-HDNNP package
----------------

To build with the ML-HDNNP package it is required to download and build the
external `n2p2 <https://github.com/CompPhysVienna/n2p2>`_ library ``v2.1.4``
(or higher). The LAMMPS build process offers an automatic download and
compilation of *n2p2* or allows you to choose the installation directory of
*n2p2* manually. Please see the boxes below for the CMake and traditional build
system for detailed information.

In case of a manual installation of *n2p2* you only need to build the *n2p2* core
library ``libnnp`` and interface library ``libnnpif``. When using GCC it should
suffice to execute ``make libnnpif`` in the *n2p2* ``src`` directory. For more
details please see ``lib/hdnnp/README`` and the `n2p2 build documentation
<https://compphysvienna.github.io/n2p2/topics/build.html>`_.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_N2P2=value    # download n2p2 for build
                                   # value = no (default) or yes
         -D N2P2_DIR=path          # n2p2 base directory
                                   # (only needed if a custom location)

      If ``DOWNLOAD_N2P2`` is set, the *n2p2* library will be downloaded and
      built inside the CMake build directory.  If the *n2p2* library is already
      on your system (in a location CMake cannot find it), set the ``N2P2_DIR``
      to path where *n2p2* is located. If *n2p2* is located directly in
      ``lib/hdnnp/n2p2`` it will be automatically found by CMake.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The ML-HDNNP package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _intel:

INTEL package
-----------------------------------

To build with this package, you must choose which hardware you want to
build for, either x86 CPUs or Intel KNLs in offload mode.  You should
also typically :ref:`install the OPENMP package <openmp>`, as it can be
used in tandem with the INTEL package to good effect, as explained
on the :doc:`Speed_intel` page.

When using Intel compilers version 16.0 or later is required.  You can
also use the GNU or Clang compilers and they will provide performance
improvements over regular styles and OPENMP styles, but less so than
with the Intel compilers.  Please also note, that some compilers have
been found to apply memory alignment constraints incompletely or
incorrectly and thus can cause segmentation faults in otherwise correct
code when using features from the INTEL package.


.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D INTEL_ARCH=value     # value = cpu (default) or knl
         -D INTEL_LRT_MODE=value # value = threads, none, or c++17

   .. tab:: Traditional make

      Choose which hardware to compile for in Makefile.machine via the
      following settings.  See ``src/MAKE/OPTIONS/Makefile.intel_cpu*``
      and ``Makefile.knl`` files for examples. and
      ``src/INTEL/README`` for additional information.

      For CPUs:

      .. code-block:: make

         OPTFLAGS =  -xHost -O2 -fp-model fast=2 -no-prec-div -qoverride-limits -qopt-zmm-usage=high
         CCFLAGS =   -g -qopenmp -DLAMMPS_MEMALIGN=64 -no-offload -fno-alias -ansi-alias -restrict $(OPTFLAGS)
         LINKFLAGS = -g -qopenmp $(OPTFLAGS)
         LIB =       -ltbbmalloc

      For KNLs:

      .. code-block:: make

         OPTFLAGS =  -xMIC-AVX512 -O2 -fp-model fast=2 -no-prec-div -qoverride-limits
         CCFLAGS =   -g -qopenmp -DLAMMPS_MEMALIGN=64 -no-offload -fno-alias -ansi-alias -restrict $(OPTFLAGS)
         LINKFLAGS = -g -qopenmp $(OPTFLAGS)
         LIB =       -ltbbmalloc

In Long-range thread mode (LRT) a modified verlet style is used, that
operates the Kspace calculation in a separate thread concurrently to
other calculations. This has to be enabled in the :doc:`package intel
<package>` command at runtime. With the setting "threads" it used the
pthreads library, while "c++17" will use the built-in thread support
of C++17 compilers. The option "none" skips compilation of this
feature. The default is to use "threads" if pthreads is available and
otherwise "none".

Best performance is achieved with Intel hardware, Intel compilers, as
well as the Intel TBB and MKL libraries. However, the code also
compiles, links, and runs with other compilers / hardware and without
TBB and MKL.

----------

.. _mdi:

MDI package
-----------

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_MDI=value    # download MDI Library for build
                                  # value = no (default) or yes

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The MDI package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _misc:

MISC package
------------

The :doc:`fix imd <fix_imd>` style in this package can be run either
synchronously (communication with IMD clients is done in the main
process) or asynchronously (the fix spawns a separate thread that can
communicate with IMD clients concurrently to the LAMMPS execution).

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D LAMMPS_ASYNC_IMD=value  # Run IMD server asynchronously
                                    # value = no (default) or yes

   .. tab:: Traditional make

      To enable asynchronous mode the ``-DLAMMPS_ASYNC_IMD`` define
      needs to be added to the ``LMP_INC`` variable in the
      ``Makefile.machine`` you are using.  For example:

      .. code-block:: make

         LMP_INC = -DLAMMPS_ASYNC_IMD -DLAMMPS_MEMALIGN=64

----------

.. _molfile:

MOLFILE package
---------------------------------------

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D MOLFILE_INCLUDE_DIR=path   # (optional) path where VMD molfile
                                       # plugin headers are installed
         -D PKG_MOLFILE=yes

      Using ``-D PKG_MOLFILE=yes`` enables the package, and setting
      ``-D MOLFILE_INCLUDE_DIR`` allows to provide a custom location for
      the molfile plugin header files. These should match the ABI of the
      plugin files used, and thus one typically sets them to include
      folder of the local VMD installation in use. LAMMPS ships with a
      couple of default header files that correspond to a popular VMD
      version, usually the latest release.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The MOLFILE package no longer supports the the traditional make
      build.  You need to build LAMMPS with CMake.

----------

.. _netcdf:

NETCDF package
-------------------------------------

To build with this package you must have the NetCDF library installed
on your system.

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_NETCDF=yes``.

      This should auto-detect the NETCDF library if it is installed on
      your system at standard locations.  Several advanced CMake NETCDF
      options exist if you need to specify where it was installed.  Use
      the ``ccmake`` (terminal window) or ``cmake-gui`` (graphical)
      tools to see these options and set them interactively from their
      user interfaces.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The NETCDF package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _openmp:

OPENMP package
-------------------------------

.. tabs::

   .. tab:: CMake build

      No additional settings are required besides ``-D
      PKG_OPENMP=yes``.  If CMake detects OpenMP compiler support, the
      OPENMP code will be compiled with multi-threading support
      enabled, otherwise as optimized serial code.

   .. tab:: Traditional make

      To enable multi-threading support in the OPENMP package (and
      other styles supporting OpenMP) the following compile and link
      flags must be added to your Makefile.machine file.  See
      ``src/MAKE/OPTIONS/Makefile.omp`` for an example.

      .. parsed-literal::

         CCFLAGS: -fopenmp               # for GNU and Clang Compilers
         CCFLAGS: -qopenmp -restrict     # for Intel compilers on Linux
         LINKFLAGS: -fopenmp             # for GNU and Clang Compilers
         LINKFLAGS: -qopenmp             # for Intel compilers on Linux

      For other platforms and compilers, please consult the
      documentation about OpenMP support for your compiler.

.. admonition:: Adding OpenMP support on macOS
   :class: note

   Apple offers the `Xcode package and IDE
   <https://developer.apple.com/xcode/>`_ for compiling software on
   macOS, so you have likely installed it to compile LAMMPS.  Their
   compiler is based on `Clang <https://clang.llvm.org/>`_, but while it
   is capable of processing OpenMP directives, the necessary header
   files and OpenMP runtime library are missing.  The `R developers
   <https://www.r-project.org/>`_ have figured out a way to build those
   in a compatible fashion. One can download them from
   `https://mac.r-project.org/openmp/
   <https://mac.r-project.org/openmp/>`_.  Simply adding those files as
   instructed enables the Xcode C++ compiler to compile LAMMPS with ``-D
   BUILD_OMP=yes``.

----------

.. _qmmm:

QMMM package
---------------------------------

For using LAMMPS to do QM/MM simulations via the QMMM package you
need to build LAMMPS as a library.  A LAMMPS executable with :doc:`fix
qmmm <fix_qmmm>` included can be built, but will not be able to do a
QM/MM simulation on as such.  You must also build a QM code - currently
only Quantum ESPRESSO (QE) is supported - and create a new executable
which links LAMMPS and the QM code together.  Details are given in the
``lib/qmmm/README`` file.  It is also recommended to read the
instructions for :doc:`linking with LAMMPS as a library <Build_link>`
for background information.  This requires compatible Quantum Espresso
and LAMMPS versions.  The current interface and makefiles have last been
verified to work in February 2020 with Quantum Espresso versions 6.3 to
6.5.

.. tabs::

   .. tab:: CMake build

      When using CMake, building a LAMMPS library is required and it is
      recommended to build a shared library, since any libraries built
      from the sources in the *lib* folder (including the essential
      libqmmm.a) are not included in the static LAMMPS library and are
      (currently) not installed, while their code is included in the
      shared LAMMPS library.  Thus a typical command to configure
      building LAMMPS for QMMM would be:

      .. code-block:: bash

         cmake -C ../cmake/presets/basic.cmake -D PKG_QMMM=yes \
             -D BUILD_LIB=yes -DBUILD_SHARED_LIBS=yes ../cmake

      After completing the LAMMPS build and also configuring and
      compiling Quantum ESPRESSO with external library support (via
      "make couple"), go back to the ``lib/qmmm`` folder and follow the
      instructions on the README file to build the combined LAMMPS/QE
      QM/MM executable (pwqmmm.x) in the ``lib/qmmm`` folder.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The QMMM package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _rheo:

RHEO package
------------

This package depends on the BPM package.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D PKG_RHEO=yes               # enable the package itself
         -D PKG_BPM=yes                # the RHEO package requires BPM
         -D USE_INTERNAL_LINALG=value  # prefer internal LAPACK if true

      Some features in the RHEO package are dependent on code in the BPM
      package so the latter one *must* be enabled as well.

      The RHEO package also requires LAPACK (and BLAS) and CMake
      can identify their locations and pass that info to the RHEO
      build script.  But on some systems this may cause problems when
      linking or the dependency is not desired.  By using the setting
      ``-D USE_INTERNAL_LINALG=yes`` when running the CMake
      configuration, you will select compiling and linking the bundled
      linear algebra library and work around the limitations.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The RHEO package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _scafacos:

SCAFACOS package
-----------------------------------------

To build with this package, you must download and build the
`ScaFaCoS Coulomb solver library <http://www.scafacos.de/>`_

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_SCAFACOS=value    # download ScaFaCoS for build, value = no (default) or yes
         -D SCAFACOS_LIBRARY=path      # ScaFaCos library file (only needed if at custom location)
         -D SCAFACOS_INCLUDE_DIR=path  # ScaFaCoS include directory (only needed if at custom location)

      If ``DOWNLOAD_SCAFACOS`` is set, the ScaFaCoS library will be
      downloaded and built inside the CMake build directory.  If the
      ScaFaCoS library is already on your system (in a location CMake
      cannot find it), ``SCAFACOS_LIBRARY`` is the filename (plus path) of
      the ScaFaCoS library file, not the directory the library file is
      in.  ``SCAFACOS_INCLUDE_DIR`` is the directory the ScaFaCoS include
      file is in.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The SCAFACOS package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.

----------

.. _vtk:

VTK package
-------------------------------

To build with this package you must have the VTK library installed on
your system.

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_VTK=yes``.

      This should auto-detect the VTK library if it is installed on your
      system at standard locations.  Several advanced VTK options exist
      if you need to specify where it was installed.  Use the ``ccmake``
      (terminal window) or ``cmake-gui`` (graphical) tools to see these
      options and set them interactively from their user interfaces.

   .. tab:: Traditional make

      .. versionchanged:: TBD

      The VTK package no longer supports the the traditional make build.
      You need to build LAMMPS with CMake.
