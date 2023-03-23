Packages with extra build options
=================================

When building with some packages, additional steps may be required,
in addition to

.. list-table::
   :align: center
   :header-rows: 1

   * - CMake build
     - Traditional make
   * - .. code-block:: bash

          cmake -D PKG_NAME=yes

     - .. code-block:: console

          make yes-name

as described on the :doc:`Build_package <Build_package>` page.

For a CMake build there may be additional optional or required
variables to set.  For a build with make, a provided library under the
lammps/lib directory may need to be built first.  Or an external
library may need to exist on your system or be downloaded and built.
You may need to tell LAMMPS where it is found on your system.

This is the list of packages that may require additional steps.

.. this list must be kept in sync with its counterpart in Build_package.rst
.. table_from_list::
   :columns: 6

   * :ref:`ADIOS <adios>`
   * :ref:`ATC <atc>`
   * :ref:`AWPMD <awpmd>`
   * :ref:`COLVARS <colvar>`
   * :ref:`COMPRESS <compress>`
   * :ref:`ELECTRODE <electrode>`
   * :ref:`GPU <gpu>`
   * :ref:`H5MD <h5md>`
   * :ref:`INTEL <intel>`
   * :ref:`KIM <kim>`
   * :ref:`KOKKOS <kokkos>`
   * :ref:`LATTE <latte>`
   * :ref:`LEPTON <lepton>`
   * :ref:`MACHDYN <machdyn>`
   * :ref:`MDI <mdi>`
   * :ref:`ML-HDNNP <ml-hdnnp>`
   * :ref:`ML-IAP <mliap>`
   * :ref:`ML-PACE <ml-pace>`
   * :ref:`ML-POD <ml-pod>`
   * :ref:`ML-QUIP <ml-quip>`
   * :ref:`MOLFILE <molfile>`
   * :ref:`MSCG <mscg>`
   * :ref:`NETCDF <netcdf>`
   * :ref:`OPENMP <openmp>`
   * :ref:`OPT <opt>`
   * :ref:`PLUMED <plumed>`
   * :ref:`POEMS <poems>`
   * :ref:`PYTHON <python>`
   * :ref:`QMMM <qmmm>`
   * :ref:`SCAFACOS <scafacos>`
   * :ref:`VORONOI <voronoi>`
   * :ref:`VTK <vtk>`

----------

.. _compress:

COMPRESS package
----------------

To build with this package you must have the `zlib compression library
<https://zlib.net>`_ available on your system to build dump styles with
a '/gz' suffix.  There are also styles using the
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
      pkg-config can find it, which may require adding it to the
      ``PKG_CONFIG_PATH`` environment variable.

   .. tab:: Traditional make

      To include support for Zstandard compression, ``-DLAMMPS_ZSTD``
      must be added to the compiler flags.  If make cannot find the
      libraries, you can edit the file ``lib/compress/Makefile.lammps``
      to specify the paths and library names.  This must be done
      **before** the package is installed.

----------

.. _gpu:

GPU package
---------------------

To build with this package, you must choose options for precision and
which GPU hardware to build for. The GPU package currently supports
three different types of backends: OpenCL, CUDA and HIP.

CMake build
^^^^^^^^^^^

.. code-block:: bash

   -D GPU_API=value             # value = opencl (default) or cuda or hip
   -D GPU_PREC=value            # precision setting
                                # value = double or mixed (default) or single
   -D GPU_ARCH=value            # primary GPU hardware choice for GPU_API=cuda
                                # value = sm_XX (see below, default is sm_50)
   -D GPU_DEBUG=value           # enable debug code in the GPU package library, mostly useful for developers
                                # value = yes or no (default)
   -D HIP_PATH=value            # value = path to HIP installation. Must be set if GPU_API=HIP
   -D HIP_ARCH=value            # primary GPU hardware choice for GPU_API=hip
                                # value depends on selected HIP_PLATFORM
                                # default is 'gfx906' for HIP_PLATFORM=amd and 'sm_50' for HIP_PLATFORM=nvcc
   -D HIP_USE_DEVICE_SORT=value # enables GPU sorting
                                # value = yes (default) or no
   -D CUDPP_OPT=value           # use GPU binning on with CUDA (should be off for modern GPUs)
                                # enables CUDA Performance Primitives, must be "no" for CUDA_MPS_SUPPORT=yes
                                # value = yes or no (default)
   -D CUDA_MPS_SUPPORT=value    # enables some tweaks required to run with active nvidia-cuda-mps daemon
                                # value = yes or no (default)
   -D USE_STATIC_OPENCL_LOADER=value  # downloads/includes OpenCL ICD loader library, no local OpenCL headers/libs needed
                                      # value = yes (default) or no

:code:`GPU_ARCH` settings for different GPU hardware is as follows:

* sm_30 for Kepler (supported since CUDA 5 and until CUDA 10.x)
* sm_35 or sm_37 for Kepler (supported since CUDA 5 and until CUDA 11.x)
* sm_50 or sm_52 for Maxwell (supported since CUDA 6)
* sm_60 or sm_61 for Pascal (supported since CUDA 8)
* sm_70 for Volta (supported since CUDA 9)
* sm_75 for Turing (supported since CUDA 10)
* sm_80 or sm_86 for Ampere (supported since CUDA 11, sm_86 since CUDA 11.1)
* sm_89 for Lovelace (supported since CUDA 11.8)
* sm_90 for Hopper (supported since CUDA 12.0)

A more detailed list can be found, for example,
at `Wikipedia's CUDA article <https://en.wikipedia.org/wiki/CUDA#GPUs_supported>`_

CMake can detect which version of the CUDA toolkit is used and thus will try
to include support for **all** major GPU architectures supported by this toolkit.
Thus the GPU_ARCH setting is merely an optimization, to have code for
the preferred GPU architecture directly included rather than having to wait
for the JIT compiler of the CUDA driver to translate it.

When compiling for CUDA or HIP with CUDA, version 8.0 or later of the CUDA toolkit
is required and a GPU architecture of Kepler or later, which must *also* be
supported by the CUDA toolkit in use **and** the CUDA driver in use.
When compiling for OpenCL, OpenCL version 1.2 or later is required and the
GPU must be supported by the GPU driver and OpenCL runtime bundled with the driver.

When building with CMake, you **must NOT** build the GPU library in ``lib/gpu``
using the traditional build procedure. CMake will detect files generated by that
process and will terminate with an error and a suggestion for how to remove them.

If you are compiling for OpenCL, the default setting is to download, build, and
link with a static OpenCL ICD loader library and standard OpenCL headers.  This
way no local OpenCL development headers or library needs to be present and only
OpenCL compatible drivers need to be installed to use OpenCL.  If this is not
desired, you can set :code:`USE_STATIC_OPENCL_LOADER` to :code:`no`.

The GPU library has some multi-thread support using OpenMP.  If LAMMPS is built
with ``-D BUILD_OMP=on`` this will also be enabled.

If you are compiling with HIP, note that before running CMake you will have to
set appropriate environment variables. Some variables such as
:code:`HCC_AMDGPU_TARGET` (for ROCm <= 4.0) or :code:`CUDA_PATH` are necessary for :code:`hipcc`
and the linker to work correctly.

Using CHIP-SPV implementation of HIP is now supported. It allows one to run HIP
code on Intel GPUs via the OpenCL or Level Zero backends. To use CHIP-SPV, you must
set :code:`-DHIP_USE_DEVICE_SORT=OFF` in your CMake command line as CHIP-SPV does not
yet support hipCUB. The use of HIP for Intel GPUs is still experimental so you
should only use this option in preparations to run on Aurora system at ANL.

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

   # CUDA target (not recommended, use GPU_ARCH=cuda)
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

Traditional make
^^^^^^^^^^^^^^^^

Before building LAMMPS, you must build the GPU library in ``lib/gpu``\ .
You can do this manually if you prefer; follow the instructions in
``lib/gpu/README``.  Note that the GPU library uses MPI calls, so you must
use the same MPI library (or the STUBS library) settings as the main
LAMMPS code.  This also applies to the ``-DLAMMPS_BIGBIG``\ ,
``-DLAMMPS_SMALLBIG``\ , or ``-DLAMMPS_SMALLSMALL`` settings in whichever
Makefile you use.

You can also build the library in one step from the ``lammps/src`` dir,
using a command like these, which simply invokes the ``lib/gpu/Install.py``
script with the specified args:

.. code-block:: bash

  make lib-gpu               # print help message
  make lib-gpu args="-b"     # build GPU library with default Makefile.linux
  make lib-gpu args="-m xk7 -p single -o xk7.single"  # create new Makefile.xk7.single, altered for single-precision
  make lib-gpu args="-m mpi -a sm_60 -p mixed -b" # build GPU library with mixed precision and P100 using other settings in Makefile.mpi

Note that this procedure starts with a Makefile.machine in lib/gpu, as
specified by the "-m" switch.  For your convenience, machine makefiles
for "mpi" and "serial" are provided, which have the same settings as
the corresponding machine makefiles in the main LAMMPS source
folder. In addition you can alter 4 important settings in the
Makefile.machine you start from via the corresponding -c, -a, -p, -e
switches (as in the examples above), and also save a copy of the new
Makefile if desired:

* ``CUDA_HOME`` = where NVIDIA CUDA software is installed on your system
* ``CUDA_ARCH`` = sm_XX, what GPU hardware you have, same as CMake GPU_ARCH above
* ``CUDA_PRECISION`` = precision (double, mixed, single)
* ``EXTRAMAKE`` = which Makefile.lammps.\* file to copy to Makefile.lammps

The file Makefile.cuda is set up to include support for multiple
GPU architectures as supported by the CUDA toolkit in use. This is done
through using the "--gencode " flag, which can be used multiple times and
thus support all GPU architectures supported by your CUDA compiler.

To enable GPU binning via CUDA performance primitives set the Makefile variable
``CUDPP_OPT = -DUSE_CUDPP -Icudpp_mini``.  This should **not** be used with
most modern GPUs.

To support the CUDA multiprocessor server you can set the define
``-DCUDA_MPS_SUPPORT``.  Please note that in this case you must **not** use
the CUDA performance primitives and thus set the variable ``CUDPP_OPT``
to empty.

The GPU library has some multi-thread support using OpenMP.  You need to add
the compiler flag that enables OpenMP to the ``CUDR_OPTS`` Makefile variable.

If the library build is successful, 3 files should be created:
``lib/gpu/libgpu.a``\ , ``lib/gpu/nvc_get_devices``\ , and
``lib/gpu/Makefile.lammps``\ .  The latter has settings that enable LAMMPS
to link with CUDA libraries.  If the settings in ``Makefile.lammps`` for
your machine are not correct, the LAMMPS build will fail, and
``lib/gpu/Makefile.lammps`` may need to be edited.

.. note::

   If you re-build the GPU library in ``lib/gpu``, you should always
   uninstall the GPU package in ``lammps/src``, then re-install it and
   re-build LAMMPS.  This is because the compilation of files in the GPU
   package uses the library settings from the ``lib/gpu/Makefile.machine``
   used to build the GPU library.

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

         -D DOWNLOAD_KIM=value           # download OpenKIM API v2 for build, value = no (default) or yes
         -D LMP_DEBUG_CURL=value         # set libcurl verbose mode on/off, value = off (default) or on
         -D LMP_NO_SSL_CHECK=value       # tell libcurl to not verify the peer, value = no (default) or yes
         -D KIM_EXTRA_UNITTESTS=value    # enables extra unit tests, value = no (default) or yes

      If ``DOWNLOAD_KIM`` is set to *yes* (or *on*), the KIM API library
      will be downloaded and built inside the CMake build directory.  If
      the KIM library is already installed on your system (in a location
      where CMake cannot find it), you may need to set the
      ``PKG_CONFIG_PATH`` environment variable so that libkim-api can be
      found, or run the command ``source kim-api-activate``.

      Extra unit tests can only be available if they are explicitly requested
      (``KIM_EXTRA_UNITTESTS`` is set to *yes* (or *on*)) and the prerequisites
      are met. See :ref:`KIM Extra unit tests <kim_extra_unittests>` for
      more details on this.

   .. tab:: Traditional make

      You can download and build the KIM library manually if you prefer;
      follow the instructions in ``lib/kim/README``.  You can also do
      this in one step from the lammps/src directory, using a command like
      these, which simply invokes the ``lib/kim/Install.py`` script with
      the specified args.

      .. code-block:: bash

         make lib-kim              # print help message
         make lib-kim args="-b "   # (re-)install KIM API lib with only example models
         make lib-kim args="-b -a Glue_Ercolessi_Adams_Al__MO_324507536345_001"  # ditto plus one model
         make lib-kim args="-b -a everything"     # install KIM API lib with all models
         make lib-kim args="-n -a EAM_Dynamo_Ackland_W__MO_141627196590_002"       # add one model or model driver
         make lib-kim args="-p /usr/local" # use an existing KIM API installation at the provided location
         make lib-kim args="-p /usr/local -a EAM_Dynamo_Ackland_W__MO_141627196590_002" # ditto but add one model or driver

      When using the "-b " option, the KIM library is built using its native
      cmake build system.  The ``lib/kim/Install.py`` script supports a
      ``CMAKE`` environment variable if the cmake executable is named other
      than ``cmake`` on your system.  Additional environment variables may be
      provided on the command line for use by cmake.  For example, to use the
      ``cmake3`` executable and tell it to use the gnu version 11 compilers
      to build KIM, one could use the following command line.

      .. code-block:: bash

         CMAKE=cmake3 CXX=g++-11 CC=gcc-11 FC=gfortran-11 make lib-kim args="-b "  # (re-)install KIM API lib using cmake3 and gnu v11 compilers with only example models

      Settings for debugging OpenKIM web queries discussed below need to
      be applied by adding them to the ``LMP_INC`` variable through
      editing the ``Makefile.machine`` you are using.  For example:

      .. code-block:: make

         LMP_INC = -DLMP_NO_SSL_CHECK

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
``KIM_EXTRA_UNITTESTS`` to *yes* (or *on*).

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
* It is also necessary to install
  ``EAM_Dynamo_MendelevAckland_2007v3_Zr__MO_004835508849_000``,
  ``EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005``, and
  ``LennardJones612_UniversalShifted__MO_959249795837_003`` KIM models.
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
   error message indicating a mismatch, if that happens.

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
      - AMD 64-bit x86 CPU (AVX 1)
   *  - ZEN
      - HOST
      - AMD Zen class CPU (AVX 2)
   *  - ZEN2
      - HOST
      - AMD Zen2 class CPU (AVX 2)
   *  - ZEN3
      - HOST
      - AMD Zen3 class CPU (AVX 2)
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
   *  - WSM
      - HOST
      - Intel Westmere CPU (SSE 4.2)
   *  - SNB
      - HOST
      - Intel Sandy/Ivy Bridge CPU (AVX 1)
   *  - HSW
      - HOST
      - Intel Haswell CPU (AVX 2)
   *  - BDW
      - HOST
      - Intel Broadwell Xeon E-class CPU (AVX 2 + transactional mem)
   *  - SKL
      - HOST
      - Intel Skylake Client CPU
   *  - SKX
      - HOST
      - Intel Skylake Xeon Server CPU (AVX512)
   *  - ICL
      - HOST
      - Intel Ice Lake Client CPU (AVX512)
   *  - ICX
      - HOST
      - Intel Ice Lake Xeon Server CPU (AVX512)
   *  - SPR
      - HOST
      - Intel Sapphire Rapids Xeon Server CPU (AVX512)
   *  - KNC
      - HOST
      - Intel Knights Corner Xeon Phi
   *  - KNL
      - HOST
      - Intel Knights Landing Xeon Phi
   *  - BGQ
      - HOST
      - IBM Blue Gene/Q CPU
   *  - POWER7
      - HOST
      - IBM POWER7 CPU
   *  - POWER8
      - HOST
      - IBM POWER8 CPU
   *  - POWER9
      - HOST
      - IBM POWER9 CPU
   *  - KEPLER30
      - GPU
      - NVIDIA Kepler generation CC 3.0 GPU
   *  - KEPLER32
      - GPU
      - NVIDIA Kepler generation CC 3.2 GPU
   *  - KEPLER35
      - GPU
      - NVIDIA Kepler generation CC 3.5 GPU
   *  - KEPLER37
      - GPU
      - NVIDIA Kepler generation CC 3.7 GPU
   *  - MAXWELL50
      - GPU
      - NVIDIA Maxwell generation CC 5.0 GPU
   *  - MAXWELL52
      - GPU
      - NVIDIA Maxwell generation CC 5.2 GPU
   *  - MAXWELL53
      - GPU
      - NVIDIA Maxwell generation CC 5.3 GPU
   *  - PASCAL60
      - GPU
      - NVIDIA Pascal generation CC 6.0 GPU
   *  - PASCAL61
      - GPU
      - NVIDIA Pascal generation CC 6.1 GPU
   *  - VOLTA70
      - GPU
      - NVIDIA Volta generation CC 7.0 GPU
   *  - VOLTA72
      - GPU
      - NVIDIA Volta generation CC 7.2 GPU
   *  - TURING75
      - GPU
      - NVIDIA Turing generation CC 7.5 GPU
   *  - AMPERE80
      - GPU
      - NVIDIA Ampere generation CC 8.0 GPU
   *  - AMPERE86
      - GPU
      - NVIDIA Ampere generation CC 8.6 GPU
   *  - ADA89
      - GPU
      - NVIDIA Ada Lovelace generation CC 8.9 GPU
   *  - HOPPER90
      - GPU
      - NVIDIA Hopper generation CC 9.0 GPU
   *  - VEGA900
      - GPU
      - AMD GPU MI25 GFX900
   *  - VEGA906
      - GPU
      - AMD GPU MI50/MI60 GFX906
   *  - VEGA908
      - GPU
      - AMD GPU MI100 GFX908
   *  - VEGA90A
      - GPU
      - AMD GPU MI200 GFX90A
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

This list was last updated for version 3.7.1 of the Kokkos library.

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
         -D CMAKE_CXX_COMPILER=wrapper # wrapper = full path to Cuda nvcc wrapper

      This will also enable executing FFTs on the GPU, either via the
      internal KISSFFT library, or - by preference - with the cuFFT
      library bundled with the CUDA toolkit, depending on whether CMake
      can identify its location.  The *wrapper* value for
      ``CMAKE_CXX_COMPILER`` variable is the path to the CUDA nvcc
      compiler wrapper provided in the Kokkos library:
      ``lib/kokkos/bin/nvcc_wrapper``\ .  The setting should include the
      full path name to the wrapper, e.g.

      .. code-block:: bash

         -D CMAKE_CXX_COMPILER=${HOME}/lammps/lib/kokkos/bin/nvcc_wrapper

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

      To simplify compilation, five preset files are included in the
      ``cmake/presets`` folder, ``kokkos-serial.cmake``,
      ``kokkos-openmp.cmake``, ``kokkos-cuda.cmake``,
      ``kokkos-hip.cmake``, and ``kokkos-sycl.cmake``.  They will enable
      the KOKKOS package and enable some hardware choice.  So to compile
      with CUDA device parallelization (for GPUs with CC 5.0 and up)
      with some common packages enabled, you can do the following:

      .. code-block:: bash

         mkdir build-kokkos-cuda
         cd build-kokkos-cuda
         cmake -C ../cmake/presets/basic.cmake -C ../cmake/presets/kokkos-cuda.cmake ../cmake
         cmake --build .

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
         KOKKOS_ARCH = HOSTARCH,GPUARCH  # HOSTARCH = HOST from list above that is hosting the GPU
         KOKKOS_CUDA_OPTIONS = "enable_lambda"
                                         # GPUARCH = GPU from list above
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
         KOKKOS_ARCH = HOSTARCH,GPUARCH  # HOSTARCH = HOST from list above that is hosting the GPU
                                         # GPUARCH = GPU from list above
         FFT_INC = -DFFT_HIPFFT           # enable use of hipFFT (optional)
         FFT_LIB = -lhipfft               # link to hipFFT library

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

Please note, that the LAMMPS KOKKOS package must **always** be compiled
with the *enable_lambda* option when using GPUs.  The CMake configuration
will thus always enable it.

----------

.. _latte:

LATTE package
-------------------------

To build with this package, you must download and build the LATTE
library.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_LATTE=value      # download LATTE for build, value = no (default) or yes
         -D LATTE_LIBRARY=path        # LATTE library file (only needed if a custom location)
         -D USE_INTERNAL_LINALG=value # Use the internal linear algebra library instead of LAPACK
                                      #   value = no (default) or yes

      If ``DOWNLOAD_LATTE`` is set, the LATTE library will be downloaded
      and built inside the CMake build directory.  If the LATTE library
      is already on your system (in a location CMake cannot find it),
      ``LATTE_LIBRARY`` is the filename (plus path) of the LATTE library
      file, not the directory the library file is in.

      The LATTE library requires LAPACK (and BLAS) and CMake can identify
      their locations and pass that info to the LATTE build script. But
      on some systems this triggers a (current) limitation of CMake and
      the configuration will fail. Try enabling ``USE_INTERNAL_LINALG`` in
      those cases to use the bundled linear algebra library and work around
      the limitation.

   .. tab:: Traditional make

      You can download and build the LATTE library manually if you
      prefer; follow the instructions in ``lib/latte/README``\ .  You
      can also do it in one step from the ``lammps/src`` dir, using a
      command like these, which simply invokes the
      ``lib/latte/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-latte                        # print help message
         make lib-latte args="-b"              # download and build in lib/latte/LATTE-master
         make lib-latte args="-p $HOME/latte"  # use existing LATTE installation in $HOME/latte
         make lib-latte args="-b -m gfortran"  # download and build in lib/latte and
                                               #   copy Makefile.lammps.gfortran to Makefile.lammps

      Note that 3 symbolic (soft) links, ``includelink`` and ``liblink``
      and ``filelink.o``, are created in ``lib/latte`` to point to
      required folders and files in the LATTE home directory.  When
      LAMMPS itself is built it will use these links.  You should also
      check that the ``Makefile.lammps`` file you create is appropriate
      for the compiler you use on your system to build LATTE.

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

      Before building LAMMPS, one must build the Lepton library in lib/lepton.

      This can be done manually in the same folder by using or adapting
      one of the provided Makefiles: for example, ``Makefile.serial`` for
      the GNU C++ compiler, or ``Makefile.mpi`` for the MPI compiler wrapper.
      The Lepton library is written in C++-11 and thus the C++ compiler
      may need to be instructed to enable support for that.

      In general, it is safer to use build setting consistent with the
      rest of LAMMPS.  This is best carried out from the LAMMPS src
      directory using a command like these, which simply invokes the
      ``lib/lepton/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-lepton                      # print help message
         make lib-lepton args="-m serial"     # build with GNU g++ compiler (settings as with "make serial")
         make lib-lepton args="-m mpi"        # build with default MPI compiler (settings as with "make mpi")

      The "machine" argument of the "-m" flag is used to find a
      Makefile.machine to use as build recipe.

      The build should produce a ``build`` folder and the library ``lib/lepton/liblmplepton.a``

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

.. _mscg:

MSCG package
-----------------------

To build with this package, you must download and build the MS-CG
library.  Building the MS-CG library requires that the GSL
(GNU Scientific Library) headers and libraries are installed on your
machine.  See the ``lib/mscg/README`` and ``MSCG/Install`` files for
more details.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_MSCG=value    # download MSCG for build, value = no (default) or yes
         -D MSCG_LIBRARY=path      # MSCG library file (only needed if a custom location)
         -D MSCG_INCLUDE_DIR=path  # MSCG include directory (only needed if a custom location)

      If ``DOWNLOAD_MSCG`` is set, the MSCG library will be downloaded
      and built inside the CMake build directory.  If the MSCG library
      is already on your system (in a location CMake cannot find it),
      ``MSCG_LIBRARY`` is the filename (plus path) of the MSCG library
      file, not the directory the library file is in.
      ``MSCG_INCLUDE_DIR`` is the directory the MSCG include file is in.

   .. tab:: Traditional make

      You can download and build the MS-CG library manually if you
      prefer; follow the instructions in ``lib/mscg/README``\ .  You can
      also do it in one step from the ``lammps/src`` dir, using a
      command like these, which simply invokes the
      ``lib/mscg/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-mscg             # print help message
         make lib-mscg args="-b -m serial"   # download and build in lib/mscg/MSCG-release-master
                                             # with the settings compatible with "make serial"
         make lib-mscg args="-b -m mpi"      # download and build in lib/mscg/MSCG-release-master
                                             # with the settings compatible with "make mpi"
         make lib-mscg args="-p /usr/local/mscg-release" # use the existing MS-CG installation in /usr/local/mscg-release

      Note that 2 symbolic (soft) links, ``includelink`` and ``liblink``,
      will be created in ``lib/mscg`` to point to the MS-CG
      ``src/installation`` dir.  When LAMMPS is built in src it will use
      these links.  You should not need to edit the
      ``lib/mscg/Makefile.lammps`` file.

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
      the :code:`CCFLAGS` line of your ``Makefile.machine``.  See
      ``src/MAKE/OPTIONS/Makefile.opt`` for an example.

----------

.. _poems:

POEMS package
-------------------------

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_OPT=yes``

   .. tab:: Traditional make

      Before building LAMMPS, you must build the POEMS library in
      ``lib/poems``\ .  You can do this manually if you prefer; follow
      the instructions in ``lib/poems/README``\ .  You can also do it in
      one step from the ``lammps/src`` dir, using a command like these,
      which simply invokes the ``lib/poems/Install.py`` script with the
      specified args:

      .. code-block:: bash

         make lib-poems                   # print help message
         make lib-poems args="-m serial"  # build with GNU g++ compiler (settings as with "make serial")
         make lib-poems args="-m mpi"     # build with default MPI C++ compiler (settings as with "make mpi")
         make lib-poems args="-m icc"     # build with Intel icc compiler

      The build should produce two files: ``lib/poems/libpoems.a`` and
      ``lib/poems/Makefile.lammps``.  The latter is copied from an
      existing ``Makefile.lammps.*`` and has settings needed to build
      LAMMPS with the POEMS library (though typically the settings are
      just blank).  If necessary, you can edit/create a new
      ``lib/poems/Makefile.machine`` file for your system, which should
      define an ``EXTRAMAKE`` variable to specify a corresponding
      ``Makefile.lammps.machine`` file.

----------

.. _python:

PYTHON package
---------------------------

Building with the PYTHON package requires you have a the Python development
headers and library available on your system, which needs to be a Python 2.7
version or a Python 3.x version.  Since support for Python 2.x has ended,
using Python 3.x is strongly recommended. See ``lib/python/README`` for
additional details.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D PYTHON_EXECUTABLE=path   # path to Python executable to use

      Without this setting, CMake will guess the default Python version
      on your system.  To use a different Python version, you can either
      create a virtualenv, activate it and then run cmake.  Or you can
      set the PYTHON_EXECUTABLE variable to specify which Python
      interpreter should be used.  Note note that you will also need to
      have the development headers installed for this version,
      e.g. python2-devel.

   .. tab:: Traditional make

      The build uses the ``lib/python/Makefile.lammps`` file in the
      compile/link process to find Python.  You should only need to
      create a new ``Makefile.lammps.*`` file (and copy it to
      ``Makefile.lammps``) if the LAMMPS build fails.

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

         -D DOWNLOAD_VORO=value    # download Voro++ for build, value = no (default) or yes
         -D VORO_LIBRARY=path      # Voro++ library file (only needed if at custom location)
         -D VORO_INCLUDE_DIR=path  # Voro++ include directory (only needed if at custom location)

      If ``DOWNLOAD_VORO`` is set, the Voro++ library will be downloaded
      and built inside the CMake build directory.  If the Voro++ library
      is already on your system (in a location CMake cannot find it),
      ``VORO_LIBRARY`` is the filename (plus path) of the Voro++ library
      file, not the directory the library file is in.
      ``VORO_INCLUDE_DIR`` is the directory the Voro++ include file is
      in.

   .. tab:: Traditional make

      You can download and build the Voro++ library manually if you
      prefer; follow the instructions in ``lib/voronoi/README``.  You
      can also do it in one step from the ``lammps/src`` dir, using a
      command like these, which simply invokes the
      ``lib/voronoi/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-voronoi                          # print help message
         make lib-voronoi args="-b"                # download and build the default version in lib/voronoi/voro++-<version>
         make lib-voronoi args="-p $HOME/voro++"   # use existing Voro++ installation in $HOME/voro++
         make lib-voronoi args="-b -v voro++0.4.6" # download and build the 0.4.6 version in lib/voronoi/voro++-0.4.6

      Note that 2 symbolic (soft) links, ``includelink`` and
      ``liblink``, are created in lib/voronoi to point to the Voro++
      source dir.  When LAMMPS builds in ``src`` it will use these
      links.  You should not need to edit the
      ``lib/voronoi/Makefile.lammps`` file.

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

.. _atc:

ATC package
-------------------------------

The ATC package requires the MANYBODY package also be installed.

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_ATC=yes``
      and ``-D PKG_MANYBODY=yes``.

   .. tab:: Traditional make

      Before building LAMMPS, you must build the ATC library in
      ``lib/atc``.  You can do this manually if you prefer; follow the
      instructions in ``lib/atc/README``.  You can also do it in one
      step from the ``lammps/src`` dir, using a command like these,
      which simply invokes the ``lib/atc/Install.py`` script with the
      specified args:

      .. code-block:: bash

         make lib-atc                      # print help message
         make lib-atc args="-m serial"     # build with GNU g++ compiler and MPI STUBS (settings as with "make serial")
         make lib-atc args="-m mpi"        # build with default MPI compiler (settings as with "make mpi")
         make lib-atc args="-m icc"        # build with Intel icc compiler

      The build should produce two files: ``lib/atc/libatc.a`` and
      ``lib/atc/Makefile.lammps``.  The latter is copied from an
      existing ``Makefile.lammps.*`` and has settings needed to build
      LAMMPS with the ATC library.  If necessary, you can edit/create a
      new ``lib/atc/Makefile.machine`` file for your system, which
      should define an ``EXTRAMAKE`` variable to specify a corresponding
      ``Makefile.lammps.<machine>`` file.

      Note that the Makefile.lammps file has settings for the BLAS and
      LAPACK linear algebra libraries.  As explained in
      ``lib/atc/README`` these can either exist on your system, or you
      can use the files provided in ``lib/linalg``.  In the latter case
      you also need to build the library in ``lib/linalg`` with a
      command like these:

      .. code-block:: bash

         make lib-linalg                   # print help message
         make lib-linalg args="-m serial"  # build with GNU C++ compiler (settings as with "make serial")
         make lib-linalg args="-m mpi"     # build with default MPI C++ compiler (settings as with "make mpi")
         make lib-linalg args="-m g++"     # build with GNU Fortran compiler

----------

.. _awpmd:

AWPMD package
-------------

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_AQPMD=yes``.

   .. tab:: Traditional make

      Before building LAMMPS, you must build the AWPMD library in
      ``lib/awpmd``.  You can do this manually if you prefer; follow the
      instructions in ``lib/awpmd/README``.  You can also do it in one
      step from the ``lammps/src`` dir, using a command like these,
      which simply invokes the ``lib/awpmd/Install.py`` script with the
      specified args:

      .. code-block:: bash

         make lib-awpmd                   # print help message
         make lib-awpmd args="-m serial"  # build with GNU g++ compiler and MPI STUBS (settings as with "make serial")
         make lib-awpmd args="-m mpi"     # build with default MPI compiler (settings as with "make mpi")
         make lib-awpmd args="-m icc"     # build with Intel icc compiler

      The build should produce two files: ``lib/awpmd/libawpmd.a`` and
      ``lib/awpmd/Makefile.lammps``.  The latter is copied from an
      existing ``Makefile.lammps.*`` and has settings needed to build
      LAMMPS with the AWPMD library.  If necessary, you can edit/create
      a new ``lib/awpmd/Makefile.machine`` file for your system, which
      should define an ``EXTRAMAKE`` variable to specify a corresponding
      ``Makefile.lammps.<machine>`` file.

      Note that the ``Makefile.lammps`` file has settings for the BLAS
      and LAPACK linear algebra libraries.  As explained in
      ``lib/awpmd/README`` these can either exist on your system, or you
      can use the files provided in ``lib/linalg``.  In the latter case
      you also need to build the library in ``lib/linalg`` with a
      command like these:

      .. code-block:: bash

         make lib-linalg                   # print help message
         make lib-linalg args="-m serial"  # build with GNU C++ compiler (settings as with "make serial")
         make lib-linalg args="-m mpi"     # build with default MPI C++ compiler (settings as with "make mpi")
         make lib-linalg args="-m g++"     # build with GNU C++ compiler

----------

.. _colvar:

COLVARS package
---------------

This package enables the use of the `Colvars <https://colvars.github.io/>`_
module included in the LAMMPS source distribution.


.. tabs::

   .. tab:: CMake build

      This is the recommended build procedure for using Colvars in
      LAMMPS. No additional settings are normally needed besides
      ``-D PKG_COLVARS=yes``.

   .. tab:: Traditional make

      As with other libraries distributed with LAMMPS, the Colvars library
      needs to be built before building the LAMMPS program with the COLVARS
      package enabled.

      From the LAMMPS ``src`` directory, this is most easily and safely done
      via one of the following commands, which implicitly rely on the
      ``lib/colvars/Install.py`` script with optional arguments:

      .. code-block:: bash

         make lib-colvars                      # print help message
         make lib-colvars args="-m serial"     # build with GNU g++ compiler (settings as with "make serial")
         make lib-colvars args="-m mpi"        # build with default MPI compiler (settings as with "make mpi")
         make lib-colvars args="-m g++-debug"  # build with GNU g++ compiler and colvars debugging enabled

      The "machine" argument of the "-m" flag is used to find a
      ``Makefile.machine`` file to use as build recipe.  If such recipe does
      not already exist in ``lib/colvars``, suitable settings will be
      auto-generated consistent with those used in the core LAMMPS makefiles.


      .. versionchanged:: 8Feb2023

      Please note that Colvars uses the Lepton library, which is now
      included with the LEPTON package; if you use anything other than
      the ``make lib-colvars`` command, please make sure to :ref:`build
      Lepton beforehand <lepton>`.

      Optional flags may be specified as environment variables:

      .. code-block:: bash

         COLVARS_DEBUG=yes make lib-colvars args="-m machine"  # Build with debug code (much slower)
         COLVARS_LEPTON=no make lib-colvars args="-m machine"  # Build without Lepton (included otherwise)

      The build should produce two files: the library
      ``lib/colvars/libcolvars.a`` and the specification file
      ``lib/colvars/Makefile.lammps``.  The latter is auto-generated,
      and normally does not need to be edited.

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
      can identify their locations and pass that info to the LATTE build
      script.  But on some systems this may cause problems when linking
      or the dependency is not desired.  Try enabling
      ``USE_INTERNAL_LINALG`` in those cases to use the bundled linear
      algebra library and work around the limitation.

   .. tab:: Traditional make

      Before building LAMMPS, you must configure the ELECTRODE support
      libraries and settings in ``lib/electrode``.  You can do this
      manually, if you prefer, or do it in one step from the
      ``lammps/src`` dir, using a command like these, which simply
      invokes the ``lib/electrode/Install.py`` script with the specified
      args:

      .. code-block:: bash

         make lib-electrode                   # print help message
         make lib-electrode args="-m serial"  # build with GNU g++ compiler and MPI STUBS (settings as with "make serial")
         make lib-electrode args="-m mpi"     # build with default MPI compiler (settings as with "make mpi")


      Note that the ``Makefile.lammps`` file has settings for the BLAS
      and LAPACK linear algebra libraries.  These can either exist on
      your system, or you can use the files provided in ``lib/linalg``.
      In the latter case you also need to build the library in
      ``lib/linalg`` with a command like these:

      .. code-block:: bash

         make lib-linalg                   # print help message
         make lib-linalg args="-m serial"  # build with GNU C++ compiler (settings as with "make serial")
         make lib-linalg args="-m mpi"     # build with default MPI C++ compiler (settings as with "make mpi")
         make lib-linalg args="-m g++"     # build with GNU C++ compiler

      The package itself is activated with ``make yes-KSPACE`` and
      ``make yes-ELECTRODE``

----------

.. _ml-pace:

ML-PACE package
-----------------------------

This package requires a library that can be downloaded and built
in lib/pace or somewhere else, which must be done before building
LAMMPS with this package. The code for the library can be found
at: `https://github.com/ICAMS/lammps-user-pace/ <https://github.com/ICAMS/lammps-user-pace/>`_

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

      You can download and build the ML-PACE library
      in one step from the ``lammps/src`` dir, using these commands,
      which invoke the ``lib/pace/Install.py`` script.

      .. code-block:: bash

         make lib-pace                          # print help message
         make lib-pace args="-b"                # download and build the default version in lib/pace

      You should not need to edit the ``lib/pace/Makefile.lammps`` file.

----------

.. _ml-pod:

ML-POD package
-----------------------------

.. tabs::

   .. tab:: CMake build

      No additional settings are needed besides ``-D PKG_ML-POD=yes``.

   .. tab:: Traditional make

      Before building LAMMPS, you must configure the ML-POD support
      settings in ``lib/mlpod``.  You can do this manually, if you
      prefer, or do it in one step from the ``lammps/src`` dir, using a
      command like the following, which simply invoke the
      ``lib/mlpod/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-mlpod                   # print help message
         make lib-mlpod args="-m serial"  # build with GNU g++ compiler and MPI STUBS (settings as with "make serial")
         make lib-mlpod args="-m mpi"     # build with default MPI compiler (settings as with "make mpi")
         make lib-mlpod args="-m mpi -e linalg"   # same as above but use the bundled linalg lib

      Note that the ``Makefile.lammps`` file has settings to use the BLAS
      and LAPACK linear algebra libraries.  These can either exist on
      your system, or you can use the files provided in ``lib/linalg``.
      In the latter case you also need to build the library in
      ``lib/linalg`` with a command like these:

      .. code-block:: bash

         make lib-linalg                   # print help message
         make lib-linalg args="-m serial"  # build with GNU C++ compiler (settings as with "make serial")
         make lib-linalg args="-m mpi"     # build with default MPI C++ compiler (settings as with "make mpi")
         make lib-linalg args="-m g++"     # build with GNU C++ compiler

      The package itself is activated with ``make yes-ML-POD``.

----------

.. _plumed:

PLUMED package
-------------------------------------

.. _plumedinstall: https://plumed.github.io/doc-master/user-doc/html/_installation.html

Before building LAMMPS with this package, you must first build PLUMED.
PLUMED can be built as part of the LAMMPS build or installed separately
from LAMMPS using the generic `PLUMED installation instructions <plumedinstall_>`_.
The PLUMED package has been tested to work with Plumed versions
2.4.x, 2.5.x, and 2.6.x and will error out, when trying to run calculations
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

.. tabs::

   .. tab:: CMake build

      When the ``-D PKG_PLUMED=yes`` flag is included in the cmake
      command you must ensure that GSL is installed in locations that
      are specified in your environment.  There are then two additional
      variables that control the manner in which PLUMED is obtained and
      linked into LAMMPS.

      .. code-block:: bash

         -D DOWNLOAD_PLUMED=value   # download PLUMED for build, value = no (default) or yes
         -D PLUMED_MODE=value       # Linkage mode for PLUMED, value = static (default), shared, or runtime

      If DOWNLOAD_PLUMED is set to "yes", the PLUMED library will be
      downloaded (the version of PLUMED that will be downloaded is
      hard-coded to a vetted version of PLUMED, usually a recent stable
      release version) and built inside the CMake build directory.  If
      ``DOWNLOAD_PLUMED`` is set to "no" (the default), CMake will try
      to detect and link to an installed version of PLUMED.  For this to
      work, the PLUMED library has to be installed into a location where
      the ``pkg-config`` tool can find it or the PKG_CONFIG_PATH
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

      PLUMED needs to be installed before the PLUMED package is
      installed so that LAMMPS can find the right settings when
      compiling and linking the LAMMPS executable.  You can either
      download and build PLUMED inside the LAMMPS plumed library folder
      or use a previously installed PLUMED library and point LAMMPS to
      its location. You also have to choose the linkage mode: "static"
      (default), "shared" or "runtime".  For a discussion of PLUMED
      linkage modes, please see above.

      Download/compilation/configuration of the plumed library can be done
      from the src folder through the following make args:

      .. code-block:: bash

         make lib-plumed                         # print help message
         make lib-plumed args="-b"               # download and build PLUMED in lib/plumed/plumed2
         make lib-plumed args="-p $HOME/.local"  # use existing PLUMED installation in $HOME/.local
         make lib-plumed args="-p /usr/local -m shared"  # use existing PLUMED installation in
                                                           # /usr/local and use shared linkage mode

      Note that 2 symbolic (soft) links, ``includelink`` and ``liblink``
      are created in lib/plumed that point to the location of the PLUMED
      build to use. A new file ``lib/plumed/Makefile.lammps`` is also
      created with settings suitable for LAMMPS to compile and link
      PLUMED using the desired linkage mode. After this step is
      completed, you can install the PLUMED package and compile
      LAMMPS in the usual manner:

      .. code-block:: bash

         make yes-plumed
         make machine

      Once this compilation completes you should be able to run LAMMPS
      in the usual way.  For shared linkage mode, libplumed.so must be
      found by the LAMMPS executable, which on many operating systems
      means, you have to set the LD_LIBRARY_PATH environment variable
      accordingly.

      Support for the different linkage modes in LAMMPS varies for
      different operating systems, using the static linkage is expected
      to be the most portable, and thus set to be the default.

      If you want to change the linkage mode, you have to re-run "make
      lib-plumed" with the desired settings **and** do a re-install if
      the PLUMED package with "make yes-plumed" to update the
      required makefile settings with the changes in the lib/plumed
      folder.

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

      Before building LAMMPS, you must build the CH5MD library in
      ``lib/h5md``.  You can do this manually if you prefer; follow the
      instructions in ``lib/h5md/README``.  You can also do it in one
      step from the ``lammps/src`` dir, using a command like these,
      which simply invokes the ``lib/h5md/Install.py`` script with the
      specified args:

      .. code-block:: bash

         make lib-h5md                     # print help message
         make lib-h5md args="-m h5cc"      # build with h5cc compiler

      The build should produce two files: ``lib/h5md/libch5md.a`` and
      ``lib/h5md/Makefile.lammps``.  The latter is copied from an
      existing ``Makefile.lammps.*`` and has settings needed to build
      LAMMPS with the system HDF5 library.  If necessary, you can
      edit/create a new ``lib/h5md/Makefile.machine`` file for your
      system, which should define an EXTRAMAKE variable to specify a
      corresponding ``Makefile.lammps.<machine>`` file.

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

         -D DOWNLOAD_N2P2=value    # download n2p2 for build, value = no (default) or yes
         -D N2P2_DIR=path          # n2p2 base directory (only needed if a custom location)

      If ``DOWNLOAD_N2P2`` is set, the *n2p2* library will be downloaded and
      built inside the CMake build directory.  If the *n2p2* library is already
      on your system (in a location CMake cannot find it), set the ``N2P2_DIR``
      to path where *n2p2* is located. If *n2p2* is located directly in
      ``lib/hdnnp/n2p2`` it will be automatically found by CMake.

   .. tab:: Traditional make

      You can download and build the *n2p2* library manually if you prefer;
      follow the instructions in ``lib/hdnnp/README``\ . You can also do it in
      one step from the ``lammps/src`` dir, using a command like these, which
      simply invokes the ``lib/hdnnp/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-hdnnp             # print help message
         make lib-hdnnp args="-b"   # download and build in lib/hdnnp/n2p2-...
         make lib-hdnnp args="-b -v 2.1.4" # download and build specific version
         make lib-hdnnp args="-p /usr/local/n2p2" # use the existing n2p2 installation in /usr/local/n2p2

      Note that 3 symbolic (soft) links, ``includelink``, ``liblink`` and
      ``Makefile.lammps``, will be created in ``lib/hdnnp`` to point to
      ``n2p2/include``, ``n2p2/lib`` and ``n2p2/lib/Makefile.lammps-extra``,
      respectively. When LAMMPS is built in ``src`` it will use these links.

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
         -D INTEL_LRT_MODE=value # value = threads, none, or c++11

   .. tab:: Traditional make

      Choose which hardware to compile for in Makefile.machine via the
      following settings.  See ``src/MAKE/OPTIONS/Makefile.intel_cpu*``
      and ``Makefile.knl`` files for examples. and
      ``src/INTEL/README`` for additional information.

      For CPUs:

      .. code-block:: make

         OPTFLAGS =      -xHost -O2 -fp-model fast=2 -no-prec-div -qoverride-limits -qopt-zmm-usage=high
         CCFLAGS =       -g -qopenmp -DLAMMPS_MEMALIGN=64 -no-offload -fno-alias -ansi-alias -restrict $(OPTFLAGS)
         LINKFLAGS =     -g -qopenmp $(OPTFLAGS)
         LIB =           -ltbbmalloc

      For KNLs:

      .. code-block:: make

         OPTFLAGS =      -xMIC-AVX512 -O2 -fp-model fast=2 -no-prec-div -qoverride-limits
         CCFLAGS =       -g -qopenmp -DLAMMPS_MEMALIGN=64 -no-offload -fno-alias -ansi-alias -restrict $(OPTFLAGS)
         LINKFLAGS =     -g -qopenmp $(OPTFLAGS)
         LIB =           -ltbbmalloc

In Long-range thread mode (LRT) a modified verlet style is used, that
operates the Kspace calculation in a separate thread concurrently to
other calculations. This has to be enabled in the :doc:`package intel
<package>` command at runtime. With the setting "threads" it used the
pthreads library, while "c++11" will use the built-in thread support
of C++11 compilers. The option "none" skips compilation of this
feature. The default is to use "threads" if pthreads is available and
otherwise "none".

Best performance is achieved with Intel hardware, Intel compilers, as
well as the Intel TBB and MKL libraries. However, the code also
compiles, links, and runs with other compilers / hardware and without
TBB and MKL.

----------

.. _mdi:

MDI package
-----------------------------

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D DOWNLOAD_MDI=value    # download MDI Library for build, value = no (default) or yes

   .. tab:: Traditional make

      Before building LAMMPS, you must build the MDI Library in
      ``lib/mdi``\ .  You can do this by executing a command like one
      of the following from the ``lib/mdi`` directory:

      .. code-block:: bash

         python Install.py -m gcc       # build using gcc compiler
         python Install.py -m icc       # build using icc compiler

      The build should produce two files: ``lib/mdi/includelink/mdi.h``
      and ``lib/mdi/liblink/libmdi.so``\ .

----------

.. _molfile:

MOLFILE package
---------------------------------------

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D MOLFILE_INCLUDE_DIR=path   # (optional) path where VMD molfile plugin headers are installed
         -D PKG_MOLFILE=yes

      Using ``-D PKG_MOLFILE=yes`` enables the package, and setting
      ``-D MOLFILE_INCLUDE_DIR`` allows to provide a custom location for
      the molfile plugin header files. These should match the ABI of the
      plugin files used, and thus one typically sets them to include
      folder of the local VMD installation in use. LAMMPS ships with a
      couple of default header files that correspond to a popular VMD
      version, usually the latest release.

   .. tab:: Traditional make

      The ``lib/molfile/Makefile.lammps`` file has a setting for a
      dynamic loading library libdl.a that is typically present on all
      systems.  It is required for LAMMPS to link with this package.  If
      the setting is not valid for your system, you will need to edit
      the Makefile.lammps file.  See ``lib/molfile/README`` and
      ``lib/molfile/Makefile.lammps`` for details. It is also possible
      to configure a different folder with the VMD molfile plugin header
      files. LAMMPS ships with a couple of default headers, but these
      are not compatible with all VMD versions, so it is often best to
      change this setting to the location of the same include files of
      the local VMD installation in use.

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

      The ``lib/netcdf/Makefile.lammps`` file has settings for NetCDF
      include and library files which LAMMPS needs to build with this
      package.  If the settings are not valid for your system, you will
      need to edit the ``Makefile.lammps`` file.  See
      ``lib/netcdf/README`` for details.

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
      libqmmm.a) are not included in the static LAMMPS library and
      (currently) not installed, while their code is included in the
      shared LAMMPS library.  Thus a typical command line to configure
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

      Before building LAMMPS, you must build the QMMM library in
      ``lib/qmmm``.  You can do this manually if you prefer; follow the
      first two steps explained in ``lib/qmmm/README``.  You can also do
      it in one step from the ``lammps/src`` dir, using a command like
      these, which simply invokes the ``lib/qmmm/Install.py`` script with
      the specified args:

      .. code-block:: bash

         make lib-qmmm                      # print help message
         make lib-qmmm args="-m serial"     # build with GNU Fortran compiler (settings as in "make serial")
         make lib-qmmm args="-m mpi"        # build with default MPI compiler (settings as in "make mpi")
         make lib-qmmm args="-m gfortran"   # build with GNU Fortran compiler

      The build should produce two files: ``lib/qmmm/libqmmm.a`` and
      ``lib/qmmm/Makefile.lammps``.  The latter is copied from an
      existing ``Makefile.lammps.*`` and has settings needed to build
      LAMMPS with the QMMM library (though typically the settings are
      just blank).  If necessary, you can edit/create a new
      ``lib/qmmm/Makefile.<machine>`` file for your system, which should
      define an ``EXTRAMAKE`` variable to specify a corresponding
      ``Makefile.lammps.<machine>`` file.

      You can then install QMMM package and build LAMMPS in the usual
      manner.  After completing the LAMMPS build and compiling Quantum
      ESPRESSO with external library support (via "make couple"), go
      back to the ``lib/qmmm`` folder and follow the instructions in the
      README file to build the combined LAMMPS/QE QM/MM executable
      (pwqmmm.x) in the lib/qmmm folder.

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

         -D DOWNLOAD_QUIP=value       # download QUIP library for build, value = no (default) or yes
         -D QUIP_LIBRARY=path         # path to libquip.a (only needed if a custom location)
         -D USE_INTERNAL_LINALG=value # Use the internal linear algebra library instead of LAPACK
                                      #   value = no (default) or yes

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

      The download/build procedure for the QUIP library, described in
      ``lib/quip/README`` file requires setting two environment
      variables, ``QUIP_ROOT`` and ``QUIP_ARCH``.  These are accessed by
      the ``lib/quip/Makefile.lammps`` file which is used when you
      compile and link LAMMPS with this package.  You should only need
      to edit ``Makefile.lammps`` if the LAMMPS build can not use its
      settings to successfully build on your system.

----------

.. _scafacos:

SCAFACOS package
-----------------------------------------

To build with this package, you must download and build the
`ScaFaCoS Coulomb solver library <http://www.scafacos.de>`_

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

      You can download and build the ScaFaCoS library manually if you
      prefer; follow the instructions in ``lib/scafacos/README``.  You
      can also do it in one step from the ``lammps/src`` dir, using a
      command like these, which simply invokes the
      ``lib/scafacos/Install.py`` script with the specified args:

      .. code-block:: bash

         make lib-scafacos                         # print help message
         make lib-scafacos args="-b"               # download and build in lib/scafacos/scafacos-<version>
         make lib-scafacos args="-p $HOME/scafacos  # use existing ScaFaCoS installation in $HOME/scafacos

      Note that 2 symbolic (soft) links, ``includelink`` and ``liblink``, are
      created in ``lib/scafacos`` to point to the ScaFaCoS src dir.  When LAMMPS
      builds in src it will use these links.  You should not need to edit
      the ``lib/scafacos/Makefile.lammps`` file.

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
         -D EIGEN3_INCLUDE_DIR=path    # path to Eigen library (only needed if a custom location)

      If ``DOWNLOAD_EIGEN3`` is set, the Eigen3 library will be
      downloaded and inside the CMake build directory.  If the Eigen3
      library is already on your system (in a location where CMake
      cannot find it), set ``EIGEN3_INCLUDE_DIR`` to the directory the
      ``Eigen3`` include file is in.

   .. tab:: Traditional make

      You can download the Eigen3 library manually if you prefer; follow
      the instructions in ``lib/smd/README``.  You can also do it in one
      step from the ``lammps/src`` dir, using a command like these,
      which simply invokes the ``lib/smd/Install.py`` script with the
      specified args:

      .. code-block:: bash

         make lib-smd                         # print help message
         make lib-smd args="-b"               # download to lib/smd/eigen3
         make lib-smd args="-p /usr/include/eigen3"    # use existing Eigen installation in /usr/include/eigen3

      Note that a symbolic (soft) link named ``includelink`` is created
      in ``lib/smd`` to point to the Eigen dir.  When LAMMPS builds it
      will use this link.  You should not need to edit the
      ``lib/smd/Makefile.lammps`` file.

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

      The ``lib/vtk/Makefile.lammps`` file has settings for accessing
      VTK files and its library, which LAMMPS needs to build with this
      package.  If the settings are not valid for your system, check if
      one of the other ``lib/vtk/Makefile.lammps.*`` files is compatible
      and copy it to Makefile.lammps.  If none of the provided files
      work, you will need to edit the ``Makefile.lammps`` file.  See
      ``lib/vtk/README`` for details.
