Packages with extra build options
=================================

When building with some packages, additional steps may be required,
in addition to:


.. code-block:: bash

   -D PKG_NAME=yes    # CMake
   make yes-name      # make

as described on the :doc:`Build\_package <Build_package>` doc page.

For a CMake build there may be additional optional or required
variables to set.  For a build with make, a provided library under the
lammps/lib directory may need to be built first.  Or an external
library may need to exist on your system or be downloaded and built.
You may need to tell LAMMPS where it is found on your system.

This is the list of packages that may require additional steps.

+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`COMPRESS <compress>`       | :ref:`GPU <gpu>`                 | :ref:`KIM <kim>`                   | :ref:`KOKKOS <kokkos>`       | :ref:`LATTE <latte>`           | :ref:`MESSAGE <message>`             |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`MSCG <mscg>`               | :ref:`OPT <opt>`                 | :ref:`POEMS <poems>`               | :ref:`PYTHON <python>`       | :ref:`VORONOI <voronoi>`       | :ref:`USER-ADIOS <user-adios>`       |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`USER-ATC <user-atc>`       | :ref:`USER-AWPMD <user-awpmd>`   | :ref:`USER-COLVARS <user-colvars>` | :ref:`USER-H5MD <user-h5md>` | :ref:`USER-INTEL <user-intel>` | :ref:`USER-MOLFILE <user-molfile>`   |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`USER-NETCDF <user-netcdf>` | :ref:`USER-PLUMED <user-plumed>` | :ref:`USER-OMP <user-omp>`         | :ref:`USER-QMMM <user-qmmm>` | :ref:`USER-QUIP <user-quip>`   | :ref:`USER-SCAFACOS <user-scafacos>` |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`USER-SMD <user-smd>`       | :ref:`USER-VTK <user-vtk>`       |                                    |                              |                                |                                      |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+


----------


.. _compress:

COMPRESS package
-------------------------------

To build with this package you must have the zlib compression library
available on your system.

**CMake build**\ :

If CMake cannot find the library, you can set these variables:


.. code-block:: bash

   -D ZLIB_INCLUDE_DIR=path    # path to zlib.h header file
   -D ZLIB_LIBRARIES=path      # path to libz.a (.so) file

**Traditional make**\ :

If make cannot find the library, you can edit the file
lib/compress/Makefile.lammps to specify the paths and library
name.


----------


.. _gpu:

GPU package
---------------------

To build with this package, you must choose options for precision and
which GPU hardware to build for.

**CMake build**\ :


.. code-block:: bash

   -D GPU_API=value          # value = opencl (default) or cuda
   -D GPU_PREC=value         # precision setting
                             # value = double or mixed (default) or single
   -D OCL_TUNE=value         # hardware choice for GPU_API=opencl
                             # generic (default) or intel (Intel CPU) or fermi, kepler, cypress (NVIDIA)
   -D GPU_ARCH=value         # primary GPU hardware choice for GPU_API=cuda
                             # value = sm_XX, see below
                             # default is sm_30
   -D CUDPP_OPT=value        # optimization setting for GPU_API=cuda
                             # enables CUDA Performance Primitives Optimizations
                             # value = yes (default) or no
   -D CUDA_MPS_SUPPORT=value # enables some tweaks required to run with active nvidia-cuda-mps daemon
                             # value = yes or no (default)

GPU\_ARCH settings for different GPU hardware is as follows:

* sm\_12 or sm\_13 for GT200 (supported by CUDA 3.2 until CUDA 6.5)
* sm\_20 or sm\_21 for Fermi (supported by CUDA 3.2 until CUDA 7.5)
* sm\_30 or sm\_35 or sm\_37 for Kepler (supported since CUDA 5)
* sm\_50 or sm\_52 for Maxwell (supported since CUDA 6)
* sm\_60 or sm\_61 for Pascal (supported since CUDA 8)
* sm\_70 for Volta (supported since CUDA 9)
* sm\_75 for Turing (supported since CUDA 10)

A more detailed list can be found, for example,
at `Wikipedia's CUDA article <https://en.wikipedia.org/wiki/CUDA#GPUs_supported>`_

CMake can detect which version of the CUDA toolkit is used and thus can
include support for **all** major GPU architectures supported by this toolkit.
Thus the GPU\_ARCH setting is merely an optimization, to have code for
the preferred GPU architecture directly included rather than having to wait
for the JIT compiler of the CUDA driver to translate it.

**Traditional make**\ :

Before building LAMMPS, you must build the GPU library in lib/gpu.
You can do this manually if you prefer; follow the instructions in
lib/gpu/README.  Note that the GPU library uses MPI calls, so you must
use the same MPI library (or the STUBS library) settings as the main
LAMMPS code.  This also applies to the -DLAMMPS\_BIGBIG,
-DLAMMPS\_SMALLBIG, or -DLAMMPS\_SMALLSMALL settings in whichever
Makefile you use.

You can also build the library in one step from the lammps/src dir,
using a command like these, which simply invoke the lib/gpu/Install.py
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

* CUDA\_HOME = where NVIDIA CUDA software is installed on your system
* CUDA\_ARCH = sm\_XX, what GPU hardware you have, same as CMake GPU\_ARCH above
* CUDA\_PRECISION = precision (double, mixed, single)
* EXTRAMAKE = which Makefile.lammps.\* file to copy to Makefile.lammps

The file Makefile.linux\_multi is set up to include support for multiple
GPU architectures as supported by the CUDA toolkit in use. This is done
through using the "--gencode " flag, which can be used multiple times and
thus support all GPU architectures supported by your CUDA compiler.

If the library build is successful, 3 files should be created:
lib/gpu/libgpu.a, lib/gpu/nvc\_get\_devices, and
lib/gpu/Makefile.lammps.  The latter has settings that enable LAMMPS
to link with CUDA libraries.  If the settings in Makefile.lammps for
your machine are not correct, the LAMMPS build will fail, and
lib/gpu/Makefile.lammps may need to be edited.

.. note::

   If you re-build the GPU library in lib/gpu, you should always
   un-install the GPU package in lammps/src, then re-install it and
   re-build LAMMPS.  This is because the compilation of files in the GPU
   package uses the library settings from the lib/gpu/Makefile.machine
   used to build the GPU library.


----------


.. _kim:

KIM package
---------------------

To build with this package, the KIM library with API v2 must be downloaded
and built on your system.  It must include the KIM models that you want to
use with LAMMPS. If you want to use the :doc:`kim_query <kim_commands>`
command, you also need to have libcurl installed with the matching
development headers and the curl-config tool.

See `Obtaining KIM Models <http://openkim.org/doc/usage/obtaining-models>`_ to
learn how to install a pre-build binary of the OpenKIM Repository of Models.
See the list of all KIM models here: https://openkim.org/browse/models

(Also note that when downloading and installing from source
the KIM API library with all its models, may take a long time (tens of
minutes to hours) to build.  Of course you only need to do that once.)

**CMake build**\ :


.. code-block:: bash

   -D DOWNLOAD_KIM=value           # download OpenKIM API v2 for build, value = no (default) or yes
   -D LMP_DEBUG_CURL=value         # set libcurl verbose mode on/off, value = off (default) or on
   -D LMP_NO_SSL_CHECK=value       # tell libcurl to not verify the peer, value = no (default) or yes

If DOWNLOAD\_KIM is set, the KIM library will be downloaded and built
inside the CMake build directory.  If the KIM library is already on
your system (in a location CMake cannot find it), set the PKG\_CONFIG\_PATH
environment variable so that libkim-api can be found.

*For using OpenKIM web queries in LAMMPS*\ :

If LMP\_DEBUG\_CURL is set, the libcurl verbose mode will be on, and any
libcurl calls within the KIM web query display a lot of information about
libcurl operations. You hardly ever want this set in production use, you will
almost always want this when you debug/report problems.

The libcurl performs peer SSL certificate verification by default. This
verification is done using a CA certificate store that the SSL library can
use to make sure the peer's server certificate is valid. If SSL reports an
error ("certificate verify failed") during the handshake and thus refuses
further communication with that server, you can set LMP\_NO\_SSL\_CHECK.
If LMP\_NO\_SSL\_CHECK is set, libcurl does not verify the peer and connection
succeeds regardless of the names in the certificate. This option is insecure.
As an alternative, you can specify your own CA cert path by setting the
environment variable CURL\_CA\_BUNDLE to the path of your choice. A call to the
KIM web query would get this value from the environmental variable.

**Traditional make**\ :

You can download and build the KIM library manually if you prefer;
follow the instructions in lib/kim/README.  You can also do it in one
step from the lammps/src dir, using a command like these, which simply
invoke the lib/kim/Install.py script with the specified args.


.. code-block:: bash

   make lib-kim              # print help message
   make lib-kim args="-b "   # (re-)install KIM API lib with only example models
   make lib-kim args="-b -a Glue_Ercolessi_Adams_Al__MO_324507536345_001"  # ditto plus one model
   make lib-kim args="-b -a everything"     # install KIM API lib with all models
   make lib-kim args="-n -a EAM_Dynamo_Ackland_W__MO_141627196590_002"       # add one model or model driver
   make lib-kim args="-p /usr/local" # use an existing KIM API installation at the provided location
   make lib-kim args="-p /usr/local -a EAM_Dynamo_Ackland_W__MO_141627196590_002" # ditto but add one model or driver

Settings for OpenKIM web queries discussed above need to be applied by adding
them to the LMP\_INC variable through editing the Makefile.machine you are
using.  For example:

.. code-block:: make

   LMP_INC =       -DLMP_NO_SSL_CHECK

----------


.. _kokkos:

KOKKOS package
---------------------------

To build with this package, you must choose which hardware you want to
build for, either CPUs (multi-threading via OpenMP) or KNLs (OpenMP)
or GPUs (NVIDIA Cuda).

For a CMake or make build, these are the possible choices for the
KOKKOS\_ARCH settings described below.  Note that for CMake, these are
really Kokkos variables, not LAMMPS variables.  Hence you must use
case-sensitive values, e.g. BDW, not bdw.

* ARMv80 = ARMv8.0 Compatible CPU
* ARMv81 = ARMv8.1 Compatible CPU
* ARMv8-ThunderX = ARMv8 Cavium ThunderX CPU
* BGQ = IBM Blue Gene/Q CPUs
* Power8 = IBM POWER8 CPUs
* Power9 = IBM POWER9 CPUs
* SNB = Intel Sandy/Ivy Bridge CPUs
* HSW = Intel Haswell CPUs
* BDW = Intel Broadwell Xeon E-class CPUs
* SKX = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)
* KNC = Intel Knights Corner Xeon Phi
* KNL = Intel Knights Landing Xeon Phi
* Kepler30 = NVIDIA Kepler generation CC 3.0
* Kepler32 = NVIDIA Kepler generation CC 3.2
* Kepler35 = NVIDIA Kepler generation CC 3.5
* Kepler37 = NVIDIA Kepler generation CC 3.7
* Maxwell50 = NVIDIA Maxwell generation CC 5.0
* Maxwell52 = NVIDIA Maxwell generation CC 5.2
* Maxwell53 = NVIDIA Maxwell generation CC 5.3
* Pascal60 = NVIDIA Pascal generation CC 6.0
* Pascal61 = NVIDIA Pascal generation CC 6.1
* Volta70 = NVIDIA Volta generation CC 7.0
* Volta72 = NVIDIA Volta generation CC 7.2
* Turing75 = NVIDIA Turing generation CC 7.5

**CMake build**\ :

For multicore CPUs using OpenMP, set these 2 variables.


.. code-block:: bash

   -D KOKKOS_ARCH=archCPU         # archCPU = CPU from list above
   -D KOKKOS_ENABLE_OPENMP=yes

For Intel KNLs using OpenMP, set these 2 variables:


.. code-block:: bash

   -D KOKKOS_ARCH=KNL
   -D KOKKOS_ENABLE_OPENMP=yes

For NVIDIA GPUs using CUDA, set these 4 variables:


.. code-block:: bash

   -D KOKKOS_ARCH="archCPU;archGPU"   # archCPU = CPU from list above that is hosting the GPU
                                      # archGPU = GPU from list above
   -D KOKKOS_ENABLE_CUDA=yes
   -D KOKKOS_ENABLE_OPENMP=yes
   -D CMAKE_CXX_COMPILER=wrapper      # wrapper = full path to Cuda nvcc wrapper

The wrapper value is the Cuda nvcc compiler wrapper provided in the
Kokkos library: lib/kokkos/bin/nvcc\_wrapper.  The setting should
include the full path name to the wrapper, e.g.


.. code-block:: bash

   -D CMAKE_CXX_COMPILER=/home/username/lammps/lib/kokkos/bin/nvcc_wrapper

**Traditional make**\ :

Choose which hardware to support in Makefile.machine via
KOKKOS\_DEVICES and KOKKOS\_ARCH settings.  See the
src/MAKE/OPTIONS/Makefile.kokkos\* files for examples.

For multicore CPUs using OpenMP:


.. code-block:: make

   KOKKOS_DEVICES = OpenMP
   KOKKOS_ARCH = archCPU      # archCPU = CPU from list above

For Intel KNLs using OpenMP:


.. code-block:: make

   KOKKOS_DEVICES = OpenMP
   KOKKOS_ARCH = KNL

For NVIDIA GPUs using CUDA:


.. code-block:: make

   KOKKOS_DEVICES = Cuda
   KOKKOS_ARCH = archCPU,archGPU    # archCPU = CPU from list above that is hosting the GPU
                                    # archGPU = GPU from list above
   FFT_INC = -DFFT_CUFFT            # enable use of cuFFT (optional)
   FFT_LIB = -lcufft                # link to cuFFT library

For GPUs, you also need the following 2 lines in your Makefile.machine
before the CC line is defined, in this case for use with OpenMPI mpicxx.
The 2 lines define a nvcc wrapper compiler, which will use nvcc for
compiling CUDA files and use a C++ compiler for non-Kokkos, non-CUDA
files.


.. code-block:: make

   KOKKOS_ABSOLUTE_PATH = $(shell cd $(KOKKOS_PATH); pwd)
   export OMPI_CXX = $(KOKKOS_ABSOLUTE_PATH)/config/nvcc_wrapper
   CC =            mpicxx


----------


.. _latte:

LATTE package
-------------------------

To build with this package, you must download and build the LATTE
library.

**CMake build**\ :


.. code-block:: bash

   -D DOWNLOAD_LATTE=value    # download LATTE for build, value = no (default) or yes
   -D LATTE_LIBRARY=path      # LATTE library file (only needed if a custom location)

If DOWNLOAD\_LATTE is set, the LATTE library will be downloaded and
built inside the CMake build directory.  If the LATTE library is
already on your system (in a location CMake cannot find it),
LATTE\_LIBRARY is the filename (plus path) of the LATTE library file,
not the directory the library file is in.

**Traditional make**\ :

You can download and build the LATTE library manually if you prefer;
follow the instructions in lib/latte/README.  You can also do it in
one step from the lammps/src dir, using a command like these, which
simply invokes the lib/latte/Install.py script with the specified
args:


.. code-block:: bash

   make lib-latte                          # print help message
   make lib-latte args="-b"                # download and build in lib/latte/LATTE-master
   make lib-latte args="-p $HOME/latte"    # use existing LATTE installation in $HOME/latte
   make lib-latte args="-b -m gfortran"    # download and build in lib/latte and
                                           #   copy Makefile.lammps.gfortran to Makefile.lammps

Note that 3 symbolic (soft) links, "includelink" and "liblink" and
"filelink.o", are created in lib/latte to point into the LATTE home
dir.  When LAMMPS itself is built it will use these links.  You should
also check that the Makefile.lammps file you create is appropriate for
the compiler you use on your system to build LATTE.


----------


.. _message:

MESSAGE package
-----------------------------

This package can optionally include support for messaging via sockets,
using the open-source `ZeroMQ library <http://zeromq.org>`_, which must
be installed on your system.

**CMake build**\ :


.. code-block:: bash

   -D MESSAGE_ZMQ=value    # build with ZeroMQ support, value = no (default) or yes
   -D ZMQ_LIBRARY=path     # ZMQ library file (only needed if a custom location)
   -D ZMQ_INCLUDE_DIR=path # ZMQ include directory (only needed if a custom location)

**Traditional make**\ :

Before building LAMMPS, you must build the CSlib library in
lib/message.  You can build the CSlib library manually if you prefer;
follow the instructions in lib/message/README.  You can also do it in
one step from the lammps/src dir, using a command like these, which
simply invoke the lib/message/Install.py script with the specified args:


.. code-block:: bash

   make lib-message               # print help message
   make lib-message args="-m -z"  # build with MPI and socket (ZMQ) support
   make lib-message args="-s"     # build as serial lib with no ZMQ support

The build should produce two files: lib/message/cslib/src/libmessage.a
and lib/message/Makefile.lammps.  The latter is copied from an
existing Makefile.lammps.\* and has settings to link with the ZeroMQ
library if requested in the build.


----------


.. _mscg:

MSCG package
-----------------------

To build with this package, you must download and build the MS-CG
library.  Building the MS-CG library and using it from LAMMPS requires
a C++11 compatible compiler and that the GSL (GNU Scientific Library)
headers and libraries are installed on your machine.  See the
lib/mscg/README and MSCG/Install files for more details.

**CMake build**\ :


.. code-block:: bash

   -D DOWNLOAD_MSCG=value    # download MSCG for build, value = no (default) or yes
   -D MSCG_LIBRARY=path      # MSCG library file (only needed if a custom location)
   -D MSCG_INCLUDE_DIR=path  # MSCG include directory (only needed if a custom location)

If DOWNLOAD\_MSCG is set, the MSCG library will be downloaded and built
inside the CMake build directory.  If the MSCG library is already on
your system (in a location CMake cannot find it), MSCG\_LIBRARY is the
filename (plus path) of the MSCG library file, not the directory the
library file is in.  MSCG\_INCLUDE\_DIR is the directory the MSCG
include file is in.

**Traditional make**\ :

You can download and build the MS-CG library manually if you prefer;
follow the instructions in lib/mscg/README.  You can also do it in one
step from the lammps/src dir, using a command like these, which simply
invoke the lib/mscg/Install.py script with the specified args:


.. code-block:: bash

   make lib-mscg             # print help message
   make lib-mscg args="-b -m serial"   # download and build in lib/mscg/MSCG-release-master
                                       # with the settings compatible with "make serial"
   make lib-mscg args="-b -m mpi"      # download and build in lib/mscg/MSCG-release-master
                                       # with the settings compatible with "make mpi"
   make lib-mscg args="-p /usr/local/mscg-release" # use the existing MS-CG installation in /usr/local/mscg-release

Note that 2 symbolic (soft) links, "includelink" and "liblink", will
be created in lib/mscg to point to the MS-CG src/installation dir.
When LAMMPS is built in src it will use these links.  You should not
need to edit the lib/mscg/Makefile.lammps file.


----------


.. _opt:

OPT package
---------------------

**CMake build**\ :

No additional settings are needed besides "-D PKG\_OPT=yes".

**Traditional make**\ :

The compile flag "-restrict" must be used to build LAMMPS with the OPT
package when using Intel compilers.  It should be added to the CCFLAGS
line of your Makefile.machine.  See src/MAKE/OPTIONS/Makefile.opt for
an example.


----------


.. _poems:

POEMS package
-------------------------

**CMake build**\ :

No additional settings are needed besides "-D PKG\_OPT=yes".

**Traditional make**\ :

Before building LAMMPS, you must build the POEMS library in lib/poems.
You can do this manually if you prefer; follow the instructions in
lib/poems/README.  You can also do it in one step from the lammps/src
dir, using a command like these, which simply invoke the
lib/poems/Install.py script with the specified args:


.. code-block:: bash

   make lib-poems                   # print help message
   make lib-poems args="-m serial"  # build with GNU g++ compiler (settings as with "make serial")
   make lib-poems args="-m mpi"     # build with default MPI C++ compiler (settings as with "make mpi")
   make lib-poems args="-m icc"     # build with Intel icc compiler

The build should produce two files: lib/poems/libpoems.a and
lib/poems/Makefile.lammps.  The latter is copied from an existing
Makefile.lammps.\* and has settings needed to build LAMMPS with the
POEMS library (though typically the settings are just blank).  If
necessary, you can edit/create a new lib/poems/Makefile.machine file
for your system, which should define an EXTRAMAKE variable to specify
a corresponding Makefile.lammps.machine file.


----------


.. _python:

PYTHON package
---------------------------

Building with the PYTHON package requires you have a Python shared
library available on your system, which needs to be a Python 2
version, 2.6 or later.  Python 3 is not yet supported.  See
lib/python/README for more details.

**CMake build**\ :


.. code-block:: bash

   -D PYTHON_EXECUTABLE=path   # path to Python executable to use

Without this setting, CMake will guess the default Python on your
system.  To use a different Python version, you can either create a
virtualenv, activate it and then run cmake.  Or you can set the
PYTHON\_EXECUTABLE variable to specify which Python interpreter should
be used.  Note note that you will also need to have the development
headers installed for this version, e.g. python2-devel.

**Traditional make**\ :

The build uses the lib/python/Makefile.lammps file in the compile/link
process to find Python.  You should only need to create a new
Makefile.lammps.\* file (and copy it to Makefile.lammps) if the LAMMPS
build fails.


----------


.. _voronoi:

VORONOI package
-----------------------------

To build with this package, you must download and build the `Voro++ library <voro-home_>`_.

.. _voro-home: http://math.lbl.gov/voro++



**CMake build**\ :


.. code-block:: bash

   -D DOWNLOAD_VORO=value    # download Voro++ for build, value = no (default) or yes
   -D VORO_LIBRARY=path      # Voro++ library file (only needed if at custom location)
   -D VORO_INCLUDE_DIR=path  # Voro++ include directory (only needed if at custom location)

If DOWNLOAD\_VORO is set, the Voro++ library will be downloaded and
built inside the CMake build directory.  If the Voro++ library is
already on your system (in a location CMake cannot find it),
VORO\_LIBRARY is the filename (plus path) of the Voro++ library file,
not the directory the library file is in.  VORO\_INCLUDE\_DIR is the
directory the Voro++ include file is in.

**Traditional make**\ :

You can download and build the Voro++ library manually if you prefer;
follow the instructions in lib/voronoi/README.  You can also do it in
one step from the lammps/src dir, using a command like these, which
simply invoke the lib/voronoi/Install.py script with the specified
args:


.. code-block:: bash

   make lib-voronoi                          # print help message
   make lib-voronoi args="-b"                # download and build the default version in lib/voronoi/voro++-<version>
   make lib-voronoi args="-p $HOME/voro++"   # use existing Voro++ installation in $HOME/voro++
   make lib-voronoi args="-b -v voro++0.4.6" # download and build the 0.4.6 version in lib/voronoi/voro++-0.4.6

Note that 2 symbolic (soft) links, "includelink" and "liblink", are
created in lib/voronoi to point to the Voro++ src dir.  When LAMMPS
builds in src it will use these links.  You should not need to edit
the lib/voronoi/Makefile.lammps file.


----------


.. _user-adios:

USER-ADIOS package
-----------------------------------

The USER-ADIOS package requires the `ADIOS I/O library <https://github.com/ornladios/ADIOS2>`_,
version 2.3.1 or newer. Make sure that you have ADIOS built either with or
without MPI to match if you build LAMMPS with or without MPI.
ADIOS compilation settings for LAMMPS are automatically detected, if the PATH
and LD\_LIBRARY\_PATH environment variables have been updated for the local ADIOS
installation and the instructions below are followed for the respective build systems.

**CMake build**\ :


.. code-block:: bash

   -D ADIOS2_DIR=path        # path is where ADIOS 2.x is installed
   -D PKG_USER-ADIOS=yes

**Traditional make**\ :

Turn on the USER-ADIOS package before building LAMMPS. If the ADIOS 2.x software is installed in PATH, there is nothing else to do:


.. code-block:: bash

   make yes-user-adios

otherwise, set ADIOS2\_DIR environment variable when turning on the package:


.. code-block:: bash

   ADIOS2_DIR=path make yes-user-adios   # path is where ADIOS 2.x is installed


----------


.. _user-atc:

USER-ATC package
-------------------------------

The USER-ATC package requires the MANYBODY package also be installed.

**CMake build**\ :

No additional settings are needed besides "-D PKG\_USER-ATC=yes"
and "-D PKG\_MANYBODY=yes".

**Traditional make**\ :

Before building LAMMPS, you must build the ATC library in lib/atc.
You can do this manually if you prefer; follow the instructions in
lib/atc/README.  You can also do it in one step from the lammps/src
dir, using a command like these, which simply invoke the
lib/atc/Install.py script with the specified args:


.. code-block:: bash

   make lib-atc                      # print help message
   make lib-atc args="-m serial"     # build with GNU g++ compiler and MPI STUBS (settings as with "make serial")
   make lib-atc args="-m mpi"        # build with default MPI compiler (settings as with "make mpi")
   make lib-atc args="-m icc"        # build with Intel icc compiler

The build should produce two files: lib/atc/libatc.a and
lib/atc/Makefile.lammps.  The latter is copied from an existing
Makefile.lammps.\* and has settings needed to build LAMMPS with the ATC
library.  If necessary, you can edit/create a new
lib/atc/Makefile.machine file for your system, which should define an
EXTRAMAKE variable to specify a corresponding Makefile.lammps.machine
file.

Note that the Makefile.lammps file has settings for the BLAS and
LAPACK linear algebra libraries.  As explained in lib/atc/README these
can either exist on your system, or you can use the files provided in
lib/linalg.  In the latter case you also need to build the library in
lib/linalg with a command like these:


.. code-block:: bash

   make lib-linalg                     # print help message
   make lib-linalg args="-m serial"    # build with GNU Fortran compiler (settings as with "make serial")
   make lib-linalg args="-m mpi"       # build with default MPI Fortran compiler (settings as with "make mpi")
   make lib-linalg args="-m gfortran"  # build with GNU Fortran compiler


----------


.. _user-awpmd:

USER-AWPMD package
-----------------------------------

**CMake build**\ :

No additional settings are needed besides "-D PKG\_USER-AQPMD=yes".

**Traditional make**\ :

Before building LAMMPS, you must build the AWPMD library in lib/awpmd.
You can do this manually if you prefer; follow the instructions in
lib/awpmd/README.  You can also do it in one step from the lammps/src
dir, using a command like these, which simply invoke the
lib/awpmd/Install.py script with the specified args:


.. code-block:: bash

   make lib-awpmd                   # print help message
   make lib-awpmd args="-m serial"  # build with GNU g++ compiler and MPI STUBS (settings as with "make serial")
   make lib-awpmd args="-m mpi"     # build with default MPI compiler (settings as with "make mpi")
   make lib-awpmd args="-m icc"     # build with Intel icc compiler

The build should produce two files: lib/awpmd/libawpmd.a and
lib/awpmd/Makefile.lammps.  The latter is copied from an existing
Makefile.lammps.\* and has settings needed to build LAMMPS with the
AWPMD library.  If necessary, you can edit/create a new
lib/awpmd/Makefile.machine file for your system, which should define
an EXTRAMAKE variable to specify a corresponding
Makefile.lammps.machine file.

Note that the Makefile.lammps file has settings for the BLAS and
LAPACK linear algebra libraries.  As explained in lib/awpmd/README
these can either exist on your system, or you can use the files
provided in lib/linalg.  In the latter case you also need to build the
library in lib/linalg with a command like these:


.. code-block:: bash

   make lib-linalg                     # print help message
   make lib-linalg args="-m serial"    # build with GNU Fortran compiler (settings as with "make serial")
   make lib-linalg args="-m mpi"       # build with default MPI Fortran compiler (settings as with "make mpi")
   make lib-linalg args="-m gfortran"  # build with GNU Fortran compiler


----------


.. _user-colvars:

USER-COLVARS package
---------------------------------------

This package includes into the LAMMPS distribution the Colvars library, which
can be built for the most part with all major versions of the C++ language.

A few of the most recent features require C++11 support.  In particular, the
library is optionally built together with the
`Lepton <https://simtk.org/projects/lepton>`_ library, a copy of which is also
included in the LAMMPS distribution.  Lepton implements the
`customFunction <http://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html#colvar|customFunction>`_
feature, and requires C++11 support.

See `here <https://colvars.github.io/README-c++11.html>`_ for a detailed list of
C++11-only features.

**CMake build**\ :

This is the recommended build recipe: no additional settings are normally
needed besides "-D PKG\_USER-COLVARS=yes".

Building and linking of Lepton (or other C++11-only features) is enabled
automatically when compilation is carried out with C++11 support, and disabled
otherwise.  Optionally, Lepton build may be manually controlled with the flag
"-D COLVARS\_LEPTON=yes\|no".

**Traditional make**\ :

Before building LAMMPS, one must build the Colvars library in lib/colvars.

This can be done manually in the same folder by using or adapting one of the
provided Makefiles: for example, Makefile.g++ for the GNU compiler.

In general, it is safer to use build setting consistent with the rest of
LAMMPS.  This is best carried out from the LAMMPS src directory using a
command like these, which simply invoke the lib/colvars/Install.py script with
the specified args:


.. code-block:: bash

   make lib-colvars                      # print help message
   make lib-colvars args="-m serial"     # build with GNU g++ compiler (settings as with "make serial")
   make lib-colvars args="-m mpi"        # build with default MPI compiler (settings as with "make mpi")
   make lib-colvars args="-m g++-debug"  # build with GNU g++ compiler and colvars debugging enabled

The "machine" argument of the "-m" flag is used to find a Makefile.machine to
use as build recipe.  If it does not already exist in lib/colvars, it will be
auto-generated by using compiler flags consistent with those parsed from the
core LAMMPS makefiles.

Optional flags may be specified as environment variables:

.. code-block:: bash

    COLVARS_DEBUG=yes make lib-colvars args="-m machine"  # Build with debug code (much slower)
    COLVARS_LEPTON=no make lib-colvars args="-m machine"  # Build without Lepton (included otherwise)

The build should produce two files: the library lib/colvars/libcolvars.a
(which also includes Lepton objects if enabled) and the specification file
lib/colvars/Makefile.lammps.  The latter is auto-generated, and normally does
not need to be edited.


----------


.. _user-plumed:

USER-PLUMED package
-------------------------------------

.. _plumedinstall: http://plumed.github.io/doc-master/user-doc/html/\_installation.html

Before building LAMMPS with this package, you must first build PLUMED.
PLUMED can be built as part of the LAMMPS build or installed separately
from LAMMPS using the generic `plumed installation instructions <plumedinstall_>`_.
The USER-PLUMED package has been tested to work with Plumed versions
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
which PLUMED kernel should be used at runtime by using the PLUMED\_KERNEL
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

**CMake build**\ :

When the "-D PKG\_USER-PLUMED" flag is included in the cmake command you
must ensure that GSL is installed in locations that are specified in
your environment.  There are then two additional commands that control
the manner in which PLUMED is obtained and linked into LAMMPS.


.. code-block:: bash

   -D DOWNLOAD_PLUMED=value   # download PLUMED for build, value = no (default) or yes
   -D PLUMED_MODE=value       # Linkage mode for PLUMED, value = static (default), shared, or runtime

If DOWNLOAD\_PLUMED is set to "yes", the PLUMED library will be
downloaded (the version of PLUMED that will be downloaded is hard-coded
to a vetted version of PLUMED, usually a recent stable release version)
and built inside the CMake build directory.  If DOWNLOAD\_PLUMED is set
to "no" (the default), CMake will try to detect and link to an installed
version of PLUMED.  For this to work, the PLUMED library has to be
installed into a location where the pkg-config tool can find it or the
PKG\_CONFIG\_PATH environment variable has to be set up accordingly.
PLUMED should be installed in such a location if you compile it using
the default make; make install commands.

The PLUMED\_MODE setting determines the linkage mode for the PLUMED
library.  The allowed values for this flag are "static" (default),
"shared", or "runtime".  For a discussion of PLUMED linkage modes,
please see above.  When DOWNLOAD\_PLUMED is enabled the static linkage
mode is recommended.

**Traditional make**\ :

PLUMED needs to be installed before the USER-PLUMED package is installed
so that LAMMPS can find the right settings when compiling and linking
the LAMMPS executable.  You can either download and build PLUMED inside
the LAMMPS plumed library folder or use a previously installed PLUMED
library and point LAMMPS to its location. You also have to choose the
linkage mode: "static" (default), "shared" or "runtime".  For a
discussion of PLUMED linkage modes, please see above.

Download/compilation/configuration of the plumed library can be done
from the src folder through the following make args:


.. code-block:: bash

   make lib-plumed                         # print help message
   make lib-plumed args="-b"               # download and build PLUMED in lib/plumed/plumed2
   make lib-plumed args="-p $HOME/.local"  # use existing PLUMED installation in $HOME/.local
   make lib-plumed args="-p /usr/local -m shared"  # use existing PLUMED installation in
                                                   # /usr/local and use shared linkage mode

Note that 2 symbolic (soft) links, "includelink" and "liblink" are
created in lib/plumed that point to the location of the PLUMED build to
use. A new file lib/plumed/Makefile.lammps is also created with settings
suitable for LAMMPS to compile and link PLUMED using the desired linkage
mode. After this step is completed, you can install the USER-PLUMED
package and compile LAMMPS in the usual manner:


.. code-block:: bash

   make yes-user-plumed
   make machine

Once this compilation completes you should be able to run LAMMPS in the
usual way.  For shared linkage mode, libplumed.so must be found by the
LAMMPS executable, which on many operating systems means, you have to
set the LD\_LIBRARY\_PATH environment variable accordingly.

Support for the different linkage modes in LAMMPS varies for different
operating systems, using the static linkage is expected to be the most
portable, and thus set to be the default.

If you want to change the linkage mode, you have to re-run "make
lib-plumed" with the desired settings **and** do a re-install if the
USER-PLUMED package with "make yes-user-plumed" to update the required
makefile settings with the changes in the lib/plumed folder.


----------


.. _user-h5md:

USER-H5MD package
---------------------------------

To build with this package you must have the HDF5 software package
installed on your system, which should include the h5cc compiler and
the HDF5 library.

**CMake build**\ :

No additional settings are needed besides "-D PKG\_USER-H5MD=yes".

This should auto-detect the H5MD library on your system.  Several
advanced CMake H5MD options exist if you need to specify where it is
installed.  Use the ccmake (terminal window) or cmake-gui (graphical)
tools to see these options and set them interactively from their user
interfaces.

**Traditional make**\ :

Before building LAMMPS, you must build the CH5MD library in lib/h5md.
You can do this manually if you prefer; follow the instructions in
lib/h5md/README.  You can also do it in one step from the lammps/src
dir, using a command like these, which simply invoke the
lib/h5md/Install.py script with the specified args:


.. code-block:: bash

   make lib-h5md                     # print help message
   make lib-h5md args="-m h5cc"      # build with h5cc compiler

The build should produce two files: lib/h5md/libch5md.a and
lib/h5md/Makefile.lammps.  The latter is copied from an existing
Makefile.lammps.\* and has settings needed to build LAMMPS with the
system HDF5 library.  If necessary, you can edit/create a new
lib/h5md/Makefile.machine file for your system, which should define an
EXTRAMAKE variable to specify a corresponding Makefile.lammps.machine
file.


----------


.. _user-intel:

USER-INTEL package
-----------------------------------

To build with this package, you must choose which hardware you want to
build for, either x86 CPUs or Intel KNLs in offload mode.  You should
also typically :ref:`install the USER-OMP package <user-omp>`, as it can be
used in tandem with the USER-INTEL package to good effect, as explained
on the :doc:`Speed intel <Speed_intel>` doc page.

**CMake build**\ :


.. code-block:: bash

   -D INTEL_ARCH=value     # value = cpu (default) or knl
   -D INTEL_LRT_MODE=value # value = threads, none, or c++11

In Long-range thread mode (LRT) a modified verlet style is used, that
operates the Kspace calculation in a separate thread concurrently to
other calculations. This has to be enabled in the :doc:`package intel <package>`
command at runtime. With the setting "threads" it used the pthreads
library, while c++11 will use the built-in thread support of C++11
compilers. The option "none" skips compilation of this feature. The
default is to use "threads" if pthreads is available and otherwise "none".

Best performance is achieved with Intel hardware, Intel compilers, as well as
the Intel TBB and MKL libraries. However, the code also compiles, links, and
runs with other compilers and without TBB and MKL.

**Traditional make**\ :

Choose which hardware to compile for in Makefile.machine via the
following settings.  See src/MAKE/OPTIONS/Makefile.intel\_cpu\* and
Makefile.knl files for examples. and src/USER-INTEL/README for
additional information.

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


----------


.. _user-molfile:

USER-MOLFILE package
---------------------------------------

**CMake build**\ :


.. code-block:: bash

   -D MOLFILE_INCLUDE_DIRS=path   # (optional) path where VMD molfile plugin headers are installed
   -D PKG_USER-MOLFILE=yes

Using "-D PKG\_USER-MOLFILE=yes" enables the package, and setting
"-D MOLFILE\_INCLUDE DIRS" allows to provide a custom location for
the molfile plugin header files. These should match the ABI of the
plugin files used, and thus one typically sets them to include
folder of the local VMD installation in use. LAMMPS ships with a
couple of default header files that correspond to a popular VMD
version, usually the latest release.

**Traditional make**\ :

The lib/molfile/Makefile.lammps file has a setting for a dynamic
loading library libdl.a that is typically present on all systems.  It
is required for LAMMPS to link with this package.  If the setting is
not valid for your system, you will need to edit the Makefile.lammps
file.  See lib/molfile/README and lib/molfile/Makefile.lammps for
details. It is also possible to configure a different folder with
the VMD molfile plugin header files. LAMMPS ships with a couple of
default headers, but these are not compatible with all VMD versions,
so it is often best to change this setting to the location of the
same include files of the local VMD installation in use.


----------


.. _user-netcdf:

USER-NETCDF package
-------------------------------------

To build with this package you must have the NetCDF library installed
on your system.

**CMake build**\ :

No additional settings are needed besides "-D PKG\_USER-NETCDF=yes".

This should auto-detect the NETCDF library if it is installed on your
system at standard locations.  Several advanced CMake NETCDF options
exist if you need to specify where it was installed.  Use the ccmake
(terminal window) or cmake-gui (graphical) tools to see these options
and set them interactively from their user interfaces.

**Traditional make**\ :

The lib/netcdf/Makefile.lammps file has settings for NetCDF include
and library files which LAMMPS needs to build with this package.  If
the settings are not valid for your system, you will need to edit the
Makefile.lammps file.  See lib/netcdf/README for details.


----------


.. _user-omp:

USER-OMP package
-------------------------------

**CMake build**\ :

No additional settings are required besides "-D PKG\_USER-OMP=yes".  If
CMake detects OpenMP support, the USER-OMP code will be compiled with
multi-threading support enabled, otherwise as optimized serial code.

**Traditional make**\ :

To enable multi-threading support in the USER-OMP package (and other
styles supporting OpenMP) the following compile and link flags must
be added to your Makefile.machine file.
See src/MAKE/OPTIONS/Makefile.omp for an example.


.. parsed-literal::

   CCFLAGS: -fopenmp               # for GNU and Clang Compilers 
   CCFLAGS: -qopenmp -restrict     # for Intel compilers on Linux
   LINKFLAGS: -fopenmp             # for GNU and Clang Compilers
   LINKFLAGS: -qopenmp             # for Intel compilers on Linux

For other platforms and compilers, please consult the documentation
about OpenMP support for your compiler. Please see the note about
how to address compatibility :ref:`issues with the 'default(none)' directive <default-none-issues>` of some compilers.


----------


.. _user-qmmm:

USER-QMMM package
---------------------------------

For using LAMMPS to do QM/MM simulations via the USER-QMMM package you
need to build LAMMPS as a library.  A LAMMPS executable with fix qmmm
included can be built, but will not be able to do a QM/MM simulation
on as such.  You must also build a QM code - currently only Quantum
ESPRESSO (QE) is supported - and create a new executable which links
LAMMPS and the QM code together.  Details are given in the
lib/qmmm/README file.  It is also recommended to read the instructions
for :doc:`linking with LAMMPS as a library <Build_link>` for
background information.  This requires compatible Quantum Espresso
and LAMMPS versions.  The current interface and makefiles have last
been verified to work in February 2020 with Quantum Espresso versions
6.3 to 6.5.

**CMake build**\ :

When using CMake, building a LAMMPS library is required and it is
recommended to build a shared library, since any libraries built from
the sources in the *lib* folder (including the essential libqmmm.a)
are not included in the static LAMMPS library and (currently) not
installed, while their code is included in the shared LAMMPS library.
Thus a typical command line to configure building LAMMPS for USER-QMMM
would be:

.. code-block:: bash

    cmake -C ../cmake/presets/minimal.cmake -D PKG_USER-QMMM=yes \
            -D BUILD_LIB=yes -DBUILD_SHARED_LIBS=yes ../cmake

After completing the LAMMPS build and also configuring and compiling
Quantum ESPRESSO with external library support (via "make couple"),
go back to the lib/qmmm folder and follow the instructions on the
README file to build the combined LAMMPS/QE QM/MM executable
(pwqmmm.x) in the lib/qmmm folder.  You need to make certain, that


**Traditional make**\ :

Before building LAMMPS, you must build the QMMM library in lib/qmmm.
You can do this manually if you prefer; follow the first two steps
explained in lib/qmmm/README.  You can also do it in one step from the
lammps/src dir, using a command like these, which simply invoke the
lib/qmmm/Install.py script with the specified args:


.. code-block:: bash

   make lib-qmmm                      # print help message
   make lib-qmmm args="-m serial"     # build with GNU Fortran compiler (settings as in "make serial")
   make lib-qmmm args="-m mpi"        # build with default MPI compiler (settings as in "make mpi")
   make lib-qmmm args="-m gfortran"   # build with GNU Fortran compiler

The build should produce two files: lib/qmmm/libqmmm.a and
lib/qmmm/Makefile.lammps.  The latter is copied from an existing
Makefile.lammps.\* and has settings needed to build LAMMPS with the
QMMM library (though typically the settings are just blank).  If
necessary, you can edit/create a new lib/qmmm/Makefile.machine file
for your system, which should define an EXTRAMAKE variable to specify
a corresponding Makefile.lammps.machine file.

You can then install QMMM package and build LAMMPS in the usual
manner.  After completing the LAMMPS build and compiling Quantum
ESPRESSO with external library support (via "make couple"), go back to
the lib/qmmm folder and follow the instructions in the README file to
build the combined LAMMPS/QE QM/MM executable (pwqmmm.x) in the
lib/qmmm folder.

----------


.. _user-quip:

USER-QUIP package
---------------------------------

To build with this package, you must download and build the QUIP
library.  It can be obtained from GitHub.  For support of GAP
potentials, additional files with specific licensing conditions need
to be downloaded and configured.  See step 1 and step 1.1 in the
lib/quip/README file for details on how to do this.

**CMake build**\ :


.. code-block:: bash

   -D QUIP_LIBRARY=path     # path to libquip.a (only needed if a custom location)

CMake will not download and build the QUIP library.  But once you have
done that, a CMake build of LAMMPS with "-D PKG\_USER-QUIP=yes" should
work.  Set QUIP\_LIBRARY if CMake cannot find the QUIP library.

**Traditional make**\ :

The download/build procedure for the QUIP library, described in
lib/quip/README file requires setting two environment variables,
QUIP\_ROOT and QUIP\_ARCH.  These are accessed by the
lib/quip/Makefile.lammps file which is used when you compile and link
LAMMPS with this package.  You should only need to edit
Makefile.lammps if the LAMMPS build can not use its settings to
successfully build on your system.


----------


.. _user-scafacos:

USER-SCAFACOS package
-----------------------------------------

To build with this package, you must download and build the `ScaFaCoS Coulomb solver library <scafacos-home_>`_

.. _scafacos-home: http://www.scafacos.de



**CMake build**\ :


.. code-block:: bash

   -D DOWNLOAD_SCAFACOS=value    # download ScaFaCoS for build, value = no (default) or yes
   -D SCAFACOS_LIBRARY=path      # ScaFaCos library file (only needed if at custom location)
   -D SCAFACOS_INCLUDE_DIR=path  # ScaFaCoS include directory (only needed if at custom location)

If DOWNLOAD\_SCAFACOS is set, the ScaFaCoS library will be downloaded
and built inside the CMake build directory.  If the ScaFaCoS library
is already on your system (in a location CMake cannot find it),
SCAFACOS\_LIBRARY is the filename (plus path) of the ScaFaCoS library
file, not the directory the library file is in.  SCAFACOS\_INCLUDE\_DIR
is the directory the ScaFaCoS include file is in.

**Traditional make**\ :

You can download and build the ScaFaCoS library manually if you
prefer; follow the instructions in lib/scafacos/README.  You can also
do it in one step from the lammps/src dir, using a command like these,
which simply invoke the lib/scafacos/Install.py script with the
specified args:

make lib-scafacos                         # print help message
make lib-scafacos args="-b"               # download and build in lib/scafacos/scafacos-<version>
make lib-scafacos args="-p $HOME/scafacos  # use existing ScaFaCoS installation in $HOME/scafacos

Note that 2 symbolic (soft) links, "includelink" and "liblink", are
created in lib/scafacos to point to the ScaFaCoS src dir.  When LAMMPS
builds in src it will use these links.  You should not need to edit
the lib/scafacos/Makefile.lammps file.


----------


.. _user-smd:

USER-SMD package
-------------------------------

To build with this package, you must download the Eigen3 library.
Eigen3 is a template library, so you do not need to build it.

**CMake build**\ :


.. code-block:: bash

   -D DOWNLOAD_EIGEN3            # download Eigen3, value = no (default) or yes
   -D EIGEN3_INCLUDE_DIR=path    # path to Eigen library (only needed if a custom location)

If DOWNLOAD\_EIGEN3 is set, the Eigen3 library will be downloaded and
inside the CMake build directory.  If the Eigen3 library is already on
your system (in a location CMake cannot find it), EIGEN3\_INCLUDE\_DIR
is the directory the Eigen3++ include file is in.

**Traditional make**\ :

You can download the Eigen3 library manually if you prefer; follow the
instructions in lib/smd/README.  You can also do it in one step from
the lammps/src dir, using a command like these, which simply invoke
the lib/smd/Install.py script with the specified args:


.. code-block:: bash

   make lib-smd                         # print help message
   make lib-smd args="-b"               # download to lib/smd/eigen3
   make lib-smd args="-p /usr/include/eigen3"    # use existing Eigen installation in /usr/include/eigen3

Note that a symbolic (soft) link named "includelink" is created in
lib/smd to point to the Eigen dir.  When LAMMPS builds it will use
this link.  You should not need to edit the lib/smd/Makefile.lammps
file.


----------


.. _user-vtk:

USER-VTK package
-------------------------------

To build with this package you must have the VTK library installed on
your system.

**CMake build**\ :

No additional settings are needed besides "-D PKG\_USER-VTK=yes".

This should auto-detect the VTK library if it is installed on your
system at standard locations.  Several advanced VTK options exist if
you need to specify where it was installed.  Use the ccmake (terminal
window) or cmake-gui (graphical) tools to see these options and set
them interactively from their user interfaces.

**Traditional make**\ :

The lib/vtk/Makefile.lammps file has settings for accessing VTK files
and its library, which LAMMPS needs to build with this package.  If
the settings are not valid for your system, check if one of the other
lib/vtk/Makefile.lammps.\* files is compatible and copy it to
Makefile.lammps.  If none of the provided files work, you will need to
edit the Makefile.lammps file.  See lib/vtk/README for details.
