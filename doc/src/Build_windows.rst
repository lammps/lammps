Notes for building LAMMPS on Windows
------------------------------------

* :ref:`General remarks <generic>`
* :ref:`Running Linux on Windows <linux>`
* :ref:`Using GNU GCC ported to Windows <gnu>`
* :ref:`Using Visual Studio <msvc>`
* :ref:`Using Intel oneAPI compilers and libraries <oneapi>`
* :ref:`Using a cross-compiler <cross>`

----------

.. _generic:

General remarks
^^^^^^^^^^^^^^^

LAMMPS is developed and tested primarily on Linux machines.  The vast
majority of HPC clusters and supercomputers today run on Linux as well.
While portability to other platforms is desired, it is not always
achieved.  That is sometimes due to non-portable code in LAMMPS itself,
but more often due to portability limitations of external libraries and
tools required to build a specific feature or package.  The LAMMPS
developers are dependent on LAMMPS users giving feedback and providing
assistance in resolving portability issues.  This is particularly true
for compiling LAMMPS on Windows, since this platform has significant
differences in some low-level functionality.  As of LAMMPS version 14
December 2021, large parts of LAMMPS can be compiled natively with the
Microsoft Visual C++ Compilers.  As of LAMMPS version 31 May 2022, also
the Intel oneAPI compilers can compile large parts of LAMMPS natively on
Windows.  This is mostly facilitated by using the
:doc:`Developer_platform` in the ``platform`` namespace and CMake.

Before trying to build LAMMPS on Windows yourself, please consider the
`pre-compiled Windows installer packages <https://packages.lammps.org/windows.html>`_
and see if they are sufficient for your needs.

.. _linux:

Running Linux on Windows
^^^^^^^^^^^^^^^^^^^^^^^^

If it is necessary for you to compile LAMMPS on a Windows machine
(e.g. because it is your main desktop), please also consider using a
virtual machine software and compile and run LAMMPS in a Linux virtual
machine, or - if you have a sufficiently up-to-date Windows 10 or
Windows 11 installation - consider using the Windows subsystem for
Linux.  This optional Windows feature allows you to run the bash shell
of a Linux system (Ubuntu by default) from within Windows and from there
on, you can pretty much use that shell like you are running on a regular
Ubuntu Linux machine (e.g. installing software via apt-get and more).
For more details on that, please see :doc:`this tutorial <Howto_wsl>`.

.. _gnu:

Using a GNU GCC ported to Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One option for compiling LAMMPS on Windows natively is to install a Bash
shell, Unix shell utilities, Perl, Python, GNU make, and a GNU compiler
ported to Windows.  The Cygwin package provides a unix/linux interface
to low-level Windows functions, so LAMMPS can be compiled on Windows.
The necessary (minor) modifications to LAMMPS are included, but may not
always up-to-date for recently added functionality and the corresponding
new code.  A machine makefile for using cygwin for the old build system
is provided.  Using CMake for this mode of compilation is untested and
not likely to work.

When compiling for Windows do **not** set the ``-DLAMMPS_MEMALIGN``
define in the LMP_INC makefile variable and add ``-lwsock32 -lpsapi`` to
the linker flags in LIB makefile variable. Try adding ``-static-libgcc``
or ``-static`` or both to the linker flags when your resulting LAMMPS
Windows executable complains about missing .dll files. The CMake
configuration should set this up automatically, but is untested.

In case of problems, you are recommended to contact somebody with
experience in using Cygwin.  If you do come across portability problems
requiring changes to the LAMMPS source code, or figure out corrections
yourself, please report them on the
`LAMMPS forum at MatSci <https://matsci.org/c/lammps/lammps-development/>`_,
or file them as an issue or pull request on the LAMMPS GitHub project.

.. _msvc:

Using Microsoft Visual Studio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following the integration of the :doc:`platform namespace
<Developer_platform>` into the LAMMPS code base, portability of LAMMPS
for native compilation on Windows using Visual Studio has been
significantly improved.  This has been tested with Visual Studio 2019
(aka version 16) and Visual Studio 2022 (aka version 17).  We strongly
recommend using Visual Studio 2022 version 17.1 or later.  Not all
features and packages in LAMMPS are currently supported out of the box,
but a preset ``cmake/presets/windows.cmake`` is provided that contains
the packages that have been compiled successfully so far.  You **must**
use the CMake based build procedure, since there is no support for GNU
make or the Unix shell utilities required for the GNU make build
procedure.

It is possible to use both the integrated CMake support of the Visual
Studio IDE or use an external CMake installation (e.g. downloaded from
cmake.org) to create build files and compile LAMMPS from the command line.

Compilation via command line and unit tests are checked automatically
for the LAMMPS development branch through
`GitHub Actions <https://github.com/lammps/lammps/actions/workflows/compile-msvc.yml>`_.

.. note::

   Versions of Visual Studio before version 17.1 may scan the entire
   LAMMPS source tree and likely miss the correct master
   ``CMakeLists.txt`` and get confused since there are multiple files
   of that name in different folders but none in top level folder.

Please note, that for either approach CMake will create a so-called
:ref:`"multi-configuration" build environment <cmake_multiconfig>`, and
the command lines for building and testing LAMMPS must be adjusted
accordingly.

The LAMMPS cmake folder contains a ``CMakeSettings.json`` file with
build configurations for MSVC compilers and the MS provided Clang
compiler package in Debug and Release mode.

To support running in parallel you can compile with OpenMP enabled using
the OPENMP package or install Microsoft MPI (including the SDK) and compile
LAMMPS with MPI enabled.

.. note::

   This is work in progress and you should contact the LAMMPS developers
   via GitHub or the `LAMMPS forum at MatSci <https://matsci.org/c/lammps/lammps-development/>`_,
   if you have questions or LAMMPS specific problems.

.. _oneapi:

Using Intel oneAPI Compilers and Libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 31May2022

After installing the `Intel oneAPI
<https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html>`_
base toolkit and the HPC toolkit, it is also possible to compile large
parts of LAMMPS natively on Windows using Intel compilers.  The HPC
toolkit provides two sets of C/C++ and Fortran compilers: the so-called
"classic" compilers (``icl.exe`` and ``ifort.exe``) and newer, LLVM
based compilers (``icx.exe`` and ``ifx.exe``).  In addition to the
compilers and their dependent modules, also the thread building blocks
(TBB) and the math kernel library (MKL) need to be installed.  Two
presets (``cmake/presets/windows-intel-llvm.cmake`` and
``cmake/presets/windows-intel-classic.cmake``) are provided for
selecting the LLVM based or classic compilers, respectively. The preset
``cmake/presets/windows.cmake`` enables compatible packages that are not
dependent on additional features or libraries.  You **must** use the
CMake based build procedure and use Ninja as build tool.  For compiling
from the command prompt, thus both `CMake <https://cmake.org>`_ and
`Ninja-build <https://ninja-build.org>`_ binaries must be installed.  It
is also possible to use Visual Studio, if it is started (``devenv.exe``)
from a command prompt that has the Intel oneAPI compilers enabled.  The
Visual Studio settings file in the ``cmake`` folder contains
configurations for both compiler variants in debug and release settings.
Those will use the CMake and Ninja binaries bundled with Visual Studio,
thus a separate installation is not required.

.. admonition:: Known Limitations
   :class: note

   In addition to portability issues with several packages and external
   libraries, the classic Intel compilers are currently not able to
   compile the googletest libraries and thus enabling the ``-DENABLE_TESTING``
   option will result in compilation failure.  The LLVM based compilers
   are compatible.

.. note::

   This is work in progress and you should contact the LAMMPS developers
   via GitHub or the `LAMMPS forum at MatSci <https://matsci.org/c/lammps/lammps-development/>`_,
   if you have questions or LAMMPS specific problems.


.. _cross:

Using a cross-compiler
^^^^^^^^^^^^^^^^^^^^^^

If you need to provide custom LAMMPS binaries for Windows, but do not
need to do the compilation on Windows, please consider using a Linux to
Windows cross-compiler.  This is how currently the Windows binary
packages are created by the LAMMPS developers.  Because of that, this is
probably the currently best tested and supported way to build LAMMPS
executables for Windows.  A CMake preset selecting all packages
compatible with this cross-compilation build is provided.  The GPU
package can only be compiled with OpenCL support.  To compile with MPI
support, a pre-compiled library and the corresponding header files are
required.  When building with CMake the matching package will be
downloaded automatically, but MPI support has to be explicitly enabled
with ``-DBUILD_MPI=on``.

Please keep in mind, though, that this only applies to **compiling** LAMMPS.
Whether the resulting binaries do work correctly is rarely tested by the
LAMMPS developers.  We instead rely on the feedback of the users
of these pre-compiled LAMMPS packages for Windows.  We will try to resolve
issues to the best of our abilities if we become aware of them. However
this is subject to time constraints and focus on HPC platforms.
