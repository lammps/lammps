Notes for building LAMMPS on Windows
====================================

* :ref:`General remarks <generic>`
* :ref:`Running Linux on Windows <linux>`
* :ref:`Using GNU GCC ported to Windows <gnu>`
* :ref:`Using a cross-compiler <cross>`

----------

.. _generic:

General remarks
-----------------------------

LAMMPS is developed and tested primarily on Linux machines.  The vast
majority of HPC clusters and supercomputers today runs on Linux as well.
While portability to other platforms is desired, it is not always achieved.
The LAMMPS developers strongly rely on LAMMPS users giving feedback and
providing assistance in resolving portability issues. This is particularly
true for compiling LAMMPS on Windows, since this platform has significant
differences with some low-level functionality.

.. _linux:

Running Linux on Windows
------------------------------------

So before trying to build LAMMPS on Windows, please consider if using
the pre-compiled Windows binary packages are sufficient for your needs
(as an aside, those packages themselves are build on a Linux machine
using cross-compilers).  If it is necessary for you to compile LAMMPS
on a Windows machine (e.g. because it is your main desktop), please also
consider using a virtual machine software and compile and run LAMMPS in
a Linux virtual machine, or - if you have a recently updated Windows 10
installation - consider using the Windows subsystem for Linux.  This
optional Windows feature allows you to run the bash shell from Ubuntu
from within Windows and from there on, you can pretty much use that
shell like you are running on an Ubuntu Linux machine (e.g. installing
software via apt-get and more). For more details on that, please
see :doc:`this tutorial <Howto_bash>`

.. _gnu:

Using a GNU GCC ported to Windows
-----------------------------------------

One option for compiling LAMMPS on Windows natively, that has been known
to work in the past is to install a bash shell, unix shell utilities,
perl, GNU make, and a GNU compiler ported to Windows. The Cygwin package
provides a unix/linux interface to low-level Windows functions, so LAMMPS
can be compiled on Windows. The necessary (minor) modifications to LAMMPS
are included, but may not always up-to-date for recently added functionality
and the corresponding new code. A machine makefile for using cygwin for
the old build system is provided. Using CMake for this mode of compilation
is untested and not likely to work.

When compiling for Windows do **not** set the -DLAMMPS\_MEMALIGN define
in the LMP\_INC makefile variable and add -lwsock32 -lpsapi to the linker
flags in LIB makefile variable. Try adding -static-libgcc or -static or
both to the linker flags when your resulting LAMMPS Windows executable
complains about missing .dll files. The CMake configuration should set
this up automatically, but is untested.

In case of problems, you are recommended to contact somebody with
experience in using cygwin.  If you do come across portability problems
requiring changes to the LAMMPS source code, or figure out corrections
yourself, please report them on the lammps-users mailing list, or file
them as an issue or pull request on the LAMMPS GitHub project.

.. _cross:

Using a cross-compiler
----------------------------------

If you need to provide custom LAMMPS binaries for Windows, but do not
need to do the compilation on Windows, please consider using a Linux
to Windows cross-compiler. This is how currently the Windows binary
packages are created by the LAMMPS developers. Because of that, this is
probably the currently best tested and supported way to build LAMMPS
executables for Windows.  There are makefiles provided for the
traditional build system, but CMake has also been successfully tested
using the mingw32-cmake and mingw64-cmake wrappers that are bundled
with the cross-compiler environment on Fedora machines. A CMake preset
selecting all packages compatible with this cross-compilation build
is provided. You will likely need to disable the GPU package unless you
download and install the contents of the pre-compiled `OpenCL ICD loader library <https://download.lammps.org/thirdparty/opencl-win-devel.tar.gz>`_
into your MinGW64 cross-compiler environment. The cross-compilation
currently will only produce non-MPI serial binaries.

Please keep in mind, though, that this only applies to **compiling** LAMMPS.
Whether the resulting binaries do work correctly is not tested by the
LAMMPS developers.  We instead rely on the feedback of the users
of these pre-compiled LAMMPS packages for Windows.  We will try to resolve
issues to the best of our abilities if we become aware of them. However
this is subject to time constraints and focus on HPC platforms.

.. _native:

Native Visual C++ support
--------------------------------------

Support for the Visual C++ compilers is currently not available. The
CMake build system is capable of creating suitable a Visual Studio
style build environment, but the LAMMPS code itself is not fully ported
to support Visual C++. Volunteers to take on this task are welcome.
