LAMMPS portability and compatibility
------------------------------------

The primary form of distributing LAMMPS is through highly portable
source code.  But also several ways of obtaining LAMMPS as :doc:`precompiled
packages or through automated build mechanisms <Install>` exist.  Most
of LAMMPS is written in C++, some support tools are written in Fortran
or Python or MATLAB.


Programming language standards
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 10Sep2025

The C++ code in LAMMPS currently requires a compiler that is compatible
with the C++17 standard. The Kokkos library used for the KOKKOS
package currently also requires at least C++17.  If your compilers are not
compatible and you cannot upgrade, please use LAMMPS version 22 July 2025,
which requires only C++11 as the minimum C++ standard.

Most of the Python code in LAMMPS is written to be compatible with Python
3.6 and later.

.. deprecated:: 2Apr2025

Python 2.x is no longer supported and trying to use it, e.g. for the
LAMMPS Python module should result in an error.  If you come across
some part of the LAMMPS distribution that is not (yet) compatible with
Python 3, please notify the LAMMPS developers.

Build systems
^^^^^^^^^^^^^

LAMMPS can be compiled from source code using the cross-platform CMake
system.  CMake must be at least version 3.20.  Alternatively, using a
(traditional) build system based on shell scripts, a few shell utilities
(grep, sed, cat, tr) and the GNU make program.  This requires running
within a Bourne shell (``/bin/sh`` or ``/bin/bash``).

.. versionchanged:: 10Sep2025

The traditional GNU make based build system no longer supports all
packages.  Details can be found in the :doc:`package specific build
instructions <Build_extras>`.

Operating systems
^^^^^^^^^^^^^^^^^

The primary development platform for LAMMPS is Linux.  Thus, the chances
for LAMMPS to compile without problems are the best on Linux machines.
Also, compilation and correct execution on macOS and Windows (using
Microsoft Visual C++) is checked automatically for the largest part of
the source code.  Some (optional) features are not compatible with all
operating systems, either through limitations of the corresponding
LAMMPS source code or through incompatibilities or build system
limitations of required external libraries or packages.

Executables for Windows may be created either natively using Cygwin,
MinGW, Intel, Clang, or Microsoft Visual C++ compilers, or with a Linux
to Windows MinGW cross-compiler.  Native compilation is supported using
Microsoft Visual Studio or a terminal window (using the CMake build
system).

Executables for macOS may be created either using Xcode or GNU compilers
installed with Homebrew.  In the latter case, building of LAMMPS through
Homebrew instead of a manual compile is also possible.

Additionally, FreeBSD and Solaris have been tested successfully to
run LAMMPS and produce results consistent with those on Linux.

Compilers
^^^^^^^^^

The most commonly used compilers are the GNU compilers, but also Clang
and the Intel compilers have been successfully used on Linux, macOS, and
Windows.  Also, the Nvidia HPC SDK (formerly PGI compilers) will compile
LAMMPS (tested on Linux).

.. versionchanged:: 10Sep2025

The GNU compilers *before* version 9.3 have known problems with supporting
C++17 and thus are **not** recommended to build LAMMPS.

CPU architectures
^^^^^^^^^^^^^^^^^

The primary CPU architecture for running LAMMPS is 64-bit x86, but also
64-bit ARM is currently regularly tested.  Further architectures are tested
by Linux distributions that bundle LAMMPS.

Portability compliance
^^^^^^^^^^^^^^^^^^^^^^

Only a subset of the LAMMPS source code is *fully* compliant to *all*
of the above mentioned standards.  This is rather typical for projects
like LAMMPS that largely depend on contributions from the user community.
Not all contributors are trained as programmers and not all of them have
access to multiple platforms for testing.  As part of the continuous
integration process, however, all contributions are automatically tested
to compile, link, and pass some runtime tests on a selection of Linux
flavors, macOS, and Windows, and on Linux with different compilers.
Thus portability issues are often found before a pull request is merged.
Other platforms may be checked occasionally or when portability bugs are
reported.
