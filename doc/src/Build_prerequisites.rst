Prerequisites
-------------

Which software you need to compile and use LAMMPS strongly depends on
which :doc:`features and settings <Build_settings>` and which
:doc:`optional packages <Packages>` you are trying to include.
Common to all is that you need a C++ and C compiler, where the C++
compiler has to support at least the C++17 standard (note that some
compilers require a command-line flag to activate C++17 support).

CMake build system
==================

If you are building with CMake, you need at least CMake version 3.20 and
a compatible build tool (e.g. GNU make or ninja-build on Linux).  The
CMake scripting includes tests for required software and will
auto-detect and auto-enable available tools and libraries for optional
features (of course those can be disabled, if desired).  If required
software is missing, CMake will try to download the compile the missing
parts automatically, if possible or stop with an error.

GNU build system
================

If you are building LAMMPS with the legacy GNU make based build system
you need GNU make (other make variants are not going to work since the
build system uses features unique to GNU make) and a Unix-like build
environment with a Bourne shell, and shell tools like "sed", "grep",
"touch", "test", "tr", "cp", "mv", "rm", "ln", "diff" and so on. Parts
of LAMMPS interface with or use Python version 3.6 or later.

Compiler and OS compatibility
=============================

The LAMMPS developers aim to keep LAMMPS very portable and usable -
at least in parts - on most operating systems commonly used for
running MD simulations.  Please see the :doc:`section on portablility
<Intro_portability>` for more details.

.. admonition:: Warning: LLVM based Intel Compilers
   :class: warning

   In recent years, Intel switched their compilers from a proprietary
   code base to compilers based on LLVM (which is also the foundation of
   the Clang compilers).  It has since stopped updating their legacy
   compilers.  Some versions of those LLVM compilers creates incorrect
   binary code that results in incorrect behavior of LAMMPS; this is
   particularly common when using high optimization levels and
   vectorization.  There have been several cases where people reported
   bugs in LAMMPS which turned out to be caused by bugs in the Intel
   LLVM compilers.

   Unfortunately there is no simple way to detect whether a binary is
   working correctly outside of running the unit and regression tests,
   but those do not cover all of LAMMPS and and would be reliable only
   for no or moderate optimization anyway.  For most of LAMMPS there is
   not much of a benefit (if any) to use the Intel compilers over the
   GCC or Clang compilers, except for the INTEL package (which *can* be
   compiled with other compilers, but most vectorization directives are
   inactive for those) or KOKKOS with SYCL.

   Our recommendation is thus to compile with GCC or Clang unless you
   are using features that *require* the Intel compilers.  In that case,
   it is recommended to have a backup GCC/Clang compiled binary
   available to validate correctness of the simulation.  If you have
   access to multiple Intel compiler versions, please try the latest
   version first and avoid "20##.0" versions which seem to have
   problems most often.
