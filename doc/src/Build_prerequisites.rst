Prerequisites
-------------

Which software you need to compile and use LAMMPS strongly depends on
which :doc:`features and settings <Build_settings>` and which
:doc:`optional packages <Packages>` you are trying to include.
Common to all is that you need a C++ and C compiler, where the C++
compiler has to support at least the C++11 standard (note that some
compilers require command-line flag to activate C++11 support).
Furthermore, if you are building with CMake, you need at least CMake
version 3.20 and a compatible build tool (make or ninja-build); if you
are building the the legacy GNU make based build system you need GNU
make (other make variants are not going to work since the build system
uses features unique to GNU make) and a Unix-like build environment with
a Bourne shell, and shell tools like "sed", "grep", "touch", "test",
"tr", "cp", "mv", "rm", "ln", "diff" and so on. Parts of LAMMPS
interface with or use Python version 3.6 or later.

The LAMMPS developers aim to keep LAMMPS very portable and usable -
at least in parts - on most operating systems commonly used for
running MD simulations.  Please see the :doc:`section on portablility
<Intro_portability>` for more details.
