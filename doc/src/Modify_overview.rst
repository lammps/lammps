Overview
========

The best way to add a new feature to LAMMPS is to find a similar
feature and look at the corresponding source and header files to figure
out what it does.  You will need some knowledge of C++ to be able to
understand the high-level structure of LAMMPS and its class
organization, but functions (class methods) that do actual
computations are written in vanilla C-style code and operate on simple
C-style data structures (vectors and arrays).

Most of the new features described on the :doc:`Modify <Modify>` doc
page require you to write a new C++ derived class (except for exceptions
described below, where you can make small edits to existing files).
Creating a new class requires 2 files, a source code file (\*.cpp) and a
header file (\*.h).  The derived class must provide certain methods to
work as a new option.  Depending on how different your new feature is
compared to existing features, you can either derive from the base class
itself, or from a derived class that already exists.  Enabling LAMMPS to
invoke the new class is as simple as putting the two source files in the
src directory and re-building LAMMPS.

The advantage of C++ and its object-orientation is that all the code
and variables needed to define the new feature are in the 2 files you
write, and thus should not make the rest of LAMMPS more complex or
cause side-effect bugs.

Here is a concrete example.  Suppose you write 2 files pair_foo.cpp
and pair_foo.h that define a new class PairFoo that computes pairwise
potentials described in the classic 1997 :ref:`paper <Foo>` by Foo, et al.
If you wish to invoke those potentials in a LAMMPS input script with a
command like

.. code-block:: LAMMPS

   pair_style foo 0.1 3.5

then your pair_foo.h file should be structured as follows:

.. code-block:: c++

   #ifdef PAIR_CLASS
   PairStyle(foo,PairFoo)
   #else
   ...
   (class definition for PairFoo)
   ...
   #endif

where "foo" is the style keyword in the pair_style command, and
PairFoo is the class name defined in your pair_foo.cpp and pair_foo.h
files.

When you re-build LAMMPS, your new pairwise potential becomes part of
the executable and can be invoked with a pair_style command like the
example above.  Arguments like 0.1 and 3.5 can be defined and
processed by your new class.

.. note:

  With the traditional make process, simply adding the new files to the
  src folder and compiling LAMMPS again for the desired configuration
  with "make machine" is sufficient.  When using CMake, you need to
  re-run CMake with "cmake ." in the build folder to have it recognize
  the added files and include them into the build system.

As illustrated by this example pair style, many kinds of options are
referred to in the LAMMPS documentation as the "style" of a particular
command.

The :doc:`Modify page <Modify>` lists all the common styles in LAMMPS,
and discusses the header file for the base class that these styles are
derived from.  Public variables in that file are ones used and set by
the derived classes which are also used by the base class.  Sometimes
they are also used by the rest of LAMMPS.  Pure functions, which means
functions declared as virtual in the base class header file which are
also set to 0, are functions you **must** implement in your new derived
class to give it the functionality LAMMPS expects. Virtual functions
that are not set to 0 are functions you may override or not.  Those
are usually defined with an empty function body.

Additionally, new output options can be added directly to the
thermo.cpp, dump_custom.cpp, and variable.cpp files.  These are also
listed on the :doc:`Modify page <Modify>`.

Here are additional guidelines for modifying LAMMPS and adding new
functionality:

* Think about whether what you want to do would be better as a pre- or
  post-processing step.  Many computations are more easily and more
  quickly done that way.
* Do not try to do anything within the timestepping of a run that is not
  parallel.  For example do not accumulate a bunch of data on a single
  processor and analyze it.  You run the risk of seriously degrading
  the parallel efficiency this way.
* If your new feature reads arguments or writes output, make sure you
  follow the unit conventions discussed by the :doc:`units <units>`
  command.

----------

.. _Foo:

**(Foo)** Foo, Morefoo, and Maxfoo, J of Classic Potentials, 75, 345 (1997).
