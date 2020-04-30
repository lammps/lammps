LAMMPS Developer Guide
**********************

This section describes the internal structure and basic algorithms
of the LAMMPS code. This is a work in progress and additional
information will be added incrementally depending on availability
of time and requests from the LAMMPS user community.


LAMMPS source files
===================

The source files of the LAMMPS code are distributed across two
directories of the distribution.  The core of the code is located in the
``src`` folder and its sub-directories. Almost all of those are C++ files
(implementation files have a ``.cpp`` extension and and headers a
``.h``).  A sizable number of these files are in the ``src`` directory
itself, but there are plenty of :doc:`packages <Packages>`, which can be
included or excluded when LAMMPS is built.  See the :doc:`Include
packages in build <Build_package>` section of the manual for more
information about that part of the build process.  LAMMPS currently
supports building with :doc:`conventional makefiles <Build_make>` and
through :doc:`CMake <Build_cmake>` which differ in how packages are
enabled or disabled for a LAMMPS binary.  The source files for each
package are in all-uppercase sub-directories of the ``src`` folder, for
example ``src/MOLECULE`` or ``src/USER-MISC``.  The ``src/STUBS``
sub-directory is not a package but contains a dummy MPI library, that is
used when building a serial version of the code. the ``src/MAKE``
directory contains makefiles with settings and flags for a variety of
configuration and machines for the build process with traditional
makefiles.

The ``lib`` directory contains the source code for several supporting
libraries or files with configuration settings to use globally installed
libraries, that are required by some of the optional packages.
Each sub-directory, like ``lib/poems`` or ``lib/gpu``, contains the
source files, some of which are in different languages such as Fortran
or CUDA. These libraries are linked to during a LAMMPS build, if the
corresponding package is installed.

LAMMPS C++ source files almost always come in pairs, such as
``src/run.cpp`` and ``src/run.h``.  The pair of files defines a C++
class, for example the :cpp:class:`LAMMPS_NS::Run` class which contains
the code invoked by the :doc:`run <run>` command in a LAMMPS input script.
As this example illustrates, source file and class names often have a
one-to-one correspondence with a command used in a LAMMPS input script.
Some source files and classes do not have a corresponding input script
command, e.g. ``src/force.cpp`` and the :cpp:class:`LAMMPS_NS::Force`
class.  They are discussed in the next section.

Overview of LAMMPS class topology
=================================

Though LAMMPS has a lot of source files and classes, its class topology
is relative flat, as outlined in the :ref:`class-topology` figure.  Each
name refers to a class and has a pair of associated source files in the
``src`` folder, for example the class :cpp:class:`LAMMPS_NS::Memory`
corresponds to the files ``memory.cpp`` and ``memory.h``, or the class
:cpp:class:`LAMMPS_NS::AtomVec` corresponds to the files
``atom_vec.cpp`` and ``atom_vec.h``.  Full lines in the figure represent
compositing: that is the class to the left holds a pointer to an
instance of the class to the right.  Dashed lines instead represent
inheritance: the class to the right is derived from the class on the
left. Classes with a red boundary are not instantiated directly, but
they represent the base classes for "styles".  Those "styles" make up
the bulk of the LAMMPS code and only a few typical examples are included
in the figure for demonstration purposes.

.. _class-topology:
.. figure:: JPG/lammps-classes.png

   LAMMPS class topology

   This figure shows some of the relations of the base classes of the
   LAMMPS simulation package.  Full lines indicate that a class holds an
   instance of the class it is pointing to; dashed lines point to
   derived classes that are given as examples of what classes may be
   instantiated during a LAMMPS run based on the input commands and
   accessed through the API define by their respective base classes.  At
   the core is the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class, which
   holds pointers to class instances with specific purposes.  Those may
   hold instances of other classes, sometimes directly, or only
   temporarily, sometimes as derived classes or derived classes or
   derived classes, which may also hold instances of other classes.

The :cpp:class:`LAMMPS_NS::LAMMPS` class is the topmost class and
represents what is referred to an "instance" of LAMMPS.  It is a
composite holding references to instances of other core classes
providing the core functionality of the MD engine in LAMMPS and through
them abstractions of the required operations.  The constructor of the
LAMMPS class will instantiate those instances, process the command line
flags, initialize MPI (if not already done) and set up file pointers for
input and output. The destructor will shut everything down and free all
associated memory.  Thus code for the standalone LAMMPS executable in
``main.cpp`` simply initializes MPI, instantiates a single instance of
LAMMPS, and passes it the command line flags and input script. It
deletes the LAMMPS instance after the method reading the input returns
and shuts down the MPI environment before it exits the executable.

The :cpp:class:`LAMMPS_NS::Pointers` is not shown in the
:ref:`class-topology` figure, it holds references to members of the
`LAMMPS_NS::LAMMPS`, so that all classes derived from
:cpp:class:`LAMMPS_NS::Pointers` have direct access to those reference.
From the class topology all classes with blue boundary are referenced in
this class and all classes in the second and third columns, that are not
listed as derived classes are instead derived from
:cpp:class:`LAMMPS_NS::Pointers`.

Since all storage is encapsulated, the LAMMPS class can also be
instantiated multiple times by a calling code, and that can be either
simultaneously or consecutively.  When running in parallel with MPI,
care has to be taken, that suitable communicators are used to not
create conflicts between different instances.

The LAMMPS class currently holds instances of 19 classes representing
different core functionalities 
There are a handful of virtual parent classes in LAMMPS that define
what LAMMPS calls ``styles``.  They are shaded red in Fig
\ref{fig:classes}.  Each of these are parents of a number of child
classes that implement the interface defined by the parent class.

-----------------

.. include:: pg_utils.rst
