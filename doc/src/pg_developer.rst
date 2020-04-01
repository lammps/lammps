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
class, for example the :ref:`Run <lammps_ns_run>` class which contains
the code invoked by the :doc:`run <run>` command in a LAMMPS input script.
As this example illustrates, source file and class names often have a
one-to-one correspondence with a command used in a LAMMPS input script.
Some source files and classes do not have a corresponding input script
command, e.g. ``src/force.cpp`` and the Force class.  They are discussed
in the next section.
