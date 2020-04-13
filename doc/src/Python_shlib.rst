Build LAMMPS as a shared library
================================

Build LAMMPS as a shared library using make
-------------------------------------------

Instructions on how to build LAMMPS as a shared library are given on
the :doc:`Build_basics <Build_basics>` doc page.  A shared library is
one that is dynamically loadable, which is what Python requires to
wrap LAMMPS.  On Linux this is a library file that ends in ".so", not
".a".

From the src directory, type

.. code-block:: bash

   make foo mode=shlib

where foo is the machine target name, such as mpi or serial.
This should create the file liblammps_foo.so in the src directory, as
well as a soft link liblammps.so, which is what the Python wrapper will
load by default.  Note that if you are building multiple machine
versions of the shared library, the soft link is always set to the
most recently built version.

.. note::

   If you are building LAMMPS with an MPI or FFT library or other
   auxiliary libraries (used by various packages), then all of these
   extra libraries must also be shared libraries.  If the LAMMPS
   shared-library build fails with an error complaining about this, see
   the :doc:`Build_basics <Build_basics>` doc page.

Build LAMMPS as a shared library using CMake
--------------------------------------------

When using CMake the following two options are necessary to generate the LAMMPS
shared library:

.. code-block:: bash

   -D BUILD_LIB=on            # enable building LAMMPS as a library
   -D BUILD_SHARED_LIBS=on    # enable building of LAMMPS shared library (both options are needed!)

What this does is create a liblammps.so which contains the majority of LAMMPS
code. The generated lmp binary also dynamically links to this library. This
means that either this liblammps.so file has to be in the same directory, a system
library path (e.g. /usr/lib64/) or in the LD_LIBRARY_PATH.

If you want to use the shared library with Python the recommended way is to create a virtualenv and use it as
CMAKE_INSTALL_PREFIX.

.. code-block:: bash

   # create virtualenv
   virtualenv --python=$(which python3) myenv3
   source myenv3/bin/activate

   # build library
   mkdir build
   cd build
   cmake -D PKG_PYTHON=on -D BUILD_LIB=on -D BUILD_SHARED_LIBS=on -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../cmake
   make -j 4

   # install into prefix
   make install

This will also install the Python module into your virtualenv. Since virtualenv
does not change your LD_LIBRARY_PATH, you still need to add its lib64 folder to
it, which contains the installed liblammps.so.

.. code-block:: bash

   export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib64:$LD_LIBRARY_PATH

Starting Python outside (!) of your build directory, but with the virtualenv
enabled and with the LD_LIBRARY_PATH set gives you access to LAMMPS via Python.
