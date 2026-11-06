Build LAMMPS with make
======================

.. versionremoved:: TBD

The legacy build system to configure and build LAMMPS based on makefiles
for GNU make is no longer available.  Please use :doc:`CMake
<Build_cmake>` instead.  There is also a detailed discussion of
:doc:`how to use CMake with LAMMPS <Howto_cmake>` and the `CMake online
documentation <https://cmake.org/cmake/help/v3.27/>`_.  For the
impatient, there is also a quick transition guide below.  Contact the
LAMMPS developers at developers@lammps.org if you have any concerns
or questions about this.

Make to CMake quick-start guide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 17 38 45

   * - Action
     - Traditional ``make`` command in ``src/``
     - CMake equivalent in main folder
   * - Initial Setup
     - None or Edit ``Makefile.<machine>``
     - ``cmake -S cmake -B build``
   * - Build
     - ``make <machine>``
     - ``cmake --build build``
   * - Parallel build
     - ``make -j 8 <machine>``
     - ``cmake --build build --parallel 8``
   * - Enable package
     - ``make yes-<package name>``
     - ``cmake -S cmake -B build -D PKG_<PACKAGE_NAME>=on``
   * - Disable package
     - ``make no-<package name>``
     - ``cmake -S cmake -B build -D PKG_<PACKAGE_NAME>=off``
   * - Load a preset configuration
     - Choose a ``Makefile.<machine>`` from ``MAKE`` folders
     - ``cmake -S cmake -B build -C cmake/presets/<preset>.cmake``
   * - List enabled packages
     - ``make ps`` (or ``make package-status``)
     - Check CMake output or use ``ccmake build`` or ``lmp -h``
   * - Change configuration settings
     - Copy ``Makefile.<machine>`` to MAKE/MINE and edit
     - Rerun CMake with changed options
   * - Clean build files
     - ``make clean-<machine>``
     - ``cmake --build build --target clean``
   * - Full reset
     - ``make no-all clean-all purge``
     - ``rm -rf build``
   * - Enable or disable MPI
     - Use ``make mpi`` or ``make serial``
     - Auto-detected, override with ``-D BUILD_MPI=off``
   * - Enable or disable OpenMP
     - Edit ``Makefile.<machine>``
     - Auto-detected, override with ``-D BUILD_OMP=off``
   * - Build LAMMPS shared library
     - ``make mode=shared <machine>``
     - Use ``-D BUILD_SHARED_LIBS=on``
   * - Install LAMMPS Python module
     - ``make install-python``
     - Use ``cmake --build build --target install-python`` (requires ``-D BUILD_SHARED_LIBS=on``)

Some general CMake tips:

* You can have multiple independent build folders with different configurations (e.g. MPI on/off or different package selections). The convention is that the name starts with "build" (e.g. ``build-serial``).
* To build a LAMMPS executable with defaults and most packages that do not need manual configurations: ``cmake -S cmake -B build -C cmake/presets/most.cmake; cmake --build build``
* To make multiple configuration changes interactively use either ``ccmake build`` or edit ``build/CMakeCache.txt`` in a text editor and run ``cmake -S cmake -B build``
