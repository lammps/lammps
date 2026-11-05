Include packages in build
=========================

In LAMMPS, a package is a group of files that enable a specific set of
features.  For example, force fields for molecular systems or
rigid-body constraints are in packages.  In the src directory, each
package is a subdirectory with the package name in capital letters.

An overview of packages is given on the :doc:`Packages <Packages>` doc
page.  Brief overviews of each package are on the :doc:`Packages details
<Packages_details>` page.

When building LAMMPS, you can choose to include or exclude each
package.  Generally, there is no need to include a package if you
never plan to use its features.

If you get a run-time error that a LAMMPS command or style is
"unknown", it is often because the command is contained in a package,
and your build did not include that package.  If the command or style
*is* available in a package included in the LAMMPS distribution,
the error message will indicate which package would be needed.
Running LAMMPS with the :doc:`-h command-line switch <Run_options>`
will print *all* optional commands and packages that were enabled
when building that executable.

For the majority of packages, if you follow the single step below to
include it, you can then build LAMMPS exactly the same as you would
without any packages installed.  A few packages may require additional
steps, as explained on the :doc:`Build extras <Build_extras>` page.

These links take you to the extra instructions for those select
packages:

.. this list must be kept in sync with its counterpart in Build_extras.rst
.. table_from_list::
   :columns: 6

   * :ref:`ADIOS <adios>`
   * :ref:`APIP <apip>`
   * :ref:`COLVARS <colvar>`
   * :ref:`COMPRESS <compress>`
   * :ref:`ELECTRODE <electrode>`
   * :ref:`GPU <gpu>`
   * :ref:`H5MD <h5md>`
   * :ref:`INTEL <intel>`
   * :ref:`KIM <kim>`
   * :ref:`KOKKOS <kokkos>`
   * :ref:`LEPTON <lepton>`
   * :ref:`MACHDYN <machdyn>`
   * :ref:`MDI <mdi>`
   * :ref:`MISC <misc>`
   * :ref:`ML-HDNNP <ml-hdnnp>`
   * :ref:`ML-IAP <mliap>`
   * :ref:`ML-PACE <ml-pace>`
   * :ref:`ML-POD <ml-pod>`
   * :ref:`ML-QUIP <ml-quip>`
   * :ref:`ML-RUNNER <ml-runner>`
   * :ref:`MOLFILE <molfile>`
   * :ref:`NETCDF <netcdf>`
   * :ref:`OPENMP <openmp>`
   * :ref:`OPT <opt>`
   * :ref:`PLUMED <plumed>`
   * :ref:`PYTHON <python>`
   * :ref:`QMMM <qmmm>`
   * :ref:`QMMM-XTB <qmmm-xtb>`
   * :ref:`RHEO <rheo>`
   * :ref:`SCAFACOS <scafacos>`
   * :ref:`VORONOI <voronoi>`

Selecting packages
------------------

The mechanism for including packages is simple:

.. code-block:: csh

   -D PKG_NAME=value          # on or off (default)

Examples:

.. code-block:: csh

   -D PKG_MANYBODY=on
   -D PKG_INTEL=on

All packages are included the same way.  By default **no** packages are
included.  See the shortcut section below for how to install many
packages at once.

Almost all packages can be included or excluded in a LAMMPS build,
independent of the other packages.  However, some packages include files
derived from files in other packages.  LAMMPS checks for this and
(usually) does the right thing automatically: Individual files are only
included if their dependencies are already included.  Likewise, if a
package is excluded, other files dependent on that package are also
excluded.  If an entire package depends on some other package that
additional package will either be included automatically or CMake will
stop with an error.  Information about this will be printed to the
screen when running the CMake configuration.

----------

.. _cmake_presets:

Presets for installing many packages
------------------------------------

Instead of specifying all the CMake options via the command-line,
CMake allows initializing its settings cache using script files.
These are regular CMake files which can manipulate and set CMake
variables (which represent selected options), and can also contain
control flow constructs for more complex operations.

LAMMPS includes several of these files to define configuration
"presets".  Using these files, you can enable/disable portions of the
available packages in LAMMPS. If you need a custom preset, you can
make a copy of one of them and modify it to suit your needs.

.. code-block:: bash

    # enable just a few core packages (MOLECULE, KSPACE, MANYBODY, RIGID, GRAPHICS)
    cmake -S cmake -B build -C cmake/presets/basic.cmake    [OPTIONS]

    # enable most packages
    cmake -S cmake -B build -C cmake/presets/most.cmake     [OPTIONS]

    # enable packages which download sources or potential files
    cmake -S cmake -B build -C cmake/presets/download.cmake [OPTIONS]

    # disable packages that do require extra libraries or tools
    cmake -S cmake -B build -C cmake/presets/nolib.cmake    [OPTIONS]

    # change settings to use the Clang compilers by default
    cmake -S cmake -B build -C cmake/presets/clang.cmake    [OPTIONS]

    # change settings to use the GNU compilers by default
    cmake -S cmake -B build -C cmake/presets/gcc.cmake      [OPTIONS]

    # like gcc.cmake but with additional compiler warnings enabled
    # (useful during development)
    cmake -S cmake -B build -C cmake/presets/pedantic.cmake [OPTIONS]

    # change settings to use the Intel oneAPI compilers by default
    cmake -S cmake -B build -C cmake/presets/oneapi.cmake   [OPTIONS]

    # change settings to use the PGI compilers by default
    cmake -S cmake -B build -C cmake/presets/pgi.cmake      [OPTIONS]

    # change settings to use the NVIDIA HPC compilers by default
    cmake -S cmake -B build -C cmake/presets/nvhpc.cmake    [OPTIONS]

    # enable all packages
    cmake -S cmake -B build -C cmake/presets/all_on.cmake   [OPTIONS]

    # disable all packages
    cmake -S cmake -B build -C cmake/presets/all_off.cmake  [OPTIONS]

    # enable all packages that contain styles with KOKKOS support
    cmake -S cmake -B build -C cmake/presets/kokkos-packages.cmake [OPTIONS]

    # enable all packages that contain styles with GPU package support
    cmake -S cmake -B build -C cmake/presets/gpu-packages.cmake [OPTIONS]

    #  compile with MinGW cross-compilers
    mingw64-cmake -S cmake -B build -C cmake/presets/mingw-cross.cmake [OPTIONS]

    # compile serial multi-arch binaries on macOS
    cmake -S cmake -B build -C cmake/presets/macos-multiarch.cmake [OPTIONS]

Presets that have names starting with "windows" are specifically for
compiling LAMMPS :doc:`natively on Windows <Build_windows>`, presets
that have names starting with "kokkos" are specifically for selecting
configurations for compiling LAMMPS with :ref:`KOKKOS <kokkos>`, and
presets that have names starting with "gpu" are specifically for
compiling LAMMPS with the :ref:`GPU package <gpu>`.  The
``kokkos-packages.cmake`` and ``gpu-packages.cmake`` presets only
select packages; they do *not* enable the KOKKOS or GPU package itself
and are meant to be combined with other presets or settings that do.

.. note::

   Running cmake this way manipulates the CMake settings cache in the
   build folder.  You can combine multiple presets and options in a
   single cmake run, or change settings incrementally by running cmake
   with new flags.  If you use a preset for selecting a set of
   compilers, it will reset all settings from previous CMake runs.

Example
"""""""

.. code-block:: bash

   # build LAMMPS with most commonly used packages, but then remove
   # those requiring additional library or tools, but still enable
   # GPU package and configure it for using CUDA. You can run:
   cmake -S cmake -B build -C cmake/presets/most.cmake -C cmake/presets/nolib.cmake \
         -D PKG_GPU=on -D GPU_API=cuda

   # to add another package, say BODY to the previous configuration you can run:
   cmake -S cmake -B build -D PKG_BODY=on

   # to reset the package selection from above to the default of no packages
   # but leaving all other settings untouched. You can run:
   cmake -S cmake -B build -C cmake/presets/all_off.cmake
