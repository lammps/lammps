Include packages in build
=========================

In LAMMPS, a package is a group of files that enable a specific set of
features.  For example, force fields for molecular systems or
rigid-body constraints are in packages.  In the src directory, each
package is a sub-directory with the package name in capital letters.

An overview of packages is given on the :doc:`Packages <Packages>` doc
page.  Brief overviews of each package are on the :doc:`Packages details <Packages_details>` doc page.

When building LAMMPS, you can choose to include or exclude each
package.  In general there is no need to include a package if you
never plan to use its features.

If you get a run-time error that a LAMMPS command or style is
"Unknown", it is often because the command is contained in a package,
and your build did not include that package.  Running LAMMPS with the
:doc:`-h command-line switch <Run_options>` will print all the included
packages and commands for that executable.

For the majority of packages, if you follow the single step below to
include it, you can then build LAMMPS exactly the same as you would
without any packages installed.  A few packages may require additional
steps, as explained on the :doc:`Build extras <Build_extras>` doc page.

These links take you to the extra instructions for those select
packages:

+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`COMPRESS <compress>`       | :ref:`GPU <gpu>`                 | :ref:`KIM <kim>`                   | :ref:`KOKKOS <kokkos>`       | :ref:`LATTE <latte>`           | :ref:`MESSAGE <message>`             |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`MSCG <mscg>`               | :ref:`OPT <opt>`                 | :ref:`POEMS <poems>`               | :ref:`PYTHON <python>`       | :ref:`VORONOI <voronoi>`       | :ref:`USER-ADIOS <user-adios>`       |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`USER-ATC <user-atc>`       | :ref:`USER-AWPMD <user-awpmd>`   | :ref:`USER-COLVARS <user-colvars>` | :ref:`USER-H5MD <user-h5md>` | :ref:`USER-INTEL <user-intel>` | :ref:`USER-MOLFILE <user-molfile>`   |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`USER-NETCDF <user-netcdf>` | :ref:`USER-PLUMED <user-plumed>` | :ref:`USER-OMP <user-omp>`         | :ref:`USER-QMMM <user-qmmm>` | :ref:`USER-QUIP <user-quip>`   | :ref:`USER-SCAFACOS <user-scafacos>` |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+
| :ref:`USER-SMD <user-smd>`       | :ref:`USER-VTK <user-vtk>`       |                                    |                              |                                |                                      |
+----------------------------------+----------------------------------+------------------------------------+------------------------------+--------------------------------+--------------------------------------+

The mechanism for including packages is simple but different for CMake
versus make.

**CMake variables**\ :


.. parsed-literal::

   -D PKG_NAME=value          # yes or no (default)

Examples:


.. parsed-literal::

   -D PKG_MANYBODY=yes
   -D PKG_USER-INTEL=yes

All standard and user packages are included the same way.  Note that
USER packages have a hyphen between USER and the rest of the package
name, not an underscore.

See the shortcut section below for how to install many packages at
once with CMake.

.. note::

   If you toggle back and forth between building with CMake vs
   make, no packages in the src directory can be installed when you
   invoke cmake.  CMake will give an error if that is not the case,
   indicating how you can un-install all packages in the src dir.

**Traditional make**\ :


.. parsed-literal::

   cd lammps/src
   make ps                    # check which packages are currently installed
   make yes-name              # install a package with name
   make no-name               # un-install a package with name
   make mpi                   # build LAMMPS with whatever packages are now installed

Examples:


.. parsed-literal::

   make no-rigid
   make yes-user-intel

All standard and user packages are included the same way.

See the shortcut section below for how to install many packages at
once with make.

.. note::

   You must always re-build LAMMPS (via make) after installing or
   un-installing a package, for the action to take effect.

.. note::

   You cannot install or un-install packages and build LAMMPS in a
   single make command with multiple targets, e.g. make yes-colloid mpi.
   This is because the make procedure creates a list of source files that
   will be out-of-date for the build if the package configuration changes
   within the same command.  You can include or exclude multiple packages
   in a single make command, e.g. make yes-colloid no-manybody.

**CMake and make info**\ :

Any package can be included or excluded in a LAMMPS build, independent
of all other packages.  However, some packages include files derived
from files in other packages.  LAMMPS checks for this and does the
right thing.  Individual files are only included if their dependencies
are already included.  Likewise, if a package is excluded, other files
dependent on that package are also excluded.

When you download a LAMMPS tarball or download LAMMPS source files
from the Git or SVN repositories, no packages are pre-installed in the
src directory.

.. note::

   Prior to Aug 2018, if you downloaded a tarball, 3 packages
   (KSPACE, MANYBODY, MOLECULE) were pre-installed in the src directory.
   That is no longer the case, so that CMake will build as-is without the
   need to un-install those packages.


----------


**CMake shortcuts for installing many packages**\ :

Instead of specifying all the CMake options via the command-line,
CMake allows initializing the variable cache using script files. These
are regular CMake files which can manipulate and set variables, and
can also contain control flow constructs.

LAMMPS includes several of these files to define configuration
"presets", similar to the options that exist for the Make based
system. Using these files you can enable/disable portions of the
available packages in LAMMPS. If you need a custom preset you can take
one of them as a starting point and customize it to your needs.

+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/all\_on.cmake  [OPTIONS] ../cmake | enable all packages                                       |
+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/all\_off.cmake [OPTIONS] ../cmake | disable all packages                                      |
+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/minimal.cmake [OPTIONS] ../cmake  | enable just a few core packages                           |
+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/most.cmake    [OPTIONS] ../cmake  | enable most common packages                               |
+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/nolib.cmake   [OPTIONS] ../cmake  | disable packages that do require extra libraries or tools |
+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/clang.cmake   [OPTIONS] ../cmake  | change settings to use the Clang compilers by default     |
+-------------------------------------------------------------+-----------------------------------------------------------+
| cmake -C ../cmake/presets/mingw.cmake [OPTIONS] ../cmake    | enable all packages compatible with MinGW compilers       |
+-------------------------------------------------------------+-----------------------------------------------------------+

.. note::

   Running cmake this way manipulates the variable cache in your
   current build directory. You can combine multiple presets and options
   in a single cmake run, or change settings incrementally by running
   cmake with new flags.

**Example:**


.. parsed-literal::

   # build LAMMPS with most commonly used packages, but then remove
   # those requiring additional library or tools, but still enable
   # GPU package and configure it for using CUDA. You can run.
   mkdir build
   cd build
   cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake -D PKG_GPU=on -D GPU_API=cuda ../cmake

   # to add another package, say BODY to the previous configuration you can run:
   cmake -D PKG_BODY=on .

   # to reset the package selection from above to the default of no packages
   # but leaving all other settings untouched. You can run:
   cmake -C ../cmake/presets/no_all.cmake .


----------


**Make shortcuts for installing many packages**\ :

The following commands are useful for managing package source files
and their installation when building LAMMPS via traditional make.
Just type "make" in lammps/src to see a one-line summary.

These commands install/un-install sets of packages:

+-----------------------------------+-----------------------------------------------------+
| make yes-all                      | install all packages                                |
+-----------------------------------+-----------------------------------------------------+
| make no-all                       | un-install all packages                             |
+-----------------------------------+-----------------------------------------------------+
| make yes-standard or make yes-std | install standard packages                           |
+-----------------------------------+-----------------------------------------------------+
| make no-standard or make no-std   | un-install standard packages                        |
+-----------------------------------+-----------------------------------------------------+
| make yes-user                     | install user packages                               |
+-----------------------------------+-----------------------------------------------------+
| make no-user                      | un-install user packages                            |
+-----------------------------------+-----------------------------------------------------+
| make yes-lib                      | install packages that require extra libraries       |
+-----------------------------------+-----------------------------------------------------+
| make no-lib                       | un-install packages that require extra libraries    |
+-----------------------------------+-----------------------------------------------------+
| make yes-ext                      | install packages that require external libraries    |
+-----------------------------------+-----------------------------------------------------+
| make no-ext                       | un-install packages that require external libraries |
+-----------------------------------+-----------------------------------------------------+

which install/un-install various sets of packages.  Typing "make
package" will list all the these commands.

.. note::

   Installing or un-installing a package works by simply copying
   files back and forth between the main src directory and
   sub-directories with the package name (e.g. src/KSPACE, src/USER-ATC),
   so that the files are included or excluded when LAMMPS is built.

The following make commands help manage files that exist in both the
src directory and in package sub-directories.  You do not normally
need to use these commands unless you are editing LAMMPS files or are
:doc:`installing a patch <Install_patch>` downloaded from the LAMMPS web
site.

Type "make package-status" or "make ps" to show which packages are
currently installed.  For those that are installed, it will list any
files that are different in the src directory and package
sub-directory.

Type "make package-installed" or "make pi" to show which packages are
currently installed, without listing the status of packages that are
not installed.

Type "make package-update" or "make pu" to overwrite src files with
files from the package sub-directories if the package is installed.
It should be used after a :doc:`patch has been applied <Install_patch>`,
since patches only update the files in the package sub-directory, but
not the src files.

Type "make package-overwrite" to overwrite files in the package
sub-directories with src files.

Type "make package-diff" to list all differences between pairs of
files in both the src dir and a package dir.
