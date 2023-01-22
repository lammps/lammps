Install LAMMPS
==============

You can download LAMMPS as an executable or as source code.

When downloading the LAMMPS source code, you also have to :doc:`build
LAMMPS <Build>`.  But you have more flexibility as to what features to
include or exclude in the build.  When you download and install
pre-compiled LAMMPS executables, you are limited to install which
version of LAMMPS is available and which features are included of these
builds.  If you plan to :doc:`modify or extend LAMMPS <Modify>`, then
you **must** build LAMMPS from the source code.

.. note::

   If you have questions about the pre-compiled LAMMPS executables, you
   need to contact the people preparing those executables.  The LAMMPS
   developers have no control over their choices of how they configure
   and build their packages and when they update them.

----

.. toctree::
   :maxdepth: 1

   Install_linux
   Install_mac
   Install_windows
   Install_conda

   Install_tarball
   Install_git

----

These are the files and subdirectories in the LAMMPS distribution:

+------------+---------------------------------------------+
| README     | Short description of the LAMMPS package     |
+------------+---------------------------------------------+
| LICENSE    | GNU General Public License (GPL)            |
+------------+---------------------------------------------+
| SECURITY.md| Security policy for the LAMMPS package      |
+------------+---------------------------------------------+
| bench      | benchmark inputs                            |
+------------+---------------------------------------------+
| cmake      | CMake build files                           |
+------------+---------------------------------------------+
| doc        | documentation and tools to build the manual |
+------------+---------------------------------------------+
| examples   | example input files                         |
+------------+---------------------------------------------+
| fortran    | Fortran module for LAMMPS library interface |
+------------+---------------------------------------------+
| lib        | additional provided or external libraries   |
+------------+---------------------------------------------+
| potentials | selected interatomic potential files        |
+------------+---------------------------------------------+
| python     | Python module for LAMMPS library interface  |
+------------+---------------------------------------------+
| src        | LAMMPS source files                         |
+------------+---------------------------------------------+
| tools      | pre- and post-processing tools              |
+------------+---------------------------------------------+
| unittest   | source code and inputs for testing LAMMPS   |
+------------+---------------------------------------------+

You will have all of these if you downloaded the LAMMPS source code.
You will have only some of them if you downloaded executables, as
explained on the pages listed above.
