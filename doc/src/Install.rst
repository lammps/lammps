Install LAMMPS
**************

You can download LAMMPS as an executable or as source code.

With source code, you also have to :doc:`build LAMMPS <Build>`.  But you
have more flexibility as to what features to include or exclude in the
build.  If you plan to :doc:`modify or extend LAMMPS <Modify>`, then you
need the source code.


.. toctree::
   :maxdepth: 1

   Install_linux
   Install_mac
   Install_windows
   Install_conda

   Install_tarball
   Install_git
   Install_patch

These are the files and sub-directories in the LAMMPS distribution:

+------------+-------------------------------------------+
| README     | text file                                 |
+------------+-------------------------------------------+
| LICENSE    | GNU General Public License (GPL)          |
+------------+-------------------------------------------+
| bench      | benchmark problems                        |
+------------+-------------------------------------------+
| cmake      | CMake build files                         |
+------------+-------------------------------------------+
| doc        | documentation                             |
+------------+-------------------------------------------+
| examples   | simple test problems                      |
+------------+-------------------------------------------+
| lib        | additional provided or external libraries |
+------------+-------------------------------------------+
| potentials | interatomic potential files               |
+------------+-------------------------------------------+
| python     | Python wrapper on LAMMPS                  |
+------------+-------------------------------------------+
| src        | source files                              |
+------------+-------------------------------------------+
| tools      | pre- and post-processing tools            |
+------------+-------------------------------------------+

You will have all of these if you download source.  You will only have
some of them if you download executables, as explained on the pages
listed above.
