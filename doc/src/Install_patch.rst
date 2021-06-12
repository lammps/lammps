Applying patches
----------------

It is easy to stay current with the most recent LAMMPS patch releases
if you use git to track the LAMMPS development.  Instructions for
how to stay current are on the
:doc:`Download the LAMMPS source with git <Install_git>` page.

If you prefer to download a tarball, as described on the
:doc:`tarball download <Install_tarball>` page, you can stay current by
downloading "patch files" when new patch releases are made.  A link to
a patch file is posted on the
`bug fixes and new feature page <https://www.lammps.org/bug.html>`_
of the LAMMPS website, along
with a list of changed files and details about what is in the new patch
release.  This page explains how to apply the patch file to your local
LAMMPS directory.

.. note::

   You should not apply patch files to a local git checkout of
   LAMMPS, only to an unpacked tarball.  Use git commands to
   update such a version of the LAMMPS source code.

Here are the steps to apply a patch file.  Note that if your version
of LAMMPS is several patch releases behind, you need to apply all the
intervening patch files in succession to bring your version of LAMMPS
up to date.

* Download the patch file.  You may have to shift-click in your browser
  to download the file instead of display it.  Patch files have names
  like patch.12Dec16.
* Put the patch file in your top-level LAMMPS directory, where the
  LICENSE and README files are.
* Apply the patch by typing the following command from your top-level
  LAMMPS directory, where the redirected file is the name of the patch
  file.

  .. code-block:: bash

     $ patch -bp1 < patch.12Dec16

* A list of updated files print out to the screen.  The -b switch
  creates backup files of your originals (e.g. src/force.cpp.orig), so
  you can manually undo the patch if something goes wrong.

* Once you have updated your local files you need to re-build LAMMPS.
  If you are applying several patches successively, you only need to
  do the rebuild once at the end. How to do it depends on the build
  system you are using.

  .. tabs::

     .. tab:: CMake build

        Change to your build folder and type:

        .. code-block:: bash

           cmake . --build

        CMake should auto-detect whether it needs to re-run the CMake
        configuration step and otherwise redo the build for all files
        that have been changed or files that depend on changed files.
        In case some build options have been changed or renamed, you
        may have to update those by running:

        .. code-block:: bash

           cmake .

        and then rebuild.

     .. tab:: Traditional make

        Switch to the src directory and type:

        .. code-block:: bash

           $ make purge             # remove any deprecated src files
           $ make package-update    # sync package files with src files
           $ make foo               # re-build for your machine (mpi, serial, etc)

        to enforce consistency of the source between the src folder
        and package directories.  This is OK to do even if you don't
        use any packages. The "make purge" command removes any deprecated
        src files if they were removed by the patch from a package
        sub-directory.

        .. warning::

           If you wish to edit/change a src file that is from a package,
           you should edit the version of the file inside the package
           sub-directory with src, then re-install the package.  The
           version in the source directory is merely a copy and will be
           wiped out if you type "make package-update".
