Applying patches
================

It is easy to stay current with the most recent LAMMPS patch releases
if you use Git or SVN to track LAMMPS development.  Instructions for
how to stay current are on the :doc:`Install git <Install_git>` and
:doc:`Install svn <Install_svn>` doc pages.

If you prefer to download a tarball, as described on the :doc:`Install git <Install_tarball>` doc page, you can stay current by
downloading "patch files" when new patch releases are made.  A link to
a patch file is posted on the `bug and feature page <http://lammps.sandia.gov/bug.html>`_ of the LAMMPS website, along
with a list of changed files and details about what is in the new patch
release.  This page explains how to apply the patch file to your local
LAMMPS directory.

.. note::

   You should not apply patch files to a local Git or SVN repo of
   LAMMPS, only to an unpacked tarball.  Use Git and SVN commands to
   update repo versions of LAMMPS.

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
  
  .. parsed-literal::
  
     patch -bp1 < patch.12Dec16

* A list of updated files print out to the screen.  The -b switch
  creates backup files of your originals (e.g. src/force.cpp.orig), so
  you can manually undo the patch if something goes wrong.
* Type the following from the src directory, to enforce consistency
  between the src and package directories.  This is OK to do even if you
  don't use one or more packages.  If you are applying several patches
  successively, you only need to type this once at the end. The purge
  command removes deprecated src files if any were removed by the patch
  from package sub-directories.
  
  .. parsed-literal::
  
     make purge
     make package-update

* Re-build LAMMPS via the "make" command.

.. warning::

   If you wish to edit/change a source file that is part of a package,
   you should edit the version of the file inside the package folder in
   src, and then re-install or update the package.  The version in the
   src directory is merely a copy and will be wiped out when you type
   "make package-update".
