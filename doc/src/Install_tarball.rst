Download source and documentation as a tarball
==============================================

You can download a current LAMMPS tarball from the `download page <download_>`_
of the `LAMMPS website <lws_>`_.

.. _download: http://lammps.sandia.gov/download.html
.. _bug: http://lammps.sandia.gov/bug.html
.. _older: http://lammps.sandia.gov/tars
.. _lws: http://lammps.sandia.gov

You have two choices of tarballs, either the most recent stable
release or the most current patch release.  Stable releases occur a
few times per year, and undergo more testing before release.  Patch
releases occur a couple times per month.  The new contents in all
releases are listed on the `bug and feature page <bug_>`_ of the website.

Both tarballs include LAMMPS documentation (HTML and PDF files)
corresponding to that version.  The download page also has an option
to download the current-version LAMMPS documentation by itself.

Older versions of LAMMPS can also be downloaded from `this page <older_>`_.

Once you have a tarball, unzip and untar it with the following
command:

.. code-block:: bash

   $ tar -xzvf lammps\*.tar.gz

This will create a LAMMPS directory with the version date
in its name, e.g. lammps-23Jun18.

----------

You can also download a zip file via the "Clone or download" button on
the `LAMMPS GitHub site <git_>`_.  The file name will be lammps-master.zip
which can be unzipped with the following command, to create
a lammps-master dir:

.. code-block:: bash

   $ unzip lammps\*.zip

This version is the most up-to-date LAMMPS development version.  It
will have the date of the most recent patch release (see the file
src/version.h).  But it will also include any new bug-fixes or
features added since the last patch release.  They will be included in
the next patch release tarball.

.. _git: https://github.com/lammps/lammps

----------

If you download a current LAMMPS tarball, one way to stay current as
new patch tarballs are released, is to download a patch file which you
can apply to your local directory to update it for each new patch
release.  (Or of course you could just download the newest tarball
periodically.)

The patch files are posted on the `bug and feature page <bug_>`_ of the
website, along with a list of changed files and details about what is
in the new patch release.  Instructions for applying a patch file are
on the :doc:`Install patch <Install_patch>` doc page.
