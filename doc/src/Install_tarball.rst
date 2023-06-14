Download source and documentation as a tarball
----------------------------------------------

You can download a current LAMMPS tarball from the `download page <download_>`_
of the `LAMMPS website <lws_>`_ or from GitHub (see below).

.. _download: https://www.lammps.org/download.html
.. _bug: https://www.lammps.org/bug.html
.. _older: https://download.lammps.org/tars/
.. _lws: https://www.lammps.org

You have two choices of tarballs, either the most recent stable release
or the most recent feature release.  Stable releases occur a few times
per year, and undergo more testing before release.  Also, between stable
releases bug fixes from the feature releases are back-ported and the
tarball occasionally updated.  Feature releases occur every 4 to 8
weeks.  The new contents in all feature releases are listed on the `bug
and feature page <bug_>`_ of the LAMMPS homepage.

Tarballs of older LAMMPS versions can also be downloaded from `this page
<older_>`_.

Tarballs downloaded from the LAMMPS homepage include the pre-translated
LAMMPS documentation (HTML and PDF files) corresponding to that version.

Once you have a tarball, uncompress and untar it with the following
command:

.. code-block:: bash

   tar -xzvf lammps*.tar.gz

This will create a LAMMPS directory with the version date in its name,
e.g. lammps-28Mar23.

----------

You can also download a compressed tar or zip archives from the
"Assets" sections of the `LAMMPS GitHub releases site <git_>`_.
The file name will be lammps-<version>.zip which can be unzipped
with the following command, to create a lammps-<version> directory:

.. code-block:: bash

   unzip lammps*.zip

This version corresponds to the selected LAMMPS feature or stable
release (as indicated by the matching git tag) and will only contain the
source code and no pre-built documentation.

.. _git: https://github.com/lammps/lammps/releases

