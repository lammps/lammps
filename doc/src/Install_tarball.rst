Download source and documentation as a tarball
----------------------------------------------

You can download a current LAMMPS tarball from the `download page <download_>`_
of the `LAMMPS website <lws_>`_.

.. _download: https://www.lammps.org/download.html
.. _bug: https://www.lammps.org/bug.html
.. _older: https://download.lammps.org/tars/
.. _lws: https://www.lammps.org

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

   tar -xzvf lammps*.tar.gz

This will create a LAMMPS directory with the version date
in its name, e.g. lammps-23Jun18.

----------

You can also download a compressed tar or zip archives from the
"Assets" sections of the `LAMMPS GitHub releases site <git_>`_.
The file name will be lammps-<version>.zip which can be unzipped
with the following command, to create a lammps-<version> dir:

.. code-block:: bash

   unzip lammps*.zip

This version corresponds to the selected LAMMPS patch or stable
release.

.. _git: https://github.com/lammps/lammps/releases

