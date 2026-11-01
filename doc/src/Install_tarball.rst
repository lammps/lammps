Download source and documentation as a tarball
----------------------------------------------

You can download a current LAMMPS tarball from the `download page <download_>`_
of the `LAMMPS website <lws_>`_ or from GitHub (see below).

.. _download: https://www.lammps.org/download/
.. _older: https://download.lammps.org/tars/
.. _lws: https://www.lammps.org
.. _git: https://github.com/lammps/lammps/releases

You have two choices of tarballs, either the most recent stable release
or the most recent feature release.  Stable releases occur a few times
per year, and undergo more testing before release.  Also, between stable
releases bug fixes from the feature releases are back-ported and the
tarball occasionally updated.  Feature releases occur every 4 to 8
weeks.  The new contents in all feature releases are listed in the
`release notes <git_>`_ of the LAMMPS GitHub page.

Tarballs of older LAMMPS versions can also be downloaded from `this page
<older_>`_.

Tarballs downloaded from the LAMMPS homepage include the pre-translated
LAMMPS documentation (HTML and PDF files) corresponding to that version.

Once you have a tarball, uncompress and untar it with the following
command:

.. code-block:: bash

   tar -xzvf lammps*.tar.gz

This will create a LAMMPS directory with the version date in its name,
e.g. ``lammps-28Mar23``.

----------

You can also download the same compressed tar archives from the
"Assets" sections of the `LAMMPS GitHub releases site <git_>`_.

----------

.. _verify_download:

Check downloaded files
^^^^^^^^^^^^^^^^^^^^^^

You can check that a downloaded file is complete and has not been
modified by comparing its SHA-256 checksum with the published checksum.
For files from the LAMMPS download server, download the file
``SHA256SUMS`` (in some folders ``SHA256SUM``) from the same web page
(for example from `this page <older_>`_) into the folder with the
downloaded files and type:

.. code-block:: bash

   sha256sum --ignore-missing -c SHA256SUMS

This prints the name of each downloaded file that is in the list
followed by either "OK" or "FAILED".  On macOS, please use ``shasum -a
256`` instead of ``sha256sum``.  On Windows, you can print the checksum
of a file in a PowerShell window with a command like the following and
compare it with the corresponding line of the ``SHA256SUMS`` file
(Windows prints the checksum with upper case letters):

.. code-block:: text

   Get-FileHash LAMMPS-Win10-x64-GUI-30Sep2026_update1.exe

On the `LAMMPS GitHub releases page <git_>`_, the SHA-256 checksum of
each file is shown next to its name (after ``sha256:``), so that you can
compare it directly with the checksum of a single downloaded file.
LAMMPS releases on GitHub that were published after the stable release
of 30 September 2026 also contain a file ``SHA256SUMS`` with the
checksums of all files of the release, and in addition a file
``SHA256SUMS.asc`` with a digital signature for it.  With this signature
you can confirm that the list of checksums was created by the LAMMPS
developer who published the release:

.. code-block:: bash

   gpg --fetch-keys https://github.com/akohlmey.gpg
   gpg --verify SHA256SUMS.asc SHA256SUMS

The first command is needed only once.  It imports the public key of
Axel Kohlmeyer, who currently signs the LAMMPS releases.  The second
command must report a "Good signature" and the following "Primary key
fingerprint":

.. parsed-literal::

   EEA1 0376 4C6C 633E DC8A  C428 D9B4 4E93 BF0C 375A

The ``gpg`` command will also print a warning that the key is "not
certified with a trusted signature".  This is expected, since it only
means that you have not told ``gpg`` that you trust this key, and it is
the reason for comparing the fingerprint.

The ``gpg`` program is part of `GnuPG <https://www.gnupg.org/>`_.  On
Linux it is usually installed already or available as a package named
``gnupg`` or ``gnupg2``.  For macOS, it is available from the `GPG Suite
<https://gpgtools.org/>`_ or from Homebrew (``brew install gnupg``), and
for Windows from `Gpg4win <https://www.gpg4win.org/>`_.  How to check a
signature and how to read the output of ``gpg`` is also explained on the
`integrity check page <https://www.gnupg.org/download/integrity_check.html>`_
of the GnuPG website, and the `GNU Privacy Handbook
<https://www.gnupg.org/gph/en/manual.html>`_ is a general introduction
to GnuPG.

All LAMMPS releases that were published on GitHub since 10 September
2025 are "immutable releases".  Once such a release is published,
neither its tag nor any of its files can be changed or replaced anymore.
This includes the files ``SHA256SUMS`` and ``SHA256SUMS.asc``, which are
uploaded before the release is published.  In addition, GitHub records
the tag and the checksums of all files of the release in a "release
attestation" that is signed by GitHub.  Releases and individual files
can be checked against it with the `GitHub command-line program
<https://cli.github.com/>`_, for example:

.. code-block:: bash

   gh release verify --repo lammps/lammps stable_30Sep2026
   gh release verify-asset --repo lammps/lammps stable_30Sep2026 lammps-src-30Sep2026.tar.gz

The checksums confirm that a file is complete and identical to the
published file, the signature of ``SHA256SUMS`` confirms that the list of
checksums was created by the LAMMPS developer who published the release
(which is also useful for copies of the files that were obtained
elsewhere), and the release attestation confirms that a file belongs to
that release on GitHub.
