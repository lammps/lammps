Security considerations
=======================

.. raw:: html

   <a href="https://www.bestpractices.dev/projects/15199"><img src="https://www.bestpractices.dev/projects/15199/badge"></a>

This page first describes what a LAMMPS input file may do on the
machine where LAMMPS is run.  It then describes which files from outside
of the LAMMPS distribution may be used when building LAMMPS, how the
LAMMPS build process checks those files, and where these checks end.  It
also explains what applies to LAMMPS when it is provided by other
projects and to add-on packages for LAMMPS that are maintained
elsewhere.  It is meant for users that have to follow security rules of
their institution or computing facility, or that want to know what
exactly they are compiling and running.

General information about the security of LAMMPS, including how to
report problems, is in the file `SECURITY.md
<https://github.com/lammps/lammps/blob/develop/SECURITY.md>`_ of the
LAMMPS repository on GitHub.  How to check that the LAMMPS source code
itself was downloaded completely and without modifications is explained
in the sections about :ref:`downloading tarballs <verify_download>` and
about :doc:`downloading with git <Install_git>`.

What a LAMMPS input file can do
-------------------------------

LAMMPS is controlled by the commands in its input.  An input file should
therefore be regarded as a program and not as data: running an input
file that you have received from somebody else is the same as running a
program or a script from that person.  LAMMPS has no means to limit what
the commands in an input may do, so an input can do everything that you
are allowed to do on the machine where LAMMPS is run.  In particular:

- **Files:** An input can read, create, overwrite, and delete any files
  and folders that you have access to, for example with the :doc:`shell
  <shell>`, :doc:`print <print>`, :doc:`write_data <write_data>`, or
  :doc:`dump <dump>` commands.
- **Other programs:** The :doc:`shell <shell>` command hands its
  arguments to the command-line interpreter of the operating system and
  can thus run any program.  LAMMPS also runs other programs on its own:
  ``gzip`` and similar programs for :ref:`reading and writing compressed
  files <gzip>`, and ``ffmpeg`` for creating movies with the :doc:`dump
  movie <dump_image>` command.
- **Program code:** The :doc:`python <python>` command and other
  features of the PYTHON package run Python code.  This includes machine
  learning models that the ML-IAP package loads through Python, because
  loading such a file can run Python code that is stored in it.  The
  :doc:`plugin <plugin>` command, the :doc:`dump molfile <dump_molfile>`
  command, and the KIM package load compiled code from files at run
  time.
- **Network:** The :doc:`geturl <geturl>` command downloads files from
  the internet, and the :doc:`kim query <kim_commands>` command sends
  requests to a server of the OpenKIM project.  The :doc:`fix imd
  <fix_imd>` command waits for other programs to connect to LAMMPS over
  the network, and the :doc:`fix ipi <fix_ipi>` command and the commands
  of the :doc:`MDI package <mdi>` exchange data with other programs,
  which may run on other machines.

This has the following consequences:

- Only run input files, Python code, plugins, and model files from
  sources that you trust, and look at input files from other people
  before running them.
- Never run LAMMPS as superuser or administrator.  The LAMMPS executable
  prints a :ref:`warning <err0040>` when it is started that way, and
  CMake prints a warning when it is used to configure LAMMPS as
  superuser.
- Most of the features mentioned above belong to optional packages
  (EXTRA-COMMAND, GRAPHICS, KIM, MDI, MISC, ML-IAP, MOLFILE, PLUGIN,
  PYTHON) and are not available when those packages are not included
  when compiling LAMMPS.  The use of ``gzip`` and ``ffmpeg`` and the
  support for downloads with the :doc:`geturl <geturl>` command can be
  turned off with the CMake settings ``-D WITH_GZIP=off``, ``-D
  WITH_FFMPEG=off``, and ``-D WITH_CURL=off``.  The :doc:`shell <shell>`
  command is always available.  Running LAMMPS with the :doc:`"-h" flag
  <Run_options>` shows which packages are included in an executable.
- If you need to restrict what LAMMPS can do beyond that, this has to be
  done with the means of the operating system, for example by running
  LAMMPS from a separate account with limited permissions or inside of a
  container.

What is included in the LAMMPS distribution
-------------------------------------------

All changes to the LAMMPS source code are submitted as pull requests on
GitHub.  They are reviewed by LAMMPS developers and have to pass
automated tests before they are included.

The LAMMPS distribution also contains copies of source code from other
projects, for example for the KOKKOS, COLVARS, and LEPTON packages.
These copies are in the ``lib``, ``src``, and ``third_party`` folders.
They are updated by LAMMPS developers from releases of those projects,
and each update is reviewed and tested like any other change to LAMMPS.
Since these copies are part of the LAMMPS source code, they are covered
by the checksums and signatures of the LAMMPS releases.

If you compile LAMMPS without any of the packages and features that are
listed :ref:`further below <security_list>`, nothing is downloaded
during the build.

How downloaded files are checked
--------------------------------

Some optional packages and features need source code from other projects
that is not included in the LAMMPS distribution.  When building LAMMPS
with CMake, this source code can be downloaded automatically.  For these
downloads the following applies:

- The LAMMPS developers select a specific version of the external source
  code and test LAMMPS with it.  The location of the corresponding
  archive file and the SHA-256 checksum of that file are recorded in the
  LAMMPS source code.
- All downloads use encrypted connections (``https``).
- The checksum of a downloaded archive is compared with the recorded
  checksum *before* the archive is unpacked.  If the checksums differ,
  the archive is not used.
- For many archives there is a copy on the LAMMPS download server
  (``download.lammps.org``).  This copy is used when the download from
  the original location fails or when its checksum does not match.  The
  copy must have the same recorded checksum.  If neither download can be
  used, the build stops with an error.
- In most cases nothing is downloaded when a suitable version of the
  external library is already installed on your machine and CMake finds
  it.  That version is then used instead.

.. note::

   A matching checksum confirms that you are using exactly the same
   archive that the LAMMPS developers have selected and tested LAMMPS
   with.  It does **not** mean that the LAMMPS developers have reviewed
   the content of that archive.

The legacy build system using GNU make does not download source code
from other projects.  It only downloads the potential file for the
MESONT package mentioned below, and checks it in the same way.

What the LAMMPS developers do not control
-----------------------------------------

The external projects are developed and managed independently of LAMMPS.
The LAMMPS developers have no control over who may change the source
code of those projects, how changes are reviewed, how releases are made,
or how well the accounts and servers of those projects are protected.
The checks described above can therefore not rule out that the selected
version of an external project contains mistakes or malicious code.

In addition, not all external sources can be checked with a checksum:

- Some sources are obtained from a git repository instead of an archive.
  When a specific commit is requested, git itself ensures that the
  content is exactly that of the requested commit.  When a tag or a
  branch is requested instead, the owners of the repository can usually
  change at any time what will be downloaded.  An exception are tags of
  releases that GitHub protects against changes after the release was
  published ("immutable releases"), as used for LAMMPS and LAMMPS-GUI.
- Python packages that are installed with ``pip`` are downloaded from
  the Python Package Index (PyPI) in their most recent compatible
  version.  Neither the versions nor the content of those packages are
  checked by the LAMMPS build process.
- Compilers, libraries, and tools that are already installed on your
  machine are used as they are.

.. _security_list:

Packages and features that use external sources
-----------------------------------------------

The following packages can download source code or data when building
LAMMPS with CMake.  Unless noted otherwise, the download only happens
when the corresponding library is not found on your machine.

.. list-table::
   :header-rows: 1
   :widths: 16 30 36 18

   * - Package
     - External source
     - Obtained as
     - Checked by
   * - :ref:`GPU <gpu>`
     - OpenCL loader (always downloaded when using OpenCL, except on
       macOS or when ``-D USE_STATIC_OPENCL_LOADER=off`` is used)
     - copy on the LAMMPS download server
     - SHA-256 checksum
   * - :ref:`GPU <gpu>`
     - CUB (when using HIP for NVIDIA GPUs)
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - :ref:`KIM <kim>`
     - KIM-API
     - release archive from the server of the OpenKIM project
     - SHA-256 checksum
   * - :ref:`KOKKOS <kokkos>`
     - Kokkos (only when requested with ``-D DOWNLOAD_KOKKOS=on``;
       otherwise the copy included in LAMMPS is used)
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - KSPACE
     - heFFTe (only with ``-D FFT_USE_HEFFTE=on``, see :ref:`FFT
       settings <fft>`)
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - :ref:`MACHDYN <machdyn>`
     - Eigen
     - copy on the LAMMPS download server
     - SHA-256 checksum
   * - :ref:`MBX <mbx>`
     - MBX
     - release archive from GitHub
     - SHA-256 checksum
   * - :ref:`MDI <mdi>`
     - MDI Library
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - MESONT
     - potential file that is too large to be included in LAMMPS
       (unless ``-D DOWNLOAD_POTENTIALS=off`` is used)
     - file on the LAMMPS download server
     - SHA-256 checksum
   * - :ref:`ML-HDNNP <ml-hdnnp>`
     - n2p2
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - :ref:`ML-PACE <ml-pace>`
     - PACE evaluator library
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - :ref:`ML-QUIP <ml-quip>`
     - QUIP and two further repositories that QUIP refers to
     - specific commit from git repositories on GitHub
     - git (commit is fixed)
   * - :ref:`ML-RUNNER <ml-runner>`
     - RuNNer (always downloaded unless ``-D DOWNLOAD_RUNNER=off``
       is used)
     - tag from a git repository on GitLab
     - not checked
   * - :ref:`PLUMED <plumed>`
     - PLUMED
     - release archive from GitHub
     - SHA-256 checksum
   * - :ref:`SCAFACOS <scafacos>`
     - ScaFaCoS
     - release archive from GitHub
     - SHA-256 checksum
   * - :ref:`VORONOI <voronoi>`
     - Voro++
     - copy on the LAMMPS download server
     - SHA-256 checksum

A "release archive" is a file that was published by the external project
for a release.  An "archive of a tagged version" is created by GitHub
from the state of the git repository that the external project has
marked with a version tag.  A "copy on the LAMMPS download server" is an
archive that the LAMMPS developers have obtained from the external
project and stored on ``download.lammps.org``.

The following other parts of the LAMMPS distribution also use external
sources when they are enabled or built.

.. list-table::
   :header-rows: 1
   :widths: 24 26 32 18

   * - Feature
     - External source
     - Obtained as
     - Checked by
   * - Unit tests (``-D ENABLE_TESTING=on``, see :ref:`testing`)
     - GoogleTest (always downloaded)
     - release archive from GitHub
     - SHA-256 checksum
   * - Unit tests
     - ``libyaml`` (when it is not found on your machine)
     - release archive from the server of the PyYAML project
     - SHA-256 checksum
   * - Tools (``-D BUILD_TOOLS=on``)
     - ``spglib`` for the :ref:`phana tool <phonon>` (always downloaded
       unless ``-D USE_SPGLIB=off`` is used)
     - archive of a tagged version from GitHub
     - SHA-256 checksum
   * - :ref:`LAMMPS-GUI <lammps_gui>` (``-D BUILD_LAMMPS_GUI=on``)
     - LAMMPS-GUI source code
     - release tag from the `LAMMPS-GUI repository
       <https://github.com/lammps/lammps-gui>`_ of the LAMMPS project on
       GitHub
     - tag of an immutable release (cannot be changed)
   * - LAMMPS-GUI
     - WHAM (unless ``-D BUILD_WHAM=off`` is used)
     - copy on the LAMMPS download server
     - SHA-256 checksum
   * - :doc:`Manual <Build_manual>`
     - Sphinx and other Python packages
     - installed with ``pip`` from PyPI; one package from a branch of a
       git repository on GitHub
     - not checked
   * - Manual
     - MathJax
     - with CMake: archive of a tagged version from GitHub; with
       ``make`` in the ``doc`` folder: tag from a git repository on
       GitHub
     - archive: SHA-256 checksum; git repository: not checked
   * - :doc:`LAMMPS Python module <Python_install>` (``install-python``)
     - Python packages for building the module
     - installed with ``pip`` from PyPI
     - not checked
   * - Packages for Windows
     - MS-MPI files for compiling; ``ffmpeg`` and ``gzip`` programs
       that are included in the installer packages
     - copies on the LAMMPS download server
     - SHA-256 checksum

Avoiding or inspecting downloads
--------------------------------

- Enable only the packages and features that you need.  Most LAMMPS
  packages do not need any external source code.
- Install the external library from a source that you trust, for example
  with the package manager of your Linux distribution, before
  configuring LAMMPS.  When CMake finds the library, it is used instead
  of a download.  Most of the packages listed above have a setting of
  the form ``-D DOWNLOAD_<NAME>=off`` to make CMake stop with an error
  instead of downloading when the library is not found.  The settings
  for each package are described in :doc:`Build_extras`.
- Use ``-D DOWNLOAD_POTENTIALS=off`` to not download potential files.
- The locations and checksums of the archives for a configured build
  folder can be listed with the command below.  Sources that are
  obtained with git or ``pip`` are not included in this list.

  .. code-block:: bash

     cmake -N -LA build | grep -E '_(URL|SHA256):'

- To use an archive that you have downloaded and inspected yourself, set
  the corresponding ``<NAME>_URL`` variable to the location of your file
  and, if it is a different version, the ``<NAME>_SHA256`` variable to
  its checksum.  Please see :ref:`this explanation <err0039>` for how
  these settings are stored in the build folder.

LAMMPS provided by other projects
---------------------------------

LAMMPS can also be obtained from other projects that prepare software
for installation.  Examples are Linux distributions like Debian, Ubuntu,
Fedora, or the EPEL repository for Red Hat Enterprise Linux, and package
managers like Homebrew, Conda, or `Spack <https://spack.io>`_.  Many
computing facilities provide their own installations of LAMMPS as well.

None of these are managed by the LAMMPS developers.  Each of those
projects follows its own conventions, security rules, and requirements.
They decide which version of LAMMPS they provide, which packages are
included, how LAMMPS is configured and compiled, whether they modify the
source code, from where the external libraries are taken, how their own
packages are checked and signed, and when they are updated.  What is
described on this page thus applies to them only as far as they use the
LAMMPS build process unchanged.  For questions about how LAMMPS from one
of those projects was built and checked, please contact the people that
have prepared it.  The :doc:`Install` section of the manual has links to
several of those projects.

The source code archives and the pre-compiled packages that are
available from the `LAMMPS website <https://www.lammps.org/download/>`_
and from the `LAMMPS releases page on GitHub
<https://github.com/lammps/lammps/releases>`_ are prepared by LAMMPS
developers.

Add-on packages maintained outside of LAMMPS
--------------------------------------------

There are many add-on packages for LAMMPS that are not part of the
LAMMPS distribution and that are developed and distributed by other
projects.  This is particularly common for software for machine learning
interatomic potentials.  There often are good reasons for this: many of
those projects provide interfaces to several simulation programs, and
managing all of those interfaces in one place is much simpler for them.
Some of these add-on packages are listed on the `External LAMMPS
packages and tools <https://www.lammps.org/ecosystem/external/>`_ page
of the LAMMPS website.

Add-on packages come in different forms: as source files that have to be
copied into the LAMMPS source code, as patches that modify the LAMMPS
source code, as :doc:`plugins <plugin>` that are loaded at run time, or
as a modified copy of the complete LAMMPS source code.  For all of them
the following needs to be considered:

- They are not reviewed or tested by the LAMMPS developers, and what is
  described on this page does not apply to them.
- They are not always in sync with the LAMMPS repository.  An add-on
  package may have been written for an older version of LAMMPS, and a
  modified copy of LAMMPS may not contain corrections that were made to
  LAMMPS after the copy was created.
- They may contain changes to the core of the LAMMPS source code, that
  is to parts that are used by all simulations and not only by the added
  features.  Such changes have not been vetted by the LAMMPS developers.

It is thus recommended to find out for which version of LAMMPS an add-on
package was written and which files of the LAMMPS distribution it
changes.  When LAMMPS was compiled from a git checkout, the output of
``lmp -h`` contains a line starting with "Git info".  If that line
contains the text ``-modified``, then files of the LAMMPS distribution
were changed before compiling.  Questions about add-on packages and
reports of problems when using them should be sent to the developers of
those packages (see also the note about external packages on the
:doc:`Packages` page).
