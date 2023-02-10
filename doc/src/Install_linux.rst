Download an executable for Linux
--------------------------------

Binaries are available for different versions of Linux:

- :ref:`Pre-built static Linux x86_64 executables <static>`
- :ref:`Pre-built Ubuntu and Debian Linux executables <ubuntu>`
- :ref:`Pre-built Fedora Linux executables <fedora>`
- :ref:`Pre-built EPEL Linux executables (RHEL, CentOS) <epel>`
- :ref:`Pre-built OpenSuse Linux executables <opensuse>`
- :ref:`Gentoo Linux executable <gentoo>`
- :ref:`Arch Linux build-script <arch>`

.. note::

   If you have questions about these pre-compiled LAMMPS executables,
   you need to contact the people preparing those packages.  The LAMMPS
   developers have no control over how they configure and build their
   packages and when they update them.  They may only provide packages
   for stable release versions and not always update the packages in a
   timely fashion after a new LAMMPS release is made.

----------

.. _static:

Pre-built static Linux x86_64 executables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-built LAMMPS executables for Linux, that are statically linked and
compiled for 64-bit x86 CPUs (x86_64 or AMD64) are available for download
at `https://download.lammps.org/static/ <https://download.lammps.org/static/>`_.
Because of that static linkage (and unlike the Linux distribution specific
packages listed below), they do not depend on any installed software and
thus should run on *any* 64-bit x86 machine with *any* Linux version.

These executable include most of the available packages and multi-thread
parallelization (via INTEL, KOKKOS, or OPENMP package).  They are **not**
compatible with MPI.  Several of the LAMMPS tools executables (e.g. ``msi2lmp``)
and the ``lammps-shell`` program are included as well.  Because of the
static linkage, there is no ``liblammps.so`` library file and thus also the
LAMMPS python module, which depends on it, is not included.

The compressed tar archives available for download have names following
the pattern `lammps-linux-x86_64-<version>.tar.gz` and will all unpack
into a ``lammps-static`` folder.  The executables are then in the
``lammps-static/bin/`` folder.  Since they do not depend on any other
software, they may be freely moved or copied around.

----------

.. _ubuntu:

Pre-built Ubuntu and Debian Linux executables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A pre-built LAMMPS executable, suitable for running on the latest Ubuntu
and Debian Linux versions, can be downloaded as a Debian package.  This
allows you to install LAMMPS with a single command, and stay (mostly)
up-to-date with the current stable version of LAMMPS by simply updating
your operating system.

To install LAMMPS do the following once:

.. code-block:: bash

   sudo apt-get install lammps

This downloads an executable named ``lmp`` to your box and multiple
packages with supporting data, examples and libraries as well as any
missing dependencies.  For example, the LAMMPS binary in this package is
built with the :ref:`KIM package <kim>` enabled, which results in the
above command also installing the ``kim-api`` binaries when LAMMPS is
installed, unless they were installed already.  In order to use
potentials from `openkim.org <openkim_>`_, you can also install the
``openkim-models`` package:

.. code-block:: bash

   sudo apt-get install openkim-models

Or use the `KIM-API commands <https://openkim.org/doc/usage/obtaining-models/#installing_api>`_
to download and install individual models.

This LAMMPS executable can then be used in the usual way to run input
scripts:

.. code-block:: bash

   lmp -in in.lj

To update LAMMPS to the latest packaged version, do the following:

.. code-block:: bash

   sudo apt-get update

This will also update other packages on your system.

To uninstall LAMMPS, do the following:

.. code-block:: bash

   sudo apt-get remove lammps

Please use ``lmp -help`` to see which compilation options, packages,
and styles are included in the binary.

Thanks to Anton Gladky (gladky.anton at gmail.com) for setting up this
Ubuntu package capability.

----------

.. _fedora:

Pre-built Fedora Linux executables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-built `LAMMPS packages for stable releases
<https://packages.fedoraproject.org/pkgs/lammps/>`_ are available in the
Fedora Linux distribution since Fedora version 28. The packages can be
installed via the dnf package manager. There are 3 basic varieties
(lammps = no MPI, lammps-mpich = MPICH MPI library, lammps-openmpi =
OpenMPI MPI library) and for each support for linking to the C library
interface (lammps-devel, lammps-mpich-devel, lammps-openmpi-devel), the
header for compiling programs using the C library interface
(lammps-headers), and the LAMMPS python module for Python 3. All
packages can be installed at the same time and the name of the LAMMPS
executable is ``lmp`` and ``lmp_openmpi`` or ``lmp_mpich`` respectively.
By default, ``lmp`` will refer to the serial executable, unless one of
the MPI environment modules is loaded (``module load mpi/mpich-x86_64``
or ``module load mpi/openmpi-x86_64``).  Then the corresponding parallel
LAMMPS executable can be used.  The same mechanism applies when loading
the LAMMPS python module.

To install LAMMPS with OpenMPI and run an input ``in.lj`` with 2 CPUs do:

.. code-block:: bash

   dnf install lammps-openmpi
   module load mpi/openmpi-x86_64
   mpirun -np 2 lmp -in in.lj

The ``dnf install`` command is needed only once.  In case of a new LAMMPS
stable release, ``dnf update`` will automatically update to the newer
version as soon as the RPM files are built and uploaded to the download
mirrors. The ``module load`` command is needed once per (shell) session
or shell terminal instance, unless it is automatically loaded from the
shell profile.

The LAMMPS binary is built with the :ref:`KIM package <kim>` which
results in the above command also installing the `kim-api` binaries when LAMMPS
is installed.  In order to use potentials from `openkim.org <openkim_>`_, you
can install the `openkim-models` package

.. code-block:: bash

   dnf install openkim-models

Please use ``lmp -help`` to see which compilation options, packages,
and styles are included in the binary.

Thanks to Christoph Junghans (LANL) for making LAMMPS available in Fedora.

.. _openkim: https://openkim.org

----------

.. _epel:

Pre-built EPEL Linux executable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-built LAMMPS (and KIM) packages for stable releases are available
in the `Extra Packages for Enterprise Linux (EPEL) repository <https://fedoraproject.org/wiki/EPEL>`_
for use with Red Hat Enterprise Linux (RHEL) or CentOS version 7.x
and compatible Linux distributions. Names of packages, executable,
and content are the same as described above for Fedora Linux.
But RHEL/CentOS 7.x uses the ``yum`` package manager instead of ``dnf``
in Fedora 28.

Please use ``lmp -help`` to see which compilation options, packages,
and styles are included in the binary.

Thanks to Christoph Junghans (LANL) for making LAMMPS available in EPEL.

----------

.. _opensuse:

Pre-built OpenSuse Linux executable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A pre-built LAMMPS package for stable releases is available
in OpenSuse as of Leap 15.0. You can install the package with:

.. code-block:: bash

   zypper install lammps

This includes support for OpenMPI. The name of the LAMMPS executable
is ``lmp``. To run an input in parallel on 2 CPUs you would do:

.. code-block:: bash

   mpirun -np 2 lmp -in in.lj

Please use ``lmp -help`` to see which compilation options, packages,
and styles are included in the binary.

The LAMMPS binary is built with the :ref:`KIM package <kim>` which
results in the above command also installing the `kim-api` binaries when LAMMPS
is installed.  In order to use potentials from `openkim.org <openkim_>`_, you
can install the `openkim-models` package

.. code-block:: bash

   zypper install openkim-models

Thanks to Christoph Junghans (LANL) for making LAMMPS available in OpenSuse.

----------

.. _gentoo:

Gentoo Linux executable
^^^^^^^^^^^^^^^^^^^^^^^

LAMMPS is part of `Gentoo's main package tree
<https://packages.gentoo.org/packages/sci-physics/lammps>`_ and can be
installed by typing:

.. code-block:: bash

   emerge --ask lammps

Note that in Gentoo the LAMMPS source code is downloaded and the package is
then compiled and installed on your machine.

Certain LAMMPS packages can be enabled via USE flags, type

.. code-block:: bash

   equery uses lammps

for details.

Thanks to Nicolas Bock and Christoph Junghans (LANL) for setting up
this Gentoo capability.

----------

.. _arch:

Archlinux build-script
^^^^^^^^^^^^^^^^^^^^^^

LAMMPS is available via Arch's unofficial Arch User repository (AUR).
There are three scripts available, named `lammps
<https://aur.archlinux.org/packages/lammps>`_, `lammps-beta
<https://aur.archlinux.org/packages/lammps>`_ and `lammps-git
<https://aur.archlinux.org/packages/lammps>`_.  They respectively
package the stable, feature, and git releases.

To install, you will need to have the git package installed. You may use
any of the above names in-place of lammps.

.. code-block:: bash

   git clone https://aur.archlinux.org/lammps.git
   cd lammps
   makepkg -s
   makepkg -i

To update LAMMPS, you may repeat the above, or change into the cloned
directory, and execute the following, after which, if there are any
changes, you may use makepkg as above.

.. code-block:: bash

   git pull

Alternatively, you may use an AUR helper to install these packages.

Note that the AUR provides build-scripts that download the source code
and then build and install the package on your machine.
