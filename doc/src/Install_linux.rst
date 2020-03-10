Download an executable for Linux
================================

Binaries are available for different versions of Linux:

| :ref:`Pre-built Ubuntu Linux executables <ubuntu>`
| :ref:`Pre-built Fedora Linux executables <fedora>`
| :ref:`Pre-built EPEL Linux executables (RHEL, CentOS) <epel>`
| :ref:`Pre-built OpenSuse Linux executables <opensuse>`
| :ref:`Gentoo Linux executable <gentoo>`
| :ref:`Arch Linux build-script <arch>`
|

----------

.. _ubuntu:

Pre-built Ubuntu Linux executables
-----------------------------------------------

A pre-built LAMMPS executable suitable for running on the latest
Ubuntu Linux versions, can be downloaded as a Debian package.  This
allows you to install LAMMPS with a single command, and stay
up-to-date with the current version of LAMMPS by simply updating your
operating system.

To install the appropriate personal-package archive (PPA), do the
following once:

.. code-block:: bash

   $ sudo add-apt-repository ppa:gladky-anton/lammps
   $ sudo apt-get update

To install LAMMPS do the following once:

.. code-block:: bash

   $ sudo apt-get install lammps-daily

This downloads an executable named "lmp\_daily" to your box, which
can then be used in the usual way to run input scripts:

.. code-block:: bash

   $ lmp_daily -in in.lj

To update LAMMPS to the most current version, do the following:

.. code-block:: bash

   $ sudo apt-get update

which will also update other packages on your system.

To get a copy of the current documentation and examples:

.. code-block:: bash

   $ sudo apt-get install lammps-daily-doc

which will download the doc files in
/usr/share/doc/lammps-daily-doc/doc and example problems in
/usr/share/doc/lammps-doc/examples.

Note that you may still wish to download the tarball to get potential
files and auxiliary tools.

To un-install LAMMPS, do the following:

.. code-block:: bash

   $ sudo apt-get remove lammps-daily

Note that the lammps-daily executable is built with the following
sequence of make commands, as if you had done the same with the
unpacked tarball files in the src directory:

.. code-block:: bash

    $ make yes-all
    $ make no-lib
    $ make mpi

Thus it builds with FFTW3 and OpenMPI.

Thanks to Anton Gladky (gladky.anton at gmail.com) for setting up this
Ubuntu package capability.

----------

.. _fedora:

Pre-built Fedora Linux executables
-----------------------------------------------

Pre-built LAMMPS packages for stable releases are available
in the Fedora Linux distribution as of version 28. The packages
can be installed via the dnf package manager. There are 3 basic
varieties (lammps = no MPI, lammps-mpich = MPICH MPI library,
lammps-openmpi = OpenMPI MPI library) and for each support for
linking to the C library interface (lammps-devel, lammps-mpich-devel,
lammps-openmpi-devel), the header for compiling programs using
the C library interface (lammps-headers), and the LAMMPS python
module for Python 3. All packages can be installed at the same
time and the name of the LAMMPS executable is *lmp* and *lmp\_openmpi*
or *lmp\_mpich* respectively.  By default, *lmp* will refer to the
serial executable, unless one of the MPI environment modules is loaded
("module load mpi/mpich-x86\_64" or "module load mpi/openmpi-x86\_64").
Then the corresponding parallel LAMMPS executable can be used.
The same mechanism applies when loading the LAMMPS python module.

To install LAMMPS with OpenMPI and run an input in.lj with 2 CPUs do:

.. code-block:: bash

   $ dnf install lammps-openmpi
   $ module load mpi/openmpi-x86_64
   $ mpirun -np 2 lmp -in in.lj

The "dnf install" command is needed only once. In case of a new LAMMPS
stable release, "dnf update" will automatically update to the newer
version as soon at the RPM files are built and uploaded to the download
mirrors. The "module load" command is needed once per (shell) session
or shell terminal instance, unless it is automatically loaded from the
shell profile.

The LAMMPS binary is built with the :ref:`KIM package <kim>` which
results in the above command also installing the `kim-api` binaries when LAMMPS
is installed.  In order to use potentials from `openkim.org <openkim_>`_, you
can install the `openkim-models` package

.. code-block:: bash

   $ dnf install openkim-models

Please use "lmp -help" to see which compilation options, packages,
and styles are included in the binary.

Thanks to Christoph Junghans (LANL) for making LAMMPS available in Fedora.

.. _openkim: https://openkim.org

----------

.. _epel:

Pre-built EPEL Linux executable
------------------------------------------

Pre-built LAMMPS (and KIM) packages for stable releases are available
in the `Extra Packages for Enterprise Linux (EPEL) repository <https://fedoraproject.org/wiki/EPEL>`_
for use with Red Hat Enterprise Linux (RHEL) or CentOS version 7.x
and compatible Linux distributions. Names of packages, executable,
and content are the same as described above for Fedora Linux.
But RHEL/CentOS 7.x uses the "yum" package manager instead of "dnf"
in Fedora 28.

Please use "lmp -help" to see which compilation options, packages,
and styles are included in the binary.

Thanks to Christoph Junghans (LANL) for making LAMMPS available in EPEL.

----------

.. _opensuse:

Pre-built OpenSuse Linux executable
--------------------------------------------------

A pre-built LAMMPS package for stable releases is available
in OpenSuse as of Leap 15.0. You can install the package with:

.. code-block:: bash

   $ zypper install lammps

This includes support for OpenMPI. The name of the LAMMPS executable
is *lmp*\ . Thus to run an input in parallel on 2 CPUs you would do:

.. code-block:: bash

   $ mpirun -np 2 lmp -in in.lj

Please use "lmp -help" to see which compilation options, packages,
and styles are included in the binary.

The LAMMPS binary is built with the :ref:`KIM package <kim>` which
results in the above command also installing the `kim-api` binaries when LAMMPS
is installed.  In order to use potentials from `openkim.org <openkim_>`_, you
can install the `openkim-models` package

.. code-block:: bash

   $ zypper install openkim-models

Thanks to Christoph Junghans (LANL) for making LAMMPS available in OpenSuse.

----------

.. _gentoo:

Gentoo Linux executable
------------------------------------

LAMMPS is part of Gentoo's main package tree and can be installed by
typing:

.. code-block:: bash

   % emerge --ask lammps

Note that in Gentoo the LAMMPS source is downloaded and the package is
built on the your machine.

Certain LAMMPS packages can be enable via USE flags, type

.. code-block:: bash

   % equery uses lammps

for details.

Thanks to Nicolas Bock and Christoph Junghans (LANL) for setting up
this Gentoo capability.

----------

.. _arch:

Archlinux build-script
---------------------------------

LAMMPS is available via Arch's unofficial Arch User repository (AUR).
There are three scripts available, named lammps, lammps-beta and lammps-git.
They respectively package the stable, patch and git releases.

To install, you will need to have the git package installed. You may use
any of the above names in-place of lammps.

.. code-block:: bash

   $ git clone https://aur.archlinux.org/lammps.git
   $ cd lammps
   $ makepkg -s
   $ makepkg -i

To update, you may repeat the above, or change into the cloned directory,
and execute the following, after which, if there are any changes, you may
use makepkg as above.

.. code-block:: bash

   $ git pull

Alternatively, you may use an AUR helper to install these packages.

Note that the AUR provides build-scripts that download the source and
the build the package on your machine.
