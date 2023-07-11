Download an executable for macOS
--------------------------------

LAMMPS can be downloaded, built, and configured for macOS with `Homebrew
<homebrew_>`_.  (Alternatively, see the installation instructions for
:doc:`downloading an executable via Conda <Install_conda>`.)  The
following LAMMPS packages are unavailable at this time because of
additional requirements not yet met: GPU, KOKKOS, MSCG, MPIIO, POEMS,
VORONOI.

After installing Homebrew, you can install LAMMPS on your system with
the following commands:

.. code-block:: bash

   brew install lammps

This will install the executables "lammps_serial" and "lammps_mpi", as
well as the LAMMPS "doc", "potentials", "tools", "bench", and "examples"
directories.

Once LAMMPS is installed, you can test the installation with the
Lennard-Jones benchmark file:

.. code-block:: bash

   brew test lammps -v

The LAMMPS binary is built with the :ref:`KIM package <kim>`, which
results in Homebrew also installing the `kim-api` binaries when LAMMPS
is installed.  In order to use potentials from `openkim.org
<openkim_>`_, you can install the `openkim-models` package

.. code-block:: bash

   brew install openkim-models

If you have problems with the installation, you can post issues to
`this link <https://github.com/Homebrew/homebrew-core/issues>`_.

.. _homebrew: https://brew.sh

Thanks to Derek Thomas (derekt at cello.t.u-tokyo.ac.jp) for setting
up the Homebrew capability.

.. _openkim: https://openkim.org
