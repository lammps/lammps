Download an executable for Mac
==============================

LAMMPS can be downloaded, built, and configured for OS X on a Mac with
`Homebrew <homebrew_>`_.  The following LAMMPS packages are unavailable at this
time because of additional needs not yet met: GPU, KOKKOS, LATTE, MSCG,
MESSAGE, MPIIO POEMS VORONOI.

After installing Homebrew, you can install LAMMPS on your system with
the following commands:


.. parsed-literal::

   % brew install lammps

This will install the executables "lammps\_serial" and "lammps\_mpi", as well as
the LAMMPS "doc", "potentials", "tools", "bench", and "examples" directories.

Once LAMMPS is installed, you can test the installation with the
Lennard-Jones benchmark file:


.. parsed-literal::

   % brew test lammps -v

The LAMMPS binary is built with the :ref:`KIM package <kim>` which
results in Homebrew also installing the `kim-api` binaries when LAMMPS is
installed.  In order to use potentials from `openkim.org <openkim_>`_, you can
install the `openkim-models` package


.. parsed-literal::

   % brew install openkim-models

If you have problems with the installation you can post issues to
`this link <homebrew_>`_.

.. _homebrew: https://github.com/Homebrew/homebrew-core/issues

Thanks to Derek Thomas (derekt at cello.t.u-tokyo.ac.jp) for setting
up the Homebrew capability.


.. _openkim: https://openkim.org




.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
