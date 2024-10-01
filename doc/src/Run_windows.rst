Running LAMMPS on Windows
=========================

To run a serial (non-MPI) executable, follow these steps:

* Install a LAMMPS installer package from https://packages.lammps.org/windows.html
* Open the "Command Prompt" or "Terminal" app.
* Change to the directory where you have your input script,
  (e.g. by typing: cd "Documents").
* At the command prompt, type "lmp -in in.file.lmp", where
  ``in.file.lmp`` is the name of your LAMMPS input script.

Note that the serial executable includes support for multi-threading
parallelization from the styles in the OPENMP and KOKKOS packages.
To run with 4 threads, you can type this:

.. code-block:: bash

   lmp -in in.lj.lmp -pk omp 4 -sf omp
   lmp -in in.lj.lmp -k on t 4 -sf kk

Alternately, you can also install a package with LAMMPS-GUI included and
open the LAMMPS-GUI app (the package includes the command line version
of LAMMPS as well) and open the input file in the GUI and run it from
there.  For details on LAMMPS-GUI, see :doc:`Howto_lammps_gui`.

----------

For the MS-MPI executables, which allow you to run LAMMPS under Windows
in parallel using MPI rather than multi-threading, follow these steps.

Download and install the MS-MPI runtime package ``msmpisetup.exe`` from
https://www.microsoft.com/en-us/download/details.aspx?id=105289 (Note
that the ``msmpisdk.msi`` is **only** required for **compilation** of
LAMMPS from source on Windows using Microsoft Visual Studio). After
installation of MS-MPI perform a reboot.

Then you can run the executable in serial like in the example above
or in parallel using MPI with one of the following commands:

.. code-block:: bash

   mpiexec -localonly 4 lmp -in in.file.lmp
   mpiexec -np 4 lmp -in in.file.lmp

where ``in.file.lmp`` is the name of your LAMMPS input script. For the
latter case, you may be prompted to enter the password that you set
during installation of the MPI library software.

In this mode, output may not immediately show up on the screen, so if
your input script takes a long time to execute, you may need to be
patient before the output shows up.

Note that the parallel executable also includes OpenMP multi-threading
through both the OPENMP and the KOKKOS package, which can be combined
with MPI using something like:

.. code-block:: bash

   mpiexec -localonly 2 lmp -in in.lj.lmp -pk omp 2 -sf omp
   mpiexec -localonly 2 lmp -in in.lj.lmp -kokkos on t 2 -sf kk

-------------

MPI parallelization will work for *all* functionality in LAMMPS and in
many cases the MPI parallelization is more efficient than
multi-threading since LAMMPS was designed from ground up for MPI
parallelization using domain decomposition.  Multi-threading is only
available for selected styles and implemented on top of the MPI
parallelization.  Multi-threading is most useful for systems with large
load imbalances when using domain decomposition and a smaller number
of threads (<= 8).
