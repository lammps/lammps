Running LAMMPS on Windows
=========================

To run a serial (non-MPI) executable, follow these steps:

* Get a command prompt by going to Start->Run... ,
  then typing "cmd".
* Move to the directory where you have your input script,
  (e.g. by typing: cd "Documents").
* At the command prompt, type "lmp -in in.file", where
  in.file is the name of your LAMMPS input script.

Note that the serial executable includes support for multi-threading
parallelization from the styles in the OPENMP packages.  To run with
4 threads, you can type this:

.. code-block:: bash

   lmp -in in.lj -pk omp 4 -sf omp

----------

For the MPI executable, which allows you to run LAMMPS under Windows
in parallel, follow these steps.

Download and install a compatible MPI library binary package:

* for 32-bit Windows: `mpich2-1.4.1p1-win-ia32.msi <http://download.lammps.org/thirdparty/mpich2-1.4.1p1-win-ia32.msi>`_
* for 64-bit Windows: `mpich2-1.4.1p1-win-x86-64.msi <http://download.lammps.org/thirdparty/mpich2-1.4.1p1-win-x86-64.msi>`_

The LAMMPS Windows installer packages will automatically adjust your
path for the default location of this MPI package. After the
installation of the MPICH2 software, it needs to be integrated into
the system.  For this you need to start a Command Prompt in
*Administrator Mode* (right click on the icon and select it). Change
into the MPICH2 installation directory, then into the sub-directory
**bin** and execute **smpd.exe -install**\ . Exit the command window.

* Get a new, regular command prompt by going to Start->Run... ,
  then typing "cmd".
* Move to the directory where you have your input file
  (e.g. by typing: cd "Documents").

Then you can run the executable in serial like in the example above
or in parallel using MPI with one of the following commands:

.. code-block:: bash

   mpiexec -localonly 4 lmp -in in.file
   mpiexec -np 4 lmp -in in.file

where in.file is the name of your LAMMPS input script. For the latter
case, you may be prompted to enter the password that you set during
installation of the MPI library software.

In this mode, output may not immediately show up on the screen, so if
your input script takes a long time to execute, you may need to be
patient before the output shows up.

The parallel executable can also run on a single processor by typing
something like this:

.. code-block:: bash

   lmp -in in.lj

Note that the parallel executable also includes OpenMP
multi-threading, which can be combined with MPI using something like:

.. code-block:: bash

   mpiexec -localonly 2 lmp -in in.lj -pk omp 2 -sf omp

