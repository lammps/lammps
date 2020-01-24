Extending Python to run in parallel
===================================

If you wish to run LAMMPS in parallel from Python, you need to extend
your Python with an interface to MPI.  This also allows you to
make MPI calls directly from Python in your script, if you desire.

We recommend use of mpi4py:

* `mpi4py <https://pypi.org/project/mpi4py/>`_

As of version 2.0.0 it allows passing a custom MPI communicator to
the LAMMPS constructor, which means one can easily run one or more
LAMMPS instances on subsets of the total MPI ranks.

To install mpi4py (version mpi4py-2.0.0 as of Oct 2015), unpack it
and from its main directory, type


.. parsed-literal::

   python setup.py build
   sudo python setup.py install

Again, the "sudo" is only needed if required to copy mpi4py files into
your Python distribution's site-packages directory. To install with
user privilege into the user local directory type


.. parsed-literal::

   python setup.py install --user

If you have successfully installed mpi4py, you should be able to run
Python and type


.. parsed-literal::

   from mpi4py import MPI

without error.  You should also be able to run python in parallel
on a simple test script


.. parsed-literal::

   % mpirun -np 4 python test.py

where test.py contains the lines


.. parsed-literal::

   from mpi4py import MPI
   comm = MPI.COMM_WORLD
   print "Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size())

and see one line of output for each processor you run on.

.. note::

   To use mpi4py and LAMMPS in parallel from Python, you must
   insure both are using the same version of MPI.  If you only have one
   MPI installed on your system, this is not an issue, but it can be if
   you have multiple MPIs.  Your LAMMPS build is explicit about which MPI
   it is using, since you specify the details in your lo-level
   src/MAKE/Makefile.foo file.  Mpi4py uses the "mpicc" command to find
   information about the MPI it uses to build against.  And it tries to
   load "libmpi.so" from the LD\_LIBRARY\_PATH.  This may or may not find
   the MPI library that LAMMPS is using.  If you have problems running
   both mpi4py and LAMMPS together, this is an issue you may need to
   address, e.g. by moving other MPI installations so that mpi4py finds
   the right one.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
