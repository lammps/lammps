Test the Python/LAMMPS interface
================================

To test if LAMMPS is callable from Python, launch Python interactively
and type:


.. parsed-literal::

   >>> from lammps import lammps
   >>> lmp = lammps()

If you get no errors, you're ready to use LAMMPS from Python.  If the
2nd command fails, the most common error to see is


.. parsed-literal::

   OSError: Could not load LAMMPS dynamic library

which means Python was unable to load the LAMMPS shared library.  This
typically occurs if the system can't find the LAMMPS shared library or
one of the auxiliary shared libraries it depends on, or if something
about the library is incompatible with your Python.  The error message
should give you an indication of what went wrong.

You can also test the load directly in Python as follows, without
first importing from the lammps.py file:


.. parsed-literal::

   >>> from ctypes import CDLL
   >>> CDLL("liblammps.so")

If an error occurs, carefully go through the steps on the
:doc:`Build\_basics <Build_basics>` doc page about building a shared
library and the :doc:`Python\_install <Python_install>` doc page about
insuring Python can find the necessary two files it needs.

**Test LAMMPS and Python in serial:**
-------------------------------------

To run a LAMMPS test in serial, type these lines into Python
interactively from the bench directory:


.. parsed-literal::

   >>> from lammps import lammps
   >>> lmp = lammps()
   >>> lmp.file("in.lj")

Or put the same lines in the file test.py and run it as


.. parsed-literal::

   % python test.py

Either way, you should see the results of running the in.lj benchmark
on a single processor appear on the screen, the same as if you had
typed something like:


.. parsed-literal::

   lmp_g++ -in in.lj

**Test LAMMPS and Python in parallel:**
---------------------------------------

To run LAMMPS in parallel, assuming you have installed the
`PyPar <https://github.com/daleroberts/pypar>`_ package as discussed
above, create a test.py file containing these lines:


.. parsed-literal::

   import pypar
   from lammps import lammps
   lmp = lammps()
   lmp.file("in.lj")
   print "Proc %d out of %d procs has" % (pypar.rank(),pypar.size()),lmp
   pypar.finalize()

To run LAMMPS in parallel, assuming you have installed the
`mpi4py <https://bitbucket.org/mpi4py/mpi4py>`_ package as discussed
above, create a test.py file containing these lines:


.. parsed-literal::

   from mpi4py import MPI
   from lammps import lammps
   lmp = lammps()
   lmp.file("in.lj")
   me = MPI.COMM_WORLD.Get_rank()
   nprocs = MPI.COMM_WORLD.Get_size()
   print "Proc %d out of %d procs has" % (me,nprocs),lmp
   MPI.Finalize()

You can either script in parallel as:


.. parsed-literal::

   % mpirun -np 4 python test.py

and you should see the same output as if you had typed


.. parsed-literal::

   % mpirun -np 4 lmp_g++ -in in.lj

Note that if you leave out the 3 lines from test.py that specify PyPar
commands you will instantiate and run LAMMPS independently on each of
the P processors specified in the mpirun command.  In this case you
should get 4 sets of output, each showing that a LAMMPS run was made
on a single processor, instead of one set of output showing that
LAMMPS ran on 4 processors.  If the 1-processor outputs occur, it
means that PyPar is not working correctly.

Also note that once you import the PyPar module, PyPar initializes MPI
for you, and you can use MPI calls directly in your Python script, as
described in the PyPar documentation.  The last line of your Python
script should be pypar.finalize(), to insure MPI is shut down
correctly.

**Running Python scripts:**
---------------------------

Note that any Python script (not just for LAMMPS) can be invoked in
one of several ways:


.. parsed-literal::

   % python foo.script
   % python -i foo.script
   % foo.script

The last command requires that the first line of the script be
something like this:


.. parsed-literal::

   #!/usr/local/bin/python
   #!/usr/local/bin/python -i

where the path points to where you have Python installed, and that you
have made the script file executable:


.. parsed-literal::

   % chmod +x foo.script

Without the "-i" flag, Python will exit when the script finishes.
With the "-i" flag, you will be left in the Python interpreter when
the script finishes, so you can type subsequent commands.  As
mentioned above, you can only run Python interactively when running
Python on a single processor, not in parallel.


