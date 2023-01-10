Running LAMMPS and Python in serial
-----------------------------------

To run a LAMMPS in serial, type these lines into Python
interactively from the ``bench`` directory:

.. code-block:: python

   from lammps import lammps
   lmp = lammps()
   lmp.file("in.lj")

Or put the same lines in the file ``test.py`` and run it as

.. code-block:: bash

   python3 test.py

Either way, you should see the results of running the ``in.lj`` benchmark
on a single processor appear on the screen, the same as if you had
typed something like:

.. code-block:: bash

   lmp_serial -in in.lj

Running LAMMPS and Python in parallel with MPI
----------------------------------------------

To run LAMMPS in parallel, assuming you have installed the
`mpi4py <https://mpi4py.readthedocs.io>`_ package as discussed
:ref:`python_install_mpi4py`, create a ``test.py`` file containing these lines:

.. code-block:: python

   from mpi4py import MPI
   from lammps import lammps
   lmp = lammps()
   lmp.file("in.lj")
   me = MPI.COMM_WORLD.Get_rank()
   nprocs = MPI.COMM_WORLD.Get_size()
   print("Proc %d out of %d procs has" % (me,nprocs),lmp)
   MPI.Finalize()

You can run the script in parallel as:

.. code-block:: bash

   mpirun -np 4 python3 test.py

and you should see the same output as if you had typed

.. code-block:: bash

   mpirun -np 4 lmp_mpi -in in.lj

Note that without the mpi4py specific lines from ``test.py``

.. code-block:: python

   from lammps import lammps
   lmp = lammps()
   lmp.file("in.lj")

running the script with ``mpirun`` on :math:`P` processors would lead to
:math:`P` independent simulations to run parallel, each with a single
processor. Therefore, if you use the mpi4py lines and you see multiple LAMMPS
single processor outputs, mpi4py is not working correctly.

Also note that once you import the mpi4py module, mpi4py initializes MPI
for you, and you can use MPI calls directly in your Python script, as
described in the mpi4py documentation.  The last line of your Python
script should be ``MPI.finalize()``, to insure MPI is shut down
correctly.


Running Python scripts
----------------------

Note that any Python script (not just for LAMMPS) can be invoked in
one of several ways:

.. code-block:: bash

   python script.py
   python -i script.py
   ./script.py

The last command requires that the first line of the script be
something like this:

.. code-block:: bash

   #!/usr/bin/python

or

.. code-block:: bash

   #!/usr/bin/env python

where the path in the first case needs to point to where you have Python
installed (the second option is workaround for when this may change),
and that you have made the script file executable:

.. code-block:: bash

   chmod +x script.py

Without the ``-i`` flag, Python will exit when the script finishes.
With the ``-i`` flag, you will be left in the Python interpreter when
the script finishes, so you can type subsequent commands.  As mentioned
above, you can only run Python interactively when running Python on a
single processor, not in parallel.
