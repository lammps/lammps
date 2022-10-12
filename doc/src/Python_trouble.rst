Troubleshooting
***************

Testing if Python can launch LAMMPS
===================================

To test if LAMMPS is callable from Python, launch Python interactively
and type:

.. code-block:: python

   >>> from lammps import lammps
   >>> lmp = lammps()

If you get no errors, you're ready to use LAMMPS from Python.  If the
second command fails, the most common error to see is

.. code-block:: bash

   OSError: Could not load LAMMPS dynamic library

which means Python was unable to load the LAMMPS shared library.  This
typically occurs if the system can't find the LAMMPS shared library or
one of the auxiliary shared libraries it depends on, or if something
about the library is incompatible with your Python.  The error message
should give you an indication of what went wrong.

If your shared library uses a suffix, such as ``liblammps_mpi.so``, change
the constructor call as follows (see :ref:`python_create_lammps` for more details):

.. code-block:: python

   >>> lmp = lammps(name='mpi')

You can also test the load directly in Python as follows, without
first importing from the ``lammps`` module:

.. code-block:: python

   >>> from ctypes import CDLL
   >>> CDLL("liblammps.so")

If an error occurs, carefully go through the steps in :ref:`python_install_guides` and on the
:doc:`Build_basics <Build_basics>` page about building a shared library.
