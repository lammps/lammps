Overview
========

The LAMMPS distribution includes a ``python`` directory with the Python
code needed to run LAMMPS from Python.  The ``python/lammps`` package
contains :doc:`the "lammps" Python module <Python_module>` that wraps the
LAMMPS C-library interface.  This module makes it is possible to do the
following either from a Python script, or interactively from a Python
prompt:

- create one or more instances of LAMMPS
- invoke LAMMPS commands or read them from an input script
- run LAMMPS incrementally
- extract LAMMPS results
- and modify internal LAMMPS data structures.

From a Python script you can do this in serial or in parallel.  Running
Python interactively in parallel does not generally work, unless you
have a version of Python that extends Python to enable multiple
instances of Python to read what you type.

To do all of this, you must build LAMMPS in :ref:`"shared" mode <exe>`
and make certain that your Python interpreter can find the ``lammps``
Python package and the LAMMPS shared library file.

.. _ctypes: https://docs.python.org/3/library/ctypes.html

The Python wrapper for LAMMPS uses the `ctypes <ctypes_>`_ package in
Python, which auto-generates the interface code needed between Python
and a set of C-style library functions.  Ctypes has been part of the
standard Python distribution since version 2.5.  You can check which
version of Python you have by simply typing "python" at a shell prompt.
Below is an example output for Python version 3.8.5.

.. code-block:: console

   $ python
   Python 3.8.5 (default, Aug 12 2020, 00:00:00)
   [GCC 10.2.1 20200723 (Red Hat 10.2.1-1)] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>>


.. warning::

   The options described in this section of the manual for using Python
   with LAMMPS support only Python 3.6 or later.  For use with Python
   2.x you will need to use an older LAMMPS version like 29 Aug 2024
   or older.  If you notice Python code in the LAMMPS distribution that
   is not compatible with Python 3, please contact the LAMMPS developers
   or submit `and issue on GitHub <https://github.com/lammps/lammps/issues>`_

---------

LAMMPS can work together with Python in two ways.  First, Python can
wrap LAMMPS through the its :doc:`library interface <Library>`, so
that a Python script can create one or more instances of LAMMPS and
launch one or more simulations.  In Python terms, this is referred to as
"extending" Python with a LAMMPS module.

.. figure:: JPG/python-invoke-lammps.png
   :figclass: align-center

   Launching LAMMPS via Python

Second, LAMMPS can use the Python interpreter, so that a LAMMPS input
script or styles can invoke Python code directly, and pass information
back-and-forth between the input script and Python functions you write.
This Python code can also call back to LAMMPS to query or change its
attributes through the LAMMPS Python module mentioned above.  In Python
terms, this is called "embedding" Python into LAMMPS.  When used in this
mode, Python can perform script operations that the simple LAMMPS input
script syntax can not.

.. figure:: JPG/lammps-invoke-python.png
   :figclass: align-center

   Calling Python code from LAMMPS
