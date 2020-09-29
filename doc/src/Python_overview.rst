Overview
========

The LAMMPS distribution includes a python directory with all you need
to run LAMMPS from Python.  The python/lammps.py file wraps the LAMMPS
library interface, with one wrapper function per LAMMPS library
function.  This file makes it is possible to do the following either
from a Python script, or interactively from a Python prompt: create
one or more instances of LAMMPS, invoke LAMMPS commands or give it an
input script, run LAMMPS incrementally, extract LAMMPS results, an
modify internal LAMMPS variables.  From a Python script you can do
this in serial or parallel.  Running Python interactively in parallel
does not generally work, unless you have a version of Python that
extends Python to enable multiple instances of Python to read what you
type.

To do all of this, you must first build LAMMPS as a shared library,
then insure that your Python can find the python/lammps.py file and
the shared library.

The Python wrapper for LAMMPS uses the "ctypes" package in Python,
which auto-generates the interface code needed between Python and a
set of C-style library functions.  Ctypes is part of standard Python
for versions 2.5 and later.  You can check which version of Python you
have by simply typing "python" at a shell prompt.

.. warning:: Python 2 support is deprecated

   While the LAMMPS Python module was originally developed to support both
   Python 2 and 3, Python 2 is no longer maintained as of `January 1, 2020 <https://www.python.org/doc/sunset-python-2/>`_.
   Therefore, we will no longer backport any new features to Python 2 and
   highly recommend using Python versions 3.6+.

---------

LAMMPS can work together with Python in three ways.  First, Python can
wrap LAMMPS through the its :doc:`library interface <Howto_library>`, so
that a Python script can create one or more instances of LAMMPS and
launch one or more simulations.  In Python lingo, this is called
"extending" Python with a LAMMPS module.

.. figure:: JPG/python-invoke-lammps.png
   :figclass: align-center

   Launching LAMMPS via Python


Second, the lower-level Python interface can be used indirectly through
the provided :code`PyLammps` and :code:`IPyLammps` wrapper classes, written in Python.
These wrappers try to simplify the usage of LAMMPS in Python by
providing an object-based interface to common LAMMPS functionality.
They also reduces the amount of code necessary to parameterize LAMMPS
scripts through Python and make variables and computes directly
accessible.

.. figure:: JPG/pylammps-invoke-lammps.png
   :figclass: align-center

   Using the PyLammps / IPyLammps wrappers

Third, LAMMPS can use the Python interpreter, so that a LAMMPS
input script or styles can invoke Python code directly, and pass
information back-and-forth between the input script and Python
functions you write.  This Python code can also callback to LAMMPS
to query or change its attributes through the LAMMPS Python module
mentioned above.  In Python lingo, this is "embedding" Python in
LAMMPS.  When used in this mode, Python can perform script operations
that the simple LAMMPS input script syntax can not.

.. figure:: JPG/lammps-invoke-python.png
   :figclass: align-center

   Calling Python code from LAMMPS
