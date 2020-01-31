Overview of Python and LAMMPS
=============================

LAMMPS can work together with Python in three ways.  First, Python can
wrap LAMMPS through the its :doc:`library interface <Howto_library>`, so
that a Python script can create one or more instances of LAMMPS and
launch one or more simulations.  In Python lingo, this is called
"extending" Python with a LAMMPS module.

Second, a lower-level Python interface can be used indirectly through
the provided PyLammps and IPyLammps wrapper classes, written in Python.
These wrappers try to simplify the usage of LAMMPS in Python by
providing an object-based interface to common LAMMPS functionality.
They also reduces the amount of code necessary to parameterize LAMMPS
scripts through Python and make variables and computes directly
accessible.

Third, LAMMPS can use the Python interpreter, so that a LAMMPS
input script or styles can invoke Python code directly, and pass
information back-and-forth between the input script and Python
functions you write.  This Python code can also callback to LAMMPS
to query or change its attributes through the LAMMPS Python module
mentioned above.  In Python lingo, this is "embedding" Python in
LAMMPS.  When used in this mode, Python can perform script operations
that the simple LAMMPS input script syntax can not.


