LAMMPS internal utility functions
=================================

The ``utils`` sub-namespace inside the ``LAMMPS_NS`` namespace provides
a collection of convenience functions and utilities that perform common
tasks that are required repeatedly throughout the LAMMPS code like
reading or writing to files with error checking or translation of
strings into specific types of numbers with checking for validity.  This
reduces redundant implementations and encourages consistent behavior.

I/O functions with validity check
---------------------------------

These are wrappers around the corresponding C library calls like
``fgets()`` or ``fread()``.  They will check if there were errors
on reading or an unexpected end-of-file state was reached.  In that
case, the functions will stop the calculation with an error message,
indicating the name of the problematic file, if possible.

.. doxygenfunction:: sfgets
   :project: progguide

.. doxygenfunction:: sfread
   :project: progguide

String to number conversions with validity check
------------------------------------------------

These are functions equivalent to those in the ``Force`` class that
were implemented with the aim to replace and supersede those.  Unlike
the versions in ``Force``, these can be used in cases where only
a single MPI rank is trying to convert strings to numbers, as you
can select through an argument, whether ``Error->all()`` or ``Error->one()``
should be called on improper strings.

These functions are preferred over C library calls like ``atoi()`` or
``atof()`` since they check if the **entire** provided string is a valid
(floating-point or integer) number, and will error out instead of silently
return the result of a partial conversion or zero in cases where the
string is not a valid number.  This behavior allows to more easily detect
typos or issues when processing input files.

.. doxygenfunction:: numeric
   :project: progguide

.. doxygenfunction:: inumeric
   :project: progguide

.. doxygenfunction:: bnumeric
   :project: progguide

.. doxygenfunction:: tnumeric
   :project: progguide


Convenience functions
---------------------

.. doxygenfunction:: strmatch
   :project: progguide

.. doxygenfunction:: check_packages_for_style
   :project: progguide
