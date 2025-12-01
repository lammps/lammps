Error Class
***********

The Error class provides centralized error handling for LAMMPS.  It supports
fatal errors that terminate execution on all or single MPI processes, as well
as warnings that are logged but allow execution to continue.

Error messages include the source file name and line number where the error
occurred.  When an error originates from parsing a LAMMPS input script command,
the class can optionally indicate which argument caused the error.

The class tracks the number of warnings issued and supports suppression after
a maximum count is reached.  The last error message is stored for retrieval,
which is useful for the library interface and testing.

.. doxygenclass:: LAMMPS_NS::Error
   :project: progguide
   :members:
