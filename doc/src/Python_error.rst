LAMMPS error handling in Python
*******************************

Compiling the shared library with :ref:`C++ exception support <exceptions>` provides a better error
handling experience. Without exceptions the LAMMPS code will terminate the
current Python process with an error message.  C++ exceptions allow capturing
them on the C++ side and rethrowing them on the Python side. This way
LAMMPS errors can be handled through the Python exception handling mechanism.

.. warning::

   Capturing a LAMMPS exception in Python can still mean that the
   current LAMMPS process is in an illegal state and must be terminated. It is
   advised to save your data and terminate the Python instance as quickly as
   possible.
