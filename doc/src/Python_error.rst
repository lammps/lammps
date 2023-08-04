Handling LAMMPS errors
**********************

The shared library is compiled with :ref:`C++ exception support
<exceptions>` to provide a better error handling experience.  C++
exceptions allow capturing errors on the C++ side and rethrowing them on
the Python side.  This way LAMMPS errors can be handled through the
Python exception handling mechanism.

.. code-block:: python

   from lammps import lammps, MPIAbortException

   lmp = lammps()

   try:
      # LAMMPS will normally terminate itself and the running process if an error
      # occurs. This would kill the Python interpreter. To avoid this, make sure to
      # compile with LAMMPS_EXCEPTIONS enabled. This ensures the library API calls
      # will not terminate the parent process. Instead, the library wrapper will
      # detect that an error has occured and throw a Python exception

      lmp.command('unknown')
   except MPIAbortException as ae:
      # Single MPI process got killed. This would normally be handled by an MPI abort
      pass
   except Exception as e:
      # All (MPI) processes have reached this error
      pass

.. warning::

   Capturing a LAMMPS exception in Python can still mean that the
   current LAMMPS process is in an illegal state and must be
   terminated. It is advised to save your data and terminate the Python
   instance as quickly as possible when running in parallel with MPI.
