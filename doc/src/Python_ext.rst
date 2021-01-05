Extending the Python interface
******************************************

As noted previously, most of the :py:class:`lammps <lammps.lammps>`
Python class methods correspond one-to-one with the functions in the
LAMMPS library interface in ``src/library.cpp`` and ``library.h``.
This means you can extend the Python wrapper by following these steps:

* Add a new interface function to ``src/library.cpp`` and
  ``src/library.h``.
* Rebuild LAMMPS as a shared library.
* Add a wrapper method to ``python/lammps/core.py`` for this interface
  function.
* Define the corresponding ``argtypes`` list and ``restype``
  in the ``lammps.__init__()`` function.
* Re-install the shared library and the python module, if needed
* You should now be able to invoke the new interface function from a
  Python script.


