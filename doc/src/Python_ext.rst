Extending the library and Python interface
******************************************

As noted above, these Python class methods correspond one-to-one with
the functions in the LAMMPS library interface in src/library.cpp and
library.h.  This means you can extend the Python wrapper via the
following steps:

* Add a new interface function to src/library.cpp and
  src/library.h.
* Rebuild LAMMPS as a shared library.
* Add a wrapper method to python/lammps.py for this interface
  function.
* You should now be able to invoke the new interface function from a
  Python script.

