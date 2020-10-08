Extending the C API
===================

The functionality of the LAMMPS library interface has historically
been motivated by the needs of its users.  Functions have been added
or expanded as they were needed and used.  Contributions to the
interface are always welcome.  However with a refactoring of the
library interface and its documentation that started in Spring 2020,
there are now a few requirements for including new changes or extensions.

  - New functions should be orthogonal to existing ones and not
    implement functionality that can already be achieved with the
    existing APIs.
  - All changes and additions should be documented with
    `Doxygen <https://doxygen.nl>`_ style comments and references
    to those functions added to the corresponding files in the
    ``doc/src`` folder.
  - If possible, new unit tests to test those new features should
    be added.
  - New features should also be implemented and documented not just
    for the C interface, but also the Python and Fortran interfaces.
  - All additions should work and be compatible with ``-DLAMMPS_BIGBIG``,
    ``-DLAMMPS_SMALLBIG``, ``-DLAMMPS_SMALLSMALL`` as well as when
    compiling with and without MPI support.
  - The ``library.h`` file should be kept compatible to C code at
    a level similar to C89. Its interfaces may not reference any
    custom data types (e.g. ``bigint``, ``tagint``, and so on) that
    are only known inside of LAMMPS; instead ``int`` and ``int64_t``
    should be used.
  - only use C style comments, not C++ style

Please note that these are not **strict** requirements, but the LAMMPS
developers very much appreciate, if they are followed and can assist
with implementing what is missing.  It helps maintaining the code base
and keeping it consistent.

