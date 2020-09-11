Adding code to the Library interface
====================================

The functionality of the LAMMPS library interface has historically
always been motivated by the needs of its users and functions were
added or expanded as they were needed and used.  Contributions to
the interface are always welcome.  However with a refactoring of
the library interface and its documentation that started in 2020,
there are now a few requirements for inclusion of changes.

  - New functions should be orthogonal to existing ones and not
    implement functionality that can already be achieved with the
    existing APIs.
  - All changes and additions should be documented with
    `Doxygen <https://doxygen.nl>`_ style comments and references
    to those functions added to the corresponding files in the
    ``doc/src`` folder.
  - If possible, new unit tests to test those new features should
    be added.
  - The new feature should also be implemented and documented for
    the Python and Fortran modules.
  - All additions should work and be compatible with ``-DLAMMPS_BIGBIG``,
    ``-DLAMMPS_SMALLBIG``, ``-DLAMMPS_SMALLSMALL`` and compiling
    with and without MPI support.
  - The ``library.h`` file should be kept compatible to C code at
    a level similar to C89. Its interfaces may not reference any
    custom data types (e.g. ``bigint``, ``tagint``, and so on) only
    known inside of LAMMPS.
  - only C style comments, not C++ style

Please note, that these are *not* *strict* requirements, but the
LAMMPS developers appreciate if they are followed closely and will
assist with implementing what is missing.
