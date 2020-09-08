LAMMPS Class
************

The LAMMPS class is encapsulating an MD simulation state and thus it is
the class that needs to be created when starting a new simulation system
state.  The LAMMPS executable essentially creates one instance of this
class and passes the command line flags and tells it to process the
provided input (a file or ``stdin``).  It shuts the class down when
control is returned to it and then exits.  When using LAMMPS as a
library from another code it is required to create an instance of this
class, either directly from C++ with ``new LAMMPS()`` or through one
of the library interface functions like :cpp:func:`lammps_open` of the
C-library interface, or the :py:class:`lammps.lammps` class constructor
of the Python module, or the :f:func:`lammps` constructor of the Fortran
module.

--------------------

.. doxygenclass:: LAMMPS_NS::LAMMPS
   :project: progguide
   :members:

.. doxygenfile:: lammps.cpp
   :project: progguide
