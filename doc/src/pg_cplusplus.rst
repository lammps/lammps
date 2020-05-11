LAMMPS C++ API Usage
********************

When using the LAMMPS library interfaces, the core
task is to create an instance of the :cpp:class:`LAMMPS_NS::LAMMPS` class.
In C++ this can be done directly through the ``new`` operator.
All further operations are then initiated through calling member functions
of some of the components of the LAMMPS class.
