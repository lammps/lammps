LAMMPS C++ base classes
***********************

LAMMPS was designed to be used as a C++ class library
where one can set up and drive a simulation with some
class instance and then some abstract operations and
commands.  These are also interfaced to the
:doc:`C library API <pg_library>`, which is the basis
for calling LAMMPS from Python or Fortran.  Behind the
scenes this is implemented through inheritance and
polymorphism where base classes define the abstract
interface and derived classes provide the specialized
implementation for specific models or optimizations or
ports to accelerator platforms.  This document will
provide an outline of the fundamental class hierarchy
and some selected examples for derived classes of
specific models.

.. note:: Thread-safety

   Please see the :ref:`note about thread-safety <thread-safety>`
   in the library Howto doc page.

.. toctree::
   :caption: Individual Base Classes
   :name: lammpsbase

   pg_region
