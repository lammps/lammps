C++ Class documentation
=======================

LAMMPS is designed to be used as a C++ class library where one can set
up and drive a simulation through creating a class instance and then
calling some abstract operations or commands on that class or its member
class instances.  These are interfaced to the :doc:`C library API
<Library>`, which providing an additional level of abstraction
simplification for common operations. The C API is also the basis for
calling LAMMPS from Python or Fortran.

When used from a C++ program, most of the symbols and functions in
LAMMPS are wrapped into the ``LAMMPS_NS`` namespace so they will not
collide with your own classes or other libraries. This, however, does
not extend to the additional libraries bundled with LAMMPS in the lib
folder and some of the low-level code of some packages.

Behind the scenes this is implemented through inheritance and
polymorphism where base classes define the abstract interface and
derived classes provide the specialized implementation for specific
models or optimizations or ports to accelerator platforms.  This
document will provide an outline of the fundamental classes and their
purposes and important member functions.

A more high-level overview and a graphical representation is on the
:doc:`page discussing the LAMMPS class hierarchy <Developer_org>`.

.. note::

   Please see the :ref:`note about thread-safety <thread-safety>`
   in the library Howto doc page.

-----------------------------------

.. toctree::
   :caption: Core LAMMPS Classes
   :name: lammpscore
   :maxdepth: 1

   Classes_lammps

-----------------------------------

The following classes provide the core infrastructure of LAMMPS and
pointers to instances of them are managed by the :cpp:class:`LAMMPS
class <LAMMPS_NS::LAMMPS>`.  Like most classes in LAMMPS they are
derived from from the :cpp:class:`Pointers class <LAMMPS_NS::Pointers>`
so they have convenient access to each other's data and methods.

.. toctree::
   :caption: Infrastructure Classes
   :name: lammpsinfra
   :maxdepth: 1

   Classes_atom
   Classes_comm
   Classes_domain
   Classes_error
   Classes_force
   Classes_group
   Classes_input
   Classes_modify
   Classes_output
   Classes_timer
   Classes_universe
   Classes_update

-----------------------------------

.. toctree::
   :caption: Utility Classes
   :name: lammpsutils
   :maxdepth: 1

   Classes_cite

-----------------------------------

.. toctree::
   :caption: Style Base Classes
   :name: lammpsbase
   :maxdepth: 1


