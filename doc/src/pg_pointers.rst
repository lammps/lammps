LAMMPS Pointers Base Class
**************************

The Pointers class is the direct or indirect base class for most classes
in LAMMPS.  It contains references to most of the members of the LAMMPS
class so they are accessible for all classes derived from Pointers.  If
you want to add a new class to LAMMPS that is not derived from one of
the base classes for styles, you probably want to derive it from
Pointers.  Because the file is so ubiquituos it also contains some
defines and macros that are intended to be globally available in LAMMPS
classes.

--------------------

.. doxygenclass:: LAMMPS_NS::Pointers
   :project: progguide
   :members:
   :protected-members:


