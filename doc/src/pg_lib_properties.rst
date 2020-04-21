Retrieving or setting LAMMPS system properties
==============================================

The library interface allows to extract all kinds of information
about the active simulation instance and also modify it.  This
allows to combine MD simulation steps with other processing and
simulation methods computed in the calling code or another code
that is coupled to LAMMPS via the library interface.

-----------------------

.. doxygenfunction:: lammps_get_natoms
   :project: progguide

-------------------

.. doxygenfunction:: lammps_extract_box
   :project: progguide

-------------------

.. doxygenfunction:: lammps_reset_box
   :project: progguide

-------------------

.. doxygenfunction:: lammps_get_thermo
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_setting
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_atom
   :project: progguide
