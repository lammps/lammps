Library interface utility functions
===================================

To simplify some of the tasks, the library interface contains
some utility functions that are not directly calling LAMMPS:

- :cpp:func:`lammps_encode_image_flags`
- :cpp:func:`lammps_decode_image_flags`
- :cpp:func:`lammps_set_fix_external_callback`
- :cpp:func:`lammps_fix_external_set_energy_global`
- :cpp:func:`lammps_fix_external_set_virial_global`
- :cpp:func:`lammps_has_error`
- :cpp:func:`lammps_get_last_error_message`

-----------------------

.. doxygenfunction:: lammps_encode_image_flags
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_decode_image_flags(int image, int *flags)
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*)
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_energy_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_virial_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_has_error
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_last_error_message
   :project: progguide
