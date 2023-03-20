Utility functions
=================

To simplify some tasks, the library interface contains these utility
functions.  They do not directly call the LAMMPS library.

- :cpp:func:`lammps_encode_image_flags`
- :cpp:func:`lammps_decode_image_flags`
- :cpp:func:`lammps_set_fix_external_callback`
- :cpp:func:`lammps_fix_external_get_force`
- :cpp:func:`lammps_fix_external_set_energy_global`
- :cpp:func:`lammps_fix_external_set_energy_peratom`
- :cpp:func:`lammps_fix_external_set_virial_global`
- :cpp:func:`lammps_fix_external_set_virial_peratom`
- :cpp:func:`lammps_fix_external_set_vector_length`
- :cpp:func:`lammps_fix_external_set_vector`
- :cpp:func:`lammps_flush_buffers`
- :cpp:func:`lammps_free`
- :cpp:func:`lammps_is_running`
- :cpp:func:`lammps_force_timeout`
- :cpp:func:`lammps_has_error`
- :cpp:func:`lammps_get_last_error_message`
- :cpp:func:`lammps_python_api_version`

The :cpp:func:`lammps_free` function is a clean-up function to free
memory that the library had allocated previously via other function
calls.  Look for notes in the descriptions of the individual commands
where such memory buffers were allocated that require the use of
:cpp:func:`lammps_free`.

-----------------------

.. doxygenfunction:: lammps_encode_image_flags
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_decode_image_flags(int image, int *flags)
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_set_fix_external_callback(void *, const char *, FixExternalFnPtr, void*)
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_get_force
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_energy_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_energy_peratom
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_virial_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_virial_peratom
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_vector_length
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_vector
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_flush_buffers
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_free
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_is_running
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_force_timeout
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_has_error
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_last_error_message
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_python_api_version
   :project: progguide

