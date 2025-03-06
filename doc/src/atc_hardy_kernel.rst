.. index:: fix_modify AtC kernel

fix_modify AtC kernel command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> kernel <type> <parameters>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* kernel = name of the AtC sub-command
* type = *step* or *cell* or *cubic_bar* or *cubic_cylinder* or
  *cubic_sphere* or *quartic_bar* or *quartic_cylinder* or
  *quartic_sphere*
* the following parameter(s) are required for each kernel:

  - *step* : <radius>
  - *cell* : <hx> <hy> <hz> or <h>
  - *cubic_bar* : <half_width>
  - *cubic_cylinder* : <radius>
  - *cubic_sphere* : <radius>
  - *quartic_bar* : <half_width>
  - *quartic_cylinder* : <radius>
  - *quartic_sphere* : <radius>



Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC kernel cell 1.0 1.0 1.0
   fix_modify AtC kernel quartic_sphere 10.0


Description
"""""""""""

Sets the localization kernel type and parameters for :doc:`fix atc hardy <fix_atc>`.

Restrictions
""""""""""""

Must be used with :doc:`fix atc hardy <fix_atc>`.  For bar kernel types,
half-width oriented along x-direction.  For cylinder kernel types,
cylindrical axis is assumed to be in z-direction.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC fields <atc_hardy_fields>`
- :doc:`fix_modify AtC gradients <atc_hardy_gradients>`
- :doc:`fix_modify AtC rates <atc_hardy_rates>`
- :doc:`fix_modify AtC computes <atc_hardy_computes>`

Default
"""""""

None.
