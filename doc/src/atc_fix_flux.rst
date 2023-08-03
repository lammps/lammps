.. index:: fix_modify AtC fix_flux

fix_modify AtC fix_flux command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> fix_flux <field> <face_set> <value|function>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* fix_flux = name of the AtC sub-command
* field = field kind name valid for type of physics: *temperature* or *electron_temperature*
* face_set = name of set of element faces
* *value* or *function* = value or name of function followed by its parameters

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC fix_flux temperature faceSet 10.0

Description
"""""""""""

Command for fixing normal fluxes e.g. heat_flux. This command only
prescribes the normal component of the physical flux, e.g. heat (energy)
flux. The units are in AtC units, i.e. derived from the LAMMPS length,
time, and mass scales.

Restrictions
""""""""""""

Only normal fluxes (Neumann data) can be prescribed.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC unfix_flux <atc_unfix_flux>`

Default
"""""""

None.
