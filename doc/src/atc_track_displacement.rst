.. index:: fix_modify AtC track_displacement

fix_modify AtC track_displacement command
=========================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> track_displacement <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* track_displacement = name of the AtC sub-command
* *on* or *off* = (undocumented)

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC track_displacement on

Description
"""""""""""

Determines whether displacement is tracked or not. For solids problems
this is a useful quantity, but for fluids it is not relevant.

Restrictions
""""""""""""

Some constitutive models require the displacement field.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

*on*
