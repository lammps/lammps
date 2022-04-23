.. index:: compute temp/uef

compute temp/uef command
========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID temp/uef

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/uef = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all temp/uef
   compute 2 sel temp/uef

Description
"""""""""""

This command is used to compute the kinetic energy tensor in
the reference frame of the applied flow field when
:doc:`fix nvt/uef <fix_nh_uef>` or
:doc:`fix npt/uef <fix_nh_uef>` is used.
It is not necessary to use this command to compute the scalar
value of the temperature. A :doc:`compute temp <compute_temp>`
may be used for that purpose.

Output information for this command can be found in the
documentation for :doc:`compute temp <compute_temp>`.

Restrictions
""""""""""""

This fix is part of the UEF package. It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package <Build_package>` page for more info.

This command can only be used when :doc:`fix nvt/uef <fix_nh_uef>`
or :doc:`fix npt/uef <fix_nh_uef>` is active.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`,
:doc:`fix nvt/uef <fix_nh_uef>`,
:doc:`compute pressure/uef <compute_pressure_uef>`

Default
"""""""

none
