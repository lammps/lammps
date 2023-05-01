.. index:: fix rheo/thermal

fix rheo/thermal command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/thermal keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/viscosity = style name of this fix command
* one or more attributes may be appended
* attribute = *conductivity* or *specific/heat* or *Tfreeze*

  .. parsed-literal::

       *conductivity* args = style param
         style = *constant* or *type*
           *constant* arg = conductivity (power/temperature)
           *type* args = list of conductivity values, one per type (power/temperature)
       *specific/heat* args = style param
         style = *constant* or *type*
           *constant* arg = specific heat (energy/(mass*temperature))
           *type* args = list of specific heat values, one per atom type (energy/(mass*temperature))
       *Tfreeze* args = style param
         style = *constant* or *type*
           *constant* arg = freezing temperature (temperature)
           *type* args = list of freezing temperature values, one per type (temperature)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/thermal conductivity constant 1.0 specific/heat constant 1.0 Tfreeze constant 1.0
   fix 1 all rheo/pressure conductivity constant 1.0 specific/heat type 1.0 2.0

Description
"""""""""""

This fix...

While the *Tfreeze* keyword is optional, the *conducitivity* and
*specific/heat* keywords are mandatory.

Multiple instances of this fix may be defined to apply different
properties to different groups. However, the union of fix groups
across all instances of fix rheo/thermal must cover all atoms.
If there are multiple instances of this fix, any intersections in
the fix groups will lead to incorrect thermal integration.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes temperature,
heatflow, and conductivity such as atom_tyle rheo/thermal This fix
must be used in conjuction with :doc:`fix rheo <fix_rheo>` with the
*thermal* setting.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none
