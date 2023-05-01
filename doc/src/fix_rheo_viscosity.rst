.. index:: fix rheo/viscosity

fix rheo/viscosity command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/viscosity vstyle args

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/viscosity = style name of this fix command
* vstyle = *constant* or *type* or *power*

  .. parsed-literal::

       *constant* arg = viscosity (mass/(length*time))
       *type* args = list of viscosity values, one per atom type (mass/(length*time))
       *power* args = *eta* *gd0* *K* *npow* *tau0*
         *eta* = (units)
         *gd0* = (units)
         *K* = (units)
         *npow* = (units)
         *tau0* = (units)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/viscosity constant 1.0
   fix 1 all rheo/viscosity power 0.1 1e-2 0.5 0.01

Description
"""""""""""

This fix...

Multiple instances of this fix may be defined to apply different
properties to different groups. However, the union of fix groups
across all instances of fix rheo/viscosity must cover all atoms.
If there are multiple instances of this fix, any intersection
between fix groups will cause the viscosity for the affected atoms
to be calculated multiple times. Any such affected atoms will enabled
up with a viscosity calculated by the latest defined fix.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes viscosity
such as atom_style rheo or rheo/thermal. This fix must be used in
conjuction with :doc:`fix rheo <fix_rheo>`.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none
