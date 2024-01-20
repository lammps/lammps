.. index:: fix rheo/pressure

fix rheo/pressure command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/pressure style args

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/pressure = style name of this fix command
* types = lists of types (see below)
* style = *linear* or *taitwater* or *cubic*

  .. parsed-literal::

       *linear* args = none
       *taitwater* args = none
       *cubic* args = cubic term prefactor :math:`A_3` (pressure/density\^2)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/pressure * linear
   fix 1 all rheo/pressure 1 linear 2 cubic 10.0

Description
"""""""""""

This fix...

Only one instance of fix rheo/pressure can be defined and the fix group must be set to all.


Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes density
such as atom_style rheo or rheo/thermal. This fix must be used in
conjuction with :doc:`fix rheo <fix_rheo>`. The fix group must be
set to all.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none
