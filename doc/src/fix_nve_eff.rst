.. index:: fix nve/eff

fix nve/eff command
===================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nve/eff

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/eff = style name of this fix command

Examples
""""""""


.. parsed-literal::

   fix 1 all nve/eff

Description
"""""""""""

Perform constant NVE integration to update position and velocity for
nuclei and electrons in the group for the :doc:`electron force field <pair_eff>` model.  V is volume; E is energy.  This creates a
system trajectory consistent with the microcanonical ensemble.

The operation of this fix is exactly like that described by the :doc:`fix nve <fix_nve>` command, except that the radius and radial velocity
of electrons are also updated.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the USER-EFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nvt/eff <fix_nh_eff>`, :doc:`fix npt/eff <fix_nh_eff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
