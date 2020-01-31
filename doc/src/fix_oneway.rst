.. index:: fix oneway

fix oneway command
==================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID oneway N region-ID direction

* ID, group-ID are documented in :doc:`fix <fix>` command
* oneway = style name of this fix command
* N = apply this fix every this many timesteps
* region-ID = ID of region where fix is active
* direction = *x* or *-x* or *y* or *-y* or *z* or *-z* = coordinate and direction of the oneway constraint


Examples
""""""""


.. parsed-literal::

   fix ions oneway 10 semi -x
   fix all oneway 1 left -z
   fix all oneway 1 right z

Description
"""""""""""

Enforce that particles in the group and in a given region can only
move in one direction.  This is done by reversing a particle's
velocity component, if it has the wrong sign in the specified
dimension.  The effect is that the particle moves in one direction
only.

This can be used, for example, as a simple model of a semi-permeable
membrane, or as an implementation of Maxwell's demon.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix wall/reflect <fix_wall_reflect>` command

**Default:** none


----------



.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
