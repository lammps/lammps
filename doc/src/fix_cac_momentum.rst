.. index:: fix cac/momentum

fix cac/momentum command
========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID momentum N keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* cac/momentum = style name of this fix command
* N = adjust the momentum every this many timesteps
  one or more keyword/value pairs may be appended
* keyword = *linear* or *angular*
  
  .. parsed-literal::
  
       *linear* values = xflag yflag zflag
         xflag,yflag,zflag = 0/1 to exclude/include each dimension
       *angular* values = none


Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all cac/momentum 1 linear 1 1 0
   fix 1 boundary cac/momentum 100 linear 1 0 1
   fix 1 all cac/momentum 100 linear 1 1 1 angular

Description
"""""""""""

Zero the linear and/or angular momentum of the group of atoms and CAC
elements every N timesteps by adjusting their velocities.  One (or both) of
the *linear* or *angular* keywords must be specified.

If the *linear* keyword is used, the linear momentum is zeroed by
subtracting the center-of-mass velocity of the group. One
or more dimensions can be excluded from this operation by setting the
corresponding flag to 0.

If the *angular* keyword is used, the angular momentum is zeroed by
subtracting a rotational component from the group.

This command can be used to insure the entire collection of atoms and CAC
elements (or a subset of them) does not drift or rotate during the simulation
due to random perturbations (e.g. :doc:`fix langevin <fix_langevin>`
thermostatting).

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

:doc:`fix momentum <fix_momentum>`

**Default:** none
