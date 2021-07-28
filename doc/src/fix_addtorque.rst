.. index:: fix addtorque

fix addtorque command
=====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID addtorque Tx Ty Tz

* ID, group-ID are documented in :doc:`fix <fix>` command
* addtorque = style name of this fix command
* Tx,Ty,Tz = torque component values (torque units)
* any of Tx,Ty,Tz can be a variable (see below)

Examples
""""""""

.. code-block:: LAMMPS

   fix kick bead addtorque 2.0 3.0 5.0
   fix kick bead addtorque 0.0 0.0 v_oscillate

Description
"""""""""""

Add a set of forces to each atom in
the group such that:

* the components of the total torque applied on the group (around its
  center of mass) are Tx,Ty,Tz
* the group would move as a rigid body in the absence of other
  forces.

This command can be used to drive a group of atoms into rotation.

Any of the 3 quantities defining the torque components can be specified
as an equal-style :doc:`variable <variable>`, namely *Tx*,
*Ty*, *Tz*\ .  If the value is a variable, it should be specified as
v_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
torque component.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent torque.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the potential "energy" inferred by the added torques
to the global potential energy of the system as part of
:doc:`thermodynamic output <thermo_style>`.  The default setting for
this fix is :doc:`fix_modify energy no <fix_modify>`.  Note that this
is a fictitious quantity but is needed so that the :doc:`minimize
<minimize>` command can include the forces added by this fix in a
consistent manner.  I.e. there is a decrease in potential energy when
atoms move in the direction of the added forces.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by
this fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its torque. Default is the outermost level.

This fix computes a global scalar and a global 3-vector, which can be
accessed by various :doc:`output commands <Howto_output>`.  The scalar
is the potential energy discussed above.  The vector is the total
torque on the group of atoms before the forces on individual atoms are
changed by the fix.  The scalar and vector values calculated by this
fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the fictitious potential energy associated with the
   added forces to be included in the total potential energy of the
   system (the quantity being minimized), you MUST enable the
   :doc:`fix_modify <fix_modify>` *energy* option for this fix.

.. note::

   You should not specify force components with a variable that has
   time-dependence for use with a minimizer, since the minimizer
   increments the timestep as the iteration count during the
   minimization.

Restrictions
""""""""""""

This fix is part of the MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix addforce <fix_addforce>`

Default
"""""""

none
