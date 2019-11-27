.. index:: fix temp/rescale

fix temp/rescale command
========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID temp/rescale N Tstart Tstop window fraction

* ID, group-ID are documented in :doc:`fix <fix>` command
* temp/rescale = style name of this fix command
* N = perform rescaling every N steps
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
  
  .. parsed-literal::
  
       Tstart can be a variable (see below)

* window = only rescale if temperature is outside this window (temperature units)
* fraction = rescale to target temperature by this fraction


Examples
""""""""


.. parsed-literal::

   fix 3 flow temp/rescale 100 1.0 1.1 0.02 0.5
   fix 3 boundary temp/rescale 1 1.0 1.5 0.05 1.0
   fix 3 boundary temp/rescale 1 1.0 1.5 0.05 1.0

Description
"""""""""""

Reset the temperature of a group of atoms by explicitly rescaling
their velocities.

The rescaling is applied to only the translational degrees of freedom
for the particles, which is an important consideration if finite-size
particles which have rotational degrees of freedom are being
thermostatted with this fix.  The translational degrees of freedom can
also have a bias velocity removed from them before thermostatting
takes place; see the description below.

Rescaling is performed every N timesteps.  The target temperature is a
ramped value between the *Tstart* and *Tstop* temperatures at the
beginning and end of the run.

.. note::

   This thermostat will generate an error if the current
   temperature is zero at the end of a timestep it is invoked on.  It
   cannot rescale a zero temperature.

*Tstart* can be specified as an equal-style :doc:`variable <variable>`.
In this case, the *Tstop* setting is ignored.  If the value is a
variable, it should be specified as v\_name, where name is the variable
name.  In this case, the variable will be evaluated each timestep, and
its value used to determine the target temperature.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo\_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent temperature.

Rescaling is only performed if the difference between the current and
desired temperatures is greater than the *window* value.  The amount
of rescaling that is applied is a *fraction* (from 0.0 to 1.0) of the
difference between the actual and desired temperature.  E.g. if
*fraction* = 1.0, the temperature is reset to exactly the desired
value.

.. note::

   Unlike the :doc:`fix nvt <fix_nh>` command which performs
   Nose/Hoover thermostatting AND time integration, this fix does NOT
   perform time integration.  It only modifies velocities to effect
   thermostatting.  Thus you must use a separate time integration fix,
   like :doc:`fix nve <fix_nve>` to actually update the positions of atoms
   using the modified velocities.  Likewise, this fix should not normally
   be used on atoms that also have their temperature controlled by
   another fix - e.g. by :doc:`fix nvt <fix_nh>` or :doc:`fix langevin <fix_langevin>` commands.

See the :doc:`Howto thermostat <Howto_thermostat>` doc page for a
discussion of different ways to compute temperature and perform
thermostatting.

This fix computes a temperature each timestep.  To do this, the fix
creates its own compute of style "temp", as if one of this command had
been issued:


.. parsed-literal::

   compute fix-ID_temp group-ID temp

See the :doc:`compute temp <compute_temp>` for details.  Note that the
ID of the new compute is the fix-ID + underscore + "temp", and the
group for the new compute is the same as the fix group.

Note that this is NOT the compute used by thermodynamic output (see
the :doc:`thermo\_style <thermo_style>` command) with ID = *thermo\_temp*.
This means you can change the attributes of this fix's temperature
(e.g. its degrees-of-freedom) via the
:doc:`compute\_modify <compute_modify>` command or print this temperature
during thermodynamic output via the :doc:`thermo\_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo\_temp* will have no
effect on this fix.

Like other fixes that perform thermostatting, this fix can be used
with :doc:`compute commands <compute>` that calculate a temperature
after removing a "bias" from the atom velocities.  E.g. removing the
center-of-mass velocity from a group of atoms or only calculating
temperature on the x-component of velocity or only calculating
temperature for atoms in a geometric region.  This is not done by
default, but only if the :doc:`fix\_modify <fix_modify>` command is used
to assign a temperature compute to this fix that includes such a bias
term.  See the doc pages for individual :doc:`compute commands <compute>` to determine which ones include a bias.  In
this case, the thermostat works in the following manner: the current
temperature is calculated taking the bias into account, bias is
removed from each atom, thermostatting is performed on the remaining
thermal degrees of freedom, and the bias is added back in.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix\_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a temperature :doc:`compute <compute>`
you have defined to this fix which will be used in its thermostatting
procedure, as described above.  For consistency, the group used by
this fix and by the compute should be the same.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change implied by a velocity rescaling to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix langevin <fix_langevin>`, :doc:`fix nvt <fix_nh>`,
:doc:`fix\_modify <fix_modify>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
