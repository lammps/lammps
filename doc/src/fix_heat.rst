.. index:: fix heat

fix heat command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID heat N eflux

* ID, group-ID are documented in :doc:`fix <fix>` command
* heat = style name of this fix command
* N = add/subtract heat every this many timesteps
* eflux = rate of heat addition or subtraction (energy/time units)
* eflux can be a variable (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region atoms must be in to have added force

Examples
""""""""

.. parsed-literal::

   fix 3 qin heat 1 1.0
   fix 3 qin heat 10 v_flux
   fix 4 qout heat 1 -1.0 region top

Description
"""""""""""

Add non-translational kinetic energy (heat) to a group of atoms in a
manner that conserves their aggregate momentum.  Two of these fixes
can be used to establish a temperature gradient across a simulation
domain by adding heat (energy) to one group of atoms (hot reservoir)
and subtracting heat from another (cold reservoir).  E.g. a simulation
sampling from the McDLT ensemble.

If the *region* keyword is used, the atom must be in both the group
and the specified geometric :doc:`region <region>` in order to have
energy added or subtracted to it.  If not specified, then the atoms in
the group are affected wherever they may move to.

Heat addition/subtraction is performed every N timesteps.  The *eflux*
parameter can be specified as a numeric constant or as a variable (see
below).  If it is a numeric constant or equal-style variable which
evaluates to a scalar value, then the *eflux* determines the change in
aggregate energy of the entire group of atoms per unit time, e.g. in
eV/psec for :doc:`metal units <units>`.  In this case it is an
"extensive" quantity, meaning its magnitude should be scaled with the
number of atoms in the group.  Note that since *eflux* has per-time
units (i.e. it is a flux), this means that a larger value of N will
add/subtract a larger amount of energy each time the fix is invoked.

.. note::

   The heat-exchange (HEX) algorithm implemented by this fix is
   known to exhibit a pronounced energy drift. An improved algorithm
   (eHEX) is available as a :doc:`fix ehex <fix_ehex>` command and might be
   preferable if energy conservation is important.

If *eflux* is specified as an atom-style variable (see below), then
the variable computes one value per atom.  In this case, each value is
the energy flux for a single atom, again in units of energy per unit
time.  In this case, each value is an "intensive" quantity, which need
not be scaled with the number of atoms in the group.

As mentioned above, the *eflux* parameter can be specified as an
equal-style or atom\_style :doc:`variable <variable>`.  If the value is a
variable, it should be specified as v\_name, where name is the variable
name.  In this case, the variable will be evaluated each timestep, and
its value(s) used to determine the flux.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent flux.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent flux
with optional time-dependence as well.

.. note::

   If heat is subtracted from the system too aggressively so that
   the group's kinetic energy would go to zero, or any individual atom's
   kinetic energy would go to zero for the case where *eflux* is an
   atom-style variable, then LAMMPS will halt with an error message.

Fix heat is different from a thermostat such as :doc:`fix nvt <fix_nh>`
or :doc:`fix temp/rescale <fix_temp_rescale>` in that energy is
added/subtracted continually.  Thus if there isn't another mechanism
in place to counterbalance this effect, the entire system will heat or
cool continuously.  You can use multiple heat fixes so that the net
energy change is 0.0 or use :doc:`fix viscous <fix_viscous>` to drain
energy from the system.

This fix does not change the coordinates of its atoms; it only scales
their velocities.  Thus you must still use an integration fix
(e.g. :doc:`fix nve <fix_nve>`) on the affected atoms.  This fix should
not normally be used on atoms that have their temperature controlled
by another fix - e.g. :doc:`fix nvt <fix_nh>` or :doc:`fix langevin <fix_langevin>` fix.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  This scalar is the most recent
value by which velocities were scaled.  The scalar value calculated by
this fix is "intensive".  If *eflux* is specified as an atom-style
variable, this fix computes the average value by which the velocities
were scaled for all of the atoms that had their velocities scaled.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix ehex <fix_ehex>`, :doc:`compute temp <compute_temp>`, :doc:`compute temp/region <compute_temp_region>`

**Default:** none
