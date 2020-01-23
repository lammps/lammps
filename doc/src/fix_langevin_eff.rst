.. index:: fix langevin/eff

fix langevin/eff command
========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID langevin/eff Tstart Tstop damp seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin/eff = style name of this fix command
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* damp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)
* zero or more keyword/value pairs may be appended
  
  .. parsed-literal::
  
     keyword = *scale* or *tally* or *zero*
       *scale* values = type ratio
         type = atom type (1-N)
         ratio = factor by which to scale the damping coefficient
       *tally* values = *no* or *yes*
         *no* = do not tally the energy added/subtracted to atoms
         *yes* = do tally the energy added/subtracted to atoms

  
  .. parsed-literal::
  
       *zero* value = *no* or *yes*
         *no* = do not set total random force to zero
         *yes* = set total random force to zero



Examples
""""""""


.. parsed-literal::

   fix 3 boundary langevin/eff 1.0 1.0 10.0 699483
   fix 1 all langevin/eff 1.0 1.1 10.0 48279 scale 3 1.5

Description
"""""""""""

Apply a Langevin thermostat as described in :ref:`(Schneider) <Schneider2>`
to a group of nuclei and electrons in the :doc:`electron force field <pair_eff>` model.  Used with :doc:`fix nve/eff <fix_nve_eff>`,
this command performs Brownian dynamics (BD), since the total force on
each atom will have the form:


.. parsed-literal::

   F = Fc + Ff + Fr
   Ff = - (m / damp) v
   Fr is proportional to sqrt(Kb T m / (dt damp))

Fc is the conservative force computed via the usual inter-particle
interactions (:doc:`pair_style <pair_style>`).

The Ff and Fr terms are added by this fix on a per-particle basis.

The operation of this fix is exactly like that described by the :doc:`fix langevin <fix_langevin>` command, except that the thermostatting
is also applied to the radial electron velocity for electron
particles.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  Because the state of the random number generator
is not saved in restart files, this means you cannot do "exact"
restarts with this fix, where the simulation continues on the same as
if no restart had taken place.  However, in a statistical sense, a
restarted simulation should produce the same behavior.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a temperature :doc:`compute <compute>`
you have defined to this fix which will be used in its thermostatting
procedure, as described above.  For consistency, the group used by
this fix and by the compute should be the same.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Langevin thermostatting to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.  Note that use of this option requires
setting the *tally* keyword to *yes*\ .

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".  Note that calculation of this quantity requires
setting the *tally* keyword to *yes*\ .

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

This fix is part of the USER-EFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix langevin <fix_langevin>`

Default
"""""""

The option defaults are scale = 1.0 for all types and tally = no.


----------


.. _Dunweg2:



**(Dunweg)** Dunweg and Paul, Int J of Modern Physics C, 2, 817-27 (1991).

.. _Schneider2:



**(Schneider)** Schneider and Stoll, Phys Rev B, 17, 1302 (1978).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
