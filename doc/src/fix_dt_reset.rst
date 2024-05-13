.. index:: fix dt/reset
.. index:: fix dt/reset/kk

fix dt/reset command
====================

Accelerator Variants: *dt/reset/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID dt/reset N Tmin Tmax Xmax keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* dt/reset = style name of this fix command
* N = re-compute dt every N timesteps
* Tmin = minimum dt allowed which can be NULL (time units)
* Tmax = maximum dt allowed which can be NULL (time units)
* Xmax = maximum distance for an atom to move in one timestep (distance units)
* zero or more keyword/value pairs may be appended
* keyword = *emax* or *units*

.. parsed-literal::

     *emax* value = Emax
       Emax = maximum kinetic energy change for an atom in one timestep (energy units)
     *units* value = *lattice* or *box*
       lattice = Xmax is defined in lattice units
       box = Xmax is defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix 5 all dt/reset 10 1.0e-5 0.01 0.1
   fix 5 all dt/reset 10 0.01 2.0 0.2 units box
   fix 5 all dt/reset 5 NULL 0.001 0.5 emax 30 units box

Description
"""""""""""

Reset the timestep size every N steps during a run, so that no atom
moves further than the specified *Xmax* distance, based on current
atom velocities and forces.  Optionally an additional criterion is
imposed by the *emax* keyword, so that no atom's kinetic energy
changes by more than the specified *Emax*\ .

This can be useful when starting from a configuration with overlapping
atoms, where forces will be large.  Or it can be useful when running
an impact simulation where one or more high-energy atoms collide with
a solid, causing a damage cascade.

This fix overrides the timestep size setting made by the
:doc:`timestep <timestep>` command.  The new timestep size *dt* is
computed in the following manner.

For each atom, the timestep is computed that would cause it to
displace *Xmax* on the next integration step, as a function of its
current velocity and force.  Since performing this calculation exactly
would require the solution to a quartic equation, a cheaper estimate
is generated.  The estimate is conservative in that the atom's
displacement is guaranteed not to exceed *Xmax*, though it may be
smaller.

In addition if the *emax* keyword is used, the specified *Emax* value
is enforced as a limit on how much an atom's kinetic energy can
change.  If the timestep required is even smaller than for the *Xmax*
displacement, then the smaller timestep is used.

Given this putative timestep for each atom, the minimum timestep value
across all atoms is computed.  Then the *Tmin* and *Tmax* bounds are
applied, if specified.  If one (or both) is specified as NULL, it is
not applied.

When the :doc:`run style <run_style>` is *respa*, this fix resets the
outer loop (largest) timestep, which is the same timestep that the
:doc:`timestep <timestep>` command sets.

Note that the cumulative simulation time (in time units), which
accounts for changes in the timestep size as a simulation proceeds,
can be accessed by the :doc:`thermo_style time <thermo_style>`
keyword.

Also note that the :doc:`dump_modify every/time <dump_modify>` option
allows dump files to be written at intervals specified by simulation
time, rather than by timesteps.  Simulation time is in time units;
see the :doc:`units <units>` doc page for details.

----------

.. include:: accel_styles.rst

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar stores the last
timestep on which the timestep was reset to a new value.

The scalar value calculated by this fix is "intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`timestep <timestep>`, :doc:`dump_modify every/time <dump_modify>`

Default
"""""""

The option defaults are units = lattice, and no emax kinetic energy
limit.
