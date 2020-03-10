.. index:: pair_style dsmc

pair_style dsmc command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style dsmc max_cell_size seed weighting Tref Nrecompute Nsample

* max\_cell\_size = global maximum cell size for DSMC interactions (distance units)
* seed = random # seed (positive integer)
* weighting = macroparticle weighting
* Tref = reference temperature (temperature units)
* Nrecompute = re-compute v\*sigma\_max every this many timesteps (timesteps)
* Nsample = sample this many times in recomputing v\*sigma\_max

Examples
""""""""

.. code-block:: LAMMPS

   pair_style dsmc 2.5 34387 10 1.0 100 20
   pair_coeff * * 1.0
   pair_coeff 1 1 1.0

Description
"""""""""""

Style *dsmc* computes collisions between pairs of particles for a
direct simulation Monte Carlo (DSMC) model following the exposition in
:ref:`(Bird) <Bird>`.  Each collision resets the velocities of the two
particles involved.  The number of pairwise collisions for each pair
or particle types and the length scale within which they occur are
determined by the parameters of the pair\_style and pair\_coeff
commands.

Stochastic collisions are performed using the variable hard sphere
(VHS) approach, with the user-defined *max\_cell\_size* value used as
the maximum DSMC cell size, and reference cross-sections for
collisions given using the pair\_coeff command.

There is no pairwise energy or virial contributions associated with
this pair style.

The following coefficient must be defined for each pair of atoms types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples above,
or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* sigma (area units, i.e. distance-squared)

The global DSMC *max\_cell\_size* determines the maximum cell length
used in the DSMC calculation.  A structured mesh is overlayed on the
simulation box such that an integer number of cells are created in
each direction for each processor's sub-domain.  Cell lengths are
adjusted up to the user-specified maximum cell size.

----------

To perform a DSMC simulation with LAMMPS, several additional options
should be set in your input script, though LAMMPS does not check for
these settings.

Since this pair style does not compute particle forces, you should use
the "fix nve/noforce" time integration fix for the DSMC particles,
e.g.

.. code-block:: LAMMPS

   fix 1 all nve/noforce

This pair style assumes that all particles will communicated to
neighboring processors every timestep as they move.  This makes it
possible to perform all collisions between pairs of particles that are
on the same processor.  To ensure this occurs, you should use
these commands:

.. code-block:: LAMMPS

   neighbor 0.0 bin
   neigh_modify every 1 delay 0 check no
   atom_modify sort 0 0.0
   communicate single cutoff 0.0

These commands ensure that LAMMPS communicates particles to
neighboring processors every timestep and that no ghost atoms are
created.  The output statistics for a simulation run should indicate
there are no ghost particles or neighbors.

In order to get correct DSMC collision statistics, users should
specify a Gaussian velocity distribution when populating the
simulation domain. Note that the default velocity distribution is
uniform, which will not give good DSMC collision rates. Specify
"dist gaussian" when using the :doc:`velocity <velocity>` command
as in the following:

.. code-block:: LAMMPS

   velocity all create 594.6 87287 loop geom dist gaussian

----------

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.  Note
that the user-specified random number seed is stored in the restart
file, so when a simulation is restarted, each processor will
re-initialize its random number generator the same way it did
initially.  This means the random forces will be random, but will not
be the same as they would have been if the original simulation had
continued past the restart time.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This style is part of the MC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix nve/noforce <fix_nve_noforce>`,
:doc:`neigh_modify <neigh_modify>`, :doc:`neighbor <neighbor>`,
:doc:`comm_modify <comm_modify>`

**Default:** none

----------

.. _Bird:

**(Bird)** G. A. Bird, "Molecular Gas Dynamics and the Direct Simulation
of Gas Flows" (1994).
