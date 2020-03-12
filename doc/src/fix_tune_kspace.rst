.. index:: fix tune/kspace

fix tune/kspace command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID tune/kspace N

* ID, group-ID are documented in :doc:`fix <fix>` command
* tune/kspace = style name of this fix command
* N = invoke this fix every N steps

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all tune/kspace 100

Description
"""""""""""

This fix tests each kspace style (Ewald, PPPM, and MSM), and
automatically selects the fastest style to use for the remainder
of the run. If the fastest style is Ewald or PPPM, the fix also
adjusts the Coulombic cutoff towards optimal speed. Future versions
of this fix will automatically select other kspace parameters
to use for maximum simulation speed. The kspace parameters may
include the style, cutoff, grid points in each direction, order,
Ewald parameter, MSM parallelization cut-point, MPI tasks to use, etc.

The rationale for this fix is to provide the user with
as-fast-as-possible simulations that include long-range electrostatics
(kspace) while meeting the user-prescribed accuracy requirement. A
simple heuristic could never capture the optimal combination of
parameters for every possible run-time scenario. But by performing
short tests of various kspace parameter sets, this fix allows
parameters to be tailored specifically to the user's machine, MPI
ranks, use of threading or accelerators, the simulated system, and the
simulation details. In addition, it is possible that parameters could
be evolved with the simulation on-the-fly, which is useful for systems
that are dynamically evolving (e.g. changes in box size/shape or
number of particles).

When this fix is invoked, LAMMPS will perform short timed tests of
various parameter sets to determine the optimal parameters. Tests are
performed on-the-fly, with a new test initialized every N steps. N should
be chosen large enough so that adequate CPU time lapses between tests,
thereby providing statistically significant timings. But N should not be
chosen to be so large that an unfortunate parameter set test takes an
inordinate amount of wall time to complete. An N of 100 for most problems
seems reasonable. Once an optimal parameter set is found, that set is
used for the remainder of the run.

This fix uses heuristics to guide it's selection of parameter sets to test,
but the actual timed results will be used to decide which set to use in the
simulation.

It is not necessary to discard trajectories produced using sub-optimal
parameter sets, or a mix of various parameter sets, since the user-prescribed
accuracy will have been maintained throughout. However, some users may prefer
to use this fix only to discover the optimal parameter set for a given setup
that can then be used on subsequent production runs.

This fix starts with kspace parameters that are set by the user with the
:doc:`kspace_style <kspace_style>` and :doc:`kspace_modify <kspace_modify>`
commands. The prescribed accuracy will be maintained by this fix throughout
the simulation.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the KSPACE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Do not set "neigh\_modify once yes" or else this fix will never be
called.  Reneighboring is required.

This fix is not compatible with a hybrid pair style, long-range dispersion,
TIP4P water support, or long-range point dipole support.

Related commands
""""""""""""""""

:doc:`kspace_style <kspace_style>`, :doc:`boundary <boundary>`
:doc:`kspace_modify <kspace_modify>`, :doc:`pair_style lj/cut/coul/long <pair_lj>`, :doc:`pair_style lj/charmm/coul/long <pair_charmm>`, :doc:`pair_style lj/long <pair_lj_long>`, :doc:`pair_style lj/long/coul/long <pair_lj_long>`,
:doc:`pair_style buck/coul/long <pair_buck>`

Default
"""""""
