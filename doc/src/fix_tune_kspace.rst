.. index:: fix tune/kspace

fix tune/kspace command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID tune/kspace N keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* tune/kspace = style name of this fix command
* N = invoke this fix every N steps
* zero or more keyword/value pairs may be appended
* keyword = *msm* or *esp*

  .. parsed-literal::

       *msm* value = *yes* or *no*
         yes = also test kspace style *msm*
         no = do not test kspace style *msm*
       *esp* value = *yes* or *no*
         yes = also test kspace style *esp*
         no = do not test kspace style *esp*

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all tune/kspace 100
   fix 2 all tune/kspace 200 msm yes esp yes

Description
"""""""""""

.. versionchanged:: TBD

This fix tests several kspace styles (Ewald and PPPM, and optionally
MSM and ESP), and automatically selects the fastest style to use for
the remainder of the run. If the fastest style is Ewald, PPPM, or ESP,
the fix also adjusts the Coulombic cutoff towards optimal speed.
Previous versions of LAMMPS always tested the MSM style, if a matching
pair style was available.  Future versions
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

This fix uses heuristics to guide its selection of parameter sets to
test, but the actual timed results will be used to decide which set to
use in the simulation.

It is not necessary to discard trajectories produced using sub-optimal
parameter sets, or a mix of various parameter sets, since the user-prescribed
accuracy will have been maintained throughout. However, some users may prefer
to use this fix only to discover the optimal parameter set for a given setup
that can then be used on subsequent production runs.

This fix starts with kspace parameters that are set by the user with the
:doc:`kspace_style <kspace_style>` and :doc:`kspace_modify <kspace_modify>`
commands. The prescribed accuracy will be maintained by this fix throughout
the simulation.

.. versionadded:: TBD

The *msm* and *esp* keywords select whether the :doc:`kspace styles
<kspace_style>` *msm* and *esp* are tested, too.  Each kspace style is
only tested if a matching pair style exists, e.g. for a simulation with
pair style *lj/cut/coul/long* the kspace style *msm* is tested with pair
style *lj/cut/coul/msm*, and the kspace style *esp* with pair style
*lj/cut/coul/esp*.  Currently, matching pair styles for the *esp* kspace
style exist only for *coul/long* and *lj/cut/coul/long*.  By default,
the kspace style *msm* computes only the scalar pressure (see the
*pressure/scalar* keyword of the :doc:`kspace_modify <kspace_modify>`
command), so switching to it would change the computed pressure tensor.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the KSPACE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Do not set "neigh_modify once yes" or else this fix will never be
called.  Reneighboring is required.

This fix is not compatible with a hybrid pair style, long-range dispersion,
TIP4P water support, or long-range point dipole support.

The *msm yes* and *esp yes* settings are not (yet) supported in
combination with the OPENMP package.

Related commands
""""""""""""""""

:doc:`kspace_style <kspace_style>`, :doc:`boundary <boundary>`
:doc:`kspace_modify <kspace_modify>`, :doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`, :doc:`pair_style lj/charmm/coul/long <pair_charmm>`, :doc:`pair_style lj/long <pair_lj_long>`, :doc:`pair_style lj/long/coul/long <pair_lj_long>`,
:doc:`pair_style buck/coul/long <pair_buck>`

Default
"""""""

The keyword defaults are msm = no and esp = no.
