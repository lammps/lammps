.. index:: fix shardlow

fix shardlow command
====================

fix shardlow/kk command
=======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID shardlow

* ID, group-ID are documented in :doc:`fix <fix>` command
* shardlow = style name of this fix command

Examples
""""""""


.. parsed-literal::

   fix 1 all shardlow

Description
"""""""""""

Specifies that the Shardlow splitting algorithm (SSA) is to be used to
integrate the DPD equations of motion.  The SSA splits the integration
into a stochastic and deterministic integration step.  The fix
*shardlow* performs the stochastic integration step and must be used
in conjunction with a deterministic integrator (e.g. :doc:`fix nve <fix_nve>` or :doc:`fix nph <fix_nh>`).  The stochastic
integration of the dissipative and random forces is performed prior to
the deterministic integration of the conservative force. Further
details regarding the method are provided in :ref:`(Lisal) <Lisal>` and
:ref:`(Larentzos1) <Larentzos1sh>`.

The fix *shardlow* must be used with the :doc:`pair\_style dpd/fdt <pair_style>` or :doc:`pair\_style dpd/fdt/energy <pair_style>` command to properly initialize the
fluctuation-dissipation theorem parameter(s) sigma (and kappa, if
necessary).

Note that numerous variants of DPD can be specified by choosing an
appropriate combination of the integrator and :doc:`pair\_style dpd/fdt <pair_style>` command.  DPD under isothermal conditions can
be specified by using fix *shardlow*\ , fix *nve* and pair\_style
*dpd/fdt*\ .  DPD under isoenergetic conditions can be specified by
using fix *shardlow*\ , fix *nve* and pair\_style *dpd/fdt/energy*\ .  DPD
under isobaric conditions can be specified by using fix shardlow, fix
*nph* and pair\_style *dpd/fdt*\ .  DPD under isoenthalpic conditions can
be specified by using fix shardlow, fix *nph* and pair\_style
*dpd/fdt/energy*\ .  Examples of each DPD variant are provided in the
examples/USER/dpd directory.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This command is part of the USER-DPD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix is currently limited to orthogonal simulation cell
geometries.

This fix must be used with an additional fix that specifies time
integration, e.g. :doc:`fix nve <fix_nve>` or :doc:`fix nph <fix_nh>`.

The Shardlow splitting algorithm requires the sizes of the sub-domain
lengths to be larger than twice the cutoff+skin.  Generally, the
domain decomposition is dependent on the number of processors
requested.

Related commands
""""""""""""""""

:doc:`pair\_style dpd/fdt <pair_dpd_fdt>`, :doc:`fix eos/cv <fix_eos_cv>`

**Default:** none


----------


.. _Lisal:



**(Lisal)** M. Lisal, J.K. Brennan, J. Bonet Avalos, "Dissipative
particle dynamics as isothermal, isobaric, isoenergetic, and
isoenthalpic conditions using Shardlow-like splitting algorithms.",
J. Chem. Phys., 135, 204105 (2011).

.. _Larentzos1sh:



**(Larentzos1)** J.P. Larentzos, J.K. Brennan, J.D. Moore, M. Lisal and
W.D. Mattson, "Parallel Implementation of Isothermal and Isoenergetic
Dissipative Particle Dynamics Using Shardlow-Like Splitting
Algorithms", Comput. Phys. Commun., 185, 1987-1998 (2014).

.. _Larentzos2sh:



**(Larentzos2)** J.P. Larentzos, J.K. Brennan, J.D. Moore, and
W.D. Mattson, "LAMMPS Implementation of Constant Energy Dissipative
Particle Dynamics (DPD-E)", ARL-TR-6863, U.S. Army Research
Laboratory, Aberdeen Proving Ground, MD (2014).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
