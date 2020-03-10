.. index:: fix tfmc

fix tfmc command
================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID tfmc Delta Temp seed keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command
* tfmc = style name of this fix command
* Delta = maximal displacement length (distance units)
* Temp = imposed temperature of the system
* seed = random number seed (positive integer)
* zero or more keyword/arg pairs may be appended
* keyword = *com* or *rot*

  .. parsed-literal::

       *com* args = xflag yflag zflag
         xflag,yflag,zflag = 0/1 to exclude/include each dimension
       *rot* args = none



Examples
""""""""


.. parsed-literal::

   fix 1 all tfmc 0.1 1000.0 159345
   fix 1 all tfmc 0.05 600.0 658943 com 1 1 0
   fix 1 all tfmc 0.1 750.0 387068 com 1 1 1 rot

Description
"""""""""""

Perform uniform-acceptance force-bias Monte Carlo (fbMC) simulations,
using the time-stamped force-bias Monte Carlo (tfMC) algorithm
described in :ref:`(Mees) <Mees>` and :ref:`(Bal) <Bal>`.

One successful use case of force-bias Monte Carlo methods is that they
can be used to extend the time scale of atomistic simulations, in
particular when long time scale relaxation effects must be considered;
some interesting examples are given in the review by :ref:`(Neyts) <Neyts>`.
An example of a typical use case would be the modelling of chemical
vapor deposition (CVD) processes on a surface, in which impacts by
gas-phase species can be performed using MD, but subsequent relaxation
of the surface is too slow to be done using MD only. Using tfMC can
allow for a much faster relaxation of the surface, so that higher
fluxes can be used, effectively extending the time scale of the
simulation. (Such an alternating simulation approach could be set up
using a :doc:`loop <jump>`.)

The initial version of tfMC algorithm in :ref:`(Mees) <Mees>` contained an
estimation of the effective time scale of such a simulation, but it
was later shown that the speed-up one can gain from a tfMC simulation
is system- and process-dependent, ranging from none to several orders
of magnitude. In general, solid-state processes such as
(re)crystallization or growth can be accelerated by up to two or three
orders of magnitude, whereas diffusion in the liquid phase is not
accelerated at all. The observed pseudodynamics when using the tfMC
method is not the actual dynamics one would obtain using MD, but the
relative importance of processes can match the actual relative
dynamics of the system quite well, provided *Delta* is chosen with
care. Thus, the system's equilibrium is reached faster than in MD,
along a path that is generally roughly similar to a typical MD
simulation (but not necessarily so). See :ref:`(Bal) <Bal>` for details.

Each step, all atoms in the selected group are displaced using the
stochastic tfMC algorithm, which is designed to sample the canonical
(NVT) ensemble at the temperature *Temp*\ . Although tfMC is a Monte
Carlo algorithm and thus strictly speaking does not perform time
integration, it is similar in the sense that it uses the forces on all
atoms in order to update their positions. Therefore, it is implemented
as a time integration fix, and no other fixes of this type (such as
:doc:`fix nve <fix_nve>`) should be used at the same time. Because
velocities do not play a role in this kind of Monte Carlo simulations,
instantaneous temperatures as calculated by :doc:`temperature computes <compute_temp>` or :doc:`thermodynamic output <thermo_style>` have no meaning: the only relevant
temperature is the sampling temperature *Temp*\ .  Similarly, performing
tfMC simulations does not require setting a :doc:`timestep <timestep>`
and the :doc:`simulated time <thermo_style>` as calculated by LAMMPS is
meaningless.

The critical parameter determining the success of a tfMC simulation is
*Delta*\ , the maximal displacement length of the lightest element in
the system: the larger it is, the longer the effective time scale of
the simulation will be (there is an approximately quadratic
dependence). However, *Delta* must also be chosen sufficiently small
in order to comply with detailed balance; in general values between 5
and 10 % of the nearest neighbor distance are found to be a good
choice. For a more extensive discussion with specific examples, please
refer to :ref:`(Bal) <Bal>`, which also describes how the code calculates
element-specific maximal displacements from *Delta*\ , based on the
fourth root of their mass.

Because of the uncorrelated movements of the atoms, the center-of-mass
of the fix group will not necessarily be stationary, just like its
orientation. When the *com* keyword is used, all atom positions will
be shifted (after every tfMC iteration) in order to fix the position
of the center-of-mass along the included directions, by setting the
corresponding flag to 1. The *rot* keyword does the same for the
rotational component of the tfMC displacements after every iteration.

.. note::

   the *com* and *rot* keywords should not be used if an external
   force is acting on the specified fix group, along the included
   directions. This can be either a true external force (e.g.  through
   :doc:`fix wall <fix_wall>`) or forces due to the interaction with atoms
   not included in the fix group. This is because in such cases,
   translations or rotations of the fix group could be induced by these
   external forces, and removing them will lead to a violation of
   detailed balance.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

This fix is not compatible with :doc:`fix shake <fix_shake>`.

Related commands
""""""""""""""""

:doc:`fix gcmc <fix_gcmc>`, :doc:`fix nvt <fix_nh>`

Default
"""""""

The option default is com = 0 0 0


----------


.. _Bal:



**(Bal)** K. M Bal and E. C. Neyts, J. Chem. Phys. 141, 204104 (2014).

.. _Mees:



**(Mees)** M. J. Mees, G. Pourtois, E. C. Neyts, B. J. Thijsse, and
A. Stesmans, Phys. Rev. B 85, 134301 (2012).

.. _Neyts:



**(Neyts)** E. C. Neyts and A. Bogaerts, Theor. Chem. Acc. 132, 1320
(2013).
