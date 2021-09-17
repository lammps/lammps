.. index:: fix smd

fix smd command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID smd type values keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd  = style name of this fix command
* mode = *cvel* or *cfor* to select constant velocity or constant force SMD

  .. parsed-literal::

       *cvel* values = K vel
         K = spring constant (force/distance units)
         vel = velocity of pulling (distance/time units)
       *cfor* values = force
         force = pulling force (force units)

* keyword = *tether* or *couple*

  .. parsed-literal::

       *tether* values = x y z R0
         x,y,z = point to which spring is tethered
         R0 = distance of end of spring from tether point (distance units)
       *couple* values = group-ID2 x y z R0
         group-ID2 = 2nd group to couple to fix group with a spring
         x,y,z = direction of spring, automatically computed with 'auto'
         R0 = distance of end of spring (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix  pull    cterm smd cvel 20.0 -0.00005 tether NULL NULL 100.0 0.0
   fix  pull    cterm smd cvel 20.0 -0.0001 tether 25.0 25 25.0 0.0
   fix  stretch cterm smd cvel 20.0  0.0001 couple nterm auto auto auto 0.0
   fix  pull    cterm smd cfor  5.0 tether 25.0 25.0 25.0 0.0

Description
"""""""""""

This fix implements several options of steered MD (SMD) as reviewed in
:ref:`(Izrailev) <Izrailev>`, which allows to induce conformational
changes in systems and to compute the potential of mean force (PMF)
along the assumed reaction coordinate :ref:`(Park) <Park>` based on
Jarzynski's equality :ref:`(Jarzynski) <Jarzynski>`.  This fix borrows
a lot from :doc:`fix spring <fix_spring>` and :doc:`fix setforce
<fix_setforce>`.

You can apply a moving spring force to a group of atoms (\ *tether*
style) or between two groups of atoms (\ *couple* style).  The spring
can then be used in either constant velocity (\ *cvel*\ ) mode or in
constant force (\ *cfor*\ ) mode to induce transitions in your systems.
When running in *tether* style, you may need some way to fix some
other part of the system (e.g. via :doc:`fix spring/self <fix_spring_self>`)

The *tether* style attaches a spring between a point at a distance of
R0 away from a fixed point *x,y,z* and the center of mass of the fix
group of atoms.  A restoring force of magnitude K (R - R0) Mi / M is
applied to each atom in the group where *K* is the spring constant, Mi
is the mass of the atom, and M is the total mass of all atoms in the
group.  Note that *K* thus represents the total force on the group of
atoms, not a per-atom force.

In *cvel* mode the distance R is incremented or decremented
monotonously according to the pulling (or pushing) velocity.
In *cfor* mode a constant force is added and the actual distance
in direction of the spring is recorded.

The *couple* style links two groups of atoms together.  The first
group is the fix group; the second is specified by group-ID2.  The
groups are coupled together by a spring that is at equilibrium when
the two groups are displaced by a vector in direction *x,y,z* with
respect to each other and at a distance R0 from that displacement.
Note that *x,y,z* only provides a direction and will be internally
normalized. But since it represents the *absolute* displacement of
group-ID2 relative to the fix group, (1,1,0) is a different spring
than (-1,-1,0).  For each vector component, the displacement can be
described with the *auto* parameter. In this case the direction is
re-computed in every step, which can be useful for steering a local
process where the whole object undergoes some other change.  When the
relative positions and distance between the two groups are not in
equilibrium, the same spring force described above is applied to atoms
in each of the two groups.

For both the *tether* and *couple* styles, any of the x,y,z values can
be specified as NULL which means do not include that dimension in the
distance calculation or force application.

For constant velocity pulling (\ *cvel* mode), the running integral
over the pulling force in direction of the spring is recorded and
can then later be used to compute the potential of mean force (PMF)
by averaging over multiple independent trajectories along the same
pulling path.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The fix stores the direction of the spring, current pulling target
distance and the running PMF to :doc:`binary restart files <restart>`.
See the :doc:`read_restart <read_restart>` command for info on how to
re-specify a fix in an input script that reads a restart file, so that
the operation of the fix continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by
this fix to add the contribution due to the added forces on atoms to
both the global pressure and per-atom stress of the system via the
:doc:`compute pressure <compute_pressure>` and :doc:`compute
stress/atom <compute_stress_atom>` commands.  The former can be
accessed by :doc:`thermodynamic output <thermo_style>`.  The default
setting for this fix is :doc:`fix_modify virial no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its forces. Default is the
outermost level.

This fix computes a vector list of 7 quantities, which can be accessed
by various :doc:`output commands <Howto_output>`.  The quantities in
the vector are in this order: the x-, y-, and z-component of the
pulling force, the total force in direction of the pull, the
equilibrium distance of the spring, the distance between the two
reference points, and finally the accumulated PMF (the sum of pulling
forces times displacement).

The force is the total force on the group of atoms by the spring.  In
the case of the *couple* style, it is the force on the fix group
(group-ID) or the negative of the force on the second group
(group-ID2).  The vector values calculated by this fix are
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix drag <fix_drag>`, :doc:`fix spring <fix_spring>`,
:doc:`fix spring/self <fix_spring_self>`,
:doc:`fix spring/rg <fix_spring_rg>`,
:doc:`fix colvars <fix_colvars>`, :doc:`fix plumed <fix_plumed>`

Default
"""""""

none

----------

.. _Izrailev:

**(Izrailev)** Izrailev, Stepaniants, Isralewitz, Kosztin, Lu, Molnar,
Wriggers, Schulten. Computational Molecular Dynamics: Challenges,
Methods, Ideas, volume 4 of Lecture Notes in Computational Science and
Engineering, pp. 39-65. Springer-Verlag, Berlin, 1998.

.. _Park:

**(Park)** Park, Schulten, J. Chem. Phys. 120 (13), 5946 (2004)

.. _Jarzynski:

**(Jarzynski)** Jarzynski, Phys. Rev. Lett. 78, 2690 (1997)
