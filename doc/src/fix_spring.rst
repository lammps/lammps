.. index:: fix spring

fix spring command
==================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID spring keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* spring = style name of this fix command
* keyword = *tether* or *couple*

  .. parsed-literal::

       *tether* values = K x y z R0
         K = spring constant (force/distance units)
         x,y,z = point to which spring is tethered
         R0 = equilibrium distance from tether point (distance units)
       *couple* values = group-ID2 K x y z R0
         group-ID2 = 2nd group to couple to fix group with a spring
         K = spring constant (force/distance units)
         x,y,z = direction of spring
         R0 = equilibrium distance of spring (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix pull ligand spring tether 50.0 0.0 0.0 0.0 0.0
   fix pull ligand spring tether 50.0 0.0 0.0 0.0 5.0
   fix pull ligand spring tether 50.0 NULL NULL 2.0 3.0
   fix 5 bilayer1 spring couple bilayer2 100.0 NULL NULL 10.0 0.0
   fix longitudinal pore spring couple ion 100.0 NULL NULL -20.0 0.0
   fix radial pore spring couple ion 100.0 0.0 0.0 NULL 5.0

Description
"""""""""""

Apply a spring force to a group of atoms or between two groups of
atoms.  This is useful for applying an umbrella force to a small
molecule or lightly tethering a large group of atoms (e.g. all the
solvent or a large molecule) to the center of the simulation box so
that it does not wander away over the course of a long simulation.  It
can also be used to hold the centers of mass of two groups of atoms at
a given distance or orientation with respect to each other.

The *tether* style attaches a spring between a fixed point *x,y,z* and
the center of mass of the fix group of atoms.  The equilibrium
position of the spring is R0.  At each timestep the distance R from
the center of mass of the group of atoms to the tethering point is
computed, taking account of wrap-around in a periodic simulation box.
A restoring force of magnitude K (R - R0) Mi / M is applied to each
atom in the group where *K* is the spring constant, Mi is the mass of
the atom, and M is the total mass of all atoms in the group.  Note
that *K* thus represents the spring constant for the total force on
the group of atoms, not for a spring applied to each atom.

The *couple* style links two groups of atoms together.  The first
group is the fix group; the second is specified by group-ID2.  The
groups are coupled together by a spring that is at equilibrium when
the two groups are displaced by a vector *x,y,z* with respect to each
other and at a distance R0 from that displacement.  Note that *x,y,z*
is the equilibrium displacement of group-ID2 relative to the fix
group.  Thus (1,1,0) is a different spring than (-1,-1,0).  When the
relative positions and distance between the two groups are not in
equilibrium, the same spring force described above is applied to atoms
in each of the two groups.

For both the *tether* and *couple* styles, any of the x,y,z values can
be specified as NULL which means do not include that dimension in the
distance calculation or force application.

The first example above pulls the ligand towards the point (0,0,0).
The second example holds the ligand near the surface of a sphere of
radius 5 around the point (0,0,0).  The third example holds the ligand
a distance 3 away from the z=2 plane (on either side).

The fourth example holds 2 bilayers a distance 10 apart in z.  For the
last two examples, imagine a pore (a slab of atoms with a cylindrical
hole cut out) oriented with the pore axis along z, and an ion moving
within the pore.  The fifth example holds the ion a distance of -20
below the z = 0 center plane of the pore (umbrella sampling).  The
last example holds the ion a distance 5 away from the pore axis
(assuming the center-of-mass of the pore in x,y is the pore axis).

.. note::

   The center of mass of a group of atoms is calculated in
   "unwrapped" coordinates using atom image flags, which means that the
   group can straddle a periodic boundary.  See the :doc:`dump <dump>` doc
   page for a discussion of unwrapped coordinates.  It also means that a
   spring connecting two groups or a group and the tether point can cross
   a periodic boundary and its length be calculated correctly.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the energy stored in the spring to the global
potential energy of the system as part of :doc:`thermodynamic output
<thermo_style>`. The default setting for this fix is :doc:`fix_modify
energy no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the spring energy
= 0.5 \* K \* r\^2.

This fix also computes global 4-vector which can be accessed by
various :doc:`output commands <Howto_output>`.  The first 3 quantities
in the vector are xyz components of the total force added to the group
of atoms by the spring.  In the case of the *couple* style, it is the
force on the fix group (group-ID) or the negative of the force on the
second group (group-ID2).  The fourth quantity in the vector is the
magnitude of the force added by the spring, as a positive value if
(r-R0) > 0 and a negative value if (r-R0) < 0.  This sign convention
can be useful when using the spring force to compute a potential of
mean force (PMF).

The scalar and vector values calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the spring energy to be included in the total
   potential energy of the system (the quantity being minimized), you
   MUST enable the :doc:`fix_modify <fix_modify>` *energy* option for this
   fix.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix drag <fix_drag>`, :doc:`fix spring/self <fix_spring_self>`,
:doc:`fix spring/rg <fix_spring_rg>`, :doc:`fix smd <fix_smd>`

Default
"""""""

none
