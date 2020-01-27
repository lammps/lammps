.. index:: displace\_atoms

displace\_atoms command
=======================

Syntax
""""""


.. parsed-literal::

   displace_atoms group-ID style args keyword value ...

* group-ID = ID of group of atoms to displace
* style = *move* or *ramp* or *random* or *rotate*
  
  .. parsed-literal::
  
       *move* args = delx dely delz
         delx,dely,delz = distance to displace in each dimension (distance units)
         any of delx,dely,delz can be a variable (see below)
       *ramp* args = ddim dlo dhi dim clo chi
         ddim = *x* or *y* or *z*
         dlo,dhi = displacement distance between dlo and dhi (distance units)
         dim = *x* or *y* or *z*
         clo,chi = lower and upper bound of domain to displace (distance units)
       *random* args = dx dy dz seed
         dx,dy,dz = random displacement magnitude in each dimension (distance units)
         seed = random # seed (positive integer)
       *rotate* args = Px Py Pz Rx Ry Rz theta
         Px,Py,Pz = origin point of axis of rotation (distance units)
         Rx,Ry,Rz = axis of rotation vector
         theta = angle of rotation (degrees)

* zero or more keyword/value pairs may be appended
  
  .. parsed-literal::
  
       keyword = *units*
         value = *box* or *lattice*



Examples
""""""""


.. parsed-literal::

   displace_atoms top move 0 -5 0 units box
   displace_atoms flow ramp x 0.0 5.0 y 2.0 20.5

Description
"""""""""""

Displace a group of atoms.  This can be used to move atoms a large
distance before beginning a simulation or to randomize atoms initially
on a lattice.  For example, in a shear simulation, an initial strain
can be imposed on the system.  Or two groups of atoms can be brought
into closer proximity.

The *move* style displaces the group of atoms by the specified 3d
displacement vector.  Any of the 3 quantities defining the vector
components can be specified as an equal-style or atom-style
:doc:`variable <variable>`.  If the value is a variable, it should be
specified as v\_name, where name is the variable name.  In this case,
the variable will be evaluated, and its value(s) used for the
displacement(s).  The scale factor implied by the *units* keyword will
also be applied to the variable result.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Atom-style variables can specify the same formulas as
equal-style variables but can also include per-atom values, such as
atom coordinates or per-atom values read from a file.  Note that if
the variable references other :doc:`compute <compute>` or :doc:`fix <fix>`
commands, those values must be up-to-date for the current timestep.
See the "Variable Accuracy" section of the :doc:`variable <variable>`
doc page for more details.

The *ramp* style displaces atoms a variable amount in one dimension
depending on the atom's coordinate in a (possibly) different
dimension.  For example, the second example command displaces atoms in
the x-direction an amount between 0.0 and 5.0 distance units.  Each
atom's displacement depends on the fractional distance its y
coordinate is between 2.0 and 20.5.  Atoms with y-coordinates outside
those bounds will be moved the minimum (0.0) or maximum (5.0) amount.

The *random* style independently moves each atom in the group by a
random displacement, uniformly sampled from a value between -dx and
+dx in the x dimension, and similarly for y and z.  Random numbers are
used in such a way that the displacement of a particular atom is the
same, regardless of how many processors are being used.

The *rotate* style rotates each atom in the group by the angle *theta*
around a rotation axis *R* = (Rx,Ry,Rz) that goes through a point *P* =
(Px,Py,Pz).  The direction of rotation for the atoms around the
rotation axis is consistent with the right-hand rule: if your
right-hand thumb points along *R*\ , then your fingers wrap around the
axis in the direction of positive theta.

If the defined :doc:`atom_style <atom_style>` assigns an orientation to
each atom (:doc:`atom styles <atom_style>` ellipsoid, line, tri, body),
then that property is also updated appropriately to correspond to the
atom's rotation.

Distance units for displacements and the origin point of the *rotate*
style are determined by the setting of *box* or *lattice* for the
*units* keyword.  *Box* means distance units as defined by the
:doc:`units <units>` command - e.g. Angstroms for *real* units.
*Lattice* means distance units are in lattice spacings.  The
:doc:`lattice <lattice>` command must have been previously used to
define the lattice spacing.


----------


.. note::

   Care should be taken not to move atoms on top of other atoms.
   After the move, atoms are remapped into the periodic simulation box if
   needed, and any shrink-wrap boundary conditions (see the
   :doc:`boundary <boundary>` command) are enforced which may change the
   box size.  Other than this effect, this command does not change the
   size or shape of the simulation box.  See the
   :doc:`change_box <change_box>` command if that effect is desired.

.. note::

   Atoms can be moved arbitrarily long distances by this command.
   If the simulation box is non-periodic and shrink-wrapped (see the
   :doc:`boundary <boundary>` command), this can change its size or shape.
   This is not a problem, except that the mapping of processors to the
   simulation box is not changed by this command from its initial 3d
   configuration; see the :doc:`processors <processors>` command.  Thus, if
   the box size/shape changes dramatically, the mapping of processors to
   the simulation box may not end up as optimal as the initial mapping
   attempted to be.


----------


Restrictions
""""""""""""


For a 2d simulation, only rotations around the a vector parallel to
the z-axis are allowed.

Related commands
""""""""""""""""

:doc:`lattice <lattice>`, :doc:`change_box <change_box>`,
:doc:`fix move <fix_move>`

Default
"""""""

The option defaults are units = lattice.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
