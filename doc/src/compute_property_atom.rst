.. index:: compute property/atom

compute property/atom command
=============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID property/atom input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* property/atom = style name of this compute command
* input = one or more atom attributes

  .. parsed-literal::

       possible attributes = id, mol, proc, type, mass,
                             x, y, z, xs, ys, zs, xu, yu, zu, ix, iy, iz,
                             vx, vy, vz, fx, fy, fz,
                             q, mux, muy, muz, mu,
                             spx, spy, spz, sp, fmx, fmy, fmz,
                             nbonds,
                             radius, diameter, omegax, omegay, omegaz,
                             angmomx, angmomy, angmomz,
                             shapex,shapey, shapez,
                             quatw, quati, quatj, quatk, tqx, tqy, tqz,
                             end1x, end1y, end1z, end2x, end2y, end2z,
                             corner1x, corner1y, corner1z,
                             corner2x, corner2y, corner2z,
                             corner3x, corner3y, corner3z,
                             i_name, d_name, i2_name[I], d2_name[I],
                             vfrac, s0, spin, eradius, ervel, erforce,
                             rho, drho, e, de, cv, buckling,

  .. parsed-literal::

           *id* = atom ID
           *mol* = molecule ID
           *proc* = ID of processor that owns atom
           *type* = atom type
           *mass* = atom mass
           *x,y,z* = unscaled atom coordinates
           *xs,ys,zs* = scaled atom coordinates
           *xu,yu,zu* = unwrapped atom coordinates
           *ix,iy,iz* = box image that the atom is in
           *vx,vy,vz* = atom velocities
           *fx,fy,fz* = forces on atoms
           *q* = atom charge
           *mux,muy,muz* = orientation of dipole moment of atom
           *mu* = magnitude of dipole moment of atom
           *spx, spy, spz* = direction of the atomic magnetic spin
           *sp* = magintude of atomic magnetic spin moment
           *fmx, fmy, fmz* = magnetic force
           *nbonds* = number of bonds assigned to an atom
           *radius,diameter* = radius,diameter of spherical particle
           *omegax,omegay,omegaz* = angular velocity of spherical particle
           *angmomx,angmomy,angmomz* = angular momentum of aspherical particle
           *shapex,shapey,shapez* = 3 diameters of aspherical particle
           *quatw,quati,quatj,quatk* = quaternion components for aspherical or body particles
           *tqx,tqy,tqz* = torque on finite-size particles
           *end12x, end12y, end12z* = end points of line segment
           *corner123x, corner123y, corner123z* = corner points of triangle
           *i_name* = custom integer vector with name
           *d_name* = custom floating point vector with name
           *i2_name[I]* = Ith column of custom integer array with name
           *d2_name[I]* = Ith column of custom floating-point array with name

  .. parsed-literal::

           PERI package per-atom properties:
           vfrac = volume fraction
           s0 = max stretch of any bond a particle is part of

  .. parsed-literal::

           EFF and AWPMD package per-atom properties:
           spin = electron spin
           eradius = electron radius
           ervel = electron radial velocity
           erforce = electron radial force

  .. parsed-literal::

           SPH package per-atom properties:
           rho = density of SPH particles
           drho = change in density
           e = energy
           de = change in thermal energy
           cv = heat capacity

  .. parsed-literal::

           MESONT package per-atom properties:
           buckling = buckling flag used in mesoscopic simulation of nanotubes

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all property/atom xs vx fx mux
   compute 2 all property/atom type
   compute 1 all property/atom ix iy iz
   compute 3 all property/atom sp spx spy spz
   compute 1 all property/atom i_myFlag d_Sxyz[1] d_Sxyz[3]

Description
"""""""""""

Define a computation that simply stores atom attributes for each atom
in the group.  This is useful so that the values can be used by other
:doc:`output commands <Howto_output>` that take computes as inputs.
See for example, the :doc:`compute reduce <compute_reduce>`, :doc:`fix
ave/atom <fix_ave_atom>`, :doc:`fix ave/histo <fix_ave_histo>`,
:doc:`fix ave/chunk <fix_ave_chunk>`, and :doc:`atom-style variable
<variable>` commands.

The list of possible attributes is essentially the same as that used
by the :doc:`dump custom <dump>` command, which describes their
meaning, with some additional quantities that are only defined for
certain :doc:`atom styles <atom_style>`.  The goal of this augmented
list gives an input script access to any per-atom quantity stored by
LAMMPS.

The values are stored in a per-atom vector or array as discussed
below.  Zeroes are stored for atoms not in the specified group or for
quantities that are not defined for a particular particle in the group
(e.g. *shapex* if the particle is not an ellipsoid).

Attributes *i_name*, *d_name*, *i2_name*, *d2_name* refer to custom
per-atom integer and floating-point vectors or arrays that have been
added via the :doc:`fix property/atom <fix_property_atom>` command.
When that command is used specific names are given to each attribute
which are the "name" portion of these attributes.  For arrays *i2_name*
and *d2_name*, the column of the array must also be included following
the name in brackets: e.g. d2_xyz[2], i2_mySpin[3].

The additional quantities only accessible via this command, and not
directly via the :doc:`dump custom <dump>` command, are as follows.

*Nbonds* is available for all molecular atom styles and refers to the
number of explicit bonds assigned to an atom.  Note that if the
:doc:`newton bond <newton>` command is set to *on*\ , which is the
default, then every bond in the system is assigned to only one of the
two atoms in the bond.  Thus a bond between atoms I,J may be tallied
for either atom I or atom J.  If :doc:`newton bond off <newton>` is
set, it will be tallied with both atom I and atom J.

*Shapex*, *shapey*, and *shapez* are defined for ellipsoidal particles
and define the 3d shape of each particle.

*Quatw*, *quati*, *quatj*, and *quatk* are defined for ellipsoidal
particles and body particles and store the 4-vector quaternion
representing the orientation of each particle.  See the :doc:`set <set>`
command for an explanation of the quaternion vector.

*End1x*, *end1y*, *end1z*, *end2x*, *end2y*, *end2z*, are
defined for line segment particles and define the end points of each
line segment.

*Corner1x*, *corner1y*, *corner1z*, *corner2x*, *corner2y*,
*corner2z*, *corner3x*, *corner3y*, *corner3z*, are defined for
triangular particles and define the corner points of each triangle.

In addition, the various per-atom quantities listed above for specific
packages are only accessible by this command.

Output info
"""""""""""

This compute calculates a per-atom vector or per-atom array depending
on the number of input values.  If a single input is specified, a
per-atom vector is produced.  If two or more inputs are specified, a
per-atom array is produced where the number of columns = the number of
inputs.  The vector or array can be accessed by any command that uses
per-atom values from a compute as input.  See the :doc:`Howto output
<Howto_output>` page for an overview of LAMMPS output options.

The vector or array values will be in whatever :doc:`units <units>` the
corresponding attribute is in, e.g. velocity units for vx, charge
units for q, etc.

For the spin quantities, sp is in the units of the Bohr magneton, spx,
spy, and spz are unitless quantities, and fmx, fmy and fmz are given
in rad/THz.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dump custom <dump>`, :doc:`compute reduce <compute_reduce>`,
:doc::doc:`fix ave/atom <fix_ave_atom>`, :doc:`fix ave/chunk
:doc:<fix_ave_chunk>`, `fix property/atom <fix_property_atom>`

Default
"""""""

none
