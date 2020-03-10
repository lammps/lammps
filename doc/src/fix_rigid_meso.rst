.. index:: fix rigid/meso

fix rigid/meso command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rigid/meso bodystyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* rigid/meso = style name of this fix command
* bodystyle = *single* or *molecule* or *group*

  .. parsed-literal::

       *single* args = none
       *molecule* args = none
       *custom* args = *i_propname* or *v_varname*
         i_propname = an integer property defined via fix property/atom
         v_varname  = an atom-style or atomfile-style variable
       *group* args = N groupID1 groupID2 ...
         N = # of groups
         groupID1, groupID2, ... = list of N group IDs

* zero or more keyword/value pairs may be appended
* keyword = *reinit* or *force* or *torque* or *infile*

  .. parsed-literal::

       *reinit* = *yes* or *no*
       *force* values = M xflag yflag zflag
         M = which rigid body from 1-Nbody (see asterisk form below)
         xflag,yflag,zflag = off/on if component of center-of-mass force is active
       *torque* values = M xflag yflag zflag
         M = which rigid body from 1-Nbody (see asterisk form below)
         xflag,yflag,zflag = off/on if component of center-of-mass torque is active
       *infile* filename
         filename = file with per-body values of mass, center-of-mass, moments of inertia

Examples
""""""""

.. parsed-literal::

   fix 1 ellipsoid rigid/meso single
   fix 1 rods      rigid/meso molecule
   fix 1 spheres   rigid/meso single force 1 off off on
   fix 1 particles rigid/meso molecule force 1\*5 off off off force 6\*10 off off on
   fix 2 spheres   rigid/meso group 3 sphere1 sphere2 sphere3 torque \* off off off

Description
"""""""""""

Treat one or more sets of mesoscopic SPH/SDPD particles as independent
rigid bodies.  This means that each timestep the total force and torque
on each rigid body is computed as the sum of the forces and torques on
its constituent particles.  The coordinates and velocities of the
particles in each body are then updated so that the body moves and
rotates as a single entity using the methods described in the paper by
:ref:`(Miller) <Miller>`. Density and internal energy of the particles will
also be updated. This is implemented by creating internal data structures
for each rigid body and performing time integration on these data
structures.  Positions and velocities of the constituent particles are
regenerated from the rigid body data structures in every time step. This
restricts which operations and fixes can be applied to rigid bodies. See
below for a detailed discussion.

The operation of this fix is exactly like that described by the
:doc:`fix rigid/nve <fix_rigid>` command, except that particles' density,
internal energy and extrapolated velocity are also updated.

.. note::

   You should not update the particles in rigid bodies via other
   time-integration fixes (e.g. :doc:`fix meso <fix_meso>`,
   :doc:`fix meso/stationary <fix_meso_stationary>`), or you will have conflicting
   updates to positions and velocities resulting in unphysical behavior in most
   cases. When performing a hybrid simulation with some atoms in rigid bodies,
   and some not, a separate time integration fix like :doc:`fix meso <fix_meso>`
   should be used for the non-rigid particles.

.. note::

   These fixes are overkill if you simply want to hold a collection
   of particles stationary or have them move with a constant velocity. To
   hold particles stationary use :doc:`fix meso/stationary <fix_meso_stationary>` instead. If you would like to
   move particles with a constant velocity use :doc:`fix meso/move <fix_meso_move>`.

.. warning::

   The aggregate properties of each rigid body are
   calculated at the start of a simulation run and are maintained in
   internal data structures. The properties include the position and
   velocity of the center-of-mass of the body, its moments of inertia, and
   its angular momentum.  This is done using the properties of the
   constituent particles of the body at that point in time (or see the *infile*
   keyword option).  Thereafter, changing these properties of individual
   particles in the body will have no effect on a rigid body's dynamics, unless
   they effect any computation of per-particle forces or torques. If the
   keyword *reinit* is set to *yes* (the default), the rigid body data
   structures will be recreated at the beginning of each *run* command;
   if the keyword *reinit* is set to *no*\ , the rigid body data structures
   will be built only at the very first *run* command and maintained for
   as long as the rigid fix is defined. For example, you might think you
   could displace the particles in a body or add a large velocity to each particle
   in a body to make it move in a desired direction before a 2nd run is
   performed, using the :doc:`set <set>` or
   :doc:`displace_atoms <displace_atoms>` or :doc:`velocity <velocity>`
   commands.  But these commands will not affect the internal attributes
   of the body unless *reinit* is set to *yes*\ . With *reinit* set to *no*
   (or using the *infile* option, which implies *reinit* *no*\ ) the position
   and velocity of individual particles in the body will be reset when time
   integration starts again.

----------

Each rigid body must have two or more particles.  A particle can belong
to at most one rigid body.  Which particles are in which bodies can be
defined via several options.

For bodystyle *single* the entire fix group of particles is treated as
one rigid body.

For bodystyle *molecule*\ , particles are grouped into rigid bodies by their
respective molecule IDs: each set of particles in the fix group with the
same molecule ID is treated as a different rigid body.  Note that particles
with a molecule ID = 0 will be treated as a single rigid body. For a
system with solvent (typically this is particles with molecule ID = 0)
surrounding rigid bodies, this may not be what you want.  Thus you
should be careful to use a fix group that only includes particles you
want to be part of rigid bodies.

Bodystyle *custom* is similar to bodystyle *molecule* except that it
is more flexible in using other per-atom properties to define the sets
of particles that form rigid bodies.  An integer vector defined by the
:doc:`fix property/atom <fix_property_atom>` command can be used.  Or an
:doc:`atom-style or atomfile-style variable <variable>` can be used; the
floating-point value produced by the variable is rounded to an
integer.  As with bodystyle *molecule*\ , each set of particles in the fix
groups with the same integer value is treated as a different rigid
body.  Since fix property/atom vectors and atom-style variables
produce values for all particles, you should be careful to use a fix group
that only includes particles you want to be part of rigid bodies.

For bodystyle *group*\ , each of the listed groups is treated as a
separate rigid body.  Only particles that are also in the fix group are
included in each rigid body.

.. note::

   To compute the initial center-of-mass position and other
   properties of each rigid body, the image flags for each particle in the
   body are used to "unwrap" the particle coordinates.  Thus you must
   insure that these image flags are consistent so that the unwrapping
   creates a valid rigid body (one where the particles are close together)
   , particularly if the particles in a single rigid body straddle a
   periodic boundary.  This means the input data file or restart file must
   define the image flags for each particle consistently or that you have
   used the :doc:`set <set>` command to specify them correctly.  If a
   dimension is non-periodic then the image flag of each particle must be
   0 in that dimension, else an error is generated.

By default, each rigid body is acted on by other particles which induce
an external force and torque on its center of mass, causing it to
translate and rotate.  Components of the external center-of-mass force
and torque can be turned off by the *force* and *torque* keywords.
This may be useful if you wish a body to rotate but not translate, or
vice versa, or if you wish it to rotate or translate continuously
unaffected by interactions with other particles.  Note that if you
expect a rigid body not to move or rotate by using these keywords, you
must insure its initial center-of-mass translational or angular
velocity is 0.0. Otherwise the initial translational or angular
momentum, the body has, will persist.

An xflag, yflag, or zflag set to *off* means turn off the component of
force or torque in that dimension.  A setting of *on* means turn on
the component, which is the default.  Which rigid body(s) the settings
apply to is determined by the first argument of the *force* and
*torque* keywords.  It can be an integer M from 1 to Nbody, where
Nbody is the number of rigid bodies defined.  A wild-card asterisk can
be used in place of, or in conjunction with, the M argument to set the
flags for multiple rigid bodies.  This takes the form "\*" or "\*n" or
"n\*" or "m\*n".  If N = the number of rigid bodies, then an asterisk
with no numeric values means all bodies from 1 to N.  A leading
asterisk means all bodies from 1 to n (inclusive).  A trailing
asterisk means all bodies from n to N (inclusive).  A middle asterisk
means all bodies from m to n (inclusive).  Note that you can use the
*force* or *torque* keywords as many times as you like.  If a
particular rigid body has its component flags set multiple times, the
settings from the final keyword are used.

For computational efficiency, you should typically define one fix
rigid/meso command which includes all the desired rigid bodies. LAMMPS
will allow multiple rigid/meso fixes to be defined, but it is more
expensive.

----------

The keyword/value option pairs are used in the following ways.

The *reinit* keyword determines, whether the rigid body properties
are re-initialized between run commands. With the option *yes* (the
default) this is done, with the option *no* this is not done. Turning
off the re-initialization can be helpful to protect rigid bodies against
unphysical manipulations between runs or when properties cannot be
easily re-computed (e.g. when read from a file). When using the *infile*
keyword, the *reinit* option is automatically set to *no*\ .

----------

The *infile* keyword allows a file of rigid body attributes to be read
in from a file, rather then having LAMMPS compute them.  There are 5
such attributes: the total mass of the rigid body, its center-of-mass
position, its 6 moments of inertia, its center-of-mass velocity, and
the 3 image flags of the center-of-mass position.  For rigid bodies
consisting of point particles or non-overlapping finite-size
particles, LAMMPS can compute these values accurately.  However, for
rigid bodies consisting of finite-size particles which overlap each
other, LAMMPS will ignore the overlaps when computing these 4
attributes.  The amount of error this induces depends on the amount of
overlap.  To avoid this issue, the values can be pre-computed
(e.g. using Monte Carlo integration).

The format of the file is as follows.  Note that the file does not
have to list attributes for every rigid body integrated by fix rigid.
Only bodies which the file specifies will have their computed
attributes overridden.  The file can contain initial blank lines or
comment lines starting with "#" which are ignored.  The first
non-blank, non-comment line should list N = the number of lines to
follow.  The N successive lines contain the following information:

.. parsed-literal::

   ID1 masstotal xcm ycm zcm ixx iyy izz ixy ixz iyz vxcm vycm vzcm lx ly lz ixcm iycm izcm
   ID2 masstotal xcm ycm zcm ixx iyy izz ixy ixz iyz vxcm vycm vzcm lx ly lz ixcm iycm izcm
   ...
   IDN masstotal xcm ycm zcm ixx iyy izz ixy ixz iyz vxcm vycm vzcm lx ly lz ixcm iycm izcm

The rigid body IDs are all positive integers.  For the *single*
bodystyle, only an ID of 1 can be used.  For the *group* bodystyle,
IDs from 1 to Ng can be used where Ng is the number of specified
groups.  For the *molecule* bodystyle, use the molecule ID for the
atoms in a specific rigid body as the rigid body ID.

The masstotal and center-of-mass coordinates (xcm,ycm,zcm) are
self-explanatory.  The center-of-mass should be consistent with what
is calculated for the position of the rigid body with all its atoms
unwrapped by their respective image flags.  If this produces a
center-of-mass that is outside the simulation box, LAMMPS wraps it
back into the box.

The 6 moments of inertia (ixx,iyy,izz,ixy,ixz,iyz) should be the
values consistent with the current orientation of the rigid body
around its center of mass.  The values are with respect to the
simulation box XYZ axes, not with respect to the principal axes of the
rigid body itself.  LAMMPS performs the latter calculation internally.

The (vxcm,vycm,vzcm) values are the velocity of the center of mass.
The (lx,ly,lz) values are the angular momentum of the body.  The
(vxcm,vycm,vzcm) and (lx,ly,lz) values can simply be set to 0 if you
wish the body to have no initial motion.

The (ixcm,iycm,izcm) values are the image flags of the center of mass
of the body.  For periodic dimensions, they specify which image of the
simulation box the body is considered to be in.  An image of 0 means
it is inside the box as defined.  A value of 2 means add 2 box lengths
to get the true value.  A value of -1 means subtract 1 box length to
get the true value.  LAMMPS updates these flags as the rigid bodies
cross periodic boundaries during the simulation.

.. note::

   If you use the *infile* keyword and write restart
   files during a simulation, then each time a restart file is written,
   the fix also write an auxiliary restart file with the name
   rfile.rigid, where "rfile" is the name of the restart file,
   e.g. tmp.restart.10000 and tmp.restart.10000.rigid.  This auxiliary
   file is in the same format described above.  Thus it can be used in a
   new input script that restarts the run and re-specifies a rigid fix
   using an *infile* keyword and the appropriate filename.  Note that the
   auxiliary file will contain one line for every rigid body, even if the
   original file only listed a subset of the rigid bodies.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information is written to :doc:`binary restart files <restart>`.
If the *infile* keyword is used, an auxiliary file is written out
with rigid body information each time a restart file is written, as
explained above for the *infile* keyword.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes a global array of values which can be accessed by
various :doc:`output commands <Howto_output>`.

The number of rows in the array is equal to the number of rigid
bodies.  The number of columns is 28.  Thus for each rigid body, 28
values are stored: the xyz coords of the center of mass (COM), the xyz
components of the COM velocity, the xyz components of the force acting
on the COM, the components of the 4-vector quaternion representing the
orientation of the rigid body, the xyz components of the angular velocity
of the body around its COM, the xyz components of the torque acting on the
COM, the 3 principal components of the moment of inertia, the xyz components
of the angular momentum of the body around its COM, and the xyz image
flags of the COM.

The center of mass (COM) for each body is similar to unwrapped
coordinates written to a dump file.  It will always be inside (or
slightly outside) the simulation box.  The image flags have the same
meaning as image flags for particle positions (see the "dump" command).
This means you can calculate the unwrapped COM by applying the image
flags to the COM, the same as when unwrapped coordinates are written
to a dump file.

The force and torque values in the array are not affected by the
*force* and *torque* keywords in the fix rigid command; they reflect
values before any changes are made by those keywords.

The ordering of the rigid bodies (by row in the array) is as follows.
For the *single* keyword there is just one rigid body.  For the
*molecule* keyword, the bodies are ordered by ascending molecule ID.
For the *group* keyword, the list of group IDs determines the ordering
of bodies.

The array values calculated by this fix are "intensive", meaning they
are independent of the number of particles in the simulation.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

This fix is part of the USER-SDPD package and also depends on the RIGID
package.  It is only enabled if LAMMPS was built with both packages. See
the :doc:`Build package <Build_package>` doc page for more info.

This fix requires that atoms store density and internal energy as
defined by the :doc:`atom_style meso <atom_style>` command.

All particles in the group must be mesoscopic SPH/SDPD particles.

Related commands
""""""""""""""""

:doc:`fix meso/move <fix_meso_move>`, :doc:`fix rigid <fix_rigid>`,
:doc:`neigh_modify exclude <neigh_modify>`

Default
"""""""

The option defaults are force \* on on on and torque \* on on on,
meaning all rigid bodies are acted on by center-of-mass force and
torque. Also reinit = yes.

----------

.. _Miller:

**(Miller)** Miller, Eleftheriou, Pattnaik, Ndirango, and Newns,
J Chem Phys, 116, 8649 (2002).
