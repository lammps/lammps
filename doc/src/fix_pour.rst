.. index:: fix pour

fix pour command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID pour N type seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* pour = style name of this fix command
* N = # of particles to insert
* type = atom type to assign to inserted particles (offset for molecule insertion)
* seed = random # seed (positive integer)
* one or more keyword/value pairs may be appended to args
* keyword = *region* or *diam* or *vol* or *rate* or *dens* or *vel* or *mol* or *rigid* or *shake* or *ignore*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region to use as insertion volume
       *diam* values = dstyle args
         dstyle = *one* or *range* or *poly*
           *one* args = D
             D = single diameter for inserted particles (distance units)
           *range* args = Dlo Dhi
             Dlo,Dhi = range of diameters for inserted particles (distance units)
           *poly* args = Npoly D1 P1 D2 P2 ...
             Npoly = # of (D,P) pairs
             D1,D2,... = diameter for subset of inserted particles (distance units)
             P1,P2,... = percentage of inserted particles with this diameter (0-1)
       *id* values = idflag
         idflag = *max* or *next* = how to choose IDs for inserted particles and molecules
       *vol* values = fraction Nattempt
         fraction = desired volume fraction for filling insertion volume
         Nattempt = max # of insertion attempts per particle
       *rate* value = V
         V = z velocity (3d) or y velocity (2d) at which
             insertion volume moves (velocity units)
       *dens* values = Rholo Rhohi
         Rholo,Rhohi = range of densities for inserted particles (mass/volume units)
       *vel* values (3d) = vxlo vxhi vylo vyhi vz
       *vel* values (2d) = vxlo vxhi vy
         vxlo,vxhi = range of x velocities for inserted particles (velocity units)
         vylo,vyhi = range of y velocities for inserted particles (velocity units)
         vz = z velocity (3d) assigned to inserted particles (velocity units)
         vy = y velocity (2d) assigned to inserted particles (velocity units)
       *mol* value = template-ID
         template-ID = ID of molecule template specified in a separate :doc:`molecule <molecule>` command
       *molfrac* values = f1 f2 ... fN
         f1 to fN = relative probability of creating each of N molecules in template-ID
       *rigid* value = fix-ID
         fix-ID = ID of :doc:`fix rigid/small <fix_rigid>` command
       *shake* value = fix-ID
         fix-ID = ID of :doc:`fix shake <fix_shake>` command
       *ignore* value = none
         skip any line or triangle particles when detecting possible
           overlaps with inserted particles

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all pour 1000 2 29494 region myblock
   fix 2 all pour 10000 1 19985583 region disk vol 0.33 100 rate 1.0 diam range 0.9 1.1
   fix 2 all pour 10000 1 19985583 region disk diam poly 2 0.7 0.4 1.5 0.6
   fix ins all pour 500 1 4767548 vol 0.8 10 region slab mol object rigid myRigid

Description
"""""""""""

Insert finite-size particles or molecules into the simulation box
every few timesteps within a specified region until N particles or
molecules have been inserted.  This is typically used to model the
pouring of granular particles into a container under the influence of
gravity.  For the remainder of this doc page, a single inserted atom
or molecule is referred to as a "particle".

If inserted particles are individual atoms, they are assigned the
specified atom type.  If they are molecules, the type of each atom in
the inserted molecule is specified in the file read by the
:doc:`molecule <molecule>` command, and those values are added to the
specified atom type.  E.g. if the file specifies atom types 1,2,3, and
those are the atom types you want for inserted molecules, then specify
*type* = 0.  If you specify *type* = 2, the in the inserted molecule
will have atom types 3,4,5.

All atoms in the inserted particle are assigned to two groups: the
default group "all" and the group specified in the fix pour command
(which can also be "all").

This command must use the *region* keyword to define an insertion
volume.  The specified region must have been previously defined with a
:doc:`region <region>` command.  It must be of type *block* or a z-axis
*cylinder* and must be defined with side = *in*\ .  The cylinder style
of region can only be used with 3d simulations.

Individual atoms are inserted, unless the *mol* keyword is used.  It
specifies a *template-ID* previously defined using the
:doc:`molecule <molecule>` command, which reads a file that defines the
molecule.  The coordinates, atom types, center-of-mass, moments of
inertia, etc, as well as any bond/angle/etc and special neighbor
information for the molecule can be specified in the molecule file.
See the :doc:`molecule <molecule>` command for details.  The only
settings required to be in this file are the coordinates and types of
atoms in the molecule.

If the molecule template contains more than one molecule, the relative
probability of depositing each molecule can be specified by the
*molfrac* keyword.  N relative probabilities, each from 0.0 to 1.0, are
specified, where N is the number of molecules in the template.  Each
time a molecule is inserted, a random number is used to sample from
the list of relative probabilities.  The N values must sum to 1.0.

If you wish to insert molecules via the *mol* keyword, that will be
treated as rigid bodies, use the *rigid* keyword, specifying as its
value the ID of a separate :doc:`fix rigid/small <fix_rigid>`
command which also appears in your input script.

.. note::

   If you wish the new rigid molecules (and other rigid molecules)
   to be thermostatted correctly via :doc:`fix rigid/small/nvt <fix_rigid>`
   or :doc:`fix rigid/small/npt <fix_rigid>`, then you need to use the
   "fix_modify dynamic/dof yes" command for the rigid fix.  This is to
   inform that fix that the molecule count will vary dynamically.

If you wish to insert molecules via the *mol* keyword, that will have
their bonds or angles constrained via SHAKE, use the *shake* keyword,
specifying as its value the ID of a separate :doc:`fix shake <fix_shake>` command which also appears in your input script.

Each timestep particles are inserted, they are placed randomly inside
the insertion volume so as to mimic a stream of poured particles.  If
they are molecules they are also oriented randomly.  Each atom in the
particle is tested for overlaps with existing particles, including
effects due to periodic boundary conditions if applicable.  If an
overlap is detected, another random insertion attempt is made; see the
*vol* keyword discussion below.  The larger the volume of the
insertion region, the more particles that can be inserted at any one
timestep.  Particles are inserted again after enough time has elapsed
that the previously inserted particles fall out of the insertion
volume under the influence of gravity.  Insertions continue every so
many timesteps until the desired # of particles has been inserted.

.. note::

   If you are monitoring the temperature of a system where the
   particle count is changing due to adding particles, you typically
   should use the :doc:`compute_modify dynamic yes <compute_modify>`
   command for the temperature compute you are using.

----------

All other keywords are optional with defaults as shown below.

The *diam* option is only used when inserting atoms and specifies the
diameters of inserted particles.  There are 3 styles: *one*, *range*,
or *poly*\ .  For *one*, all particles will have diameter *D*\ .  For
*range*, the diameter of each particle will be chosen randomly and
uniformly between the specified *Dlo* and *Dhi* bounds.  For *poly*, a
series of *Npoly* diameters is specified.  For each diameter a
percentage value from 0.0 to 1.0 is also specified.  The *Npoly*
percentages must sum to 1.0.  For the example shown above with "diam 2
0.7 0.4 1.5 0.6", all inserted particles will have a diameter of 0.7
or 1.5.  40% of the particles will be small; 60% will be large.

Note that for molecule insertion, the diameters of individual atoms in
the molecule can be specified in the file read by the
:doc:`molecule <molecule>` command.  If not specified, the diameter of
each atom in the molecule has a default diameter of 1.0.

The *id* option has two settings which are used to determine the atom
or molecule IDs to assign to inserted particles/molecules.  In both
cases a check is done of the current system to find the maximum
current atom and molecule ID of any existing particle.  Newly inserted
particles and molecules are assigned IDs that increment those max
values.  For the *max* setting, which is the default, this check is
done at every insertion step, which allows for particles to leave the
system, and their IDs to potentially be re-used.  For the *next*
setting this check is done only once when the fix is specified, which
can be more efficient if you are sure particles will not be added in
some other way.

The *vol* option specifies what volume fraction of the insertion
volume will be filled with particles.  For particles with a size
specified by the *diam range* keyword, they are assumed to all be of
maximum diameter *Dhi* for purposes of computing their contribution to
the volume fraction.

The higher the volume fraction value, the more particles are inserted
each timestep.  Since inserted particles cannot overlap, the maximum
volume fraction should be no higher than about 0.6.  Each timestep
particles are inserted, LAMMPS will make up to a total of M tries to
insert the new particles without overlaps, where M = # of inserted
particles \* Nattempt.  If LAMMPS is unsuccessful at completing all
insertions, it prints a warning.

The *dens* and *vel* options enable inserted particles to have a range
of densities or xy velocities.  The specific values for a particular
inserted particle will be chosen randomly and uniformly between the
specified bounds.  Internally, the density value for a particle is
converted to a mass, based on the radius (volume) of the particle.
The *vz* or *vy* value for option *vel* assigns a z-velocity (3d) or
y-velocity (2d) to each inserted particle.

The *rate* option moves the insertion volume in the z direction (3d)
or y direction (2d).  This enables pouring particles from a
successively higher height over time.

The *ignore* option is useful when running a simulation that used line
segment (2d) or triangle (3d) particles, typically to define
boundaries for spherical granular particles to interact with.  See the
:doc:`atom_style line or tri <atom_style>` command for details.  Lines
and triangles store their size, and if the size is large it may
overlap (in a spherical sense) with the insertion region, even if the
line/triangle is oriented such that there is no actual overlap.  This
can prevent particles from being inserted.  The *ignore* keyword
causes the overlap check to skip any line or triangle particles.
Obviously you should only use it if there is in fact no overlap of the
line or triangle particles with the insertion region.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  This means you must be careful when restarting a
pouring simulation, when the restart file was written in the middle of
the pouring operation.  Specifically, you should use a new fix pour
command in the input script for the restarted simulation that
continues the operation.  You will need to adjust the arguments of the
original fix pour command to do this.

Also note that because the state of the random number generator is not
saved in restart files, you cannot do "exact" restarts with this fix,
where the simulation continues on the same as if no restart had taken
place.  However, in a statistical sense, a restarted simulation should
produce the same behavior if you adjust the fix pour parameters
appropriately.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

For 3d simulations, a gravity fix in the -z direction must be defined
for use in conjunction with this fix.  For 2d simulations, gravity
must be defined in the -y direction.

The specified insertion region cannot be a "dynamic" region, as
defined by the :doc:`region <region>` command.

Related commands
""""""""""""""""

:doc:`fix deposit <fix_deposit>`, :doc:`fix gravity <fix_gravity>`,
:doc:`region <region>`

Default
"""""""

Insertions are performed for individual particles, i.e. no *mol*
setting is defined.  If the *mol* keyword is used, the default for
*molfrac* is an equal probabilities for all molecules in the template.
Additional option defaults are diam = one 1.0, dens = 1.0 1.0, vol =
0.25 50, rate = 0.0, vel = 0.0 0.0 0.0 0.0 0.0 (for 3d), vel = 0.0 0.0 0.0
(for 2d), and id = max.
