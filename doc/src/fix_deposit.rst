.. index:: fix deposit

fix deposit command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID deposit N type M seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* deposit = style name of this fix command
* N = # of atoms or molecules to insert
* type = atom type to assign to inserted atoms (offset for molecule insertion)
* M = insert a single atom or molecule every M steps
* seed = random # seed (positive integer)
* one or more keyword/value pairs may be appended to args
* keyword = *region* or *id* or *global* or *local* or *near* or *gaussian* or *attempt* or *rate* or *vx* or *vy* or *vz* or *target* or *mol* or *molfrac* or *rigid* or *shake* or *orient* or *units*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region to use as insertion volume
       *id* value = *max* or *next*
         max = atom ID for new atom(s) is max ID of all current atoms plus one
         next = atom ID for new atom(s) increments by one for every deposition
       *global* values = lo hi
         lo,hi = put new atom/molecule a distance lo-hi above all other atoms (distance units)
       *local* values = lo hi delta
         lo,hi = put new atom/molecule a distance lo-hi above any nearby atom beneath it (distance units)
         delta = lateral distance within which a neighbor is considered "nearby" (distance units)
       *near* value = R
         R = only insert atom/molecule if further than R from existing particles (distance units)
       *gaussian* values = xmid ymid zmid sigma
         xmid,ymid,zmid = center of the gaussian distribution (distance units)
         sigma = width of gaussian distribution (distance units)
       *attempt* value = Q
         Q = attempt a single insertion up to Q times
       *rate* value = V
         V = z velocity (y in 2d) at which insertion volume moves (velocity units)
       *vx* values = vxlo vxhi
         vxlo,vxhi = range of x velocities for inserted atom/molecule (velocity units)
       *vy* values = vylo vyhi
         vylo,vyhi = range of y velocities for inserted atom/molecule (velocity units)
       *vz* values = vzlo vzhi
         vzlo,vzhi = range of z velocities for inserted atom/molecule (velocity units)
       *target* values = tx ty tz
         tx,ty,tz = location of target point (distance units)
       *mol* value = template-ID
         template-ID = ID of molecule template specified in a separate :doc:`molecule <molecule>` command
       *molfrac* values = f1 f2 ... fN
         f1 to fN = relative probability of creating each of N molecules in template-ID
       *rigid* value = fix-ID
         fix-ID = ID of :doc:`fix rigid/small <fix_rigid>` command
       *shake* value = fix-ID
         fix-ID = ID of :doc:`fix shake <fix_shake>` command
       *orient* values = rx ry rz
         rx,ry,rz = vector to randomly rotate an inserted molecule around
       *units* value = *lattice* or *box*
         lattice = the geometry is defined in lattice units
         box = the geometry is defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all deposit 1000 2 100 29494 region myblock local 1.0 1.0 1.0 units box
   fix 2 newatoms deposit 10000 1 500 12345 region disk near 2.0 vz -1.0 -0.8
   fix 4 sputter deposit 1000 2 500 12235 region sphere vz -1.0 -1.0 target 5.0 5.0 0.0 units lattice
   fix 5 insert deposit 200 2 100 777 region disk gaussian 5.0 5.0 9.0 1.0 units box

Description
"""""""""""

Insert a single atom or molecule into the simulation domain every M
timesteps until N atoms or molecules have been inserted.  This is
useful for simulating deposition onto a surface.  For the remainder of
this doc page, a single inserted atom or molecule is referred to as a
"particle".

If inserted particles are individual atoms, they are assigned the
specified atom type.  If they are molecules, the type of each atom in
the inserted molecule is specified in the file read by the
:doc:`molecule <molecule>` command, and those values are added to the
specified atom type.  E.g. if the file specifies atom types 1,2,3, and
those are the atom types you want for inserted molecules, then specify
*type* = 0.  If you specify *type* = 2, the in the inserted molecule
will have atom types 3,4,5.

All atoms in the inserted particle are assigned to two groups: the
default group "all" and the group specified in the fix deposit command
(which can also be "all").

If you are computing temperature values which include inserted
particles, you will want to use the
:doc:`compute_modify <compute_modify>` dynamic option, which ensures the
current number of atoms is used as a normalizing factor each time the
temperature is computed.

Care must be taken that inserted particles are not too near existing
atoms, using the options described below.  When inserting particles
above a surface in a non-periodic box (see the
:doc:`boundary <boundary>` command), the possibility of a particle
escaping the surface and flying upward should be considered, since the
particle may be lost or the box size may grow infinitely large.  A
:doc:`fix wall/reflect <fix_wall_reflect>` command can be used to
prevent this behavior.  Note that if a shrink-wrap boundary is used,
it is OK to insert the new particle outside the box, however the box
will immediately be expanded to include the new particle. When
simulating a sputtering experiment it is probably more realistic to
ignore those atoms using the :doc:`thermo_modify <thermo_modify>`
command with the *lost ignore* option and a fixed
:doc:`boundary <boundary>`.

The fix deposit command must use the *region* keyword to define an
insertion volume.  The specified region must have been previously
defined with a :doc:`region <region>` command.  It must be defined with
side = *in*\ .

.. note::

   LAMMPS checks that the specified region is wholly inside the
   simulation box.  It can do this correctly for orthonormal simulation
   boxes.  However for :doc:`triclinic boxes <Howto_triclinic>`, it only
   tests against the larger orthonormal box that bounds the tilted
   simulation box.  If the specified region includes volume outside the
   tilted box, then an insertion will likely fail, leading to a "lost
   atoms" error.  Thus for triclinic boxes you should ensure the
   specified region is wholly inside the simulation box.

The locations of inserted particles are taken from uniform distributed
random numbers, unless the *gaussian* keyword is used. Then the
individual coordinates are taken from a gaussian distribution of
width *sigma* centered on *xmid,ymid,zmid*\ .

Individual atoms are inserted, unless the *mol* keyword is used.  It
specifies a *template-ID* previously defined using the
:doc:`molecule <molecule>` command, which reads files that define one or
more molecules.  The coordinates, atom types, charges, etc, as well as
any bond/angle/etc and special neighbor information for the molecule
can be specified in the molecule file.  See the
:doc:`molecule <molecule>` command for details.  The only settings
required to be in each file are the coordinates and types of atoms in
the molecule.

If the molecule template contains more than one molecule, the relative
probability of depositing each molecule can be specified by the
*molfrac* keyword.  N relative probabilities, each from 0.0 to 1.0, are
specified, where N is the number of molecules in the template.  Each
time a molecule is deposited, a random number is used to sample from
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

Each timestep a particle is inserted, the coordinates for its atoms
are chosen as follows.  For insertion of individual atoms, the
"position" referred to in the following description is the coordinate
of the atom.  For insertion of molecule, the "position" is the
geometric center of the molecule; see the :doc:`molecule <molecule>` doc
page for details.  A random rotation of the molecule around its center
point is performed, which determines the coordinates all the
individual atoms.

A random position within the region insertion volume is generated.  If
neither the *global* or *local* keyword is used, the random position
is the trial position.  If the *global* keyword is used, the random
x,y values are used, but the z position of the new particle is set
above the highest current atom in the simulation by a distance
randomly chosen between lo/hi.  (For a 2d simulation, this is done for
the y position.)  If the *local* keyword is used, the z position is
set a distance between lo/hi above the highest current atom in the
simulation that is "nearby" the chosen x,y position.  In this context,
"nearby" means the lateral distance (in x,y) between the new and old
particles is less than the *delta* setting.

Once a trial x,y,z position has been selected, the insertion is only
performed if no current atom in the simulation is within a distance R
of any atom in the new particle, including the effect of periodic
boundary conditions if applicable.  R is defined by the *near*
keyword.  Note that the default value for R is 0.0, which will allow
atoms to strongly overlap if you are inserting where other atoms are
present.  This distance test is performed independently for each atom
in an inserted molecule, based on the randomly rotated configuration
of the molecule.  If this test fails, a new random position within the
insertion volume is chosen and another trial is made.  Up to Q
attempts are made.  If the particle is not successfully inserted,
LAMMPS prints a warning message.

.. note::

   If you are inserting finite size particles or a molecule or
   rigid body consisting of finite-size particles, then you should
   typically set R larger than the distance at which any inserted
   particle may overlap with either a previously inserted particle or an
   existing particle.  LAMMPS will issue a warning if R is smaller than
   this value, based on the radii of existing and inserted particles.

The *rate* option moves the insertion volume in the z direction (3d)
or y direction (2d).  This enables particles to be inserted from a
successively higher height over time.  Note that this parameter is
ignored if the *global* or *local* keywords are used, since those
options choose a z-coordinate for insertion independently.

The vx, vy, and vz components of velocity for the inserted particle
are set by sampling a uniform distribution between the bounds set by
the values specified for the *vx*, *vy*, and *vz* keywords. Note that 
normally, new particles should be a assigned a negative vertical 
velocity so that they move towards the surface.  For molecules, the 
same velocity is given to every particle (no rotation or bond vibration).

If the *target* option is used, the velocity vector of the inserted
particle is changed so that it points from the insertion position
towards the specified target point.  The magnitude of the velocity is
unchanged.  This can be useful, for example, for simulating a
sputtering process.  E.g. the target point can be far away, so that
all incident particles strike the surface as if they are in an
incident beam of particles at a prescribed angle.

The *orient* keyword is only used when molecules are deposited.  By
default, each molecule is inserted at a random orientation.  If this
keyword is specified, then (rx,ry,rz) is used as an orientation
vector, and each inserted molecule is rotated around that vector with
a random value from zero to 2*PI.  For a 2d simulation, rx = ry = 0.0
is required, since rotations can only be performed around the z axis.

The *id* keyword determines how atom IDs and molecule IDs are assigned
to newly deposited particles.  Molecule IDs are only assigned if
molecules are being inserted.  For the *max* setting, the atom and
molecule IDs of all current atoms are checked.  Atoms in the new
particle are assigned IDs starting with the current maximum plus one.
If a molecule is inserted it is assigned an ID = current maximum plus
one.  This means that if particles leave the system, the new IDs may
replace the lost ones.  For the *next* setting, the maximum ID of any
atom and molecule is stored at the time the fix is defined.  Each time
a new particle is added, this value is incremented to assign IDs to
the new atom(s) or molecule.  Thus atom and molecule IDs for deposited
particles will be consecutive even if particles leave the system over
time.

The *units* keyword determines the meaning of the distance units used
for the other deposition parameters.  A *box* value selects standard
distance units as defined by the :doc:`units <units>` command,
e.g. Angstroms for units = real or metal.  A *lattice* value means the
distance units are in lattice spacings.  The :doc:`lattice <lattice>`
command must have been previously used to define the lattice spacing.
Note that the units choice affects all the keyword values that have
units of distance or velocity.

.. note::

   If you are monitoring the temperature of a system where the atom
   count is changing due to adding particles, you typically should use
   the :doc:`compute_modify dynamic yes <compute_modify>` command for the
   temperature compute you are using.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the deposition to :doc:`binary restart files <restart>`.  This includes information about how many
particles have been deposited, the random number generator seed, the
next timestep for deposition, etc.  See the
:doc:`read_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file, so that the
operation of the fix continues in an uninterrupted fashion.

.. note::

   For this to work correctly, the timestep must **not** be changed
   after reading the restart with :doc:`reset_timestep <reset_timestep>`.
   The fix will try to detect it and stop with an error.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The specified insertion region cannot be a "dynamic" region, as
defined by the :doc:`region <region>` command.

Related commands
""""""""""""""""

:doc:`fix pour <fix_pour>`, :doc:`region <region>`

Default
"""""""

Insertions are performed for individual atoms, i.e. no *mol* setting
is defined.  If the *mol* keyword is used, the default for *molfrac*
is an equal probabilities for all molecules in the template.
Additional option defaults are id = max, delta = 0.0, near = 0.0,
attempt = 10, rate = 0.0, vx = 0.0 0.0, vy = 0.0 0.0, vz = 0.0 0.0,
and units = lattice.
