.. index:: fix srd

fix srd command
===============

Syntax
""""""


.. parsed-literal::

   fix ID group-ID srd N groupbig-ID Tsrd hgrid seed keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* srd = style name of this fix command
* N = reset SRD particle velocities every this many timesteps
* groupbig-ID = ID of group of large particles that SRDs interact with
* Tsrd = temperature of SRD particles (temperature units)
* hgrid = grid spacing for SRD grouping (distance units)
* seed = random # seed (positive integer)

* zero or more keyword/value pairs may be appended
* keyword = *lamda* or *collision* or *overlap* or *inside* or *exact* or *radius* or *bounce* or *search* or *cubic* or *shift* or *tstat* or *rescale*
  
  .. parsed-literal::
  
       *lamda* value = mean free path of SRD particles (distance units)
       *collision* value = *noslip* or *slip* = collision model
       *overlap* value = *yes* or *no* = whether big particles may overlap
       *inside* value = *error* or *warn* or *ignore* = how SRD particles which end up inside a big particle are treated
       *exact* value = *yes* or *no*
       *radius* value = rfactor = scale collision radius by this factor
       *bounce* value = Nbounce = max # of collisions an SRD particle can undergo in one timestep
       *search* value = sgrid = grid spacing for collision partner searching (distance units)
       *cubic* values = style tolerance
         style = *error* or *warn*
         tolerance = fractional difference allowed (0 <= tol <= 1)
       *shift* values = flag shiftseed
         flag = *yes* or *no* or *possible* = SRD bin shifting for better statistics
           *yes* = perform bin shifting each time SRD velocities are rescaled
           *no* = no shifting
           *possible* = shift depending on mean free path and bin size
         shiftseed = random # seed (positive integer)
       *tstat* value = *yes* or *no* = thermostat SRD particles or not
       *rescale* value = *yes* or *no* or *rotate* or *collide* = rescaling of SRD velocities
         *yes* = rescale during velocity rotation and collisions
         *no* = no rescaling
         *rotate* = rescale during velocity rotation, but not collisions
         *collide* = rescale during collisions, but not velocity rotation



Examples
""""""""


.. parsed-literal::

   fix 1 srd srd 10 big 1.0 0.25 482984
   fix 1 srd srd 10 big 0.5 0.25 482984 collision slip search 0.5

Description
"""""""""""

Treat a group of particles as stochastic rotation dynamics (SRD)
particles that serve as a background solvent when interacting with big
(colloidal) particles in groupbig-ID.  The SRD formalism is described
in :ref:`(Hecht) <Hecht>`.  The key idea behind using SRD particles as a
cheap coarse-grained solvent is that SRD particles do not interact
with each other, but only with the solute particles, which in LAMMPS
can be spheroids, ellipsoids, or line segments, or triangles, or rigid
bodies containing multiple spheroids or ellipsoids or line segments
or triangles.  The collision and rotation properties of the model
imbue the SRD particles with fluid-like properties, including an
effective viscosity.  Thus simulations with large solute particles can
be run more quickly, to measure solute properties like diffusivity
and viscosity in a background fluid.  The usual LAMMPS fixes for such
simulations, such as :doc:`fix deform <fix_deform>`, :doc:`fix viscosity <fix_viscosity>`, and :doc:`fix nvt/sllod <fix_nvt_sllod>`,
can be used in conjunction with the SRD model.

For more details on how the SRD model is implemented in LAMMPS, :ref:`this paper <Petersen1>` describes the implementation and usage of pure SRD
fluids.  :ref:`This paper <Lechman>`, which is nearly complete, describes
the implementation and usage of mixture systems (solute particles in
an SRD fluid).  See the examples/srd directory for sample input
scripts using SRD particles in both settings.

This fix does 2 things:

(1) It advects the SRD particles, performing collisions between SRD
and big particles or walls every timestep, imparting force and torque
to the big particles.  Collisions also change the position and
velocity of SRD particles.

(2) It resets the velocity distribution of SRD particles via random
rotations every N timesteps.

SRD particles have a mass, temperature, characteristic timestep
dt\_SRD, and mean free path between collisions (lamda).  The
fundamental equation relating these 4 quantities is


.. parsed-literal::

   lamda = dt_SRD \* sqrt(Kboltz \* Tsrd / mass)

The mass of SRD particles is set by the :doc:`mass <mass>` command
elsewhere in the input script.  The SRD timestep dt\_SRD is N times the
step dt defined by the :doc:`timestep <timestep>` command.  Big
particles move in the normal way via a time integration :doc:`fix <fix>`
with a short timestep dt.  SRD particles advect with a large timestep
dt\_SRD >= dt.

If the *lamda* keyword is not specified, the SRD temperature
*Tsrd* is used in the above formula to compute lamda.  If the *lamda*
keyword is specified, then the *Tsrd* setting is ignored and the above
equation is used to compute the SRD temperature.

The characteristic length scale for the SRD fluid is set by *hgrid*
which is used to bin SRD particles for purposes of resetting their
velocities.  Normally hgrid is set to be 1/4 of the big particle
diameter or smaller, to adequately resolve fluid properties around the
big particles.

Lamda cannot be smaller than 0.6 \* hgrid, else an error is generated
(unless the *shift* keyword is used, see below).  The velocities of
SRD particles are bounded by Vmax, which is set so that an SRD
particle will not advect further than Dmax = 4\*lamda in dt\_SRD.  This
means that roughly speaking, Dmax should not be larger than a big
particle diameter, else SRDs may pass through big particles without
colliding.  A warning is generated if this is the case.

Collisions between SRD particles and big particles or walls are
modeled as a lightweight SRD point particle hitting a heavy big
particle of given diameter or a wall at a point on its surface and
bouncing off with a new velocity.  The collision changes the momentum
of the SRD particle.  It imparts a force and torque to the big
particle.  It imparts a force to a wall.  Static or moving SRD walls
are setup via the :doc:`fix wall/srd <fix_wall_srd>` command.  For the
remainder of this doc page, a collision of an SRD particle with a wall
can be viewed as a collision with a big particle of infinite radius
and mass.

The *collision* keyword sets the style of collisions.  The *slip*
style means that the tangential component of the SRD particle momentum
is preserved.  Thus a force is imparted to a big particle, but no
torque.  The normal component of the new SRD velocity is sampled from
a Gaussian distribution at temperature *Tsrd*\ .

For the *noslip* style, both the normal and tangential components of
the new SRD velocity are sampled from a Gaussian distribution at
temperature *Tsrd*\ .  Additionally, a new tangential direction for the
SRD velocity is chosen randomly.  This collision style imparts torque
to a big particle.  Thus a time integrator :doc:`fix <fix>` that rotates
the big particles appropriately should be used.


----------


The *overlap* keyword should be set to *yes* if two (or more) big
particles can ever overlap.  This depends on the pair potential
interaction used for big-big interactions, or could be the case if
multiple big particles are held together as rigid bodies via the :doc:`fix rigid <fix_rigid>` command.  If the *overlap* keyword is *no* and
big particles do in fact overlap, then SRD/big collisions can generate
an error if an SRD ends up inside two (or more) big particles at once.
How this error is treated is determined by the *inside* keyword.
Running with *overlap* set to *no* allows for faster collision
checking, so it should only be set to *yes* if needed.

The *inside* keyword determines how a collision is treated if the
computation determines that the timestep started with the SRD particle
already inside a big particle.  If the setting is *error* then this
generates an error message and LAMMPS stops.  If the setting is *warn*
then this generates a warning message and the code continues.  If the
setting is *ignore* then no message is generated.  One of the output
quantities logged by the fix (see below) tallies the number of such
events, so it can be monitored.  Note that once an SRD particle is
inside a big particle, it may remain there for several steps until it
drifts outside the big particle.

The *exact* keyword determines how accurately collisions are computed.
A setting of *yes* computes the time and position of each collision as
SRD and big particles move together.  A setting of *no* estimates the
position of each collision based on the end-of-timestep positions of
the SRD and big particle.  If *overlap* is set to yes, the setting of
the *exact* keyword is ignored since time-accurate collisions are
needed.

The *radius* keyword scales the effective size of big particles.  If
big particles will overlap as they undergo dynamics, then this keyword
can be used to scale down their effective collision radius by an
amount *rfactor*\ , so that SRD particle will only collide with one big
particle at a time.  For example, in a Lennard-Jones system at a
temperature of 1.0 (in reduced LJ units), the minimum separation
between two big particles is as small as about 0.88 sigma.  Thus an
*rfactor* value of 0.85 should prevent dual collisions.

The *bounce* keyword can be used to limit the maximum number of
collisions an SRD particle undergoes in a single timestep as it
bounces between nearby big particles.  Note that if the limit is
reached, the SRD can be left inside a big particle.  A setting of 0 is
the same as no limit.


----------


There are 2 kinds of bins created and maintained when running an SRD
simulation.  The first are "SRD bins" which are used to bin SRD
particles and reset their velocities, as discussed above.  The second
are "search bins" which are used to identify SRD/big particle
collisions.

The *search* keyword can be used to choose a search bin size for
identifying SRD/big particle collisions.  The default is to use the
*hgrid* parameter for SRD bins as the search bin size.  Choosing a
smaller or large value may be more efficient, depending on the
problem.  But, in a statistical sense, it should not change the
simulation results.

The *cubic* keyword can be used to generate an error or warning when
the bin size chosen by LAMMPS creates SRD bins that are non-cubic or
different than the requested value of *hgrid* by a specified
*tolerance*\ .  Note that using non-cubic SRD bins can lead to
undetermined behavior when rotating the velocities of SRD particles,
hence LAMMPS tries to protect you from this problem.

LAMMPS attempts to set the SRD bin size to exactly *hgrid*\ .  However,
there must be an integer number of bins in each dimension of the
simulation box.  Thus the actual bin size will depend on the size and
shape of the overall simulation box.  The actual bin size is printed
as part of the SRD output when a simulation begins.

If the actual bin size in non-cubic by an amount exceeding the
tolerance, an error or warning is printed, depending on the style of
the *cubic* keyword.  Likewise, if the actual bin size differs from
the requested *hgrid* value by an amount exceeding the tolerance, then
an error or warning is printed.  The *tolerance* is a fractional
difference.  E.g. a tolerance setting of 0.01 on the shape means that
if the ratio of any 2 bin dimensions exceeds (1 +/- tolerance) then an
error or warning is generated.  Similarly, if the ratio of any bin
dimension with *hgrid* exceeds (1 +/- tolerance), then an error or
warning is generated.

.. note::

   The fix srd command can be used with simulations the size and/or
   shape of the simulation box changes.  This can be due to non-periodic
   boundary conditions or the use of fixes such as the :doc:`fix deform <fix_deform>` or :doc:`fix wall/srd <fix_wall_srd>` commands
   to impose a shear on an SRD fluid or an interaction with an external
   wall.  If the box size changes then the size of SRD bins must be
   recalculated every reneighboring.  This is not necessary if only the
   box shape changes.  This re-binning is always done so as to fit an
   integer number of bins in the current box dimension, whether it be a
   fixed, shrink-wrapped, or periodic boundary, as set by the
   :doc:`boundary <boundary>` command.  If the box size or shape changes,
   then the size of the search bins must be recalculated every
   reneighboring.  Note that changing the SRD bin size may alter the
   properties of the SRD fluid, such as its viscosity.

The *shift* keyword determines whether the coordinates of SRD
particles are randomly shifted when binned for purposes of rotating
their velocities.  When no shifting is performed, SRD particles are
binned and the velocity distribution of the set of SRD particles in
each bin is adjusted via a rotation operator.  This is a statistically
valid operation if SRD particles move sufficiently far between
successive rotations.  This is determined by their mean-free path
lamda.  If lamda is less than 0.6 of the SRD bin size, then shifting
is required.  A shift means that all of the SRD particles are shifted
by a vector whose coordinates are chosen randomly in the range [-1/2
bin size, 1/2 bin size].  Note that all particles are shifted by the
same vector.  The specified random number *shiftseed* is used to
generate these vectors.  This operation sufficiently randomizes which
SRD particles are in the same bin, even if lamda is small.

If the *shift* flag is set to *no*\ , then no shifting is performed, but
bin data will be communicated if bins overlap processor boundaries.
An error will be generated if lamda < 0.6 of the SRD bin size.  If the
*shift* flag is set to *possible*\ , then shifting is performed only if
lamda < 0.6 of the SRD bin size.  A warning is generated to let you
know this is occurring.  If the *shift* flag is set to *yes* then
shifting is performed regardless of the magnitude of lamda.  Note that
the *shiftseed* is not used if the *shift* flag is set to *no*\ , but
must still be specified.

Note that shifting of SRD coordinates requires extra communication,
hence it should not normally be enabled unless required.

The *tstat* keyword will thermostat the SRD particles to the specified
*Tsrd*\ .  This is done every N timesteps, during the velocity rotation
operation, by rescaling the thermal velocity of particles in each SRD
bin to the desired temperature.  If there is a streaming velocity
associated with the system, e.g. due to use of the :doc:`fix deform <fix_deform>` command to perform a simulation undergoing
shear, then that is also accounted for.  The mean velocity of each bin
of SRD particles is set to the position-dependent streaming velocity,
based on the coordinates of the center of the SRD bin.  Note that
collisions of SRD particles with big particles or walls has a
thermostatting effect on the colliding particles, so it may not be
necessary to thermostat the SRD particles on a bin by bin basis in
that case.  Also note that for streaming simulations, if no
thermostatting is performed (the default), then it may take a long
time for the SRD fluid to come to equilibrium with a velocity profile
that matches the simulation box deformation.

The *rescale* keyword enables rescaling of an SRD particle's velocity
if it would travel more than 4 mean-free paths in an SRD timestep.  If
an SRD particle exceeds this velocity it is possible it will be lost
when migrating to other processors or that collisions with big
particles will be missed, either of which will generate errors.  Thus
the safest mode is to run with rescaling enabled.  However rescaling
removes kinetic energy from the system (the particle's velocity is
reduced).  The latter will not typically be a problem if
thermostatting is enabled via the *tstat* keyword or if SRD collisions
with big particles or walls effectively thermostat the system.  If you
wish to turn off rescaling (on is the default), e.g. for a pure SRD
system with no thermostatting so that the temperature does not decline
over time, the *rescale* keyword can be used.  The *no* value turns
rescaling off during collisions and the per-bin velocity rotation
operation.  The *collide* and *rotate* values turn it on for
one of the operations and off for the other.


----------


.. note::

   This fix is normally used for simulations with a huge number of
   SRD particles relative to the number of big particles, e.g. 100 to 1.
   In this scenario, computations that involve only big particles
   (neighbor list creation, communication, time integration) can slow
   down dramatically due to the large number of background SRD particles.

Three other input script commands will largely overcome this effect,
speeding up an SRD simulation by a significant amount.  These are the
:doc:`atom_modify first <atom_modify>`, :doc:`neigh_modify include <neigh_modify>`, and :doc:`comm_modify group <comm_modify>`
commands.  Each takes a group-ID as an argument, which in this case
should be the group-ID of the big solute particles.

Additionally, when a :doc:`pair_style <pair_style>` for big/big particle
interactions is specified, the :doc:`pair_coeff <pair_coeff>` command
should be used to turn off big/SRD interactions, e.g. by setting their
epsilon or cutoff length to 0.0.

The "delete\_atoms overlap" command may be useful in setting up an SRD
simulation to insure there are no initial overlaps between big and SRD
particles.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix tabulates several SRD statistics which are stored in a vector
of length 12, which can be accessed by various :doc:`output commands <Howto_output>`.  The vector values calculated by this fix
are "intensive", meaning they do not scale with the size of the
simulation.  Technically, the first 8 do scale with the size of the
simulation, but treating them as intensive means they are not scaled
when printed as part of thermodynamic output.

These are the 12 quantities.  All are values for the current timestep,
except for quantity 5 and the last three, each of which are
cumulative quantities since the beginning of the run.

* (1) # of SRD/big collision checks performed
* (2) # of SRDs which had a collision
* (3) # of SRD/big collisions (including multiple bounces)
* (4) # of SRD particles inside a big particle
* (5) # of SRD particles whose velocity was rescaled to be < Vmax
* (6) # of bins for collision searching
* (7) # of bins for SRD velocity rotation
* (8) # of bins in which SRD temperature was computed
* (9) SRD temperature
* (10) # of SRD particles which have undergone max # of bounces
* (11) max # of bounces any SRD particle has had in a single step
* (12) # of reneighborings due to SRD particles moving too far

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the SRD
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`fix wall/srd <fix_wall_srd>`

Default
"""""""

The option defaults are lamda inferred from Tsrd, collision = noslip,
overlap = no, inside = error, exact = yes, radius = 1.0, bounce = 0,
search = hgrid, cubic = error 0.01, shift = no, tstat = no, and
rescale = yes.


----------


.. _Hecht:



**(Hecht)** Hecht, Harting, Ihle, Herrmann, Phys Rev E, 72, 011408 (2005).

.. _Petersen1:



**(Petersen)** Petersen, Lechman, Plimpton, Grest, in' t Veld, Schunk, J
Chem Phys, 132, 174106 (2010).

.. _Lechman:



**(Lechman)** Lechman, et al, in preparation (2010).
