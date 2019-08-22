.. index:: fix langevin

fix langevin command
====================

fix langevin/kk command
=======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID langevin Tstart Tstop damp seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin = style name of this fix command
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* Tstart can be a variable (see below)
* damp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *angmom* or *omega* or *scale* or *tally* or *zero*
  
  .. parsed-literal::
  
       *angmom* value = *no* or factor
         *no* = do not thermostat rotational degrees of freedom via the angular momentum
         factor = do thermostat rotational degrees of freedom via the angular momentum and apply numeric scale factor as discussed below
       *gjf* value = *no* or *yes*
         *no* = use standard formulation
         *yes* = use Gronbech-Jensen/Farago formulation
       *omega* value = *no* or *yes*
         *no* = do not thermostat rotational degrees of freedom via the angular velocity
         *yes* = do thermostat rotational degrees of freedom via the angular velocity
       *scale* values = type ratio
         type = atom type (1-N)
         ratio = factor by which to scale the damping coefficient
       *tally* value = *no* or *yes*
         *no* = do not tally the energy added/subtracted to atoms
         *yes* = do tally the energy added/subtracted to atoms
       *zero* value = *no* or *yes*
         *no* = do not set total random force to zero
         *yes* = set total random force to zero



Examples
""""""""


.. parsed-literal::

   fix 3 boundary langevin 1.0 1.0 1000.0 699483
   fix 1 all langevin 1.0 1.1 100.0 48279 scale 3 1.5
   fix 1 all langevin 1.0 1.1 100.0 48279 angmom 3.333

Description
"""""""""""

Apply a Langevin thermostat as described in :ref:`(Schneider) <Schneider1>`
to a group of atoms which models an interaction with a background
implicit solvent.  Used with :doc:`fix nve <fix_nve>`, this command
performs Brownian dynamics (BD), since the total force on each atom
will have the form:


.. parsed-literal::

   F = Fc + Ff + Fr
   Ff = - (m / damp) v
   Fr is proportional to sqrt(Kb T m / (dt damp))

Fc is the conservative force computed via the usual inter-particle
interactions (:doc:`pair\_style <pair_style>`,
:doc:`bond\_style <bond_style>`, etc).

The Ff and Fr terms are added by this fix on a per-particle basis.
See the :doc:`pair\_style dpd/tstat <pair_dpd>` command for a
thermostatting option that adds similar terms on a pairwise basis to
pairs of interacting particles.

Ff is a frictional drag or viscous damping term proportional to the
particle's velocity.  The proportionality constant for each atom is
computed as m/damp, where m is the mass of the particle and damp is
the damping factor specified by the user.

Fr is a force due to solvent atoms at a temperature T randomly bumping
into the particle.  As derived from the fluctuation/dissipation
theorem, its magnitude as shown above is proportional to sqrt(Kb T m /
dt damp), where Kb is the Boltzmann constant, T is the desired
temperature, m is the mass of the particle, dt is the timestep size,
and damp is the damping factor.  Random numbers are used to randomize
the direction and magnitude of this force as described in
:ref:`(Dunweg) <Dunweg1>`, where a uniform random number is used (instead of
a Gaussian random number) for speed.

Note that unless you use the *omega* or *angmom* keywords, the
thermostat effect of this fix is applied to only the translational
degrees of freedom for the particles, which is an important
consideration for finite-size particles, which have rotational degrees
of freedom, are being thermostatted.  The translational degrees of
freedom can also have a bias velocity removed from them before
thermostatting takes place; see the description below.

.. note::

   Unlike the :doc:`fix nvt <fix_nh>` command which performs
   Nose/Hoover thermostatting AND time integration, this fix does NOT
   perform time integration.  It only modifies forces to effect
   thermostatting.  Thus you must use a separate time integration fix,
   like :doc:`fix nve <fix_nve>` to actually update the velocities and
   positions of atoms using the modified forces.  Likewise, this fix
   should not normally be used on atoms that also have their temperature
   controlled by another fix - e.g. by :doc:`fix nvt <fix_nh>` or :doc:`fix temp/rescale <fix_temp_rescale>` commands.

See the :doc:`Howto thermostat <Howto_thermostat>` doc page for
a discussion of different ways to compute temperature and perform
thermostatting.

The desired temperature at each timestep is a ramped value during the
run from *Tstart* to *Tstop*\ .

*Tstart* can be specified as an equal-style or atom-style
:doc:`variable <variable>`.  In this case, the *Tstop* setting is
ignored.  If the value is a variable, it should be specified as
v\_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
target temperature.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo\_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent temperature.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent
temperature with optional time-dependence as well.

Like other fixes that perform thermostatting, this fix can be used
with :doc:`compute commands <compute>` that remove a "bias" from the
atom velocities.  E.g. removing the center-of-mass velocity from a
group of atoms or removing the x-component of velocity from the
calculation.  This is not done by default, but only if the
:doc:`fix\_modify <fix_modify>` command is used to assign a temperature
compute to this fix that includes such a bias term.  See the doc pages
for individual :doc:`compute commands <compute>` to determine which ones
include a bias.  In this case, the thermostat works in the following
manner: bias is removed from each atom, thermostatting is performed on
the remaining thermal degrees of freedom, and the bias is added back
in.

The *damp* parameter is specified in time units and determines how
rapidly the temperature is relaxed.  For example, a value of 100.0
means to relax the temperature in a timespan of (roughly) 100 time
units (tau or fmsec or psec - see the :doc:`units <units>` command).
The damp factor can be thought of as inversely related to the
viscosity of the solvent.  I.e. a small relaxation time implies a
hi-viscosity solvent and vice versa.  See the discussion about gamma
and viscosity in the documentation for the :doc:`fix viscous <fix_viscous>` command for more details.

The random # *seed* must be a positive integer.  A Marsaglia random
number generator is used.  Each processor uses the input seed to
generate its own unique seed and its own stream of random numbers.
Thus the dynamics of the system will not be identical on two runs on
different numbers of processors.


----------


The keyword/value option pairs are used in the following ways.

The keyword *angmom* and *omega* keywords enable thermostatting of
rotational degrees of freedom in addition to the usual translational
degrees of freedom.  This can only be done for finite-size particles.

A simulation using atom\_style sphere defines an omega for finite-size
spheres.  A simulation using atom\_style ellipsoid defines a finite
size and shape for aspherical particles and an angular momentum.
The Langevin formulas for thermostatting the rotational degrees of
freedom are the same as those above, where force is replaced by
torque, m is replaced by the moment of inertia I, and v is replaced by
omega (which is derived from the angular momentum in the case of
aspherical particles).

The rotational temperature of the particles can be monitored by the
:doc:`compute temp/sphere <compute_temp_sphere>` and :doc:`compute temp/asphere <compute_temp_asphere>` commands with their rotate
options.

For the *omega* keyword there is also a scale factor of 10.0/3.0 that
is applied as a multiplier on the Ff (damping) term in the equation
above and of sqrt(10.0/3.0) as a multiplier on the Fr term.  This does
not affect the thermostatting behavior of the Langevin formalism but
insures that the randomized rotational diffusivity of spherical
particles is correct.

For the *angmom* keyword a similar scale factor is needed which is
10.0/3.0 for spherical particles, but is anisotropic for aspherical
particles (e.g. ellipsoids).  Currently LAMMPS only applies an
isotropic scale factor, and you can choose its magnitude as the
specified value of the *angmom* keyword.  If your aspherical particles
are (nearly) spherical than a value of 10.0/3.0 = 3.333 is a good
choice.  If they are highly aspherical, a value of 1.0 is as good a
choice as any, since the effects on rotational diffusivity of the
particles will be incorrect regardless.  Note that for any reasonable
scale factor, the thermostatting effect of the *angmom* keyword on the
rotational temperature of the aspherical particles should still be
valid.

The keyword *scale* allows the damp factor to be scaled up or down by
the specified factor for atoms of that type.  This can be useful when
different atom types have different sizes or masses.  It can be used
multiple times to adjust damp for several atom types.  Note that
specifying a ratio of 2 increases the relaxation time which is
equivalent to the solvent's viscosity acting on particles with 1/2 the
diameter.  This is the opposite effect of scale factors used by the
:doc:`fix viscous <fix_viscous>` command, since the damp factor in fix
*langevin* is inversely related to the gamma factor in fix *viscous*\ .
Also note that the damping factor in fix *langevin* includes the
particle mass in Ff, unlike fix *viscous*\ .  Thus the mass and size of
different atom types should be accounted for in the choice of ratio
values.

The keyword *tally* enables the calculation of the cumulative energy
added/subtracted to the atoms as they are thermostatted.  Effectively
it is the energy exchanged between the infinite thermal reservoir and
the particles.  As described below, this energy can then be printed
out or added to the potential energy of the system to monitor energy
conservation.

The keyword *zero* can be used to eliminate drift due to the
thermostat. Because the random forces on different atoms are
independent, they do not sum exactly to zero.  As a result, this fix
applies a small random force to the entire system, and the
center-of-mass of the system undergoes a slow random walk.  If the
keyword *zero* is set to *yes*\ , the total random force is set exactly
to zero by subtracting off an equal part of it from each atom in the
group.  As a result, the center-of-mass of a system with zero initial
momentum will not drift over time.

The keyword *gjf* can be used to run the :ref:`Gronbech-Jensen/Farago <Gronbech-Jensen>` time-discretization of the Langevin model.  As
described in the papers cited below, the purpose of this method is to
enable longer timesteps to be used (up to the numerical stability
limit of the integrator), while still producing the correct Boltzmann
distribution of atom positions.  It is implemented within LAMMPS, by
changing how the random force is applied so that it is composed of
the average of two random forces representing half-contributions from
the previous and current time intervals.

In common with all methods based on Verlet integration, the
discretized velocities generated by this method in conjunction with
velocity-Verlet time integration are not exactly conjugate to the
positions.  As a result the temperature (computed from the discretized
velocities) will be systematically lower than the target temperature,
by a small amount which grows with the timestep.  Nonetheless, the
distribution of atom positions will still be consistent with the
target temperature.

As an example of using the *gjf* keyword, for molecules containing C-H
bonds, configurational properties generated with dt = 2.5 fs and tdamp
= 100 fs are indistinguishable from dt = 0.5 fs.  Because the velocity
distribution systematically decreases with increasing timestep, the
method should not be used to generate properties that depend on the
velocity distribution, such as the velocity auto-correlation function
(VACF). In this example, the velocity distribution at dt = 2.5fs
generates an average temperature of 220 K, instead of 300 K.


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


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  Because the state of the random number generator
is not saved in restart files, this means you cannot do "exact"
restarts with this fix, where the simulation continues on the same as
if no restart had taken place.  However, in a statistical sense, a
restarted simulation should produce the same behavior.

The :doc:`fix\_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a temperature :doc:`compute <compute>`
you have defined to this fix which will be used in its thermostatting
procedure, as described above.  For consistency, the group used by
this fix and by the compute should be the same.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Langevin thermostatting to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.  Note that use of this option requires
setting the *tally* keyword to *yes*\ .

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".  Note that calculation of this quantity requires
setting the *tally* keyword to *yes*\ .

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix viscous <fix_viscous>`, :doc:`fix nvt <fix_nh>`, :doc:`pair\_style dpd/tstat <pair_dpd>`

Default
"""""""

The option defaults are angmom = no, omega = no, scale = 1.0 for all
types, tally = no, zero = no, gjf = no.


----------


.. _Dunweg1:



**(Dunweg)** Dunweg and Paul, Int J of Modern Physics C, 2, 817-27 (1991).

.. _Schneider1:



**(Schneider)** Schneider and Stoll, Phys Rev B, 17, 1302 (1978).

.. _Gronbech-Jensen:



**(Gronbech-Jensen)** Gronbech-Jensen and Farago, Mol Phys, 111, 983
(2013); Gronbech-Jensen, Hayre, and Farago, Comp Phys Comm,
185, 524 (2014)


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
