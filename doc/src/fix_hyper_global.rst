.. index:: fix hyper/global

fix hyper/global command
========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID hyper/global cutbond qfactor Vmax Tequil

* ID, group-ID are documented in :doc:`fix <fix>` command
* hyper/global = style name of this fix command
* cutbond = max distance at which a pair of atoms is considered bonded (distance units)
* qfactor = max strain at which bias potential goes to 0.0 (unitless)
* Vmax = height of bias potential (energy units)
* Tequil = equilibration temperature (temperature units)

Examples
""""""""


.. parsed-literal::

   fix 1 all hyper/global 1.0 0.3 0.8 300.0

Description
"""""""""""

This fix is meant to be used with the :doc:`hyper <hyper>` command to
perform a bond-boost global hyperdynamics (GHD) simulation.  The role
of this fix is to a select a single pair of atoms in the system at
each timestep to add a global bias potential to, which will alter the
dynamics of the system in a manner that effectively accelerates time.
This is in contrast to the :doc:`fix hyper/local <fix_hyper_local>`
command, which can be user to perform a local hyperdynamics (LHD)
simulation, by adding a local bias potential to multiple pairs of
atoms at each timestep.  GHD can time accelerate a small simulation
with up to a few 100 atoms.  For larger systems, LHD is needed to
achieve good time acceleration.

For a system that undergoes rare transition events, where one or more
atoms move over an energy barrier to a new potential energy basin, the
effect of the bias potential is to induce more rapid transitions.
This can lead to a dramatic speed-up in the rate at which events
occurs, without altering their relative frequencies, thus leading to
an overall increase in the elapsed real time of the simulation as
compared to running for the same number of timesteps with normal MD.
See the :doc:`hyper <hyper>` doc page for a more general discussion of
hyperdynamics and citations that explain both GHD and LHD.

The equations and logic used by this fix and described here to perform
GHD follow the description given in :ref:`(Voter2013) <Voter2013ghd>`.  The
bond-boost form of a bias potential for HD is due to Miron and
Fichthorn as described in :ref:`(Miron) <Mironghd>`.  In LAMMPS we use a
simplified version of bond-boost GHD where a single bond in the system
is biased at any one timestep.

Bonds are defined between each pair of I,J atoms whose R0ij distance
is less than *cutbond*\ , when the system is in a quenched state
(minimum) energy.  Note that these are not "bonds" in a covalent
sense.  A bond is simply any pair of atoms that meet the distance
criterion.  *Cutbond* is an argument to this fix; it is discussed
below.  A bond is only formed if one or both of the I.J atoms are in
the specified group.

The current strain of bond IJ (when running dynamics) is defined as


.. parsed-literal::

   Eij = (Rij - R0ij) / R0ij

where Rij is the current distance between atoms I,J, and R0ij is the
equilibrium distance in the quenched state.

The bias energy Vij of any bond IJ is defined as


.. parsed-literal::

   Vij = Vmax \* (1 - (Eij/q)\^2) for abs(Eij) < qfactor
       = 0 otherwise

where the prefactor *Vmax* and the cutoff *qfactor* are arguments to
this fix; they are discussed below.  This functional form is an
inverse parabola centered at 0.0 with height Vmax and which goes to
0.0 at +/- qfactor.

Let Emax = the maximum of abs(Eij) for all IJ bonds in the system on a
given timestep.  On that step, Vij is added as a bias potential to
only the single bond with strain Emax, call it Vij(max).  Note that
Vij(max) will be 0.0 if Emax >= qfactor on that timestep.  Also note
that Vij(max) is added to the normal interatomic potential that is
computed between all atoms in the system at every step.

The derivative of Vij(max) with respect to the position of each atom
in the Emax bond gives a bias force Fij(max) acting on the bond as


.. parsed-literal::

   Fij(max) = - dVij(max)/dEij = 2 Vmax Eij / qfactor\^2   for abs(Eij) < qfactor
            = 0 otherwise

which can be decomposed into an equal and opposite force acting on
only the two I,J atoms in the Emax bond.

The time boost factor for the system is given each timestep I by


.. parsed-literal::

   Bi = exp(beta \* Vij(max))

where beta = 1/kTequil, and *Tequil* is the temperature of the system
and an argument to this fix.  Note that Bi >= 1 at every step.

.. note::

   To run a GHD simulation, the input script must also use the :doc:`fix langevin <fix_langevin>` command to thermostat the atoms at the
   same *Tequil* as specified by this fix, so that the system is running
   constant-temperature (NVT) dynamics.  LAMMPS does not check that this
   is done.

The elapsed time t\_hyper for a GHD simulation running for *N*
timesteps is simply


.. parsed-literal::

   t_hyper = Sum (i = 1 to N) Bi \* dt

where dt is the timestep size defined by the :doc:`timestep <timestep>`
command.  The effective time acceleration due to GHD is thus t\_hyper /
N\*dt, where N\*dt is elapsed time for a normal MD run of N timesteps.

Note that in GHD, the boost factor varies from timestep to timestep.
Likewise, which bond has Emax strain and thus which pair of atoms the
bias potential is added to, will also vary from timestep to timestep.
This is in contrast to local hyperdynamics (LHD) where the boost
factor is an input parameter; see the :doc:`fix hyper/local <fix_hyper_local>` doc page for details.


----------


Here is additional information on the input parameters for GHD.

The *cutbond* argument is the cutoff distance for defining bonds
between pairs of nearby atoms.  A pair of I,J atoms in their
equilibrium, minimum-energy configuration, which are separated by a
distance Rij < *cutbond*\ , are flagged as a bonded pair.  Setting
*cubond* to be ~25% larger than the nearest-neighbor distance in a
crystalline lattice is a typical choice for solids, so that bonds
exist only between nearest neighbor pairs.

The *qfactor* argument is the limiting strain at which the bias
potential goes to 0.0.  It is dimensionless, so a value of 0.3 means a
bond distance can be up to 30% larger or 30% smaller than the
equilibrium (quenched) R0ij distance and the two atoms in the bond
could still experience a non-zero bias force.

If *qfactor* is set too large, then transitions from one energy basin
to another are affected because the bias potential is non-zero at the
transition state (e.g. saddle point).  If *qfactor* is set too small
than little boost is achieved because the Eij strain of some bond in
the system will (nearly) always exceed *qfactor*\ .  A value of 0.3 for
*qfactor* is typically reasonable.

The *Vmax* argument is the prefactor on the bias potential.  Ideally,
tt should be set to a value slightly less than the smallest barrier
height for an event to occur.  Otherwise the applied bias potential
may be large enough (when added to the interatomic potential) to
produce a local energy basin with a maxima in the center.  This can
produce artificial energy minima in the same basin that trap an atom.
Or if *Vmax* is even larger, it may induce an atom(s) to rapidly
transition to another energy basin.  Both cases are "bad dynamics"
which violate the assumptions of GHD that guarantee an accelerated
time-accurate trajectory of the system.

Note that if *Vmax* is set too small, the GHD simulation will run
correctly.  There will just be fewer events because the hyper time
(t\_hyper equation above) will be shorter.

.. note::

   If you have no physical intuition as to the smallest barrier
   height in your system, a reasonable strategy to determine the largest
   *Vmax* you can use for a GHD model, is to run a sequence of
   simulations with smaller and smaller *Vmax* values, until the event
   rate does not change (as a function of hyper time).

The *Tequil* argument is the temperature at which the system is
simulated; see the comment above about the :doc:`fix langevin <fix_langevin>` thermostatting.  It is also part of the
beta term in the exponential factor that determines how much boost is
achieved as a function of the bias potential.

In general, the lower the value of *Tequil* and the higher the value
of *Vmax*\ , the more time boost will be achievable by the GHD
algorithm.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy of the bias potential to the system's
potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes a global scalar and global vector of length 12, which
can be accessed by various :doc:`output commands <Howto_output>`.  The
scalar is the magnitude of the bias potential (energy units) applied on
the current timestep.  The vector stores the following quantities:

* 1 = boost factor on this step (unitless)
* 2 = max strain Eij of any bond on this step (absolute value, unitless)
* 3 = ID of first atom in the max-strain bond
* 4 = ID of second atom in the max-strain bond
* 5 = average # of bonds/atom on this step

* 6 = fraction of timesteps where the biased bond has bias = 0.0 during this run
* 7 = fraction of timesteps where the biased bond has negative strain during this run
* 8 = max drift distance of any atom during this run (distance units)
* 9 = max bond length during this run (distance units)

* 10 = cumulative hyper time since fix was defined (time units)
* 11 = cumulative count of event timesteps since fix was defined
* 12 = cumulative count of atoms in events since fix was defined

The first 5 quantities are for the current timestep.  Quantities 6-9
are for the current hyper run.  They are reset each time a new hyper
run is performed.  Quantities 19-12 are cumulative across multiple
runs (since the point in the input script the fix was defined).

For value 8, drift is the distance an atom moves between two quenched
states when the second quench determines an event has occurred.  Atoms
involved in an event will typically move the greatest distance since
others typically remain near their original quenched position.

For value 11, events are checked for by the :doc:`hyper <hyper>` command
once every *Nevent* timesteps.  This value is the count of those
timesteps on which one (or more) events was detected.  It is NOT the
number of distinct events, since more than one event may occur in the
same *Nevent* time window.

For value 12, each time the :doc:`hyper <hyper>` command checks for an
event, it invokes a compute to flag zero or more atoms as
participating in one or more events.  E.g. atoms that have displaced
more than some distance from the previous quench state.  Value 11 is
the cumulative count of the number of atoms participating in any of
the events that were found.

The scalar and vector values calculated by this fix are all
"intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`hyper <hyper>`, :doc:`fix hyper/local <fix_hyper_local>`

**Default:** None


----------


.. _Voter2013ghd:



**(Voter2013)** S. Y. Kim, D. Perez, A. F. Voter, J Chem Phys, 139,
144110 (2013).

.. _Mironghd:



**(Miron)** R. A. Miron and K. A. Fichthorn, J Chem Phys, 119, 6210 (2003).
