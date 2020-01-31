.. index:: hyper

hyper command
=============

Syntax
""""""


.. parsed-literal::

   hyper N Nevent fix-ID compute-ID keyword values ...

* N = # of timesteps to run
* Nevent = check for events every this many steps
* fix-ID = ID of a fix that applies a global or local bias potential, can be NULL
* compute-ID = ID of a compute that identifies when an event has occurred
* zero or more keyword/value pairs may be appended
* keyword = *min* or *dump* or *rebond*
  
  .. parsed-literal::
  
       *min* values = etol ftol maxiter maxeval
         etol = stopping tolerance for energy, used in quenching
         ftol = stopping tolerance for force, used in quenching
         maxiter = max iterations of minimize, used in quenching
         maxeval = max number of force/energy evaluations, used in quenching
       *dump* value = dump-ID
         dump-ID = ID of dump to trigger whenever an event takes place
       *rebond* value = Nrebond
         Nrebond = frequency at which to reset bonds, even if no event has occurred



Examples
""""""""


.. parsed-literal::

   compute event all event/displace 1.0
   fix HG mobile hyper/global 3.0 0.3 0.4 800.0
   hyper 5000 100 HG event min 1.0e-6 1.0e-6 100 100 dump 1 dump 5

Description
"""""""""""

Run a bond-boost hyperdynamics (HD) simulation where time is
accelerated by application of a bias potential to one or more pairs of
nearby atoms in the system.  This command can be used to run both
global and local hyperdynamics.  In global HD a single bond within the
system is biased on each timestep.  In local HD multiple bonds
(separated by a sufficient distance) can be biased simultaneously at
each timestep.  In the bond-boost hyperdynamics context, a "bond" is
not a covalent bond between a pair of atoms in a molecule.  Rather it
is simply a pair of nearby atoms as discussed below.

Both global and local HD are described in :ref:`(Voter2013) <Voter2013>` by
Art Voter and collaborators.  Similar to parallel replica dynamics
(PRD), global and local HD are methods for performing accelerated
dynamics that are suitable for infrequent-event systems that obey
first-order kinetics.  A good overview of accelerated dynamics methods
for such systems in given in :ref:`(Voter2002) <Voter2002hd>` from the same
group.  To quote from the review paper: "The dynamical evolution is
characterized by vibrational excursions within a potential basin,
punctuated by occasional transitions between basins."  The transition
probability is characterized by p(t) = k\*exp(-kt) where k is the rate
constant.  Running multiple replicas gives an effective enhancement in
the timescale spanned by the multiple simulations, while waiting for
an event to occur.

Both HD and PRD produce a time-accurate trajectory that effectively
extends the timescale over which a system can be simulated, but they
do it differently.  HD uses a single replica of the system and
accelerates time by biasing the interaction potential in a manner such
that each timestep is effectively longer.  PRD creates Nr replicas of
the system and runs dynamics on each independently with a normal
unbiased potential until an event occurs in one of the replicas.  The
time between events is reduced by a factor of Nr replicas.  For both
methods, per CPU second, more physical time elapses and more events
occur.  See the :doc:`prd <prd>` doc page for more info about PRD.

An HD run has several stages, which are repeated each time an event
occurs, as explained below.  The logic for an HD run is as follows:


.. parsed-literal::

   quench
   create initial list of bonds

   while (time remains):
     run dynamics for Nevent steps
     quench
     check for an event
     if event occurred: reset list of bonds
     restore pre-quench state

The list of bonds is the list of atom pairs of atoms that are within a
short cutoff distance of each other after the system energy is
minimized (quenched).  This list is created and reset by a :doc:`fix hyper/global <fix_hyper_global>` or :doc:`fix hyper/local <fix_hyper_local>` command specified as *fix-ID*\ .  At
every dynamics timestep, the same fix selects one of more bonds to
apply a bias potential to.

.. note::

   The style of fix associated with the specified *fix-ID*
   determines whether you are running the global versus local
   hyperdynamics algorithm.

Dynamics (with the bias potential) is run continuously, stopping every
*Nevent* steps to check if a transition event has occurred.  The
specified *N* for total steps must be a multiple of *Nevent*\ .  check
is performed by quenching the system and comparing the resulting atom
coordinates to the coordinates from the previous basin.

A quench is an energy minimization and is performed by whichever
algorithm has been defined by the :doc:`min_style <min_style>` command.
Minimization parameters may be set via the
:doc:`min_modify <min_modify>` command and by the *min* keyword of the
hyper command.  The latter are the settings that would be used with
the :doc:`minimize <minimize>` command.  Note that typically, you do not
need to perform a highly-converged minimization to detect a transition
event, though you may need to in order to prevent a set of atoms in
the system from relaxing to a saddle point.

The event check is performed by a compute with the specified
*compute-ID*\ .  Currently there is only one compute that works with the
hyper command, which is the :doc:`compute event/displace <compute_event_displace>` command.  Other
event-checking computes may be added.  :doc:`Compute event/displace <compute_event_displace>` checks whether any atom in
the compute group has moved further than a specified threshold
distance.  If so, an event has occurred.

If this happens, the list of bonds is reset, since some bond pairs
are likely now too far apart, and new pairs are likely close enough
to be considered a bond.  The pre-quenched state of the
system (coordinates and velocities) is restored, and dynamics continue.

At the end of the hyper run, a variety of statistics are output to the
screen and logfile.  These include info relevant to both global and
local hyperdynamics, such as the number of events and the elapsed
hyper time (accelerated time), And it includes info specific to one or
the other, depending on which style of fix was specified by *fix-ID*\ .


----------


The optional keywords operate as follows.

As explained above, the *min* keyword can be used to specify
parameters for the quench.  Their meaning is the same
as for the :doc:`minimize <minimize>` command

The *dump* keyword can be used to trigger a specific dump command with
the specified *dump-ID* to output a snapshot each time an event is
detected.  It can be specified multiple times with different *dump-ID*
values, as in the example above.  These snapshots will be for the
quenched state of the system on a timestep that is a multiple of
*Nevent*\ , i.e. a timestep after the event has occurred.  Note that any
dump command in the input script will also output snapshots at
whatever timestep interval it defines via its *N* argument; see the
:doc:`dump <dump>` command for details.  This means if you only want a
particular dump to output snapshots when events are detected, you
should specify its *N* as a value larger than the length of the
hyperdynamics run.

As in the code logic above, the bond list is normally only reset when
an event occurs.  The *rebond* keyword will force a reset of the bond
list every *Nrebond* steps, even if an event has not occurred.
*Nrebond* must be a multiple of *Nevent*\ .  This can be useful to check
if more frequent resets alter event statistics, perhaps because the
parameters chosen for defining what is a bond and what is an event are
producing bad dynamics in the presence of the bias potential.


----------


Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`fix hyper/global <fix_hyper_global>`, :doc:`fix hyper/local <fix_hyper_local>`, :doc:`compute event/displace <compute_event_displace>`, :doc:`prd <prd>`

Default
"""""""

The option defaults are min = 0.1 0.1 40 50 and time = steps.


----------


.. _Voter2013:



**(Voter2013)** S. Y. Kim, D. Perez, A. F. Voter, J Chem Phys, 139,
144110 (2013).

.. _Voter2002hd:



**(Voter2002)** Voter, Montalenti, Germann, Annual Review of Materials
Research 32, 321 (2002).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
