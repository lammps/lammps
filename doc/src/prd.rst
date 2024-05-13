.. index:: prd

prd command
===========

Syntax
""""""

.. code-block:: LAMMPS

   prd N t_event n_dephase t_dephase t_correlate compute-ID seed keyword value ...

* N = # of timesteps to run (not including dephasing/quenching)
* t_event = timestep interval between event checks
* n_dephase = number of velocity randomizations to perform in each dephase run
* t_dephase = number of timesteps to run dynamics after each velocity randomization during dephase
* t_correlate = number of timesteps within which 2 consecutive events are considered to be correlated
* compute-ID = ID of the compute used for event detection
* random_seed = random # seed (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *min* or *temp* or *vel* or *time*

  .. parsed-literal::

       *min* values = etol ftol maxiter maxeval
         etol = stopping tolerance for energy, used in quenching
         ftol = stopping tolerance for force, used in quenching
         maxiter = max iterations of minimize, used in quenching
         maxeval = max number of force/energy evaluations, used in quenching
       *temp* value = Tdephase
         Tdephase = target temperature for velocity randomization, used in dephasing
       *vel* values = loop dist
         loop = *all* or *local* or *geom*, used in dephasing
         dist = *uniform* or *gaussian*, used in dephasing
       *time* value = *steps* or *clock*
         *steps* = simulation runs for N timesteps on each replica (default)
         *clock* = simulation runs for N timesteps across all replicas

Examples
""""""""

.. code-block:: LAMMPS

   prd 5000 100 10 10 100 1 54982
   prd 5000 100 10 10 100 1 54982 min 0.1 0.1 100 200

Description
"""""""""""

Run a parallel replica dynamics (PRD) simulation using multiple
replicas of a system.  One or more replicas can be used.  The total
number of steps *N* to run can be interpreted in one of two ways; see
discussion of the *time* keyword below.

PRD is described in :ref:`(Voter1998) <Voter1998>` by Art Voter.
Similar to global or local hyperdynamics (HD), PRD is a method for
performing accelerated dynamics that is suitable for infrequent-event
systems that obey first-order kinetics.  A good overview of
accelerated dynamics methods (AMD) for such systems in given in this
review paper :ref:`(Voter2002) <Voter2002prd>` from Art's group.  To
quote from the paper: "The dynamical evolution is characterized by
vibrational excursions within a potential basin, punctuated by
occasional transitions between basins.  The transition probability is
characterized by p(t) = k\*exp(-kt) where k is the rate constant."

Both PRD and HD produce a time-accurate trajectory that effectively
extends the timescale over which a system can be simulated, but they
do it differently.  PRD creates Nr replicas of the system and runs
dynamics on each independently with a normal unbiased potential until
an event occurs in one of the replicas.  The time between events is
reduced by a factor of Nr replicas.  HD uses a single replica of the
system and accelerates time by biasing the interaction potential in a
manner such that each timestep is effectively longer.  For both
methods, per CPU second, more physical time elapses and more events
occur.  See the :doc:`hyper <hyper>` page for more info about HD.

In PRD, each replica runs on a partition of one or more processors.
Processor partitions are defined at run-time using the :doc:`-partition command-line switch <Run_options>`.  Note that if you have MPI
installed, you can run a multi-replica simulation with more replicas
(partitions) than you have physical processors, e.g you can run a
10-replica simulation on one or two processors.  However for PRD, this
makes little sense, since running a replica on virtual instead of
physical processors,offers no effective parallel speed-up in searching
for infrequent events.  See the :doc:`Howto replica <Howto_replica>` doc
page for further discussion.

When a PRD simulation is performed, it is assumed that each replica is
running the same model, though LAMMPS does not check for this.
I.e. the simulation domain, the number of atoms, the interaction
potentials, etc should be the same for every replica.

A PRD run has several stages, which are repeated each time an "event"
occurs in one of the replicas, as explained below.  The logic for a
PRD run is as follows:

.. parsed-literal::

   while (time remains):
     dephase for n_dephase\*t_dephase steps
     until (event occurs on some replica):
       run dynamics for t_event steps
       quench
       check for uncorrelated event on any replica
     until (no correlated event occurs):
       run dynamics for t_correlate steps
       quench
       check for correlated event on this replica
     event replica shares state with all replicas

Before this loop begins, the state of the system on replica 0 is
shared with all replicas, so that all replicas begin from the same
initial state. The first potential energy basin is identified by
quenching (an energy minimization, see below) the initial state and
storing the resulting coordinates for reference.

In the first stage, dephasing is performed by each replica
independently to eliminate correlations between replicas.  This is
done by choosing a random set of velocities, based on the
*random_seed* that is specified, and running *t_dephase* timesteps of
dynamics.  This is repeated *n_dephase* times.  At each of the
*n_dephase* stages, if an event occurs during the *t_dephase* steps of
dynamics for a particular replica, the replica repeats the stage until
no event occurs.

If the *temp* keyword is not specified, the target temperature for
velocity randomization for each replica is the current temperature of
that replica.  Otherwise, it is the specified *Tdephase* temperature.
The style of velocity randomization is controlled using the keyword
*vel* with arguments that have the same meaning as their counterparts
in the :doc:`velocity <velocity>` command.

In the second stage, each replica runs dynamics continuously, stopping
every *t_event* steps to check if a transition event has occurred.
This check is performed by quenching the system and comparing the
resulting atom coordinates to the coordinates from the previous basin.
The first time through the PRD loop, the "previous basin" is the set
of quenched coordinates from the initial state of the system.

A quench is an energy minimization and is performed by whichever
algorithm has been defined by the :doc:`min_style <min_style>` command.
Minimization parameters may be set via the
:doc:`min_modify <min_modify>` command and by the *min* keyword of the
PRD command.  The latter are the settings that would be used with the
:doc:`minimize <minimize>` command.  Note that typically, you do not
need to perform a highly-converged minimization to detect a transition
event, though you may need to in order to prevent a set of atoms in
the system from relaxing to a saddle point.

The event check is performed by a compute with the specified
*compute-ID*\ .  Currently there is only one compute that works with the
PRD command, which is the :doc:`compute event/displace <compute_event_displace>` command.  Other
event-checking computes may be added.  :doc:`Compute event/displace <compute_event_displace>` checks whether any atom in
the compute group has moved further than a specified threshold
distance.  If so, an "event" has occurred.

In the third stage, the replica on which the event occurred (event
replica) continues to run dynamics to search for correlated events.
This is done by running dynamics for *t_correlate* steps, quenching
every *t_event* steps, and checking if another event has occurred.

The first time no correlated event occurs, the final state of the
event replica is shared with all replicas, the new basin reference
coordinates are updated with the quenched state, and the outer loop
begins again. While the replica event is searching for correlated
events, all the other replicas also run dynamics and event checking
with the same schedule, but the final states are always overwritten by
the state of the event replica.

The outer loop of the pseudocode above continues until *N* steps of
dynamics have been performed.  Note that *N* only includes the
dynamics of stages 2 and 3, not the steps taken during dephasing or
the minimization iterations of quenching.  The specified *N* is
interpreted in one of two ways, depending on the *time* keyword.  If
the *time* value is *steps*, which is the default, then each replica
runs for *N* timesteps.  If the *time* value is *clock*, then the
simulation runs until *N* aggregate timesteps across all replicas have
elapsed.  This aggregate time is the "clock" time defined below, which
typically advances nearly M times faster than the timestepping on a
single replica, where M is the number of replicas.

----------

Four kinds of output can be generated during a PRD run: event
statistics, thermodynamic output by each replica, dump files, and
restart files.

When running with multiple partitions (each of which is a replica in
this case), the print-out to the screen and master log.lammps file is
limited to event statistics.  Note that if a PRD run is performed on
only a single replica then the event statistics will be intermixed
with the usual thermodynamic output discussed below.

The quantities printed each time an event occurs are the timestep, CPU
time, clock, event number, a correlation flag, the number of
coincident events, and the replica number of the chosen event.

The timestep is the usual LAMMPS timestep, except that time does not
advance during dephasing or quenches, but only during dynamics.  Note
that are two kinds of dynamics in the PRD loop listed above that
contribute to this timestepping.  The first is when all replicas are
performing independent dynamics, waiting for an event to occur.  The
second is when correlated events are being searched for, but only one
replica is running dynamics.

The CPU time is the total elapsed time on each processor, since the
start of the PRD run.

The clock is the same as the timestep except that it advances by M
steps per timestep during the first kind of dynamics when the M
replicas are running independently.  The clock advances by only 1 step
per timestep during the second kind of dynamics, when only a single
replica is checking for a correlated event.  Thus "clock" time
represents the aggregate time (in steps) that has effectively elapsed
during a PRD simulation on M replicas.  If most of the PRD run is
spent in the second stage of the loop above, searching for infrequent
events, then the clock will advance nearly M times faster than it
would if a single replica was running.  Note the clock time between
successive events should be drawn from p(t).

The event number is a counter that increments with each event, whether
it is uncorrelated or correlated.

The correlation flag will be 0 when an uncorrelated event occurs
during the second stage of the loop listed above, i.e. when all
replicas are running independently.  The correlation flag will be 1
when a correlated event occurs during the third stage of the loop
listed above, i.e. when only one replica is running dynamics.

When more than one replica detects an event at the end of the same
event check (every *t_event* steps) during the second stage, then
one of them is chosen at random.  The number of coincident events is
the number of replicas that detected an event.  Normally, this value
should be 1.  If it is often greater than 1, then either the number of
replicas is too large, or *t_event* is too large.

The replica number is the ID of the replica (from 0 to M-1) in which
the event occurred.

----------

When running on multiple partitions, LAMMPS produces additional log
files for each partition, e.g. log.lammps.0, log.lammps.1, etc.  For
the PRD command, these contain the thermodynamic output for each
replica.  You will see short runs and minimizations corresponding to
the dynamics and quench operations of the loop listed above.  The
timestep will be reset appropriately depending on whether the
operation advances time or not.

After the PRD command completes, timing statistics for the PRD run are
printed in each replica's log file, giving a breakdown of how much CPU
time was spent in each stage (dephasing, dynamics, quenching, etc).

----------

Any :doc:`dump files <dump>` defined in the input script, will be
written to during a PRD run at timesteps corresponding to both
uncorrelated and correlated events.  This means the requested dump
frequency in the :doc:`dump <dump>` command is ignored.  There will be
one dump file (per dump command) created for all partitions.

The atom coordinates of the dump snapshot are those of the minimum
energy configuration resulting from quenching following a transition
event.  The timesteps written into the dump files correspond to the
timestep at which the event occurred and NOT the clock.  A dump
snapshot corresponding to the initial minimum state used for event
detection is written to the dump file at the beginning of each PRD
run.

----------

If the :doc:`restart <restart>` command is used, a single restart file
for all the partitions is generated, which allows a PRD run to be
continued by a new input script in the usual manner.

The restart file is generated at the end of the loop listed above.  If
no correlated events are found, this means it contains a snapshot of
the system at time T + *t_correlate*, where T is the time at which the
uncorrelated event occurred.  If correlated events were found, then it
contains a snapshot of the system at time T + *t_correlate*, where T
is the time of the last correlated event.

The restart frequency specified in the :doc:`restart <restart>` command
is interpreted differently when performing a PRD run.  It does not
mean the timestep interval between restart files.  Instead it means an
event interval for uncorrelated events.  Thus a frequency of 1 means
write a restart file every time an uncorrelated event occurs.  A
frequency of 10 means write a restart file every 10th uncorrelated
event.

When an input script reads a restart file from a previous PRD run, the
new script can be run on a different number of replicas or processors.
However, it is assumed that *t_correlate* in the new PRD command is
the same as it was previously.  If not, the calculation of the "clock"
value for the first event in the new run will be slightly off.

----------

Restrictions
""""""""""""

This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

The *N* and *t_correlate* settings must be integer multiples of
*t_event*.

Runs restarted from restart file written during a PRD run will not
produce identical results due to changes in the random numbers used
for dephasing.

This command cannot be used when any fixes are defined that keep track
of elapsed time to perform time-dependent operations.  Examples
include the "ave" fixes such as :doc:`fix ave/chunk <fix_ave_chunk>`.
Also :doc:`fix dt/reset <fix_dt_reset>` and :doc:`fix deposit <fix_deposit>`.

Related commands
""""""""""""""""

:doc:`compute event/displace <compute_event_displace>`,
:doc:`min_modify <min_modify>`, :doc:`min_style <min_style>`,
:doc:`run_style <run_style>`, :doc:`minimize <minimize>`,
:doc:`velocity <velocity>`, :doc:`temper <temper>`, :doc:`neb <neb>`,
:doc:`tad <tad>`, :doc:`hyper <hyper>`

Default
"""""""

The option defaults are min = 0.1 0.1 40 50, no temp setting, vel =
geom gaussian, and time = steps.

----------

.. _Voter1998:

**(Voter1998)** Voter, Phys Rev B, 57, 13985 (1998).

.. _Voter2002prd:

**(Voter2002)** Voter, Montalenti, Germann, Annual Review of Materials
Research 32, 321 (2002).
