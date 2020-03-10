.. index:: tad

tad command
===========

Syntax
""""""

.. parsed-literal::

   tad N t_event T_lo T_hi delta tmax compute-ID keyword value ...

* N = # of timesteps to run (not including dephasing/quenching)
* t\_event = timestep interval between event checks
* T\_lo = temperature at which event times are desired
* T\_hi = temperature at which MD simulation is performed
* delta = desired confidence level for stopping criterion
* tmax = reciprocal of lowest expected pre-exponential factor (time units)
* compute-ID = ID of the compute used for event detection
* zero or more keyword/value pairs may be appended
* keyword = *min* or *neb* or *min\_style* or *neb\_style* or *neb\_log*

  .. parsed-literal::

       *min* values = etol ftol maxiter maxeval
         etol = stopping tolerance for energy (energy units)
         ftol = stopping tolerance for force (force units)
         maxiter = max iterations of minimize
         maxeval = max number of force/energy evaluations
       *neb* values = ftol N1 N2 Nevery
         etol = stopping tolerance for energy (energy units)
         ftol = stopping tolerance for force (force units)
         N1 = max # of iterations (timesteps) to run initial NEB
         N2 = max # of iterations (timesteps) to run barrier-climbing NEB
         Nevery = print NEB statistics every this many timesteps
       *neb_style* value = *quickmin* or *fire*
       *neb_step* value = dtneb
         dtneb = timestep for NEB damped dynamics minimization
       *neb_log* value = file where NEB statistics are printed

Examples
""""""""

.. parsed-literal::

   tad 2000 50 1800 2300 0.01 0.01 event
   tad 2000 50 1800 2300 0.01 0.01 event &
       min 1e-05 1e-05 100 100 &
       neb 0.0 0.01 200 200 20 &
       min_style cg &
       neb_style fire &
       neb_log log.neb

Description
"""""""""""

Run a temperature accelerated dynamics (TAD) simulation. This method
requires two or more partitions to perform NEB transition state
searches.

TAD is described in :ref:`this paper <Voter2000>` by Art Voter.  It is a method
that uses accelerated dynamics at an elevated temperature to generate
results at a specified lower temperature.  A good overview of
accelerated dynamics methods for such systems is given in :ref:`this review paper <Voter2002>` from the same group. In general, these methods assume
that the long-time dynamics is dominated by infrequent events i.e. the
system is confined to low energy basins for long periods,
punctuated by brief, randomly-occurring transitions to adjacent
basins.  TAD is suitable for infrequent-event systems, where in
addition, the transition kinetics are well-approximated by harmonic
transition state theory (hTST). In hTST, the temperature dependence of
transition rates follows the Arrhenius relation.  As a consequence a
set of event times generated in a high-temperature simulation can be
mapped to a set of much longer estimated times in the low-temperature
system. However, because this mapping involves the energy barrier of
the transition event, which is different for each event, the first
event at the high temperature may not be the earliest event at the low
temperature. TAD handles this by first generating a set of possible
events from the current basin. After each event, the simulation is
reflected backwards into the current basin.  This is repeated until
the stopping criterion is satisfied, at which point the event with the
earliest low-temperature occurrence time is selected.  The stopping
criterion is that the confidence measure be greater than
1-\ *delta*\ . The confidence measure is the probability that no earlier
low-temperature event will occur at some later time in the
high-temperature simulation.  hTST provides an lower bound for this
probability, based on the user-specified minimum pre-exponential
factor (reciprocal of *tmax*\ ).

In order to estimate the energy barrier for each event, the TAD method
invokes the :doc:`NEB <neb>` method. Each NEB replica runs on a
partition of processors. The current NEB implementation in LAMMPS
restricts you to having exactly one processor per replica. For more
information, see the documentation for the :doc:`neb <neb>` command.  In
the current LAMMPS implementation of TAD, all the non-NEB TAD
operations are performed on the first partition, while the other
partitions remain idle. See the :doc:`Howto replica <Howto_replica>` doc
page for further discussion of multi-replica simulations.

A TAD run has several stages, which are repeated each time an event is
performed.  The logic for a TAD run is as follows:

.. parsed-literal::

   while (time remains):
     while (time < tstop):
       until (event occurs):
         run dynamics for t_event steps
         quench
       run neb calculation using all replicas
       compute tlo from energy barrier
       update earliest event
       update tstop
       reflect back into current basin
     execute earliest event

Before this outer loop begins, the initial potential energy basin is
identified by quenching (an energy minimization, see below) the
initial state and storing the resulting coordinates for reference.

Inside the inner loop, dynamics is run continuously according to
whatever integrator has been specified by the user, stopping every
*t\_event* steps to check if a transition event has occurred.  This
check is performed by quenching the system and comparing the resulting
atom coordinates to the coordinates from the previous basin.

A quench is an energy minimization and is performed by whichever
algorithm has been defined by the :doc:`min_style <min_style>` command;
its default is the CG minimizer.  The tolerances and limits for each
quench can be set by the *min* keyword.  Note that typically, you do
not need to perform a highly-converged minimization to detect a
transition event.

The event check is performed by a compute with the specified
*compute-ID*\ .  Currently there is only one compute that works with the
TAD command, which is the :doc:`compute event/displace <compute_event_displace>` command.  Other
event-checking computes may be added.  :doc:`Compute event/displace <compute_event_displace>` checks whether any atom in
the compute group has moved further than a specified threshold
distance.  If so, an "event" has occurred.

The NEB calculation is similar to that invoked by the :doc:`neb <neb>`
command, except that the final state is generated internally, instead
of being read in from a file.  The style of minimization performed by
NEB is determined by the *neb\_style* keyword and must be a damped
dynamics minimizer.  The tolerances and limits for each NEB
calculation can be set by the *neb* keyword.  As discussed on the
:doc:`neb <neb>`, it is often advantageous to use a larger timestep for
NEB than for normal dynamics.  Since the size of the timestep set by
the :doc:`timestep <timestep>` command is used by TAD for performing
dynamics, there is a *neb\_step* keyword which can be used to set a
larger timestep for each NEB calculation if desired.

----------

A key aspect of the TAD method is setting the stopping criterion
appropriately.  If this criterion is too conservative, then many
events must be generated before one is finally executed.  Conversely,
if this criterion is too aggressive, high-entropy high-barrier events
will be over-sampled, while low-entropy low-barrier events will be
under-sampled. If the lowest pre-exponential factor is known fairly
accurately, then it can be used to estimate *tmax*\ , and the value of
*delta* can be set to the desired confidence level e.g. *delta* = 0.05
corresponds to 95% confidence. However, for systems where the dynamics
are not well characterized (the most common case), it will be
necessary to experiment with the values of *delta* and *tmax* to get a
good trade-off between accuracy and performance.

A second key aspect is the choice of *t\_hi*. A larger value greatly
increases the rate at which new events are generated.  However, too
large a value introduces errors due to anharmonicity (not accounted
for within hTST). Once again, for any given system, experimentation is
necessary to determine the best value of *t\_hi*.

----------

Five kinds of output can be generated during a TAD run: event
statistics, NEB statistics, thermodynamic output by each replica, dump
files, and restart files.

Event statistics are printed to the screen and master log.lammps file
each time an event is executed. The quantities are the timestep, CPU
time, global event number *N*\ , local event number *M*\ , event status,
energy barrier, time margin, *t\_lo* and *delt\_lo*.  The timestep is
the usual LAMMPS timestep, which corresponds to the high-temperature
time at which the event was detected, in units of timestep.  The CPU
time is the total processor time since the start of the TAD run.  The
global event number *N* is a counter that increments with each
executed event. The local event number *M* is a counter that resets to
zero upon entering each new basin.  The event status is *E* when an
event is executed, and is *D* for an event that is detected, while
*DF* is for a detected event that is also the earliest (first) event
at the low temperature.

The time margin is the ratio of the high temperature time in the
current basin to the stopping time. This last number can be used to
judge whether the stopping time is too short or too long (see above).

*t\_lo* is the low-temperature event time when the current basin was
entered, in units of timestep.  del*t\_lo* is the time of each detected
event, measured relative to *t\_lo*.  *delt\_lo* is equal to the
high-temperature time since entering the current basin, scaled by an
exponential factor that depends on the hi/lo temperature ratio and the
energy barrier for that event.

On lines for executed events, with status *E*\ , the global event number
is incremented by one,
the local event number and time margin are reset to zero,
while the global event number, energy barrier, and
*delt\_lo* match the last event with status *DF*
in the immediately preceding block of detected events.
The low-temperature event time *t\_lo* is incremented by *delt\_lo*.

NEB statistics are written to the file specified by the *neb\_log*
keyword. If the keyword value is "none", then no NEB statistics are
printed out. The statistics are written every *Nevery* timesteps.  See
the :doc:`neb <neb>` command for a full description of the NEB
statistics. When invoked from TAD, NEB statistics are never printed to
the screen.

Because the NEB calculation must run on multiple partitions, LAMMPS
produces additional screen and log files for each partition,
e.g. log.lammps.0, log.lammps.1, etc. For the TAD command, these
contain the thermodynamic output of each NEB replica. In addition, the
log file for the first partition, log.lammps.0, will contain
thermodynamic output from short runs and minimizations corresponding
to the dynamics and quench operations, as well as a line for each new
detected event, as described above.

After the TAD command completes, timing statistics for the TAD run are
printed in each replica's log file, giving a breakdown of how much CPU
time was spent in each stage (NEB, dynamics, quenching, etc).

Any :doc:`dump files <dump>` defined in the input script will be written
to during a TAD run at timesteps when an event is executed.  This
means the requested dump frequency in the :doc:`dump <dump>` command
is ignored.  There will be one dump file (per dump command) created
for all partitions.  The atom coordinates of the dump snapshot are
those of the minimum energy configuration resulting from quenching
following the executed event.  The timesteps written into the dump
files correspond to the timestep at which the event occurred and NOT
the clock.  A dump snapshot corresponding to the initial minimum state
used for event detection is written to the dump file at the beginning
of each TAD run.

If the :doc:`restart <restart>` command is used, a single restart file
for all the partitions is generated, which allows a TAD run to be
continued by a new input script in the usual manner.  The restart file
is generated after an event is executed. The restart file contains a
snapshot of the system in the new quenched state, including the event
number and the low-temperature time.  The restart frequency specified
in the :doc:`restart <restart>` command is interpreted differently when
performing a TAD run.  It does not mean the timestep interval between
restart files.  Instead it means an event interval for executed
events.  Thus a frequency of 1 means write a restart file every time
an event is executed.  A frequency of 10 means write a restart file
every 10th executed event.  When an input script reads a restart file
from a previous TAD run, the new script can be run on a different
number of replicas or processors.

Note that within a single state, the dynamics will typically
temporarily continue beyond the event that is ultimately chosen, until
the stopping criterion is satisfied.  When the event is eventually
executed, the timestep counter is reset to the value when the event
was detected. Similarly, after each quench and NEB minimization, the
timestep counter is reset to the value at the start of the
minimization. This means that the timesteps listed in the replica log
files do not always increase monotonically. However, the timestep
values printed to the master log file, dump files, and restart files
are always monotonically increasing.

----------

Restrictions
""""""""""""

This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

*N* setting must be integer multiple of *t\_event*.

Runs restarted from restart files written during a TAD run will only
produce identical results if the user-specified integrator supports
exact restarts. So :doc:`fix nvt <fix_nh>` will produce an exact
restart, but :doc:`fix langevin <fix_langevin>` will not.

This command cannot be used when any fixes are defined that keep track
of elapsed time to perform time-dependent operations.  Examples
include the "ave" fixes such as :doc:`fix ave/chunk <fix_ave_chunk>`.
Also :doc:`fix dt/reset <fix_dt_reset>` and :doc:`fix deposit <fix_deposit>`.

Related commands
""""""""""""""""

:doc:`compute event/displace <compute_event_displace>`,
:doc:`min_modify <min_modify>`, :doc:`min_style <min_style>`,
:doc:`run_style <run_style>`, :doc:`minimize <minimize>`,
:doc:`temper <temper>`, :doc:`neb <neb>`,
:doc:`prd <prd>`

Default
"""""""

The option defaults are *min* = 0.1 0.1 40 50, *neb* = 0.01 100 100
10, *neb\_style* = *quickmin*\ , *neb\_step* = the same timestep set by
the :doc:`timestep <timestep>` command, and *neb\_log* = "none".

----------

.. _Voter2000:

**(Voter2000)** Sorensen and Voter, J Chem Phys, 112, 9599 (2000)

.. _Voter2002:

**(Voter2002)** Voter, Montalenti, Germann, Annual Review of Materials
Research 32, 321 (2002).
