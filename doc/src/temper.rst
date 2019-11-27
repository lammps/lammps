.. index:: temper

temper command
==============

Syntax
""""""


.. parsed-literal::

   temper N M temp fix-ID seed1 seed2 index

* N = total # of timesteps to run
* M = attempt a tempering swap every this many steps
* temp = initial temperature for this ensemble
* fix-ID = ID of the fix that will control temperature during the run
* seed1 = random # seed used to decide on adjacent temperature to partner with
* seed2 = random # seed for Boltzmann factor in Metropolis swap
* index = which temperature (0 to N-1) I am simulating (optional)

Examples
""""""""


.. parsed-literal::

   temper 100000 100 $t tempfix 0 58728
   temper 40000 100 $t tempfix 0 32285 $w

Description
"""""""""""

Run a parallel tempering or replica exchange simulation using multiple
replicas (ensembles) of a system.  Two or more replicas must be used.

Each replica runs on a partition of one or more processors.  Processor
partitions are defined at run-time using the :doc:`-partition command-line switch <Run_options>`.  Note that if you have MPI installed, you
can run a multi-replica simulation with more replicas (partitions)
than you have physical processors, e.g you can run a 10-replica
simulation on one or two processors.  You will simply not get the
performance speed-up you would see with one or more physical
processors per replica.  See the :doc:`Howto replica <Howto_replica>`
doc page for further discussion.

Each replica's temperature is controlled at a different value by a fix
with *fix-ID* that controls temperature. Most thermostat fix styles
(with and without included time integration) are supported. The command
will print an error message and abort, if the chosen fix is unsupported.
The desired temperature is specified by *temp*\ , which is typically a
variable previously set in the input script, so that each partition is
assigned a different temperature.  See the :doc:`variable <variable>`
command for more details.  For example:


.. parsed-literal::

   variable t world 300.0 310.0 320.0 330.0
   fix myfix all nvt temp $t $t 100.0
   temper 100000 100 $t myfix 3847 58382

would define 4 temperatures, and assign one of them to the thermostat
used by each replica, and to the temper command.

As the tempering simulation runs for *N* timesteps, a temperature swap
between adjacent ensembles will be attempted every *M* timesteps.  If
*seed1* is 0, then the swap attempts will alternate between odd and
even pairings.  If *seed1* is non-zero then it is used as a seed in a
random number generator to randomly choose an odd or even pairing each
time.  Each attempted swap of temperatures is either accepted or
rejected based on a Boltzmann-weighted Metropolis criterion which uses
*seed2* in the random number generator.

As a tempering run proceeds, multiple log files and screen output
files are created, one per replica.  By default these files are named
log.lammps.M and screen.M where M is the replica number from 0 to N-1,
with N = # of replicas.  See the :doc:`-log and -screen command-line swiches <Run_options>` for info on how to change these names.

The main screen and log file (log.lammps) will list information about
which temperature is assigned to each replica at each thermodynamic
output timestep.  E.g. for a simulation with 16 replicas:


.. parsed-literal::

   Running on 16 partitions of processors
   Step T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15
   0    0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
   500  1 0 3 2 5 4 6 7 8 9 10 11 12 13 14 15
   1000 2 0 4 1 5 3 6 7 8 9 10 11 12 14 13 15
   1500 2 1 4 0 5 3 6 7 9 8 10 11 12 14 13 15
   2000 2 1 3 0 6 4 5 7 10 8 9 11 12 14 13 15
   2500 2 1 3 0 6 4 5 7 11 8 9 10 12 14 13 15
   ...

The column headings T0 to TN-1 mean which temperature is currently
assigned to the replica 0 to N-1.  Thus the columns represent replicas
and the value in each column is its temperature (also numbered 0 to
N-1).  For example, a 0 in the 4th column (column T3, step 2500) means
that the 4th replica is assigned temperature 0, i.e. the lowest
temperature.  You can verify this time sequence of temperature
assignments for the Nth replica by comparing the Nth column of screen
output to the thermodynamic data in the corresponding log.lammps.N or
screen.N files as time proceeds.

You can have each replica create its own dump file in the following
manner:


.. parsed-literal::

   variable rep world 0 1 2 3 4 5 6 7
   dump 1 all atom 1000 dump.temper.$\ *rep*

.. note::

   Each replica's dump file will contain a continuous trajectory
   for its atoms where the temperature varies over time as swaps take
   place involving that replica.  If you want a series of dump files,
   each with snapshots (from all replicas) that are all at a single
   temperature, then you will need to post-process the dump files using
   the information from the log.lammps file.  E.g. you could produce one
   dump file with snapshots at 300K (from all replicas), another with
   snapshots at 310K, etc.  Note that these new dump files will not
   contain "continuous trajectories" for individual atoms, because two
   successive snapshots (in time) may be from different replicas. The 
   reorder\_remd\_traj python script can do the reordering for you 
   (and additionally also calculated configurational log-weights of 
   trajectory snapshots in the canonical ensemble). The script can be found
   in the tools/replica directory while instructions on how to use it is
   available in doc/Tools (in brief) and as a README file in tools/replica
   (in detail).

The last argument *index* in the temper command is optional and is
used when restarting a tempering run from a set of restart files (one
for each replica) which had previously swapped to new temperatures.
The *index* value (from 0 to N-1, where N is the # of replicas)
identifies which temperature the replica was simulating on the
timestep the restart files were written.  Obviously, this argument
must be a variable so that each partition has the correct value.  Set
the variable to the *N* values listed in the log file for the previous
run for the replica temperatures at that timestep.  For example if the
log file listed the following for a simulation with 5 replicas:


.. parsed-literal::

   500000 2 4 0 1 3

then a setting of


.. parsed-literal::

   variable w world 2 4 0 1 3

would be used to restart the run with a tempering command like the
example above with $w as the last argument.


----------


Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`variable <variable>`, :doc:`prd <prd>`, :doc:`neb <neb>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
