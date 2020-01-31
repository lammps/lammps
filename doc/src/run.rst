.. index:: run

run command
===========

Syntax
""""""


.. parsed-literal::

   run N keyword values ...

* N = # of timesteps
* zero or more keyword/value pairs may be appended
* keyword = *upto* or *start* or *stop* or *pre* or *post* or *every*
  
  .. parsed-literal::
  
       *upto* value = none
       *start* value = N1
         N1 = timestep at which 1st run started
       *stop* value = N2
         N2 = timestep at which last run will end
       *pre* value = *no* or *yes*
       *post* value = *no* or *yes*
       *every* values = M c1 c2 ...
         M = break the run into M-timestep segments and invoke one or more commands between each segment
         c1,c2,...,cN = one or more LAMMPS commands, each enclosed in quotes
         c1 = NULL means no command will be invoked



Examples
""""""""


.. parsed-literal::

   run 10000
   run 1000000 upto
   run 100 start 0 stop 1000
   run 1000 pre no post yes
   run 100000 start 0 stop 1000000 every 1000 "print 'Protein Rg = $r'"
   run 100000 every 1000 NULL

Description
"""""""""""

Run or continue dynamics for a specified number of timesteps.

When the :doc:`run style <run_style>` is *respa*\ , N refers to outer
loop (largest) timesteps.

A value of N = 0 is acceptable; only the thermodynamics of the system
are computed and printed without taking a timestep.

The *upto* keyword means to perform a run starting at the current
timestep up to the specified timestep.  E.g. if the current timestep
is 10,000 and "run 100000 upto" is used, then an additional 90,000
timesteps will be run.  This can be useful for very long runs on a
machine that allocates chunks of time and terminate your job when time
is exceeded.  If you need to restart your script multiple times
(reading in the last restart file), you can keep restarting your
script with the same run command until the simulation finally
completes.

The *start* or *stop* keywords can be used if multiple runs are being
performed and you want a :doc:`fix <fix>` command that changes some
value over time (e.g. temperature) to make the change across the
entire set of runs and not just a single run.  See the doc page for
individual fixes to see which ones can be used with the *start/stop*
keywords.

For example, consider this fix followed by 10 run commands:


.. parsed-literal::

   fix         1 all nvt 200.0 300.0 1.0
   run         1000 start 0 stop 10000
   run         1000 start 0 stop 10000
   ...
   run         1000 start 0 stop 10000

The NVT fix ramps the target temperature from 200.0 to 300.0 during a
run.  If the run commands did not have the start/stop keywords (just
"run 1000"), then the temperature would ramp from 200.0 to 300.0
during the 1000 steps of each run.  With the start/stop keywords, the
ramping takes place over the 10000 steps of all runs together.

The *pre* and *post* keywords can be used to streamline the setup,
clean-up, and associated output to the screen that happens before and
after a run.  This can be useful if you wish to do many short runs in
succession (e.g. LAMMPS is being called as a library which is doing
other computations between successive short LAMMPS runs).

By default (pre and post = yes), LAMMPS creates neighbor lists,
computes forces, and imposes fix constraints before every run.  And
after every run it gathers and prints timings statistics.  If a run is
just a continuation of a previous run (i.e. no settings are changed),
the initial computation is not necessary; the old neighbor list is
still valid as are the forces.  So if *pre* is specified as "no" then
the initial setup is skipped, except for printing thermodynamic info.
Note that if *pre* is set to "no" for the very 1st run LAMMPS
performs, then it is overridden, since the initial setup computations
must be done.

.. note::

   If your input script changes the system between 2 runs, then the
   initial setup must be performed to insure the change is recognized by
   all parts of the code that are affected.  Examples are adding a
   :doc:`fix <fix>` or :doc:`dump <dump>` or :doc:`compute <compute>`, changing
   a :doc:`neighbor <neigh_modify>` list parameter, or writing restart file
   which can migrate atoms between processors.  LAMMPS has no easy way to
   check if this has happened, but it is an error to use the *pre no*
   option in this case.

If *post* is specified as "no", the full timing summary is skipped;
only a one-line summary timing is printed.

The *every* keyword provides a means of breaking a LAMMPS run into a
series of shorter runs.  Optionally, one or more LAMMPS commands (c1,
c2, ..., cN) will be executed in between the short runs.  If used, the
*every* keyword must be the last keyword, since it has a variable
number of arguments.  Each of the trailing arguments is a single
LAMMPS command, and each command should be enclosed in quotes, so that
the entire command will be treated as a single argument.  This will
also prevent any variables in the command from being evaluated until
it is executed multiple times during the run.  Note that if a command
itself needs one of its arguments quoted (e.g. the :doc:`print <print>`
command), then you can use a combination of single and double quotes,
as in the example above or below.

The *every* keyword is a means to avoid listing a long series of runs
and interleaving commands in your input script.  For example, a
:doc:`print <print>` command could be invoked or a :doc:`fix <fix>` could
be redefined, e.g. to reset a thermostat temperature.  Or this could
be useful for invoking a command you have added to LAMMPS that wraps
some other code (e.g. as a library) to perform a computation
periodically during a long LAMMPS run.  See the :doc:`Modify <Modify>`
doc page for info about how to add new commands to LAMMPS.  See the
:doc:`Howto couple <Howto_couple>` doc page for ideas about how to
couple LAMMPS to other codes.

With the *every* option, N total steps are simulated, in shorter runs
of M steps each.  After each M-length run, the specified commands are
invoked.  If only a single command is specified as NULL, then no
command is invoked.  Thus these lines:


.. parsed-literal::

   variable q equal x[100]
   run 6000 every 2000 "print 'Coord = $q'"

are the equivalent of:


.. parsed-literal::

   variable q equal x[100]
   run 2000
   print "Coord = $q"
   run 2000
   print "Coord = $q"
   run 2000
   print "Coord = $q"

which does 3 runs of 2000 steps and prints the x-coordinate of a
particular atom between runs.  Note that the variable "$q" will
be evaluated afresh each time the print command is executed.

Note that by using the line continuation character "&", the run every
command can be spread across many lines, though it is still a single
command:


.. parsed-literal::

   run 100000 every 1000 &
     "print 'Minimum value = $a'" &
     "print 'Maximum value = $b'" &
     "print 'Temp = $c'" &
     "print 'Press = $d'"

If the *pre* and *post* options are set to "no" when used with the
*every* keyword, then the 1st run will do the full setup and the last
run will print the full timing summary, but these operations will be
skipped for intermediate runs.

.. note::

   You might wish to specify a command that exits the run by
   jumping out of the loop, e.g.


.. parsed-literal::

   variable t equal temp
   run 10000 every 100 "if '$t < 300.0' then 'jump SELF afterrun'"

However, this will not work.  The run command simply executes each
command one at a time each time it pauses, then continues the run.

Instead, you should use the :doc:`fix halt <fix_halt>` command, which
has additional options for how to exit the run.

Restrictions
""""""""""""


When not using the *upto* keyword, the number of specified timesteps N
must fit in a signed 32-bit integer, so you are limited to slightly
more than 2 billion steps (2\^31) in a single run.  When using *upto*\ ,
N can be larger than a signed 32-bit integer, however the difference
between N and the current timestep must still be no larger than
2\^31 steps.

However, with or without the *upto* keyword, you can perform
successive runs to run a simulation for any number of steps (ok, up to
2\^63 total steps).  I.e. the timestep counter within LAMMPS is a
64-bit signed integer.

Related commands
""""""""""""""""

:doc:`minimize <minimize>`, :doc:`run_style <run_style>`,
:doc:`temper <temper>`, :doc:`fix halt <fix_halt>`

Default
"""""""

The option defaults are start = the current timestep, stop = current
timestep + N, pre = yes, and post = yes.
