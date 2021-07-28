.. index:: timer

timer command
=============

Syntax
""""""

.. parsed-literal::

   timer args

* *args* = one or more of *off* or *loop* or *normal* or *full* or *sync* or *nosync* or *timeout* or *every*

.. parsed-literal::

     *off* = do not collect or print any timing information
     *loop* = collect only the total time for the simulation loop
     *normal* = collect timer information broken down by sections (default)
     *full* = like *normal* but also include CPU and thread utilization
     *sync* = explicitly synchronize MPI tasks between sections
     *nosync* = do not synchronize MPI tasks between sections (default)
     *timeout* elapse = set wall time limit to *elapse*
     *every* Ncheck = perform timeout check every *Ncheck* steps

Examples
""""""""

.. code-block:: LAMMPS

   timer full sync
   timer timeout 2:00:00 every 100
   timer loop

Description
"""""""""""

Select the level of detail at which LAMMPS performs its CPU timings.
Multiple keywords can be specified with the *timer* command.  For
keywords that are mutually exclusive, the last one specified takes
precedence.

During a simulation run LAMMPS collects information about how much
time is spent in different sections of the code and thus can provide
information for determining performance and load imbalance problems.
This can be done at different levels of detail and accuracy.  For more
information about the timing output, see the :doc:`Run output <Run_output>` doc page.

The *off* setting will turn all time measurements off. The *loop*
setting will only measure the total time for a run and not collect any
detailed per section information.  With the *normal* setting, timing
information for portions of the timestep (pairwise calculations,
neighbor list construction, output, etc) are collected as well as
information about load imbalances for those sections across
processors.  The *full* setting adds information about CPU
utilization and thread utilization, when multi-threading is enabled.

With the *sync* setting, all MPI tasks are synchronized at each timer
call which measures load imbalance for each section more accurately,
though it can also slow down the simulation by prohibiting overlapping
independent computations on different MPI ranks  Using the *nosync*
setting (which is the default) turns this synchronization off.

With the *timeout* keyword a wall time limit can be imposed, that
affects the :doc:`run <run>` and :doc:`minimize <minimize>` commands.
This can be convenient when calculations have to comply with execution
time limits, e.g. when running under a batch system when you want to
maximize the utilization of the batch time slot, especially for runs
where the time per timestep varies much and thus it becomes difficult
to predict how many steps a simulation can perform for a given wall time
limit. This also applies for difficult to converge minimizations.
The timeout *elapse* value should be somewhat smaller than the maximum
wall time requested from the batch system, as there is usually
some overhead to launch jobs, and it is advisable to write
out a restart after terminating a run due to a timeout.

The timeout timer starts when the command is issued. When the time
limit is reached, the run or energy minimization will exit on the
next step or iteration that is a multiple of the *Ncheck* value
which can be set with the *every* keyword. Default is checking
every 10 steps. After the timer timeout has expired all subsequent
run or minimize commands in the input script will be skipped.
The remaining time or timer status can be accessed with the
:doc:`thermo <thermo_style>` variable *timeremain*, which will be
zero, if the timeout is inactive (default setting), it will be
negative, if the timeout time is expired and positive if there
is time remaining and in this case the value of the variable are
the number of seconds remaining.

When the *timeout* key word is used a second time, the timer is
restarted with a new time limit. The timeout *elapse* value can
be specified as *off* or *unlimited* to impose a no timeout condition
(which is the default).  The *elapse* setting can be specified as
a single number for seconds, two numbers separated by a colon (MM:SS)
for minutes and seconds, or as three numbers separated by colons for
hours, minutes, and seconds (H:MM:SS).

The *every* keyword sets how frequently during a run or energy
minimization the wall clock will be checked.  This check count applies
to the outer iterations or time steps during minimizations or :doc:`r-RESPA runs <run_style>`, respectively.  Checking for timeout too often,
can slow a calculation down.  Checking too infrequently can make the
timeout measurement less accurate, with the run being stopped later
than desired.

.. note::

   Using the *full* and *sync* options provides the most detailed
   and accurate timing information, but can also have a negative
   performance impact due to the overhead of the many required system
   calls. It is thus recommended to use these settings only when testing
   tests to identify performance bottlenecks. For calculations with few
   atoms or a very large number of processors, even the *normal* setting
   can have a measurable negative performance impact. In those cases you
   can just use the *loop* or *off* setting.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`run post no <run>`, :doc:`kspace_modify fftbench <kspace_modify>`

Default
"""""""

.. code-block:: LAMMPS

   timer normal nosync
   timer timeout off
   timer every 10
