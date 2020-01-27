.. index:: fix halt

fix halt command
================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID halt N attribute operator avalue keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* halt = style name of this fix command
* N = check halt condition every N steps
* attribute = *bondmax* or *tlimit* or v\_name
  
  .. parsed-literal::
  
       bondmax = length of longest bond in the system
       tlimit = elapsed CPU time
       v_name = name of :doc:`equal-style variable <variable>`

* operator = "<" or "<=" or ">" or ">=" or "==" or "!=" or "\|\^"
* avalue = numeric value to compare attribute to
* zero or more keyword/value pairs may be appended
* keyword = *error* or *message*
  
  .. parsed-literal::
  
       *error* value = *hard* or *soft* or *continue*
       *message* value = *yes* or *no*



Examples
""""""""


.. parsed-literal::

   fix 10 all halt 1 bondmax > 1.5
   fix 10 all print 10 v_myCheck != 0 error soft

Description
"""""""""""

Check a condition every N steps during a simulation run.  N must be >=
1.  If the condition is met, exit the run immediately.  In this
context a "run" can be dynamics or minimization iterations, as
specified by the :doc:`run <run>` or :doc:`minimize <minimize>` command.

The specified group-ID is ignored by this fix.

The specified *attribute* can be one of the options listed above,
namely *bondmax* or *tlimit*\ , or an :doc:`equal-style variable <variable>` referenced as *v\_name*, where "name" is the
name of a variable that has been defined previously in the input
script.

The *bondmax* attribute will loop over all bonds in the system,
compute their current lengths, and set *attribute* to the longest bond
distance.

The *tlimit* attribute queries the elapsed CPU time (in seconds) since
the current run began, and sets *attribute* to that value.  This is an
alternative way to limit the length of a simulation run, similar to
the :doc:`timer <timer>` timeout command.  There are two differences in
using this method versus the timer command option.  The first is that
the clock starts at the beginning of the current run (not when the
timer or fix command is specified), so that any setup time for the run
is not included in the elapsed time.  The second is that the timer
invocation and syncing across all processors (via MPI\_Allreduce) is
not performed once every *N* steps by this command.  Instead it is
performed (typically) only a small number of times and the elapsed
times are used to predict when the end-of-the-run will be.  Both of
these attributes can be useful when performing benchmark calculations
for a desired length of time with minimal overhead.  For example, if
a run is performing 1000s of timesteps/sec, the overhead for syncing
the timer frequently across a large number of processors may be
non-negligible.

Equal-style variables evaluate to a numeric value.  See the
:doc:`variable <variable>` command for a description.  They calculate
formulas which can involve mathematical operations, atom properties,
group properties, thermodynamic properties, global values calculated
by a :doc:`compute <compute>` or :doc:`fix <fix>`, or references to other
:doc:`variables <variable>`.  Thus they are a very general means of
computing some attribute of the current system.  For example, the
following "bondmax" variable will calculate the same quantity as the
hstyle = bondmax option.


.. parsed-literal::

   compute         bdist all bond/local dist
   compute         bmax all reduce max c_bdist
   variable        bondmax equal c_bmax

Thus these two versions of a fix halt command will do the same thing:


.. parsed-literal::

   fix 10 all halt 1 bondmax > 1.5
   fix 10 all halt 1 v_bondmax > 1.5

The version with "bondmax" will just run somewhat faster, due to less
overhead in computing bond lengths and not storing them in a separate
compute.

The choice of operators listed above are the usual comparison
operators.  The XOR operation (exclusive or) is also included as "\|\^".
In this context, XOR means that if either the attribute or avalue is
0.0 and the other is non-zero, then the result is "true".  Otherwise
it is "false".

The specified *avalue* must be a numeric value.


----------


The optional *error* keyword determines how the current run is halted.
If its value is *hard*\ , then LAMMPS will stop with an error message.

If its value is *soft*\ , LAMMPS will exit the current run, but continue
to execute subsequent commands in the input script.  However,
additional :doc:`run <run>` or :doc:`minimize <minimize>` commands will be
skipped.  For example, this allows a script to output the current
state of the system, e.g. via a :doc:`write_dump <write_dump>` or
:doc:`write_restart <write_restart>` command.

If its value is *continue*\ , the behavior is the same as for *soft*\ ,
except subsequent :doc:`run <run>` or :doc:`minimize <minimize>` commands
are executed.  This allows your script to remedy the condition that
triggered the halt, if necessary.  Note that you may wish use the
:doc:`unfix <unfix>` command on the fix halt ID, so that the same
condition is not immediately triggered in a subsequent run.

The optional *message* keyword determines whether a message is printed
to the screen and logfile when the halt condition is triggered.  If
*message* is set to yes, a one line message with the values that
triggered the halt is printed.  If *message* is set to no, no message
is printed; the run simply exits.  The latter may be desirable for
post-processing tools that extract thermodynamic information from log
files.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`variable <variable>`

Default
"""""""

The option defaults are error = hard and message = yes.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
