Common problems
===============

If two LAMMPS runs do not produce the exact same answer on different
machines or different numbers of processors, this is typically not a
bug.  In theory you should get identical answers on any number of
processors and on any machine.  In practice, numerical round-off can
cause slight differences and eventual divergence of molecular dynamics
phase space trajectories within a few 100s or few 1000s of timesteps.
However, the statistical properties of the two runs (e.g. average
energy or temperature) should still be the same.

If the :doc:`velocity <velocity>` command is used to set initial atom
velocities, a particular atom can be assigned a different velocity
when the problem is run on a different number of processors or on
different machines.  If this happens, the phase space trajectories of
the two simulations will rapidly diverge.  See the discussion of the
*loop* option in the :doc:`velocity <velocity>` command for details and
options that avoid this issue.

Similarly, the :doc:`create_atoms <create_atoms>` command generates a
lattice of atoms.  For the same physical system, the ordering and
numbering of atoms by atom ID may be different depending on the number
of processors.

Some commands use random number generators which may be setup to
produce different random number streams on each processor and hence
will produce different effects when run on different numbers of
processors.  A commonly-used example is the :doc:`fix langevin <fix_langevin>` command for thermostatting.

A LAMMPS simulation typically has two stages, setup and run.  Most
LAMMPS errors are detected at setup time; others like a bond
stretching too far may not occur until the middle of a run.

LAMMPS tries to flag errors and print informative error messages so
you can fix the problem.  For most errors it will also print the last
input script command that it was processing.  Of course, LAMMPS cannot
figure out your physics or numerical mistakes, like choosing too big a
timestep, specifying erroneous force field coefficients, or putting 2
atoms on top of each other!  If you run into errors that LAMMPS
does not catch that you think it should flag, please send an email to
the `developers <https://www.lammps.org/authors.html>`_.

If you get an error message about an invalid command in your input
script, you can determine what command is causing the problem by
looking in the log.lammps file or using the :doc:`echo command <echo>`
to see it on the screen.  If you get an error like "Invalid ...
style", with ... being fix, compute, pair, etc, it means that you
mistyped the style name or that the command is part of an optional
package which was not compiled into your executable.  The list of
available styles in your executable can be listed by using
:doc:`the -h command-line switch <Run_options>`.  The installation and
compilation of optional packages is explained on the
:doc:`Build packages <Build_package>` doc page.

For a given command, LAMMPS expects certain arguments in a specified
order.  If you mess this up, LAMMPS will often flag the error, but it
may also simply read a bogus argument and assign a value that is
valid, but not what you wanted.  E.g. trying to read the string "abc"
as an integer value of 0.  Careful reading of the associated doc page
for the command should allow you to fix these problems. In most cases,
where LAMMPS expects to read a number, either integer or floating point,
it performs a stringent test on whether the provided input actually
is an integer or floating-point number, respectively, and reject the
input with an error message (for instance, when an integer is required,
but a floating-point number 1.0 is provided):

.. parsed-literal::

   ERROR: Expected integer parameter instead of '1.0' in input script or data file

Some commands allow for using variable references in place of numeric
constants so that the value can be evaluated and may change over the
course of a run.  This is typically done with the syntax *v_name* for a
parameter, where name is the name of the variable. On the other hand,
immediate variable expansion with the syntax ${name} is performed while
reading the input and before parsing commands,

.. note::

   Using a variable reference (i.e. *v_name*) is only allowed if
   the documentation of the corresponding command explicitly says it is.
   Otherwise, you will receive an error message of this kind:

.. parsed-literal::

   ERROR: Expected floating point parameter instead of 'v_name' in input script or data file

Generally, LAMMPS will print a message to the screen and logfile and
exit gracefully when it encounters a fatal error.  Sometimes it will
print a WARNING to the screen and logfile and continue on; you can
decide if the WARNING is important or not.  A WARNING message that is
generated in the middle of a run is only printed to the screen, not to
the logfile, to avoid cluttering up thermodynamic output.  If LAMMPS
crashes or hangs without spitting out an error message first then it
could be a bug (see :doc:`this section <Errors_bugs>`) or one of the following
cases:

LAMMPS runs in the available memory a processor allows to be
allocated.  Most reasonable MD runs are compute limited, not memory
limited, so this should not be a bottleneck on most platforms.  Almost
all large memory allocations in the code are done via C-style malloc's
which will generate an error message if you run out of memory.
Smaller chunks of memory are allocated via C++ "new" statements.  If
you are unlucky you could run out of memory just when one of these
small requests is made, in which case the code will crash or hang (in
parallel), since LAMMPS does not trap on those errors.

Illegal arithmetic can cause LAMMPS to run slow or crash.  This is
typically due to invalid physics and numerics that your simulation is
computing.  If you see wild thermodynamic values or NaN values in your
LAMMPS output, something is wrong with your simulation.  If you
suspect this is happening, it is a good idea to print out
thermodynamic info frequently (e.g. every timestep) via the
:doc:`thermo <thermo>` so you can monitor what is happening.
Visualizing the atom movement is also a good idea to ensure your model
is behaving as you expect.

In parallel, one way LAMMPS can hang is due to how different MPI
implementations handle buffering of messages.  If the code hangs
without an error message, it may be that you need to specify an MPI
setting or two (usually via an environment variable) to enable
buffering or boost the sizes of messages that can be buffered.
