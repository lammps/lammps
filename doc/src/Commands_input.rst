LAMMPS input scripts
====================

LAMMPS executes by reading commands from a input script (text file),
one line at a time.  When the input script ends, LAMMPS exits.  Each
command causes LAMMPS to take some action.  It may set an internal
variable, read in a file, or run a simulation.  Most commands have
default settings, which means you only need to use the command if you
wish to change the default.

In many cases, the ordering of commands in an input script is not
important.  However the following rules apply:

(1) LAMMPS does not read your entire input script and then perform a
simulation with all the settings.  Rather, the input script is read
one line at a time and each command takes effect when it is read.
Thus this sequence of commands:


.. parsed-literal::

   timestep 0.5
   run      100
   run      100

does something different than this sequence:


.. parsed-literal::

   run      100
   timestep 0.5
   run      100

In the first case, the specified timestep (0.5 fs) is used for two
simulations of 100 timesteps each.  In the 2nd case, the default
timestep (1.0 fs) is used for the 1st 100 step simulation and a 0.5 fs
timestep is used for the 2nd one.

(2) Some commands are only valid when they follow other commands.  For
example you cannot set the temperature of a group of atoms until atoms
have been defined and a group command is used to define which atoms
belong to the group.

(3) Sometimes command B will use values that can be set by command A.
This means command A must precede command B in the input script if it
is to have the desired effect.  For example, the
:doc:`read\_data <read_data>` command initializes the system by setting
up the simulation box and assigning atoms to processors.  If default
values are not desired, the :doc:`processors <processors>` and
:doc:`boundary <boundary>` commands need to be used before read\_data to
tell LAMMPS how to map processors to the simulation box.

Many input script errors are detected by LAMMPS and an ERROR or
WARNING message is printed.  The :doc:`Errors <Errors>` doc page gives
more information on what errors mean.  The documentation for each
command lists restrictions on how the command can be used.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
