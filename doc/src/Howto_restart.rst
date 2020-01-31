Restart a simulation
====================

There are 3 ways to continue a long LAMMPS simulation.  Multiple
:doc:`run <run>` commands can be used in the same input script.  Each
run will continue from where the previous run left off.  Or binary
restart files can be saved to disk using the :doc:`restart <restart>`
command.  At a later time, these binary files can be read via a
:doc:`read_restart <read_restart>` command in a new script.  Or they can
be converted to text data files using the :doc:`-r command-line switch <Run_options>` and read by a :doc:`read_data <read_data>`
command in a new script.

Here we give examples of 2 scripts that read either a binary restart
file or a converted data file and then issue a new run command to
continue where the previous run left off.  They illustrate what
settings must be made in the new script.  Details are discussed in the
documentation for the :doc:`read_restart <read_restart>` and
:doc:`read_data <read_data>` commands.

Look at the *in.chain* input script provided in the *bench* directory
of the LAMMPS distribution to see the original script that these 2
scripts are based on.  If that script had the line


.. parsed-literal::

   restart         50 tmp.restart

added to it, it would produce 2 binary restart files (tmp.restart.50
and tmp.restart.100) as it ran.

This script could be used to read the 1st restart file and re-run the
last 50 timesteps:


.. parsed-literal::

   read_restart    tmp.restart.50

   neighbor        0.4 bin
   neigh_modify    every 1 delay 1

   fix             1 all nve
   fix             2 all langevin 1.0 1.0 10.0 904297

   timestep        0.012

   run             50

Note that the following commands do not need to be repeated because
their settings are included in the restart file: *units, atom\_style,
special\_bonds, pair\_style, bond\_style*.  However these commands do
need to be used, since their settings are not in the restart file:
*neighbor, fix, timestep*\ .

If you actually use this script to perform a restarted run, you will
notice that the thermodynamic data match at step 50 (if you also put a
"thermo 50" command in the original script), but do not match at step
100.  This is because the :doc:`fix langevin <fix_langevin>` command
uses random numbers in a way that does not allow for perfect restarts.

As an alternate approach, the restart file could be converted to a data
file as follows:


.. parsed-literal::

   lmp_g++ -r tmp.restart.50 tmp.restart.data

Then, this script could be used to re-run the last 50 steps:


.. parsed-literal::

   units           lj
   atom_style      bond
   pair_style      lj/cut 1.12
   pair_modify     shift yes
   bond_style      fene
   special_bonds   0.0 1.0 1.0

   read_data       tmp.restart.data

   neighbor        0.4 bin
   neigh_modify    every 1 delay 1

   fix             1 all nve
   fix             2 all langevin 1.0 1.0 10.0 904297

   timestep        0.012

   reset_timestep  50
   run             50

Note that nearly all the settings specified in the original *in.chain*
script must be repeated, except the *pair\_coeff* and *bond\_coeff*
commands since the new data file lists the force field coefficients.
Also, the :doc:`reset_timestep <reset_timestep>` command is used to tell
LAMMPS the current timestep.  This value is stored in restart files,
but not in data files.


