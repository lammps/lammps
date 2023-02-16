Run multiple simulations from one input script
==============================================

This can be done in several ways.  See the documentation for
individual commands for more details on how these examples work.

If "multiple simulations" means to continue a previous simulation for
more timesteps, then you simply use the :doc:`run <run>` command
multiple times.  For example, this script

.. code-block:: LAMMPS

   units lj
   atom_style atomic
   read_data data.lj
   run 10000
   run 10000
   run 10000
   run 10000
   run 10000

would run 5 successive simulations of the same system for a total of
50,000 timesteps.

If you wish to run totally different simulations, one after the other,
the :doc:`clear <clear>` command can be used in between them to
re-initialize LAMMPS.  For example, this script

.. code-block:: LAMMPS

   units lj
   atom_style atomic
   read_data data.lj
   run 10000
   clear
   units lj
   atom_style atomic
   read_data data.lj.new
   run 10000

would run 2 independent simulations, one after the other.

For large numbers of independent simulations, you can use
:doc:`variables <variable>` and the :doc:`next <next>` and
:doc:`jump <jump>` commands to loop over the same input script
multiple times with different settings.  For example, this
script, named in.polymer

.. code-block:: LAMMPS

   variable d index run1 run2 run3 run4 run5 run6 run7 run8
   shell cd $d
   read_data data.polymer
   run 10000
   shell cd ..
   clear
   next d
   jump in.polymer

would run 8 simulations in different directories, using a data.polymer
file in each directory.  The same concept could be used to run the
same system at 8 different temperatures, using a temperature variable
and storing the output in different log and dump files, for example

.. code-block:: LAMMPS

   variable a loop 8
   variable t index 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15
   log log.$a
   read data.polymer
   velocity all create $t 352839
   fix 1 all nvt $t $t 100.0
   dump 1 all atom 1000 dump.$a
   run 100000
   clear
   next t
   next a
   jump in.polymer

All of the above examples work whether you are running on 1 or
multiple processors, but assumed you are running LAMMPS on a single
partition of processors.  LAMMPS can be run on multiple partitions via
the :doc:`-partition command-line switch <Run_options>`.

In the last 2 examples, if LAMMPS were run on 3 partitions, the same
scripts could be used if the "index" and "loop" variables were
replaced with *universe*\ -style variables, as described in the
:doc:`variable <variable>` command.  Also, the "next t" and "next a"
commands would need to be replaced with a single "next a t" command.
With these modifications, the 8 simulations of each script would run
on the 3 partitions one after the other until all were finished.
Initially, 3 simulations would be started simultaneously, one on each
partition.  When one finished, that partition would then start
the fourth simulation, and so forth, until all 8 were completed.
