OPT package
===========

The OPT package was developed by James Fischer (High Performance
Technologies), David Richie, and Vincent Natoli (Stone Ridge
Technologies).  It contains a handful of pair styles whose compute()
methods were rewritten in C++ templated form to reduce the overhead
due to if tests and other conditional code.

**Required hardware/software:**

None.

**Building LAMMPS with the OPT package:**

See the :ref:`Build extras <opt>` doc page for instructions.

**Run with the OPT package from the command line:**

.. code-block:: bash

   lmp_mpi -sf opt -in in.script                # run in serial
   mpirun -np 4 lmp_mpi -sf opt -in in.script   # run in parallel

Use the "-sf opt" :doc:`command-line switch <Run_options>`, which will
automatically append "opt" to styles that support it.

**Or run with the OPT package by editing an input script:**

Use the :doc:`suffix opt <suffix>` command, or you can explicitly add an
"opt" suffix to individual styles in your input script, e.g.

.. code-block:: LAMMPS

   pair_style lj/cut/opt 2.5

**Speed-ups to expect:**

You should see a reduction in the "Pair time" value printed at the end
of a run.  On most machines for reasonable problem sizes, it will be a
5 to 20% savings.

**Guidelines for best performance:**

Just try out an OPT pair style to see how it performs.

Restrictions
""""""""""""

None.
