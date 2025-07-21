.. index:: fix fail

fix fail command
================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID fail keyword args

* keyword = *rank* or *timestep* or *step* or *var*

   .. parsed-literal::

      *rank* arg = fail_rank
         fail_rank = rank(s) allowed to fail. This can be "none", "all", or any
           combination of CSV and ranges
      *timestep* arg = fail_timestep
         fail_timestep = simulation timestep required for a rank to fail
      *step* arg = fail_step
         fail_step = name of the simulation step during which a rank may fail
      *var* args = fail_var fail_var_val
         fail_var = name of an equal-style variable to check when deciding to fail
         fail_var_val = value that fail_var must equal for a rank to fail
      *wait_only* args = none

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all fail rank 1 timestep 10
   fix 2 all fail rank 1,2,3 timestep 10
   fix 2 all fail rank 1-3 timestep 10
   fix 2 all fail rank 1-2,3 timestep 10
   fix 2 all fail var should_fail 1 step *
   fix 2 all fail rank ALL timestep 10 wait_only

Description
"""""""""""

This command is used to test error recovering by artificially generating failures.
You must at least specify either a variable or both a rank and timestep. If a variable
is specified, any rank or timestep information passed will further limit the possibly
generated failures.

Step can be an asterisk to indicate that failure requirements should be checked at
every step of the solver or the name of a function to use for checking for failures
(pre_force, post_neighbor, etc.). Multiple steps can be passed separated by an
ampersand (pre_force&post_neighbor).

If wait_only is given, the identified ranks will not inject any failures. Instead, they
will busy-loop on an MPI_Barrier until a failure occurs. This allows users to more easily
use their own error injection strategies.


Restrictions
""""""""""""

This fix is part of the FENIX package. It is only enabled if LAMMPS was built
with that package. See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""
:doc:`kill` :doc:`fenix`

Default
"""""""

step = end_of_step
