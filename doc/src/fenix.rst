.. index:: fenix

fenix command
=============

Syntax
""""""

.. code-block:: LAMMPS

   fenix keyword value ...

* keyword = *restart_file*, *restart_label*, *universal*, or *spares*

   .. parsed-literal::

      *restart_file* arg = file
         file = the file argument to pass to a jump command
      *restart_label* arg = label
         label = the label argument to pass to a jump command
      *spares* arg = N
         N = the number of ranks for Fenix to use as spares
      *universal* arg = none
         initialize Fenix with universal scope

Examples
""""""""

.. code-block:: LAMMPS

   fenix spares 10
   fenix spares 10 restart_label restart
   fenix universal restart_label restart

----------

Description
"""""""""""

.. versionadded:: 2Sep2026

This command initializes `Fenix <https://github.com/sandialabs/fenix>`_
for online process recovery.  From this point on, Fenix claims a number
of ranks to be used as spares.  The LAMMPS world is rebuilt without the
spare ranks.  Importantly, this means Fenix should be initialized before
you perform any steps that rely on knowing the number of ranks or
communicating.  In general, the safest approach is to only have
:doc:`variable <variable>` setup commands before calling this command.

On detecting a process failure, Fenix will automatically regain control.
It will rebuild the LAMMPS world to the same number of ranks as before
(when possible), using the spare ranks held at initialization. These
spares take the place of the failed ranks, leaving surviving processes
with the same rank number.  Fenix then tears down LAMMPS and restarts
according to whatever restart_file and restart_label are configured.
The default restart file is "SELF".  Upon restart, the `fenix_restarted`
internal variable is defined and set to 1.

If this command is called more than once, subsequent invocations will
update arguments.  Fenix is not reinitialized, so changes to the spares
argument are ignored.

By default, Fenix initializes assuming LAMMPS is run on a :doc:`single
partition <Run_options>`.  This is the best option for runs that do not
involve partitions, or runs with partitions that run independently and
thus have no communication between them.  For runs with partitions that
communicate, the *universal* flag is needed to ensure the universe's
uworld is always valid to communicate on.

When the *universal* flag is not used, communication on the universe's
uworld communicator may cause MPI to abort if any processes have failed.
When the *universal* flag is used, you must not specify a number of
spare ranks - instead, all ranks in the last partition are used as
spares.  When the *universal* flag is used, shrinking recovery
(recovering when all spare ranks have been consumed) is not currently
supported.

----------

Restrictions
""""""""""""

This command is part of the FENIX package. It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`restart <restart>`, :doc:`write restart <write_restart>`,
:doc:`read restart <read_restart>`

Default
"""""""

The default *restart_file* = SELF.
