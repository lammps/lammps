.. index:: fenix

fenix command
=============

Syntax
""""""

.. code-block:: LAMMPS
   fenix keyword value ...

* keyword = *checkpoint_every*, *restart_jump*, or *spares*

   .. parsed-literal::

      *checkpoint_every* arg = N
         N = make in-memory checkpoints on timesteps which are multiples of N (except 0)
      *restart_jump* arg = file label
         file = the file argument to pass to a jump command
         label = the label argument to pass to a jump command
      *restart_file* arg = filename
         filename = the argument to pass to a read_restart command
      *spares* arg = N
         N = the number of ranks for Fenix to use as spares

Examples
""""""""

.. code-block:: LAMMPS

   fenix checkpoint_every 50 spares 10
   fenix spares 10 restart_jump SELF restart

----------

Description
"""""""""""

This command initializes Fenix for online process recovery. From this point on,
Fenix claims a number of ranks to be used as spares. The LAMMPS world is rebuilt
without the spare ranks. Importantly, this means Fenix should be initialized
before you perform any steps that rely on knowing the number of ranks or
communicating. In general, the safest approach is to only have variable setup
commands before calling this command.

On detecting a process failure, Fenix will automatically regain control. It will
rebuild the LAMMPS world to the same number of ranks as before, using the spare
ranks held at initialization. These spares take the place of the failed ranks,
leaving surviving processes with the same rank number. Fenix then tears down
LAMMPS and restarts the file from the beginning.

When this command is reached again on subsequent passes through the file, it
will not repeat the Fenix initialization (which is still valid). Instead, this
command recovers from the checkpointed data if checkpoint_every or restart_file
were specified. This recovery follows roughly the same semantics as read_restart
- meaning some commands in your input script should not be repeated, while
others should.

.. note::

   Memory-based recovery is not precisely the same as file-based. Some steps are
   taken to remove unnecessary communication during checkpoint and restart.

   During checkpoints, the value of comm->me is set to 0, to encourage all ranks
   to save each fix's relevant checkpoint data to their local checkpoints. This
   works well for most fixes, but some fixes use communication during checkpoint
   writes. These few fixes will likely need to be modified - feel free to reach
   out to the Fenix team for help with this if needed.

   Memory-based checkpoints do not currently support the write_restart_file
   function. Again, feel free to reach out to the Fenix team if you need support
   to be implemented.

   Memory-based recovery restores the rank-to-grid assignments that existed
   before the failure - this is important to avoid the need to move atom data
   after each recovery. Since recovery always has the same number of ranks as
   checkpointed, memory-based recovery also skips rebalancing.

   Memory-based recovery uses the POSIX function fmemopen, which may not be
   portable to all operating systems (e.g. definitely not Windows).

----------

If restart_jump is specified, this command will issue a jump command when it is
reached on subsequent passes through the file. This is useful for skipping
commands that should not be repeated (e.g. commands initializing your input
state) or for running commands that should only run for recovery (e.g.
read_restart, if using file based checkpoints instead of Fenix).

----------

Restrictions
""""""""""""

This fix is part of the FENIX package. It is only enabled if LAMMPS was built
with that package. See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`fix fail`, :doc:`kill`,
:doc:`restart`, :doc:`write_restart`, :doc:`read_restart`,

Default
"""""""

none
