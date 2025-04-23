.. index:: fenix

fenix command
=============

Syntax
""""""

.. code-block:: LAMMPS
   fenix keyword value ...

* keyword = *restart_file*, *restart_label*, or *spares*

   .. parsed-literal::

      *restart_file* arg = file
         file = the file argument to pass to a jump command
      *restart_label* arg = label
         file = the label argument to pass to a jump command
      *spares* arg = N
         N = the number of ranks for Fenix to use as spares

Examples
""""""""

.. code-block:: LAMMPS

   fenix spares 10
   fenix spares 10 restart_label restart

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
rebuild the LAMMPS world to the same number of ranks as before (when possible),
using the spare ranks held at initialization. These spares take the place of
the failed ranks, leaving surviving processes with the same rank number. Fenix
then tears down LAMMPS and restarts according to whatever restart_file and
restart_label are configured. The default restart file is "SELF".

If this command is called more than once, subsequent invocations will update
arguments. Fenix is not reinitialized, so changes to the spares argument are
ignored.

----------

Restrictions
""""""""""""

This fix is part of the FENIX package. It is only enabled if LAMMPS was built
with that package. See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`restart`, :doc:`write_restart`, :doc:`read_restart`,

Default
"""""""

none
